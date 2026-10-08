#!/usr/bin/env python3
"""The band-window report, the unusual-Bandlist warning and the Exciton_cutoff default.

Bandlist is a LIST of band offsets from the Fermi level (0 = top valence, 1 = bottom conduction), not a range:
'-1 2' selects two bands and skips both frontier bands. Such a window is legal but easy to write by mistake, and it
cuts degenerate pairs (on MoS2 it gave a symmetry-violating rectification background 30-500x the shift current).
OptiX prints the resolved window and warns when it omits offset 0 or 1, has gaps, or repeats an offset.

Checks (MoS2 tb model, 6x6 mesh, linear response only, a few seconds; hBN Xatu files for the .states path):
  '-1 0 1 2' and the full window: window printed, no warning;  '-1 2': warning naming 0 and 1 and the gap;
  '-2 0 1': warning naming the missing -1;  '0 1 3': missing 2;  '0 1 1': repeated offset;
  Xatu input (bands taken from hBN_N30.states, offsets 0 1): no warning.
Exciton_cutoff is optional: without it every exciton in the Xatu files is used (the count is the third line of
the .eigval file; the first line of the .states file is the basis size, not the count). A negative value stops.
  no Exciton_cutoff: 'using all 900 excitons' on hBN_N30;  Exciton_cutoff = -5: stops with a message.

Usage: python3 tools/check_bandlist_guard.py --opticx bin/opticx --root . --workdir bin/check_bandlist_guard
"""
import argparse, os, subprocess, sys

SP = """# Periodic dimensions
2
# Wannier90_filename
{tb}
# Xatu_interface
false
# Ncells
6
# Bandlist
{bands}
# Nfermi
26
# OME_sp
linear
# Response
absorbance
# Energy_variables
1.0 4.0 0.1 31
"""
EX = """# Periodic dimensions
2
# Wannier90_filename
{tb}
# Xatu_interface
true
{eig}
{sta}
{cutoff}
# Nfermi
1
# OME_sp
linear
# OME_ex
linear
# Response
absorbance
# Energy_variables
4.0 8.0 0.1 21
"""
WARN = 'WARNING (Bandlist)'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--opticx', required=True); ap.add_argument('--root', default='.')
    ap.add_argument('--workdir', default='check_bandlist_guard')
    a = ap.parse_args()
    root, work, opticx = os.path.abspath(a.root), os.path.abspath(a.workdir), os.path.abspath(a.opticx)
    mos2 = os.path.join(root, 'MoS2_spin_wannier_07032024_tb.dat')
    hbn, eig, sta = (os.path.join(root, f) for f in ('hBN_tb.dat', 'hBN_N30.eigval', 'hBN_N30.states'))
    for f in (mos2, hbn, eig, sta):
        if not os.path.exists(f):
            print('SKIP: missing ' + f); sys.exit(0)
    os.makedirs(work, exist_ok=True)
    nfail = 0

    def check(label, ok, msg=''):
        nonlocal nfail
        print(('  PASS   ' if ok else '  FAIL   ') + label + ('   ' + msg if msg else '')); nfail += 0 if ok else 1

    def run(tag, text):
        d = os.path.join(work, tag); os.makedirs(d, exist_ok=True)
        name = f'{tag}.txt'
        open(os.path.join(d, name), 'w').write(text)
        r = subprocess.run([opticx, name], cwd=d, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log = r.stdout.decode('utf-8', 'replace'); open(os.path.join(d, 'run.log'), 'w').write(log)
        for f in os.listdir(d):
            if f.endswith('.omesp') or f.endswith('.omeex'):
                os.remove(os.path.join(d, f))
        return r.returncode, log

    print('Band-window report and unusual-Bandlist warning\n')
    full = ' '.join(str(k) for k in range(-25, 9))
    for tag, bands, expect in (('MoS2_bandlist_4', '-1 0 1 2', []), ('MoS2_bandlist_all', full, []),
                               ('MoS2_bandlist_m1_2', '-1 2', ['offset 0', 'offset 1', 'missing between the window edges: 0 1']),
                               ('MoS2_bandlist_gapv', '-2 0 1', ['missing between the window edges: -1']),
                               ('MoS2_bandlist_gapc', '0 1 3', ['missing between the window edges: 2']),
                               ('MoS2_bandlist_dup', '0 1 1', ['listed more than once'])):
        rc, log = run(tag, SP.format(tb=mos2, bands=bands))
        shown = 'Band window (offset from Nfermi -> band):' in log
        if not expect:
            check(f"Bandlist '{bands if len(bands) < 20 else '-25 ... 8'}': window printed, no warning",
                  rc == 0 and shown and WARN not in log, f'rc={rc}')
        else:
            ok = rc == 0 and shown and WARN in log and all(e in log for e in expect)
            check(f"Bandlist '{bands}': warned ({', '.join(expect)})", ok, f'rc={rc}')
    rc, log = run('hBN_N30_xatu_bands', EX.format(tb=hbn, eig=eig, sta=sta, cutoff='# Exciton_cutoff\n50'))
    check('Xatu input (bands from the .states file): window printed, no warning',
          rc == 0 and 'Band window' in log and WARN not in log, f'rc={rc}')
    rc, log = run('hBN_N30_xatu_allexc', EX.format(tb=hbn, eig=eig, sta=sta, cutoff=''))
    check('no Exciton_cutoff: every exciton in the Xatu files is used (900 on hBN_N30)',
          rc == 0 and 'Exciton_cutoff not given: using all 900 excitons' in log, f'rc={rc}')
    rc, log = run('hBN_N30_xatu_cutneg', EX.format(tb=hbn, eig=eig, sta=sta, cutoff='# Exciton_cutoff\n-5'))
    check('Exciton_cutoff = -5: stops with a message', rc != 0 and 'Exciton_cutoff must be positive' in log, f'rc={rc}')
    if nfail:
        print(f'\nFAILED: {nfail} check(s)'); sys.exit(1)
    print('\nALL CHECKS PASSED')


if __name__ == '__main__':
    main()
