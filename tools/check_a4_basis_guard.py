#!/usr/bin/env python3
"""The Eq. (A4) basis with OME_sp = none: stored in the .omesp, and refused when it is missing.

Why this needs a test of its own. get_ome_sp rotates the single-particle states into the Eq. (A4) basis and
rotate_fk_ex_to_a4_basis carries Xatu's envelopes into the same basis. With OME_sp = none the rotation is not
recomputed, so since 2026-10-06 the .omesp stores it (the linear file a marked block, the nonlinear file a tagged
tail that also holds the exciton-window states the covariant X_nm needs). If the envelopes were left in Xatu's
basis the error would be O(1): on In2Se3 (4 valence x 1 conduction) it moved the excitonic SHG by 13% of the
largest tensor component and by more than 100% of several smaller ones. With the stored basis, OME_sp = none on
In2Se3 reproduces the full calculation to 4e-16 (SHG) and 4e-16 (absorbance), 14 s -> 2 s.

hBN cannot expose the size of the error -- with a single valence-conduction pair the rotation is a phase -- so
this test asserts the mechanism: OME_sp = none with a current .omesp is ALLOWED and EXACT (linear and second
order, and the single-particle covariant arrays behind the tag are read back intact); with an .omesp written
before the change (simulated by cutting the file to the old layout) the excitonic path is REFUSED, nothing is
computed, and the message names the fix. The cache exemption still holds: a second-order cache hit never reads
the envelopes, so it stays allowed even with an old .omesp, while a cache MISS under those settings still stops.

Usage: python3 tools/check_a4_basis_guard.py --opticx bin/opticx --root . --workdir bin/check_a4
"""
import argparse, os, re, shutil, struct, subprocess, sys
import numpy as np

IN = """# Periodic dimensions
2
# Wannier90_filename
{tb}
# Xatu_interface
true
{eig}
{sta}
# Exciton_cutoff
120
# Nfermi
1
# OME_sp
{omesp}
# OME_ex
{omeex}
# Response
{resp}
# Energy_variables
2.0 8.0 0.05 60
# Cache_ome_ex
{mode}
"""
CACHE = 'ome_second_ex_hBN.omeex2'
NL = 'ome_nonlinear_sp_hBN.omesp'
LIN = 'ome_linear_sp_hBN.omesp'
OMESP_TAG = 1330464562


def old_layout_nonlinear(src, dst):
    """Cut a tagged nonlinear .omesp back to the layout written before 2026-10-06 (no tail)."""
    b = open(src, 'rb').read()
    nk, nb = struct.unpack('<ii', b[4:12])
    pos = 4 + 8 + 3*nk*8 + nk*(nb*8 + 2*3*nb*nb*16 + 9*nb*nb*8 + 9*nb*nb*16) + nk*9*nb*nb*8 + nk*9*nb*nb*16
    tag = struct.unpack('<i', b[pos:pos+4])[0] if len(b) >= pos + 4 else None
    open(dst, 'wb').write(b[:pos])
    return tag


def old_layout_linear(src, dst):
    lines = open(src).read().split('\n')
    k = lines.index('#A4W') if '#A4W' in lines else len(lines)
    open(dst, 'w').write('\n'.join(lines[:k]) + '\n')
    return k < len(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--opticx', required=True)
    ap.add_argument('--root', default='.')
    ap.add_argument('--workdir', default='check_a4')
    a = ap.parse_args()
    root, work = os.path.abspath(a.root), os.path.abspath(a.workdir)
    opticx = os.path.abspath(a.opticx)
    tb = os.path.join(root, 'hBN_tb.dat')
    eig = os.path.join(root, 'hBN_N30.eigval')
    sta = os.path.join(root, 'hBN_N30.states')
    for f in (tb, eig, sta):
        if not os.path.exists(f):
            print('SKIP: missing ' + f); sys.exit(0)
    shutil.rmtree(work, ignore_errors=True); os.makedirs(work)
    nfail = 0

    def check(label, ok, msg=''):
        nonlocal nfail
        print(('  PASS   ' if ok else '  FAIL   ') + label + ('   ' + msg if msg else ''))
        nfail += 0 if ok else 1

    def run(tag, omesp, mode, omeex='nonlinear', resp='shift_shiftvector', copy=()):
        d = os.path.join(work, tag); os.makedirs(d, exist_ok=True)
        name = f'hBN_a4_{tag}.txt'
        open(os.path.join(d, name), 'w').write(
            IN.format(tb=tb, eig=eig, sta=sta, omesp=omesp, omeex=omeex, resp=resp, mode=mode))
        for src, fname in copy:
            p = os.path.join(work, src, fname)
            if os.path.exists(p):
                shutil.copy2(p, os.path.join(d, fname))
        r = subprocess.run([opticx, name], cwd=d, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log = r.stdout.decode('utf-8', 'replace')
        open(os.path.join(d, 'run.log'), 'w').write(log)
        return r.returncode, log, d

    def spec(d, pat):
        f = [x for x in os.listdir(d) if x.startswith(pat)]
        return np.loadtxt(os.path.join(d, f[0]), comments='#') if f else None

    def rel(x, y):
        if x is None or y is None: return np.inf
        return abs(x - y).max() / max(abs(y).max(), 1e-300)

    kloop = lambda log: len(re.findall(r'OME \(ex\): k-point', log))
    refused = lambda log: 'not in the same basis' in log.lower()
    loaded = lambda log: 'basis of the exciton window read from' in log.lower()

    print('Eq. (A4) basis with OME_sp = none: stored in the .omesp, refused when missing\n')

    # -- reference runs that compute everything and write the .omesp files ---------------------------
    rc, log, d = run('nl_full', 'nonlinear', 'write')
    check('OME_sp = nonlinear          runs', rc == 0 and kloop(log) > 0, f'{kloop(log)} k-points')
    b = open(os.path.join(d, NL), 'rb').read() if os.path.exists(os.path.join(d, NL)) else b''
    tmp = os.path.join(work, 'probe.omesp')
    if b:
        open(tmp, 'wb').write(b)
    tag = old_layout_nonlinear(tmp, os.path.join(work, 'old_' + NL)) if b else None
    check('nonlinear .omesp carries the tagged tail', tag == OMESP_TAG, f'tag {tag}')
    ref_ex = spec(d, 'shift_ex_'); ref_sp = spec(d, 'shift_sp_')

    rc, log, d = run('lin_full', 'linear', 'off', omeex='linear', resp='absorbance')
    check('OME_sp = linear             runs', rc == 0)
    has = old_layout_linear(os.path.join(d, LIN), os.path.join(work, 'old_' + LIN)) if os.path.exists(os.path.join(d, LIN)) else False
    check('linear .omesp carries the A4 block', has)
    ref_lin = spec(d, 'sigma_first_ex_real_')

    # -- current .omesp: OME_sp = none is ALLOWED and EXACT --------------------------------------------
    rc, log, d = run('nl_none', 'none', 'off', copy=(('nl_full', NL),))
    check('OME_sp = none, current .omesp   allowed (2nd order)', rc == 0 and not refused(log), f'rc={rc}')
    check('OME_sp = none, current .omesp   basis read from the .omesp', loaded(log))
    e = rel(spec(d, 'shift_ex_'), ref_ex)
    check('OME_sp = none, current .omesp   excitonic result exact', e < 1e-12, f'{e:.2e}')

    rc, log, d = run('lin_none', 'none', 'off', omeex='linear', resp='absorbance', copy=(('lin_full', LIN),))
    check('OME_sp = none, current .omesp   allowed (linear)', rc == 0 and not refused(log), f'rc={rc}')
    e = rel(spec(d, 'sigma_first_ex_real_'), ref_lin)
    check('OME_sp = none, current .omesp   linear excitonic result exact', e < 1e-12, f'{e:.2e}')

    # -- the covariant single-particle arrays sit behind the tag and must still be found --------------
    rc, log, d = run('cov_full', 'nonlinear', 'off', resp='shift')
    ref_cov = spec(d, 'shift_sp_')
    rc, log, d = run('cov_none', 'none', 'off', resp='shift', copy=(('cov_full', NL),))
    e = rel(spec(d, 'shift_sp_'), ref_cov)
    check('OME_sp = none, Response = shift   covariant arrays read back exactly', rc == 0 and e < 1e-12,
          f'rc={rc} {e:.2e}')

    # -- an .omesp written before the change: refused, nothing computed, the fix is named --------------
    os.makedirs(os.path.join(work, 'old'), exist_ok=True)
    for f in (NL, LIN):
        if os.path.exists(os.path.join(work, 'old_' + f)):
            shutil.copy2(os.path.join(work, 'old_' + f), os.path.join(work, 'old', f))
    rc, log, d = run('nl_none_old', 'none', 'off', copy=(('old', NL),))
    check('OME_sp = none, old .omesp   REFUSED', rc != 0 and refused(log),
          f'rc={rc}' + ('' if refused(log) else '  -- ran with mismatched bases'))
    check('OME_sp = none, old .omesp   no k-loop was run', kloop(log) == 0, f'{kloop(log)} k-points')
    check('OME_sp = none, old .omesp   no spectrum written', spec(d, 'shift_ex_') is None)
    check('OME_sp = none, old .omesp   message names the fix',
          'regenerate the .omesp' in log.lower() and 'cache_ome_ex = read' in log.lower())
    rc, log, d = run('lin_none_old', 'none', 'off', omeex='linear', resp='absorbance', copy=(('old', LIN),))
    check('OME_sp = none, old linear .omesp   REFUSED', rc != 0 and refused(log), f'rc={rc}')

    # -- the cache exemption (old .omesp): a hit never reads the envelopes; a miss must still stop -----
    rc, log, d = run('old_cache_hit', 'none', 'read', copy=(('nl_full', CACHE), ('old', NL)))
    e = rel(spec(d, 'shift_ex_'), ref_ex)
    check('old .omesp, cache HIT   allowed and exact', rc == 0 and not refused(log) and e < 1e-12,
          f'rc={rc} {e:.2e}')
    check('old .omesp, cache HIT   k-loop skipped', kloop(log) == 0, f'{kloop(log)} k-points')
    rc, log, d = run('old_cache_miss', 'none', 'read', copy=(('old', NL),))
    check('old .omesp, cache MISS  REFUSED', rc != 0 and refused(log) and kloop(log) == 0, f'rc={rc}')

    if nfail:
        print(f'\nFAILED: {nfail} check(s)'); sys.exit(1)
    print('\nALL CHECKS PASSED')


if __name__ == '__main__':
    main()
