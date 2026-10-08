#!/usr/bin/env python3
"""Xatu HDF5 exciton archive as opticx input (Xatu_interface followed by one '<label>.h5' path).

Needs an opticx built with 'make HDF5=1', h5py, and hBN_tb.dat + hBN_N30.eigval/.states in --root.

A. Same answer as the text files. hBN_N30 is converted with tools/xatu_text_to_h5.py, which stores every number
   exactly as parsed, and a second-order run (Response = shg: single-particle and excitonic, X_nm from the
   envelope derivatives) must give the text-file spectra to 1e-12 (run-to-run noise of the OpenMP
   reductions is ~5e-15), and the same .omesp byte for byte (it once carried unset stencil energies of the
   inactive direction, so its bytes depended on the heap layout). The log must say the archive was read and
   its mesh matched opticx's.
B. Guards, each on a doctored copy of the archive: they must stop with a message naming the problem and
   compute nothing. Wrong Nfermi, Xatu submesh, shifted k mesh, basis rows out of order, archive without
   states (Xatu run without -c), Exciton_cutoff above the archive's count, Bandlist that disagrees, a file
   that is not a Xatu archive, a missing file. Bandlist equal to the archive's window and a full-BSE archive
   (warning) must run; no Exciton_cutoff must use all 900 excitons.
C. Optional (--xatu path/to/xatu, built with HDF5=1): a real Xatu run written both ways. The archive keeps full
   precision, the text files 7-8 digits, so the excitonic spectra agree to text precision, and the single-
   particle SHG obeys D3h better from the archive (its k mesh is exact; the text one is rounded to 1e-7).

Usage: python3 tools/check_xatu_h5.py --opticx bin/opticx --root . --workdir bin/check_xatu_h5 [--xatu PATH]
"""
import argparse
import os
import shutil
import subprocess
import sys

import h5py
import numpy as np

HEAD = """# Periodic dimensions
2
# Wannier90_filename
{tb}
# Xatu_interface
true
{xatu}
{extra}# Nfermi
{nf}
"""
LINEAR = """# OME_sp
linear
# OME_ex
linear
# Response
absorbance
# Energy_variables
4.0 8.0 0.1 21
"""
SHG = """# OME_sp
nonlinear
# OME_ex
nonlinear
# Response
shg
# Energy_variables
2.0 8.0 0.05 60
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--opticx', required=True)
    ap.add_argument('--root', default='.')
    ap.add_argument('--workdir', default='check_xatu_h5')
    ap.add_argument('--xatu', default='', help='Xatu binary built with HDF5=1 (enables part C)')
    a = ap.parse_args()
    root, work, opticx = os.path.abspath(a.root), os.path.abspath(a.workdir), os.path.abspath(a.opticx)
    tb, eig, sta = (os.path.join(root, f) for f in ('hBN_tb.dat', 'hBN_N30.eigval', 'hBN_N30.states'))
    for f in (tb, eig, sta):
        if not os.path.exists(f):
            print('SKIP: missing ' + f)
            sys.exit(0)
    if os.path.exists(work):
        shutil.rmtree(work)
    os.makedirs(work)
    nfail = 0

    def check(label, ok, msg=''):
        nonlocal nfail
        print(('  PASS   ' if ok else '  FAIL   ') + label + ('   ' + msg if msg else ''))
        nfail += 0 if ok else 1

    def run(tag, xatu, body, nf=1, extra=''):
        d = os.path.join(work, tag)
        os.makedirs(d, exist_ok=True)
        name = f'{tag}.txt'
        with open(os.path.join(d, name), 'w') as f:
            f.write(HEAD.format(tb=tb, xatu=xatu, extra=extra, nf=nf) + body)
        r = subprocess.run([opticx, name], cwd=d, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log = r.stdout.decode('utf-8', 'replace')
        open(os.path.join(d, 'run.log'), 'w').write(log)
        return r.returncode, log, d

    def spectra(d):
        return {f: np.loadtxt(os.path.join(d, f)) for f in sorted(os.listdir(d))
                if f.endswith('.dat') and not f.startswith('bands')}

    def max_rel(sa, sb):
        worst = 0.0
        for f, x in sa.items():
            y = sb[f]
            worst = max(worst, np.abs(x - y)[:, 1:].max() / np.abs(x[:, 1:]).max())
        return worst

    def doctored(name, edit):
        path = os.path.join(work, name)
        shutil.copy(archive, path)
        with h5py.File(path, 'r+') as f:
            edit(f)
        return path

    # ---------------------------------------------------------------------------------------------
    print('A. Archive = text files (hBN_N30, converted exactly)\n')
    archive = os.path.join(work, 'hBN_N30.h5')
    r = subprocess.run([sys.executable, os.path.join(root, 'tools', 'xatu_text_to_h5.py'), eig, sta,
                        '--nfermi', '1', '--ndim', '2', '-o', archive], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    check('converter: hBN_N30 text -> archive', r.returncode == 0, r.stdout.decode().strip().splitlines()[-1])
    if r.returncode != 0:
        sys.exit(1)
    cut = '# Exciton_cutoff\n120\n'
    rc_t, _, d_t = run('shg_text', f'{eig}\n{sta}', SHG, extra=cut)
    rc_h, log_h, d_h = run('shg_h5', archive, SHG, extra=cut)
    ok = rc_t == 0 and rc_h == 0
    check('second-order run from text and from archive', ok, f'rc text {rc_t}, archive {rc_h}')
    if ok:
        st, sh = spectra(d_t), spectra(d_h)
        same = sorted(st) == sorted(sh) and len(st) >= 2
        diff = max_rel(st, sh) if same else float('inf')
        check(f'same spectra ({", ".join(sorted(st))})', same and diff < 1e-12, f'max rel diff {diff:.1e}')
        om = 'ome_nonlinear_sp_hBN.omesp'
        a_om, b_om = (os.path.join(d, om) for d in (d_t, d_h))
        identical = os.path.exists(a_om) and os.path.exists(b_om) and open(a_om, 'rb').read() == open(b_om, 'rb').read()
        check(f'same {om}, byte for byte', identical)
    check('log: archive read, k mesh matched', 'Reading Xatu archive' in log_h and 'Xatu k mesh = opticx mesh' in log_h
          and 'Reading exciton wavefunctions from the Xatu archive' in log_h)

    rc, log, _ = run('all_excitons', archive, LINEAR)
    check('no Exciton_cutoff: all 900 excitons of the archive', rc == 0 and
          'Exciton_cutoff not given: using all 900 excitons in the Xatu archive' in log, f'rc={rc}')

    # ---------------------------------------------------------------------------------------------
    print('\nB. Guards\n')

    def stops(label, xatu, expect, nf=1, extra=''):
        rc, log, d = run(label, xatu, LINEAR, nf=nf, extra=extra)
        computed = any(f.endswith('.dat') and not f.startswith('bands') for f in os.listdir(d))
        found = all(e in log for e in expect)
        check(f'{label}: stops ({expect[0]})', rc != 0 and found and not computed, f'rc={rc}')

    stops('wrong_nfermi', archive, ['filled bands', 'Nfermi = 2'], nf=2)

    def submesh(f):
        f.attrs['submesh_factor'] = np.int64(2)
    stops('submesh', doctored('submesh.h5', submesh), ['full ncells^ndim'])

    def shift(f):
        f['kpoints'][:, 0] += 0.01
    stops('shifted_mesh', doctored('shifted.h5', shift), ['differs from the opticx mesh'])

    def swap_rows(f):
        b = f['basis'][()]
        b[0, 2], b[1, 2] = b[1, 2], b[0, 2]      # first two rows now list their k points in the wrong order
        f['basis'][...] = b
    stops('basis_order', doctored('basis_order.h5', swap_rows), ['exciton basis row'])

    def no_states(f):
        del f['basis']
    stops('no_states', doctored('no_states.h5', no_states), ['holds no exciton states', '-c -H'])

    stops('cutoff_too_large', archive, ['but the Xatu archive holds 900'], extra='# Exciton_cutoff\n1000\n')
    stops('bandlist_mismatch', archive, ['Bandlist disagrees'], extra='# Bandlist\n-1 0 1\n')

    def other_format(f):
        f.attrs['format'] = 'something-else'
    stops('not_an_archive', doctored('other.h5', other_format), ['is not a Xatu exciton archive'])
    stops('missing_file', os.path.join(work, 'nonexistent.h5'), ['Xatu archive not found'])

    rc, log, _ = run('bandlist_match', archive, LINEAR, extra='# Bandlist\n0 1\n# Exciton_cutoff\n50\n')
    check('Bandlist equal to the archive window: runs', rc == 0, f'rc={rc}')

    def full_bse(f):
        f.attrs['tamm_dancoff'] = np.bool_(False)
    rc, log, _ = run('full_bse', doctored('full_bse.h5', full_bse), LINEAR, extra='# Exciton_cutoff\n50\n')
    check('full-BSE archive: runs with a warning', rc == 0 and 'no Tamm-Dancoff' in log, f'rc={rc}')

    # ---------------------------------------------------------------------------------------------
    print('\nC. Real Xatu run written as text and as archive\n')
    if not a.xatu:
        print('  SKIP   no --xatu given')
    else:
        xd = os.path.join(work, 'xatu_hbn30')
        os.makedirs(xd)
        for mode, flags in (('text', ['-e', '-c']), ('h5', ['-c', '-H'])):
            with open(os.path.join(xd, f'hbn30_{mode}_exciton.txt'), 'w') as f:
                f.write(f'# label\nhbn30_{mode}\n# ncells\n30\n# bands\n1\n# Dielectric\n1 1 10\n#\n')
            r = subprocess.run([os.path.abspath(a.xatu), '-w', '1', tb, f'hbn30_{mode}_exciton.txt', '-n', '0'] + flags,
                               cwd=xd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            check(f'Xatu {mode} run', r.returncode == 0, f'rc={r.returncode}')
        x_eig, x_sta, x_h5 = (os.path.join(xd, f) for f in ('hbn30_text.eigval', 'hbn30_text.states', 'hbn30_h5.h5'))
        rc_t, _, d_t = run('xatu_text', f'{x_eig}\n{x_sta}', SHG)
        rc_h, _, d_h = run('xatu_h5', x_h5, SHG)
        if rc_t == 0 and rc_h == 0:
            st, sh = spectra(d_t), spectra(d_h)
            ex = {f: v for f, v in st.items() if '_ex_' in f}
            diff = max_rel(ex, {f: sh[f] for f in ex})
            check('excitonic spectra agree to text precision (< 1e-5)', diff < 1e-5, f'max rel diff {diff:.1e}')

            def d3h(s):
                x = s[:, 1] + 1j*s[:, 2]
                comp = {'xyy': 4, 'yxy': 10, 'yyx': 12}
                return max(np.abs(x + s[:, 1 + 2*i] + 1j*s[:, 2 + 2*i]).max() for i in comp.values()) / np.abs(x).max()
            f = 'shg_sp_lengthgauge_hBN.dat'
            rt, rh = d3h(st[f]), d3h(sh[f])
            check('sp SHG D3h residual: archive (exact mesh) below text (rounded mesh)', rh < rt,
                  f'archive {rh:.1e}, text {rt:.1e}')
        else:
            check('opticx on the real Xatu outputs', False, f'rc text {rc_t}, archive {rc_h}')

    if nfail:
        print(f'\nFAILED: {nfail} check(s)')
        sys.exit(1)
    print('\nALL CHECKS PASSED')


if __name__ == '__main__':
    main()
