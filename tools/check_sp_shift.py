#!/usr/bin/env python3
"""Regression check of the single-particle shift conductivity (shift_shiftvector) on hBN, with tolerances.

What it does (about 1 minute):
  1. builds synthetic non-interacting excitons (ni_excitons.py) from hBN_tb.dat + hBN_N30.states;
  2. runs opticx once with Response = shift_shiftvector: the sp code and, in the same run, the excitonic code, which for
     non-interacting excitons is the independent-particle limit of paper Eq. (10)/(11), on the same k-grid;
  3. evaluates the exact IPA expression, paper Eq. (9), independently in NumPy (ipa_shift_numpy.py);
  4. compares, component by component, and exits with status 1 if a tolerance is violated.
What it catches: the overall sign (paper Eq. 9), the normalisation, the missing amplitude-gradient term for b != c
(before the fix hBN yxy came out 0.324x the exact value), and the C3 relations xxx = -xyy = -yxy = -yyx.

Usage:  python3 tools/check_sp_shift.py --opticx bin/opticx --root . --workdir bin/check_sp_shift
Tolerances are set from the measured values (30x30 grid): sp vs exact xxx/xyy scale 1.001/0.981, yxy/yyx 0.981 (the
2-3% shortfall is the mirror line where r^y = 0 exactly: a central difference of |v| is zero there), C3-forbidden
fraction 1.6e-2, sp vs excitonic limit 0.984-1.019.
"""
import argparse, os, subprocess, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ni_excitons import make_ni
from ipa_shift_numpy import read_tb, k_from_states, ipa_shift
import plot_test_outputs

U = np.array([[1, 1j], [1, -1j]]) / np.sqrt(2)


def c3_forbidden(T):          # T[w,a,b,c] in-plane; fraction of the weight in the charge +-1 (C3-forbidden) components
    Tp = np.einsum('ia,jb,kc,wabc->wijk', U, U, U, T)
    tot = (abs(Tp) ** 2).sum(axis=(1, 2, 3)); al = abs(Tp[:, 0, 0, 0]) ** 2 + abs(Tp[:, 1, 1, 1]) ** 2
    return np.sqrt(((tot - al).clip(0)).sum() / tot.sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--opticx', required=True); ap.add_argument('--root', default='.'); ap.add_argument('--workdir', default='check_sp_shift')
    a = ap.parse_args()
    root = os.path.abspath(a.root); work = os.path.abspath(a.workdir); os.makedirs(work, exist_ok=True)
    tb = os.path.join(root, 'hBN_tb.dat'); states = os.path.join(root, 'hBN_N30.states')
    make_ni(tb, states, os.path.join(work, 'hBN_NI30'))
    inp = f"""# Periodic dimensions
2
# Wannier90_filename
{tb}
# Xatu_interface
true
{work}/hBN_NI30.eigval
{work}/hBN_NI30.states
# Exciton_cutoff
900
# Nfermi
1
# OME_sp
nonlinear
# OME_ex
nonlinear
# Response
shift_shiftvector
# Energy_variables
6.5 12.5 0.15 600
"""
    open(os.path.join(work, 'hBN_NI30_shift_check.txt'), 'w').write(inp)
    r = subprocess.run([os.path.abspath(a.opticx), 'hBN_NI30_shift_check.txt'], cwd=work, stdout=open(os.path.join(work, 'out.log'), 'w'), stderr=subprocess.STDOUT)
    if r.returncode != 0:
        print('FAIL: opticx exited with status', r.returncode, '(see', os.path.join(work, 'out.log') + ')'); sys.exit(1)
    ld = lambda fn: np.loadtxt(os.path.join(work, fn))[:, 1:].reshape(-1, 3, 3, 3)[:, :2, :2, :2]
    sp = ld('shift_sp_lengthgauge_hBN.dat'); ex = ld('shift_ex_lengthgauge_hBN.dat')
    w = 6.5 + (12.5 - 6.5) / 600 * np.arange(600)
    exact = ipa_shift(read_tb(tb), k_from_states(states), 1, w, 0.15)

    cols = ['E_eV']; data = [w]
    for name, (i, j, k) in {'xxx': (0, 0, 0), 'xyy': (0, 1, 1), 'yxy': (1, 0, 1)}.items():
        cols += [name + '_sp', name + '_exact', name + '_ex']; data += [sp[:, i, j, k], exact[:, i, j, k], ex[:, i, j, k]]
    spec = os.path.join(work, 'check_sp_shift_spectra.dat')
    np.savetxt(spec, np.array(data).T, header='kind: sp_shift\ncolumns: ' + ' '.join(cols), comments='# ')
    plot_test_outputs.plot(spec)

    nfail = 0
    def check(label, ok, msg):
        nonlocal nfail
        print(('  PASS   ' if ok else '  FAIL   ') + label + '   ' + msg); nfail += 0 if ok else 1
    def fit(x, ref):
        s = (x * ref).sum() / (ref * ref).sum(); return s, np.linalg.norm(x - s * ref) / np.linalg.norm(x)

    print(f'sp shift_shiftvector on non-interacting hBN 30x30 (eta 0.15 eV), {len(w)} frequencies')
    for name, ((i, j, k), (slo, shi, rmax)) in {'xxx': ((0, 0, 0), (0.97, 1.02, 0.03)), 'xyy': ((0, 1, 1), (0.97, 1.02, 0.03)),
                                                'yxy': ((1, 0, 1), (0.95, 1.03, 0.05)), 'yyx': ((1, 1, 0), (0.95, 1.03, 0.05))}.items():
        s, r = fit(sp[:, i, j, k], exact[:, i, j, k])
        check(f'sp vs exact IPA (Eq. 9) {name}', slo < s < shi and r < rmax, f'scale {s:+.3f} (want {slo}..{shi}), residual {r:.3f} (max {rmax})')
    for name, (i, j, k) in {'xxx': (0, 0, 0), 'xyy': (0, 1, 1), 'yxy': (1, 0, 1)}.items():
        s, r = fit(sp[:, i, j, k], ex[:, i, j, k])
        check(f'sp vs excitonic IPA limit (Eq. 10) {name}', 0.95 < s < 1.06 and r < 0.03, f'scale {s:+.3f} (want 0.95..1.06, sign +), residual {r:.3f} (max 0.03)')
    f = c3_forbidden(sp)
    check('C3: forbidden fraction of the sp tensor', f < 0.05, f'{f:.3e} (max 5e-2; exact IPA {c3_forbidden(exact):.1e})')

    # Excitonic rectification in the same non-interacting limit against the single-particle (causal, covariant)
    # rectification of the same run. Both are the whole causal sigma(0; w, -w); the excitonic one is method A
    # with term 3 on the bare exciton-exciton current, which makes the non-interacting limit exact on the mesh
    # (with Pi_nm = i(E_n-E_m) X_nm instead it converges only with the k-mesh).
    d = os.path.join(work, 'rect'); os.makedirs(d, exist_ok=True)
    name = 'hBN_NI30_exrect_causal.txt'
    open(os.path.join(d, name), 'w').write(inp.replace('shift_shiftvector', 'rectification').replace('6.5 12.5 0.15 600', '4.0 10.0 0.15 61'))
    rr = subprocess.run([os.path.abspath(a.opticx), name], cwd=d, stdout=open(os.path.join(d, 'out.log'), 'w'),
                        stderr=subprocess.STDOUT)
    if rr.returncode != 0:
        print('FAIL: opticx rectification run exited with', rr.returncode); sys.exit(1)
    lr = lambda fn: np.loadtxt(os.path.join(d, fn), comments='#')
    A = lr('second_ex_rectification_lengthgauge_hBN.dat'); B = lr('second_rectification_lengthgauge_hBN.dat')
    wr, exr, spr = A[:, 0], A[:, 2], B[:, 2]                         # Re xxx
    below = wr < 6.5
    e = abs(exr[below] - spr[below]).max() / abs(spr[below]).max()
    check('excitonic rectification = sp below the gap (NI)', e < 5e-3, f'{e:.2e} (max 5e-3)')
    ip = np.argmax(abs(spr)); e = abs(exr[ip] - spr[ip]) / abs(spr[ip])
    check('excitonic rectification = sp at the peak (NI)', e < 5e-3, f'{e:.2e} at {wr[ip]:.1f} eV (max 5e-3)')
    if nfail:
        print(f'\nFAILED: {nfail} check(s)'); sys.exit(1)
    print('\nALL CHECKS PASSED')


if __name__ == '__main__':
    main()
