#!/usr/bin/env python3
"""Single-particle shift current (Response = shift, Sp_method = covariant) against the velocity sum rule of Kameda,
Nagashio and Wakabayashi, arXiv:2611.05990, Eqs. 23-25, generalised to a full position matrix and degenerate bands
(tools/sumrule_shift.py). The sum rule shares nothing with opticx but the tight-binding file: analytic derivatives of
H(k) and of the position matrices, no k-mesh stencil, no parallel transport, an explicit sum over every band of the
model. Same Monkhorst-Pack mesh and Lorentzian on both sides, so they must agree to the finite-difference error of
opticx's derivatives.

 1. buckled hBN 30x30 (two bands, C3v, Wannier centres off the plane): every component, the out-of-plane ones (zxx,
    xxz, zzz) included, to 1e-8 of the tensor (measured 1e-11).
 2. MoS2 30x30, 2v x 2c (34-band spinor model, spin-degenerate Gamma-M lines handled by blocks on both sides): to
    2e-3 (measured 1.6e-4; ~1 min in all), and the formula as printed (orbital centres only) is > 5% off -- the
    full position matrix matters.
MoS2 is skipped (loudly) when MoS2_spin_wannier_07032024_tb.dat is not in --root.

Usage: python3 tools/check_sumrule_shift.py --opticx bin/opticx --root . --workdir bin/check_sumrule_shift
"""
import argparse, os, subprocess, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sumrule_shift as S
from ipa_shift_numpy import HA
from check_realtime_sign import C_AU

INP = """# Periodic dimensions
2
# Wannier90_filename
{tb}
# Xatu_interface
false
# Ncells
{nc}
# Bandlist
{bl}
# Nfermi
{nf}
# OME_sp
nonlinear
# Response
shift
# Energy_variables
{ev}
"""
CASES = (('buckled hBN', 'buckled-hBN_tb.dat', 1, '0 1', (5.0, 13.0, 0.1, 400)),
         ('MoS2 2v x 2c', 'MoS2_spin_wannier_07032024_tb.dat', 26, '-1 0 1 2', (0.5, 4.5, 0.05, 400)))
NC = 30


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--opticx', required=True); ap.add_argument('--root', default='.')
    ap.add_argument('--workdir', default='check_sumrule_shift')
    a = ap.parse_args()
    root, work, opticx = os.path.abspath(a.root), os.path.abspath(a.workdir), os.path.abspath(a.opticx)
    os.makedirs(work, exist_ok=True)
    nfail = 0

    def check(label, ok, msg=''):
        nonlocal nfail
        print(('  PASS   ' if ok else '  FAIL   ') + label + ('   ' + msg if msg else '')); nfail += 0 if ok else 1

    rel = lambda X, Y: np.linalg.norm(X - Y)/np.linalg.norm(Y)
    print(f'Single-particle shift current vs the velocity sum rule of arXiv:2611.05990 (generalised), {NC}x{NC}\n')
    for label, tbname, nf, bl, (e1, e2, eta, nw) in CASES:
        tb = os.path.join(root, tbname)
        if not os.path.exists(tb):
            print(f'SKIP {label}: missing {tb}'); continue
        tag = tbname.replace('_tb.dat', ''); d = os.path.join(work, tag); os.makedirs(d, exist_ok=True)
        open(os.path.join(d, f'{tag}_sp_shift.txt'), 'w').write(INP.format(tb=tb, nc=NC, bl=bl, nf=nf, ev=f'{e1} {e2} {eta} {nw}'))
        r = subprocess.run([opticx, f'{tag}_sp_shift.txt'], cwd=d, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        open(os.path.join(d, 'run.log'), 'wb').write(r.stdout)
        if r.returncode != 0:
            print(f'FAIL: opticx run {tag} exited with {r.returncode} (see {d}/run.log)'); sys.exit(1)
        O = np.loadtxt(os.path.join(d, f'shift_sp_lengthgauge_{tag}.dat'), comments='#')[:, 1:28]
        O = O.reshape(-1, 3, 3, 3).transpose(1, 2, 3, 0)/C_AU
        w = (e1 + (e2 - e1)/nw*np.arange(nw))/HA
        window = [nf - 1 + int(o) for o in bl.split()]
        G = S.spectra(S.model(tb), NC, nf, window, w, eta/HA)[0].real
        print(label)
        tol = 1e-8 if label.startswith('buckled') else 2e-3
        e = rel(O, G); check('opticx = sum rule, whole tensor', e < tol, f'{e:.1e} (tol {tol:.0e})')
        if label.startswith('buckled'):
            for c in ('zxx', 'xxz', 'zzz'):
                i, j, l = ('xyz'.index(x) for x in c)
                e = rel(O[i, j, l], G[i, j, l]); check(f'  {c}', e < tol, f'{e:.1e} (tol {tol:.0e}; max |{c}| {abs(G[i, j, l]).max():.3f} a.u.)')
        else:
            P = S.spectra(S.model(tb, centres_only=True), NC, nf, window, w, eta/HA)[0].real
            e = rel(P, O); check('as printed (orbital centres only) differs', e > 5e-2, f'{e:.3f} (> 0.05)')
    print(f'\n{"ALL PASS" if nfail == 0 else f"{nfail} FAIL"}')
    sys.exit(1 if nfail else 0)


if __name__ == '__main__':
    main()
