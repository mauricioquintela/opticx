#!/usr/bin/env python3
"""Regression check of Response = shift_covariant: the single-particle shift current from paper Eq. (9)
with the U(n)-covariant generalised derivative over blocks of (near-)degenerate bands.

hBN, 30x30 opticx mesh (no degeneracies):
  * equals the exact IPA of tools/ipa_shift_numpy.py on the same k list, every in-plane component (tol 1e-5);
  * unchanged by OPTICX_GAUGE_TEST (smooth, rapidly varying eigenvector phases; tol 1e-8).
Monolayer MoS2, the repo's MoS2_spin_wannier_07032024_tb.dat (Gamma-M lines degenerate to the model's precision),
30x30, Bandlist -1..2, Nfermi 26:
  * C3 residual of the in-plane tensor < 2% and the D3h ratios xyy/xxx, yxy/xxx within 2% of -1 (measured 0.62%,
    -1.002, -1.001; the NumPy prototype gives the same to 1%);
  * SENSITIVITY: shift_shiftvector on the same grid breaks C3 by > 5% (measured 12.2%), so the check can see the
    defect the new path removes;
  * unchanged by OPTICX_GAUGE_TEST (tol 1e-8) and by OPTICX_BLOCK_SCRAMBLE (a unitary inside every numerically
    degenerate group at the evaluated k; tol 1e-5, measured 1.2e-7), while the stored velocity inside those groups
    DOES change (sensitivity: the hook really rotates the basis);
  * the log announces both test hooks.
Covariant single-particle SHG / rectification (Sp_method = covariant, the default since 2026-10-06):
  * hBN (no degeneracies): Response = shg and Response = rectification equal their Sp_method = per_band (Eq. A3a)
    counterparts (1e-6 / 1e-4; both now in the 2018 convention with the physical sign);
  * MoS2: Response = shg C3 < 2%, D3h ratios within 2% of -1 (measured 0.62%, -1.001, -1.002); unchanged by
    OPTICX_GAUGE_TEST (1e-8) and OPTICX_BLOCK_SCRAMBLE (1e-5, measured 9.6e-8); SENSITIVITY: Sp_method = per_band
    breaks C3 by > 50% on the same grid (measured 82%).
"""
import argparse, os, shutil, subprocess, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ipa_shift_numpy import read_tb, k_from_omesp, ipa_shift      # noqa: E402

INP = """# Periodic dimensions
2
# Wannier90_filename
{tb}
# Xatu_interface
false
# Ncells
30
# Bandlist
{bands}
# Nfermi
{nf}
# OME_sp
nonlinear
# Response
{resp}
# Energy_variables
{ev}
"""
C, S = np.cos(2 * np.pi / 3), np.sin(2 * np.pi / 3)
C3 = np.array([[C, -S], [S, C]])


def c3res(T):
    rot = lambda R: np.einsum('ai,bj,ck,wijk->wabc', R, R, R, T)
    return np.linalg.norm(T - (T + rot(C3) + rot(C3 @ C3)) / 3) / np.linalg.norm(T)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--opticx', required=True); ap.add_argument('--root', default='.')
    ap.add_argument('--workdir', default='check_shift_covariant')
    a = ap.parse_args()
    root, work, opticx = os.path.abspath(a.root), os.path.abspath(a.workdir), os.path.abspath(a.opticx)
    hbn, mos2 = os.path.join(root, 'hBN_tb.dat'), os.path.join(root, 'MoS2_spin_wannier_07032024_tb.dat')
    for f in (hbn, mos2):
        if not os.path.exists(f):
            sys.exit(f'missing {f}')
    shutil.rmtree(work, ignore_errors=True); os.makedirs(work)
    nfail = 0

    def check(label, ok, msg=''):
        nonlocal nfail
        print(('  PASS   ' if ok else '  FAIL   ') + label + ('   ' + msg if msg else ''))
        nfail += 0 if ok else 1

    def run(tag, tb, bands, nf, resp, ev, env_extra=None, keep_omesp=False, extra=''):
        d = os.path.join(work, tag); os.makedirs(d)
        name = f'{tag}.txt'                                   # descriptive input name (htop)
        open(os.path.join(d, name), 'w').write(INP.format(tb=tb, bands=bands, nf=nf, resp=resp, ev=ev) + extra)
        env = dict(os.environ)
        for k in ('OPTICX_GAUGE_TEST', 'OPTICX_BLOCK_SCRAMBLE'):
            env.pop(k, None)
        env.update(env_extra or {})
        r = subprocess.run([opticx, name], cwd=d, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log = r.stdout.decode('utf-8', 'replace'); open(os.path.join(d, 'run.log'), 'w').write(log)
        if r.returncode != 0:
            check(f'{tag} ran', False, f'rc={r.returncode}')
        if not keep_omesp:
            for f in os.listdir(d):
                if f.endswith('.omesp'):
                    os.remove(os.path.join(d, f))
        return d, log

    load = lambda d: np.loadtxt([os.path.join(d, f) for f in os.listdir(d) if f.startswith('shift_sp_')][0])[:, 1:28] \
        .reshape(-1, 3, 3, 3)[:, :2, :2, :2]
    rel = lambda x, y: abs(x - y).max() / abs(y).max()

    print('shift_covariant on hBN 30x30')
    HEV = '6.5 12.5 0.15 600'
    d0, _ = run('hBN_N30_shift_covariant', hbn, '0 1', 1, 'shift_covariant', HEV, keep_omesp=True)
    T = load(d0)
    K = k_from_omesp(os.path.join(d0, 'ome_nonlinear_sp_hBN.omesp'))
    w = 6.5 + 6.0 / 600 * np.arange(600)
    ex = ipa_shift(read_tb(hbn), K, 1, w, 0.15)
    for name, (i, j, k) in {'xxx': (0, 0, 0), 'xyy': (0, 1, 1), 'yxy': (1, 0, 1), 'yyx': (1, 1, 0)}.items():
        e = rel(T[:, i, j, k], ex[:, i, j, k])
        check(f'hBN  {name} equals the exact IPA (Eq. 9)', e < 1e-5, f'{e:.2e} (tol 1e-5)')
    for f in os.listdir(d0):
        if f.endswith('.omesp'):
            os.remove(os.path.join(d0, f))
    d1, log = run('hBN_N30_shift_covariant_gauge', hbn, '0 1', 1, 'shift_covariant', HEV, {'OPTICX_GAUGE_TEST': '1'})
    check('hBN  log announces the gauge test', 'gauge test' in log.lower())
    e = rel(load(d1), T); check('hBN  unchanged by OPTICX_GAUGE_TEST', e < 1e-8, f'{e:.2e} (tol 1e-8)')

    print('\nshift_covariant on MoS2 30x30 (Gamma-M degenerate)')
    MEV = '1.8 4.5 0.025 1081'
    dm, _ = run('MoS2_N30_shift_covariant', mos2, '-1 0 1 2', 26, 'shift_covariant', MEV, keep_omesp=True)
    M = load(dm); x = M[:, 0, 0, 0]
    r = c3res(M); check('MoS2 C3 residual', r < 0.02, f'{100 * r:.2f}% (max 2%)')
    q1, q2 = (x @ M[:, 0, 1, 1]) / (x @ x), (x @ M[:, 1, 0, 1]) / (x @ x)
    check('MoS2 D3h xyy/xxx = -1', abs(q1 + 1) < 0.02, f'{q1:.4f}')
    check('MoS2 D3h yxy/xxx = -1', abs(q2 + 1) < 0.02, f'{q2:.4f}')
    ds, _ = run('MoS2_N30_shift_shiftvector', mos2, '-1 0 1 2', 26, 'shift_shiftvector', MEV)
    r = c3res(load(ds)); check('SENSITIVITY  shift_shiftvector breaks C3 on the same grid', r > 0.05,
                               f'{100 * r:.2f}% (must exceed 5%)')
    dg, log = run('MoS2_N30_shift_covariant_gauge', mos2, '-1 0 1 2', 26, 'shift_covariant', MEV, {'OPTICX_GAUGE_TEST': '1'})
    e = rel(load(dg), M); check('MoS2 unchanged by OPTICX_GAUGE_TEST', e < 1e-8, f'{e:.2e} (tol 1e-8)')
    db, log = run('MoS2_N30_shift_covariant_blockscramble', mos2, '-1 0 1 2', 26, 'shift_covariant', MEV,
                  {'OPTICX_BLOCK_SCRAMBLE': '1'}, keep_omesp=True)
    check('MoS2 log announces the block scramble', 'block scramble test' in log.lower())
    e = rel(load(db), M); check('MoS2 unchanged by OPTICX_BLOCK_SCRAMBLE', e < 1e-5, f'{e:.2e} (tol 1e-5)')
    # sensitivity of the hook: the stored velocity inside numerically degenerate groups must differ
    def vraw(d):
        f = open(os.path.join(d, 'ome_nonlinear_sp_MoS2_spin_wannier_07032024.omesp'), 'rb')
        i4 = lambda n: np.fromfile(f, '<i4', n)
        i4(1); npt, nb = i4(2); np.fromfile(f, '<f8', 3 * npt)
        rec = np.dtype([('e', '<f8', nb), ('v', '<c16', 3 * nb * nb), ('b', '<c16', 3 * nb * nb),
                        ('s', '<f8', 9 * nb * nb), ('g', '<c16', 9 * nb * nb)])
        np.fromfile(f, rec, npt); np.fromfile(f, '<f8', npt * 9 * nb * nb); np.fromfile(f, '<c16', npt * 9 * nb * nb)
        here = f.tell(); tag, flags, norb = i4(3)
        if tag == 1330464562:                         # tagged tail (since 2026-10-06): skip A4 rotation, window states
            if flags & 1: np.fromfile(f, '<c16', npt * nb * nb)
            if flags & 2: np.fromfile(f, '<c16', 2 * norb * nb * npt + npt * 3 * nb * nb)
        else:
            f.seek(here)
        np.fromfile(f, '<c16', npt * 9 * nb * nb)     # block-covariant derivative
        return np.fromfile(f, '<c16', npt * 3 * nb * nb)
    dv = abs(vraw(db) - vraw(dm)).max()
    check('SENSITIVITY  the block scramble really rotates the basis', dv > 1e-3, f'max |dv| {dv:.2e} (must exceed 1e-3)')
    for d in (dm, db):
        for f in os.listdir(d):
            if f.endswith('.omesp'):
                os.remove(os.path.join(d, f))

    print('\ncovariant single-particle SHG / rectification (Sp_method = covariant, the default)')
    SEV = '2.0 6.0 0.15 200'
    PB = '# Sp_method\nper_band\n'
    lshg = lambda d: (lambda a: (a[:, 1::2] + 1j * a[:, 2::2])[:, :27].reshape(-1, 3, 3, 3))(
        np.loadtxt([os.path.join(d, f) for f in os.listdir(d) if f.startswith('shg_sp_')][0], comments='#'))
    lsec = lambda d: (lambda a: (a[:, 2::2] + 1j * a[:, 3::2])[:, :27].reshape(-1, 3, 3, 3))(
        np.loadtxt([os.path.join(d, f) for f in os.listdir(d) if f.startswith('second_')][0], comments='#'))
    h1, _ = run('hBN_N30_shg_per_band', hbn, '0 1', 1, 'shg', SEV, extra=PB); h2, _ = run('hBN_N30_shg', hbn, '0 1', 1, 'shg', SEV)
    e = rel(lshg(h2), lshg(h1)); check('hBN  shg (covariant) equals Sp_method = per_band (Eq. A3a)', e < 1e-6, f'{e:.2e} (tol 1e-6)')
    r1, _ = run('hBN_N30_rect_per_band', hbn, '0 1', 1, 'rectification', '6.5 12.5 0.15 60', extra=PB)
    r2, _ = run('hBN_N30_rect', hbn, '0 1', 1, 'rectification', '6.5 12.5 0.15 60')
    # tolerance 1e-4, measured 3.5e-5: at the DC point method B's prefactor i w_2 = -2 eta multiplies terms of order 1/eta^2,
    # which amplifies the finite-difference rounding of both routes; both equal the real-time current (check_realtime_sign)
    e = rel(lsec(r2), lsec(r1)); check('hBN  rectification (covariant) equals Sp_method = per_band', e < 1e-4, f'{e:.2e} (tol 1e-4)')
    MSEV = '0.8 2.4 0.025 641'
    m1, _ = run('MoS2_N30_shg', mos2, '-1 0 1 2', 26, 'shg', MSEV)
    G = lshg(m1); g = G[:, :2, :2, :2]; xg = g[:, 0, 0, 0]
    r = c3res(g); check('MoS2 shg C3 residual', r < 0.02, f'{100 * r:.2f}% (max 2%)')
    q1 = (np.vdot(xg, g[:, 0, 1, 1]) / np.vdot(xg, xg)).real; q2 = (np.vdot(xg, g[:, 1, 0, 1]) / np.vdot(xg, xg)).real
    check('MoS2 shg D3h xyy/xxx = -1', abs(q1 + 1) < 0.02, f'{q1:.4f}')
    check('MoS2 shg D3h yxy/xxx = -1', abs(q2 + 1) < 0.02, f'{q2:.4f}')
    m0, _ = run('MoS2_N30_shg_per_band', mos2, '-1 0 1 2', 26, 'shg', MSEV, extra=PB)
    r = c3res(lshg(m0)[:, :2, :2, :2]); check('SENSITIVITY  Sp_method = per_band breaks C3 on the same grid', r > 0.5,
                                            f'{100 * r:.1f}% (must exceed 50%)')
    m2, _ = run('MoS2_N30_shg_gauge', mos2, '-1 0 1 2', 26, 'shg', MSEV, {'OPTICX_GAUGE_TEST': '1'})
    e = rel(lshg(m2), G); check('MoS2 shg unchanged by OPTICX_GAUGE_TEST', e < 1e-8, f'{e:.2e} (tol 1e-8)')
    m3, _ = run('MoS2_N30_shg_block', mos2, '-1 0 1 2', 26, 'shg', MSEV, {'OPTICX_BLOCK_SCRAMBLE': '1'})
    e = rel(lshg(m3), G); check('MoS2 shg unchanged by OPTICX_BLOCK_SCRAMBLE', e < 1e-5, f'{e:.2e} (tol 1e-5)')

    if nfail:
        print(f'\nFAILED: {nfail} check(s)'); sys.exit(1)
    print('\nALL CHECKS PASSED')


if __name__ == '__main__':
    main()
