#!/usr/bin/env python3
"""Independent NumPy evaluation of the single-particle (IPA) shift conductivity, paper Eq. (9):

    sigma^{abc}(w) = (i*pi*e^3/(2*hbar*V)) sum_{k,n,m} f_nm (I^{abc}_{nm} + I^{acb}_{nm}) delta(w - w_nm),
    I^{abc}_{nm}   = r^b_{nm} r^{c;a}_{mn},     r^{c;a} = d_a r^c - i (xi^a_nn - xi^a_mm) r^c

(Esteve-Paredes et al., npj Comput. Mater. 11, 13; hbar = 1 and e = -|e| = -1 in atomic units, applied through FEPS).
Works for any number of bands of a wannier90-style *_tb.dat model (in-plane a,b,c = x,y only), on a k list taken from a
Xatu .states file or an opticx .omesp file. Every derivative is a nested central finite difference in a LOCAL
parallel-transport gauge (each band's phase fixed so that <u(k)|u(k+dk)> > 0), so no global smooth gauge is needed;
degenerate bands are not handled (they are also where the code's own derivative needs clipping).

The output is decomposed as Im(I^{abc}+I^{acb}) = A + B, A = shift-vector part, B = amplitude-gradient part (with
r^b = rho_b e^{i theta_b}: B = sin(th_b - th_c)(rho_b d_a rho_c - rho_c d_a rho_b)), returned only for two-band models where both are well defined (rho > 0).

Usage as a library:  from ipa_shift_numpy import read_tb, k_from_states, k_from_omesp, ipa_shift
Usage as a script:   python3 ipa_shift_numpy.py --tb X_tb.dat --states X.states --nv 1 --eta 0.15 --emin 6.5 --emax 12.5 --nw 600 out.npz
Output units: uA nm/V^2 (same a.u.->SI factor as opticx' shift output), energies in eV.
"""
import argparse
import numpy as np

ANG = 0.52917721067121      # bohr per Angstrom
HA = 27.211385              # eV per Hartree
# a.u. -> uA nm / V^2, times e^3 with e = -|e| = -1 in a.u.: the formulas below are written with e = 1, and opticx
# applies the same e^3 = -1 in constants_math::sigma2_au_to_si (2026-10-04). ipa_shg_numpy imports this.
E_CHARGE = -1.0
FEPS = E_CHARGE**3 * (6.623618e-03) * 1.0e+06 * (27.211386 ** -2) * 5.291772e-11 * 1.0e+09


def read_tb(fn):
    """wannier90-style tb file (first line is skipped, as opticx does). Returns dict with lattice (Angstrom), H[R],
    RP[R] (3,norb,norb) Angstrom, norb, and the cell area in bohr^2 (in-plane a1 x a2)."""
    lines = open(fn).read().split('\n')
    lat = np.array([list(map(float, lines[i].split())) for i in (1, 2, 3)])
    norb, nR = int(lines[4].split()[0]), int(lines[5].split()[0])
    tok = iter(' '.join(lines[6:]).split())
    deg = np.array([int(next(tok)) for _ in range(nR)])
    H, RP = {}, {}
    for iR in range(nR):
        R = (int(next(tok)), int(next(tok)), int(next(tok))); M = np.zeros((norb, norb), complex)
        for _ in range(norb * norb):
            i, j = int(next(tok)) - 1, int(next(tok)) - 1
            M[i, j] = float(next(tok)) + 1j * float(next(tok))
        H[R] = M / deg[iR]
    for iR in range(nR):
        R = (int(next(tok)), int(next(tok)), int(next(tok))); M = np.zeros((3, norb, norb), complex)
        for _ in range(norb * norb):
            i, j = int(next(tok)) - 1, int(next(tok)) - 1
            for c in range(3):
                M[c, i, j] = float(next(tok)) + 1j * float(next(tok))
        RP[R] = M
    area = abs(lat[0, 0] * lat[1, 1] - lat[0, 1] * lat[1, 0]) / ANG ** 2
    return dict(lat=lat, H=H, RP=RP, norb=norb, area=area)


def k_from_states(fn):
    """k list (bohr^-1, in-plane) from a Xatu .states file (stored in Angstrom^-1)."""
    with open(fn) as f:
        n = int(f.readline().split()[0])
        K = np.array([[float(x) for x in f.readline().split()[:2]] for _ in range(n)])
    return K * ANG


def k_from_omesp(fn):
    """k list (bohr^-1, in-plane) from an opticx .omesp stream file."""
    with open(fn, 'rb') as f:
        np.fromfile(f, np.int32, 1); npt, nb = np.fromfile(f, np.int32, 2); npt = int(npt)
        rk = np.fromfile(f, np.float64, 3 * npt).reshape(3, npt)
    return rk[:2].T.copy()


def ipa_shift(tb, K, nv, w_ev, eta_ev, dk=1e-5, decompose=False, clip=None, val_bands=None, cond_bands=None, kres=False):
    """Exact IPA shift conductivity on the k list K (bohr^-1, shape (N,2)). Bands 0..nv-1 valence, the rest conduction.
    clip: if given (bohr), a (k, n, m) pair is dropped when max_{a,c} |r^{c;a}_{nm}| / |r^c_{nm}| exceeds it. This is the analogue of the
    clip_threshold = 50 that opticx applies to shift vectors / Berry connections: near band crossings the finite-difference
    derivatives are garbage (SI Note 7 of the code paper removes degenerate points for the same reason).
    kres: if True return a dict instead: 'pairs' (list of (n,m)), 'isym' (npairs,2,2,2,N) = Im(I^{abc}+I^{acb}) per k-point and pair
    (atomic units, before the delta function and prefactor), 'energies' (N,norb) in Hartree, 'sigma'.
    Returns sigma[w, a, b, c] (uA nm/V^2, in-plane indices) [, A, B spectra for 2-band models]."""
    norb = tb['norb']; H = tb['H']; RP = tb['RP']
    Rv = {R: (R[0] * tb['lat'][0] + R[1] * tb['lat'][1])[:2] / ANG for R in H}      # bohr
    Hk = lambda k: sum(np.exp(1j * (k @ Rv[R]))[:, None, None] * H[R] / HA for R in H)
    Ak = lambda k, c: sum(np.exp(1j * (k @ Rv[R]))[:, None, None] * RP[R][c] / ANG for R in H)
    E2 = np.eye(2)

    def eig(k):
        return np.linalg.eigh(Hk(k))

    def PT(U, Uref):
        p = np.einsum('kab,kab->kb', np.conj(Uref), U)
        return U * (np.conj(p) / abs(p))[:, None, :]

    def rmat(k, Uref):
        U = PT(eig(k)[1], Uref); r = np.zeros((2, len(k), norb, norb), complex)
        for c in range(2):
            Up = PT(eig(k + dk * E2[c])[1], U); Um = PT(eig(k - dk * E2[c])[1], U)
            r[c] = 1j * np.einsum('kan,kam->knm', np.conj(U), (Up - Um) / (2 * dk)) \
                   + np.einsum('kan,kab,kbm->knm', np.conj(U), Ak(k, c), U)
        return r

    w0, U0 = eig(K)
    r0 = rmat(K, U0)
    dr = np.zeros((2, 2, len(K), norb, norb), complex); xi = np.zeros((2, len(K), norb))
    for a in range(2):
        dr[a] = (rmat(K + dk * E2[a], U0) - rmat(K - dk * E2[a], U0)) / (2 * dk)
        Up = PT(eig(K + dk * E2[a])[1], U0); Um = PT(eig(K - dk * E2[a])[1], U0)
        xi[a] = np.real(1j * np.einsum('kan,kan->kn', np.conj(U0), (Up - Um) / (2 * dk))
                        + np.einsum('kan,kab,kbn->kn', np.conj(U0), Ak(K, a), U0))
    cond = list(range(nv, norb)) if cond_bands is None else list(cond_bands)      # optional band windows (0-based), like
    val = list(range(nv)) if val_bands is None else list(val_bands)               # opticx's Bandlist
    w_ha = np.asarray(w_ev) / HA; eta = eta_ev / HA
    sig = np.zeros((len(w_ev), 2, 2, 2)); N = len(K)
    Asp = np.zeros_like(sig); Bsp = np.zeros_like(sig)
    pairs_kres, isym_kres = [], []
    for n in cond:
        for m in val:
            om = w0[:, n] - w0[:, m]
            L = (eta / np.pi) / ((w_ha[:, None] - om[None, :]) ** 2 + eta ** 2)      # (nw, N)
            r = r0[:, :, n, m]                                                       # (c, k)
            rg = np.zeros((2, 2, N), complex)
            for a in range(2):
                for c in range(2):
                    rg[a, c] = dr[a, c][:, n, m] - 1j * (xi[a][:, n] - xi[a][:, m]) * r[c]
            keep = np.ones(N, bool)
            if clip is not None:
                with np.errstate(all='ignore'):
                    ratio = np.max([abs(rg[a, c]) / abs(r[c]) for a in range(2) for c in range(2)], axis=0)
                keep = np.nan_to_num(ratio, nan=np.inf, posinf=np.inf) <= clip
            Iabc = np.zeros((2, 2, 2, N), complex)
            for a in range(2):
                for b in range(2):
                    for c in range(2):
                        Iabc[a, b, c] = r[b] * np.conj(rg[a, c])
            Isym = np.imag(Iabc + np.swapaxes(Iabc, 1, 2)) * keep                    # Re sigma = pi/2 * Im(I+I)
            sig += np.einsum('wk,abck->wabc', L, Isym)
            if kres:
                pairs_kres.append((n, m)); isym_kres.append(np.imag(Iabc + np.swapaxes(Iabc, 1, 2)))
            if decompose and norb == 2:
                rho = abs(r); th = np.angle(r)
                for a in range(2):
                    for b in range(2):
                        for c in range(2):
                            d = dr[a, c][:, n, m]
                            with np.errstate(all='ignore'):
                                drho = np.real(np.conj(r[c]) * d) / rho[c]
                                Rt = np.imag(np.conj(r[c]) * d) / rho[c] ** 2 - (xi[a][:, n] - xi[a][:, m])
                                A = -rho[b] * rho[c] * Rt * np.cos(th[b] - th[c]); B = rho[b] * drho * np.sin(th[b] - th[c])
                                A2 = -rho[c] * rho[b] * (np.imag(np.conj(r[b]) * dr[a, b][:, n, m]) / rho[b] ** 2 - (xi[a][:, n] - xi[a][:, m])) * np.cos(th[c] - th[b])
                                B2 = rho[c] * (np.real(np.conj(r[b]) * dr[a, b][:, n, m]) / rho[b]) * np.sin(th[c] - th[b])
                            Asp[:, a, b, c] += np.einsum('wk,k->w', L, np.nan_to_num(A + A2))
                            Bsp[:, a, b, c] += np.einsum('wk,k->w', L, np.nan_to_num(B + B2))
    pref = FEPS * np.pi / (2 * N * tb['area'])
    if kres:
        return dict(pairs=pairs_kres, isym=np.array(isym_kres), energies=w0 / 1.0, sigma=pref * sig)
    if decompose and norb == 2:
        return pref * sig, pref * Asp, pref * Bsp
    return pref * sig


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--tb', required=True); ap.add_argument('--nv', type=int, required=True)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--states'); g.add_argument('--omesp')
    ap.add_argument('--eta', type=float, default=0.15); ap.add_argument('--emin', type=float, default=6.5)
    ap.add_argument('--emax', type=float, default=12.5); ap.add_argument('--nw', type=int, default=600)
    ap.add_argument('--val-bands', help='comma-separated 0-based valence band indices (default: all below --nv)')
    ap.add_argument('--cond-bands', help='comma-separated 0-based conduction band indices (default: all from --nv)')
    ap.add_argument('--clip', type=float, default=None, help='drop k-pairs with |r^{c;a}|/|r^c| > clip (bohr), e.g. 50')
    ap.add_argument('out')
    a = ap.parse_args()
    tb = read_tb(a.tb); K = k_from_states(a.states) if a.states else k_from_omesp(a.omesp)
    w = a.emin + (a.emax - a.emin) / a.nw * np.arange(a.nw)                          # same grid as opticx
    ints = lambda t: None if t is None else [int(x) for x in t.split(',')]
    res = ipa_shift(tb, K, a.nv, w, a.eta, decompose=(tb['norb'] == 2), clip=a.clip,
                    val_bands=ints(a.val_bands), cond_bands=ints(a.cond_bands))
    if isinstance(res, tuple):
        np.savez(a.out, w=w, sigma=res[0], A=res[1], B=res[2])
    else:
        np.savez(a.out, w=w, sigma=res)
    print(f"wrote {a.out}: {len(K)} k-points, {tb['norb']} bands, nv={a.nv}")


if __name__ == '__main__':
    main()
