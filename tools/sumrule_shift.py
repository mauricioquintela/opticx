#!/usr/bin/env python3
"""Single-particle shift current from the velocity sum rule of Kameda, Nagashio and Wakabayashi, arXiv:2611.05990,
Eqs. 23-25, generalised to Wannier models with a full position matrix and to degenerate bands. An evaluation that
shares nothing with opticx but the tight-binding file: no k-mesh derivatives, every derivative analytic.

Operators in the orbital basis (lattice gauge: H(k) = sum_R e^{ik.R} H(R), A_c(k) = sum_R e^{ik.R} r_c(R)), r_a = i d_a + A_a:
    V_b  = -i [r_b, H]   = d_b H + i [H, A_b]
    W_ab = -i [r_a, V_b] = d_a V_b - i [A_a, V_b]            (d_a = 0 along a non-periodic direction)
For centre-only positions (A diagonal, constant) these are the paper's d_b H and d_a d_b H, in the gauge with the
orbital positions in the Bloch phases. Band basis, blocks B(n) of exactly degenerate bands, Ra_nl = v^a_nl/E_nl off-block:
    (v^b_nm);a = W^ab_nm + sum_{l not in B(n)} v^a_nl v^b_lm / E_nl - sum_{l not in B(m)} v^b_nl v^a_lm / E_lm   [Eq. 25]
    (r^b_nm);a = [ -i (v^b_nm);a - (Vd^a r^b - r^b Vd^a)_nm ] / E_nm ,   Vd^a = block-diagonal part of v^a
(per band, Vd is the diagonal and the last term is the paper's v^j_nm (v^i_mm - v^i_nn)/E_nm). [r_a, r_b] = 0 is not assumed.
    sigma^abc(0; w, -w) = 1/(2 Nk A) sum_k sum_nm f_nm [r^b_mn (r^c_nm);a + r^c_mn (r^b_nm);a] / (E_mn - w - i eta)  [Eq. 23]
    shift: the same with 1/(E_mn - w - i eta) -> i pi delta_eta(w - E_mn)       (npj Comput. Mater. 11, 13, Eq. 9)
Atomic units, e = 1: opticx's shift_sp_lengthgauge file is C_AU times the real part. n, m run over the window, the
intermediate sum over ALL bands of the model. Dictionary to npj / Taghizadeh & Pedersen notation: the preprint's
I^{jk}_mn R^{ik}_nm = i eps_mn^2 r^b_mn r^c_{nm;a}, with R^{ik}_nm = i r^c_{nm;a}/r^c_nm = R^{a,c}_nm + i d_a ln|r^c_nm|
(R^{a,c} the shift vector of Ahn, Guo and Nagaosa, PRX 10, 041041). W, the second derivative of H, is nonzero for
every lattice model; with it the sum rule is exact for a complete model, without it only with many remote bands.

coherence_shift: the part of the shift current that the excitonic shift route (npj SI Eqs. 9-11) does not contain
without interaction. In Eq. 25 the same-type three-band pieces (valence l for n = v, conduction l for m = c) carry
1/eps_vl or 1/eps_lc; by partial fractions
    1/eps_cl = eps_cv/(eps_lv eps_cl) + 1/eps_lv          1/eps_vl = eps_cv/(eps_cl eps_vl) + 1/eps_cl
the first part is a coherence between the non-interacting excitons (v,c,k) and (l,c,k) or (v,l,k) -- term 3 of
Taghizadeh & Pedersen PRB 97, 205432, Eq. (B1a), cancelled by the shift-current convention -- and the second a
transition at the other exciton's energy, which terms 1-2 keep. Window bands only, no degenerate pairs.

The 2D SSH model of the preprint (write_ssh_tb) is a complete four-band model, gap 1.09 eV, the two valence bands
>= 2.9 eV and the two conduction bands >= 2.4 eV apart at every k (no near-degenerate pairs).

Usage: python3 tools/sumrule_shift.py model_tb.dat out.npz --nfermi 26 --bandlist '-1 0 1 2' --ev 0 5 0.05 1000
"""
import os, sys, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ipa_shift_numpy import read_tb, ANG, HA
from check_realtime_sign import mesh


def model(tbfile, centres_only=False):
    """H(R) (Ha), Hermitian position matrices A(R) (bohr) and lattice vectors R (bohr, in-plane) of a 2D model.
    centres_only keeps the diagonal of r(R = 0) alone: the paper's Eqs. as printed."""
    tb = read_tb(tbfile); lat = tb['lat']; Rs = sorted(tb['H']); norb = tb['norb']
    H = np.array([tb['H'][R] for R in Rs])/HA
    zero = np.zeros((3, norb, norb), complex)
    A = np.array([0.5*(tb['RP'].get(R, zero) + tb['RP'].get(tuple(-x for x in R), zero).conj().transpose(0, 2, 1))
                  for R in Rs])/ANG
    if centres_only:
        i0 = Rs.index((0, 0, 0)); A0 = np.array([np.diag(np.diag(A[i0, c]).real) for c in range(3)])
        A = np.zeros_like(A); A[i0] = A0
    Rv = np.array([R[0]*lat[0] + R[1]*lat[1] for R in Rs])/ANG; Rv[:, 2] = 0.0          # 2D: no d/dk_z
    return dict(H=H, A=A, Rv=Rv, lat=lat, norb=norb, area=abs(np.linalg.det(lat[:2, :2]))/ANG**2)


def comm(X, Y):
    return X @ Y - Y @ X


def kpoint_ops(M, k):
    """H(k), V_b(k), W_ab(k) in the orbital basis."""
    ph = np.exp(1j*(M['Rv'][:, :2] @ k)); iR = 1j*M['Rv']
    Hk = np.einsum('r,rij->ij', ph, M['H']); Ak = np.einsum('r,rcij->cij', ph, M['A'])
    dH = np.einsum('r,ra,rij->aij', ph, iR, M['H']); dA = np.einsum('r,ra,rcij->acij', ph, iR, M['A'])
    ddH = np.einsum('r,ra,rb,rij->abij', ph, iR, iR, M['H'])
    V = np.array([dH[b] + 1j*comm(Hk, Ak[b]) for b in range(3)])
    Wm = np.empty((3, 3) + Hk.shape, complex)
    for a in range(3):
        for b in range(3):
            dV = ddH[a, b] + 1j*comm(dH[a], Ak[b]) + 1j*comm(Hk, dA[a, b])
            Wm[a, b] = dV - 1j*comm(Ak[a], V[b])
    return 0.5*(Hk + Hk.conj().T), V, Wm


def gen_der(Hk, V, Wm, tol):
    """Band energies E, band velocities v[a]_nm, r[b]_nm (off-block) and rg[a,b]_nm = (r^b_nm);a."""
    E, U = np.linalg.eigh(Hk); n = len(E)
    lab = np.concatenate([[0], np.cumsum(np.diff(E) > tol)]); same = lab[:, None] == lab[None, :]
    dE = E[:, None] - E[None, :]; inv = np.where(same, 0.0, 1.0/np.where(same, 1.0, dE))
    v = np.array([U.conj().T @ V[a] @ U for a in range(3)])
    Ra = v*inv[None]; r = -1j*Ra; Vd = np.where(same[None], v, 0.0)
    rg = np.empty((3, 3, n, n), complex)
    for a in range(3):
        for b in range(3):
            gv = U.conj().T @ Wm[a, b] @ U + Ra[a] @ v[b] - v[b] @ Ra[a]
            rg[a, b] = (-1j*gv - (Vd[a] @ r[b] - r[b] @ Vd[a]))*inv
    return E, v, r, rg


def coherence_gen_der(E, v, val, con):
    """The coherence part of (r^c_nm);a for n = v in val, m = c in con (zero elsewhere); see the module docstring."""
    rgc = np.zeros((3, 3, len(E), len(E)), complex)
    for iv in val:
        for ic in con:
            Ecv = E[ic] - E[iv]
            for a in range(3):
                for c in range(3):
                    g = sum(v[a, iv, l]*v[c, l, ic]/(E[iv] - E[l])*Ecv/(E[ic] - E[l]) for l in val if l != iv) \
                      - sum(v[c, iv, l]*v[a, l, ic]/(E[l] - E[ic])*Ecv/(E[l] - E[iv]) for l in con if l != ic)
                    rgc[a, c, iv, ic] = -1j*g/(E[iv] - E[ic])
    return rgc


def spectra(M, ncells, nfermi, window, w, eta, tol=1e-6, scissor=0.0, coherence=False):
    """window: 0-based band indices. Returns (shift, causal) sigma^abc(w) in a.u. (e = 1). The shift part keeps the
    resonant delta(w - eps_cv) only, as opticx writes it (the antiresonant one sits at w < 0; Lorentzian tails of it
    reach 5e-3 of the SSH spectrum at eta 0.1 eV). With coherence = True, the coherence parts only."""
    K = mesh(M['lat'], ncells); nw = len(w); win = np.array(window)
    val = win[win < nfermi]; con = win[win >= nfermi]
    S_shift = np.zeros((3, 3, 3, nw), complex); S_caus = np.zeros((3, 3, 3, nw), complex)
    for k in K:
        E, v, r, rg = gen_der(*kpoint_ops(M, k), tol)
        pairs = ((val, con, 1.0), (con, val, -1.0))
        if coherence:
            rg = coherence_gen_der(E, v, val, con); pairs = pairs[:1]
        for (ns, ms, f) in pairs:                                                       # f_nm = f_n - f_m
            Emn = (E[ms][None, :] - E[ns][:, None]) + f*scissor                          # E_m - E_n, (n, m)
            rmn = r[:, ms][:, :, ns].transpose(0, 2, 1)                                  # r^b_mn as (b, n, m)
            rgnm = rg[:, :, ns][:, :, :, ms]                                             # (a, c, n, m)
            I = np.einsum('bnm,acnm->abcnm', rmn, rgnm); I = I + I.transpose(0, 2, 1, 3, 4)
            x = Emn.ravel(); Iv = I.reshape(3, 3, 3, -1)
            den = 1.0/(x[:, None] - w[None, :] - 1j*eta)                                 # (pairs, nw)
            S_caus += f*np.einsum('abcp,pw->abcw', Iv, den)
            if f > 0:
                S_shift += np.einsum('abcp,pw->abcw', Iv, 1j*den.imag)                   # i pi delta_eta = i Im 1/(x - w - i eta)
    norm = 1.0/(2*len(K)*M['area'])
    return norm*S_shift, norm*S_caus


def write_ssh_tb(path, a=4.0, g=1.7, gp=0.3, buckle=0.0):
    """2D SSH model of arXiv:2611.05990, single domain 1: sites A (0,0), B (a/2,0), C (0,a/2), D (a/2,a/2) (Angstrom),
    intracell -gamma and intercell -gamma' along x (A-B, C-D) and y (A-C, B-D), on-site 2/3, -2, 2/3, 2/3 eV;
    orbital centres only. buckle moves the centres to z = +b, -b, -b, +b (H unchanged): out-of-plane components.
    Wannier90 tb format."""
    tau = np.array([[0, 0, buckle], [a/2, 0, -buckle], [0, a/2, -buckle], [a/2, a/2, buckle]])
    Rs = [(0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)]
    H = {R: np.zeros((4, 4)) for R in Rs}
    A, B, C, D = range(4)
    for i, e in enumerate((2/3, -2.0, 2/3, 2/3)):
        H[(0, 0, 0)][i, i] = e
    for (p, q, d) in ((A, B, (1, 0, 0)), (C, D, (1, 0, 0)), (A, C, (0, 1, 0)), (B, D, (0, 1, 0))):
        H[(0, 0, 0)][p, q] = H[(0, 0, 0)][q, p] = -g                                   # p on the low side, q on the high side
        H[d][q, p] = -gp; H[tuple(-x for x in d)][p, q] = -gp                           # q (cell 0) - p (cell +d)
    with open(path, 'w') as f:
        f.write('2D SSH model, single domain 1 (arXiv:2611.05990)' + (f', centres buckled +-{buckle} A' if buckle else '') + '\n')
        f.write(f'{a:.10f} 0.0 0.0\n0.0 {a:.10f} 0.0\n0.0 0.0 20.0\n4\n{len(Rs)}\n' + ' '.join('1' for _ in Rs) + '\n')
        for R in Rs:
            f.write(f'\n{R[0]} {R[1]} {R[2]}\n')
            for j in range(4):
                for i in range(4):
                    f.write(f'{i+1} {j+1} {H[R][i, j]:.16e} 0.0\n')
        for R in Rs:
            f.write(f'\n{R[0]} {R[1]} {R[2]}\n')
            for j in range(4):
                for i in range(4):
                    x = tau[i] if (R == (0, 0, 0) and i == j) else np.zeros(3)
                    f.write(f'{i+1} {j+1} ' + ' '.join(f'{x[c]:.16e} 0.0' for c in range(3)) + '\n')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('tb'); ap.add_argument('out'); ap.add_argument('--ncells', type=int, default=30)
    ap.add_argument('--nfermi', type=int, required=True)
    ap.add_argument('--bandlist', required=True, help='offsets from Nfermi, as opticx')
    ap.add_argument('--ev', nargs=4, type=float, required=True, metavar=('E1', 'E2', 'ETA', 'NW'))
    ap.add_argument('--centres-only', action='store_true'); ap.add_argument('--tol', type=float, default=1e-6)
    ap.add_argument('--scissor', type=float, default=0.0)
    a = ap.parse_args()
    e1, e2, eta, nw = a.ev; nw = int(nw); w = (e1 + (e2 - e1)/nw*np.arange(nw))/HA
    window = [a.nfermi - 1 + int(o) for o in a.bandlist.split()]
    sh, ca = spectra(model(a.tb, a.centres_only), a.ncells, a.nfermi, window, w, eta/HA, a.tol, a.scissor/HA)
    np.savez(a.out, w=w*HA, shift=sh, causal=ca)
