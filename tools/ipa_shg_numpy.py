#!/usr/bin/env python3
"""Independent NumPy evaluation of the single-particle (IPA) SHG conductivity sigma^{abc}(2w;w,w), method
A (length gauge, direct generalized derivative), Taghizadeh, Hipolito & Pedersen, PRB 96, 195413 (2017),
Eq. (A3a), specialised to a cold intrinsic semiconductor (occupations f_n in {0,1}, k-independent), for
which the two terms containing df/dk vanish and only two of the four terms in Eq. (A3a) survive:

  sigma_lab(w) = C * sum_k { sum_{n,m,l; l!=n,l!=m}  p^l_nm (g^a_ln p^b_ml - g^a_ml p^b_ln)
                                                        / (E_ml E_ln [2hw - E_mn])                  (term 1)
                            + sum_{n,m; n!=m}  [ -p^l_nm / (2hw - E_mn) ] * (g^a_mn/E_mn)_{;k^b} }   (term 2)

  g^a_mn(w) = f_nm p^a_mn / (hw - E_mn),   f_nm = f_n - f_m,   E_mn = E_m - E_n,   hw = hbar*w + i*eta
  (O_nm)_{;k^b} = d_b O_nm - i(xi^b_nn - xi^b_mm) O_nm     [generalized derivative, Eq. (6b)]
  C = g_spin * e^3 * hbar^2 / (4 m^3 A_uc), reduced (per unit cell k-sum -> (1/Nk) sum_k) to
  C/A_uc/Nk in atomic units; g_spin is DROPPED (=1), consistent with the rest of opticx.

term 1's exclusion, printed under Eq. (A3a) as "n != l != m", states only that l must differ from BOTH n and
m (E_ml, E_ln, g_ml, g_ln would otherwise be singular); it does NOT forbid n=m (E_mn=0 there is not
singular, and p^lam_nn is then the intraband/group velocity -- a legitimate, generically nonzero
intraband-times-interband contribution; contrast Eq. A3b's 3-band term, whose paper-stated condition is the
genuinely cyclic n!=m!=l!=n and has no n=m piece). An earlier version of this file dropped n=m from term 1 by
analogy with A3b. That piece is not small: on two-band models of lower symmetry (buckled hBN, C3v) it is
44% of the tensor, all of it b<->c antisymmetric (the injection channel); on flat hBN, where injection is
symmetry forbidden, it vanishes. Because n=m survives, term 1 is nonzero already for a
2-band model (l is then just "the other band"); only the n!=m branch of term 1 needs >=3 distinct bands.

(g^a_mn/E_mn)_{;k^b} is expanded with the product rule (f_nm and E_mn are locally-constant/real scalars at
T=0, so their GD reduces to an ordinary k-derivative, Eq. (6a)/(A5b) of the paper):
  (g^a_mn/E_mn)_{;k^b} = f_nm { dphi/dE(E_mn) * (dE_mn/dk^b) * p^a_mn  +  phi(E_mn) * (p^a_mn)_{;k^b} }
  [note the index order p^a_MN, opposite to the p^lam_nm vertex of term 2 -- see the RETRACTED note below]
  phi(E) = 1/[(hw-E)E],  dphi/dE = (1/hw)[1/(hw-E)^2 - 1/E^2],  dE_mn/dk^b = p^b_mm - p^b_nn (Hellmann-Feynman)
(p^a_nm)_{;k^b} is computed via a DIRECT k-derivative (Eq. 6b), NOT opticx's own gen_der_ex_band, which is
the sum-rule form (Eq. 7) and gives an unreliable, near-zero result on small band windows -- see below; opticx's sp SHG code (get_shg_intens_sp, sigma_second_sp.f90) uses a NEW, separately
validated direct-derivative array, vme_der_pt_ex_band, not gen_der_ex_band either.

Momentum matrix elements: p_nn (diagonal) = <u_n|dH/dk|u_n> (Hellmann-Feynman, exact, gauge invariant).
p_nm (off-diagonal, n!=m) is built as -i*E_nm*r_nm from the ALREADY-VALIDATED r_nm (same construction as
ipa_shift_numpy.py's rmat, matched to opticx's X_nm previously). An initial version built p_nm directly from
a velocity kernel plus an ad hoc intra-atomic correction (mirroring opticx's get_vme_eigen_ome, ome_sp.f90);
that matched opticx's own v_nm in magnitude but its SPECIFIC per-band phase inevitably differs from opticx's
own gauge (phase_eigvec_nk) -- comparing individual off-diagonal matrix elements between two different but
individually valid gauges is not meaningful (confirmed directly). The
k-derivative needed for the GD is a finite difference, done in a LOCAL parallel-transport gauge (as in
ipa_shift_numpy.py), so no global smooth gauge is required, but note this tool's OWN internal gauge choice
still generally differs from opticx's -- fine for gauge-invariant quantities (a properly multi-band total),
but NOT meaningful for comparing individual matrix elements against a Fortran debug dump. Degenerate bands
are not handled.

RETRACTED 2026-09-24 (was: "on a strictly 2-band model term 2 alone is not gauge invariant ... the 2-band
SHG total under this formula is gauge-DEPENDENT by construction, so hBN should not be used as a validation
target"). That was wrong, and it was masking two real bugs, both now fixed:
  (1) term 2 used p^a_nm where g^a_mn requires p^a_MN (Eq. 10). The two momentum factors of term 2 must
      carry OPPOSITE band-index order so their gauge phases cancel; with the same order the term scales as
      e^{2i(phi_m-phi_n)} and sigma is gauge dependent. opticx's Fortran kernel shared this exact slip.
  (2) xi_nn (diagonal Berry connection) omitted the Wannier-centre <n|A|n> term that rmat includes.
Term 1 was always gauge invariant: p^lam_nm (g^a_ln p^b_ml - g^a_ml p^b_ln) -> e^{i(phi_m-phi_n)}
e^{i(phi_n-phi_m)} = 1. With both fixed, hBN's 2-band SHG is gauge invariant to machine precision and
satisfies D3h (sigma_xxx=-sigma_xyy=-sigma_yxy=-sigma_yyx) to 6e-6, so hBN IS a valid validation target.
Why it hid for so long: a k-INDEPENDENT gauge change multiplies a 2-band tensor by a common factor, which
cancels in any ratio-based symmetry test (C3 leakage stayed ~0.9%); only opticx's k-dependent
parallel-transport gauge exposed it (leakage 110%). Any test that does not VARY THE GAUGE is blind to this.

Usage as a library:  from ipa_shg_numpy import ipa_shg
Usage as a script:   python3 ipa_shg_numpy.py --tb X_tb.dat --states X.states --nv 1 --eta 0.15 --emin 3 --emax 6 --nw 300 out.npz
Output units: uA nm/V^2 (same convention as opticx' shg output), energies in eV, frequency axis is
hbar*omega (the fundamental, not 2*omega). This now genuinely matches opticx's shg_ex/_sp output
convention: until 2026-09-24 opticx wrote 2*hbar*omega in column 1 and this sentence was wrong, so the two
axes were off by a factor of 2.
"""
import argparse
import numpy as np
from ipa_shift_numpy import read_tb, k_from_states, k_from_omesp, ANG, HA, FEPS

E2 = np.eye(2)


def eig(Hk_fun, k):
    return np.linalg.eigh(Hk_fun(k))


def PT(U, Uref):
    p = np.einsum('kab,kab->kb', np.conj(Uref), U)
    return U * (np.conj(p) / abs(p))[:, None, :]


def ipa_shg(tb, K, nv, w_ev, eta_ev, dk=1e-5, val_bands=None, cond_bands=None, band_window=None):
    """sigma^{abc}(2w;w,w) on the k list K (bohr^-1, shape (N,2)), in-plane a,b,c = x,y. Bands 0..nv-1 filled
    (f=1), the rest empty (f=0), UNLESS band_window is given (list of 0-based band indices): then only those
    bands enter the n,m,l sums (a numerical "Bandlist", to isolate term 1 with a small band count on a
    many-band material), with occupation still set by n < nv. Returns sigma[w,a,b,c] (uA nm/V^2)."""
    norb = tb['norb']; H = tb['H']; RP = tb['RP']
    Rv = {R: (R[0] * tb['lat'][0] + R[1] * tb['lat'][1])[:2] / ANG for R in H}      # bohr

    def Hk(k):
        return sum(np.exp(1j * (k @ Rv[R]))[:, None, None] * H[R] / HA for R in H)

    def dHk(k, c):
        return sum((1j * Rv[R][c]) * np.exp(1j * (k @ Rv[R]))[:, None, None] * H[R] / HA for R in H)

    def Ak(k, c):
        return sum(np.exp(1j * (k @ Rv[R]))[:, None, None] * RP[R][c] / ANG for R in H)

    def rmat(k, Uref):
        # Same construction as ipa_shift_numpy.py's rmat (already validated there against opticx's own
        # X_nm), reused here rather than re-deriving the velocity operator from scratch.
        U = PT(eig(Hk, k)[1], Uref); r = np.zeros((2, len(k), norb, norb), complex)
        for c in range(2):
            Up = PT(eig(Hk, k + dk * E2[c])[1], U); Um = PT(eig(Hk, k - dk * E2[c])[1], U)
            r[c] = 1j * np.einsum('kan,kam->knm', np.conj(U), (Up - Um) / (2 * dk)) \
                   + np.einsum('kan,kab,kbm->knm', np.conj(U), Ak(k, c), U)
        return r

    def pmat(k, Uref):
        # p_nm built from r_nm (off-diagonal, n!=m: p_nm = -i*E_nm*r_nm) and from Hellmann-Feynman
        # (diagonal, n=m: p_nn = dE_n/dk = <n|dH/dk|n>, gauge invariant so no A-kernel correction needed).
        # An initial attempt built p directly from the velocity kernel plus an A-kernel correction term
        # meant to mirror opticx's get_vme_eigen_ome (ome_sp.f90); that matched opticx's own v_nm in
        # magnitude but not sign for n!=m (checked against a direct k-point debug dump), and
        # going through the ALREADY-VALIDATED r_nm (same object ipa_shift_numpy.py uses, matched to
        # opticx's X_nm previously) turned out to have the same sign relative to opticx's v_nm -- i.e. the
        # remaining sign is a genuine convention difference (which sign of e/r opticx's velocity operator
        # uses for the interband part), not a bug in either derivation. Flipped here to match opticx's own
        # convention, the one the SHG formula must ultimately agree with.
        ek, Uraw = eig(Hk, k)
        U = PT(Uraw, Uref)
        r = rmat(k, Uref)
        p = np.zeros((2, len(k), norb, norb), complex)
        for c in range(2):
            p[c] = np.einsum('kan,kab,kbm->knm', np.conj(U), dHk(k, c), U)   # diagonal only used, off-diag overwritten
        Emn = ek[:, :, None] - ek[:, None, :]
        off = ~np.eye(norb, dtype=bool)
        for c in range(2):
            p[c][:, off] = (-1j * Emn * r[c])[:, off]
        return p

    w0, U0 = eig(Hk, K)
    N = len(K)
    p0 = pmat(K, U0)                                     # (2,N,norb,norb)

    # Diagonal Berry connection xi_nn = r_nn, which is i<u_n|d u_n> PLUS the Wannier-centre (A-block)
    # term <n|A|n>, exactly as rmat builds it. FIXED 2026-09-24: this used only i U^dag dU. The A term is
    # what makes the velocity operator C3-covariant (dH/dk alone violates C3 by 67% on hBN), so omitting
    # it here broke the generalized derivative's covariance and hence the C3 symmetry of sigma. opticx's
    # Fortran never had this defect (berry_eigen = berry_eigen1 + berry_eigen2, ome_sp.f90).
    r0 = rmat(K, U0)
    xi = np.zeros((2, N, norb))
    for a in range(2):
        xi[a] = np.real(np.einsum('knn->kn', r0[a]))

    dp = np.zeros((2, 2, N, norb, norb), complex)         # dp[b,a] = d p^a / d k^b (raw, before gauge term)
    for b in range(2):
        dp[b] = (pmat(K + dk * E2[b], U0) - pmat(K - dk * E2[b], U0)) / (2 * dk)
    gdp = np.zeros_like(dp)                               # generalized derivative (p^a_nm)_{;k^b}
    for b in range(2):
        for a in range(2):
            gdp[b, a] = dp[b, a] - 1j * (xi[b][:, :, None] - xi[b][:, None, :]) * p0[a]

    bands = range(norb) if band_window is None else list(band_window)
    val = list(range(nv)) if val_bands is None else list(val_bands)
    cond = list(range(nv, norb)) if cond_bands is None else list(cond_bands)
    f = np.zeros(norb); f[val] = 1.0

    w_ha = np.asarray(w_ev) / HA; eta = eta_ev / HA
    hw = w_ha + 1j * eta                                  # (nw,)
    sig = np.zeros((len(w_ev), 2, 2, 2), complex)

    E = w0                                                # (N,norb)

    # ---- term 2 (2-band): n,m in the band window, n != m ----
    for n in bands:
        for m in bands:
            if m == n:
                continue
            Emn = E[:, m] - E[:, n]                       # (N,)
            fnm = f[n] - f[m]
            if fnm == 0.0:
                continue
            denom2 = 2.0 * hw[:, None] - Emn[None, :]     # (nw,N)
            phi = 1.0 / ((hw[:, None] - Emn[None, :]) * Emn[None, :])          # (nw,N)
            dphi = (1.0 / hw[:, None]) * (1.0 / (hw[:, None] - Emn[None, :]) ** 2 - 1.0 / Emn[None, :] ** 2)
            dEdk = {}
            for b in range(2):
                dEdk[b] = np.real(p0[b, :, m, m] - p0[b, :, n, n])             # (N,)
            for a in range(2):
                for b in range(2):
                    # g^a_mn = f_nm p^a_MN/(hw-E_mn) (Eq. 10): index order (m,n), OPPOSITE to the
                    # p^lam_nm vertex below, so their gauge phases cancel. Using p^a_nm in both places
                    # (as this file did until 2026-09-24) makes sigma gauge dependent -- see module docstring.
                    gd_h = fnm * (dphi * dEdk[b][None, :] * p0[a, :, m, n][None, :]
                                  + phi * gdp[b, a, :, m, n][None, :])          # (nw,N)
                    for lam in range(2):
                        term2 = -p0[lam, :, n, m][None, :] / denom2 * gd_h
                        sig[:, lam, a, b] += term2.sum(axis=1)

    # ---- term 1: l must differ from both n and m (E_ml, E_ln, g_ml, g_ln would otherwise be singular/
    # ill-defined), but n=m IS allowed (E_mn=0 there is not singular; p^lam_nn is then the intraband/group
    # velocity, giving a legitimate intraband-times-interband contribution). The paper's "n != l != m" under
    # Eq. (A3a) states exactly these two conditions (contrast Eq. A3b, whose 3-band term is over a genuinely
    # cyclic n!=m!=l!=n and needs no n=m case). Dropping n=m here silently broke the C3 (D3h) relation
    # sigma_xxx=-sigma_xyy=-sigma_yxy on hBN by ~80% -- caught by that symmetry check, not by inspection.
    for n in bands:
        for m in bands:
            Emn = E[:, m] - E[:, n]
            denom_outer = 2.0 * hw[:, None] - Emn[None, :]                     # (nw,N)
            for l in bands:
                if l == n or l == m:
                    continue
                Eml = E[:, m] - E[:, l]; Eln = E[:, l] - E[:, n]
                fln = f[n] - f[l]; fml = f[l] - f[m]
                if fln == 0.0 and fml == 0.0:
                    continue
                g_ln = fln * p0[:, None, :, l, n] / (hw[None, :, None] - Eln[None, None, :])   # (2,nw,N)
                g_ml = fml * p0[:, None, :, m, l] / (hw[None, :, None] - Eml[None, None, :])   # (2,nw,N)
                denom13 = Eml[None, :] * Eln[None, :] * denom_outer            # (nw,N)
                for a in range(2):
                    for b in range(2):
                        bracket = g_ln[a] * p0[b, :, m, l][None, :] - g_ml[a] * p0[b, :, l, n][None, :]
                        for lam in range(2):
                            term1 = p0[lam, :, n, m][None, :] * bracket / denom13
                            sig[:, lam, a, b] += term1.sum(axis=1)

    Auc = tb['area']                                       # bohr^2 (in-plane cell area)
    # C_ee = C_ie = 1/4 in these units (e = hbar = m = 1), i.e. the 1/4 of C = e^3 hbar^2/(4 m^3 A_uc)
    # quoted in the module docstring, combined with the usual 1/(Nk*A_uc) BZ-sum discretisation.
    # ADDED 2026-09-24: this factor was missing, making the tool 4x larger than opticx. It could not be
    # noticed earlier because both codes were gauge dependent (see the RETRACTED note above),
    # so their absolute values were not comparable at all. With both fixed they now agree to ~1e-7.
    # x (-4) since 2026-10-06: opticx-wide convention J(2) = 1/4 sum sigma E E (Taghizadeh & Pedersen
    # 2018 / npj; Eq. A3a is written for J(2) = sum sigma E E) and the physical sign (Eq. A3a already contains the
    # electron charge, FEPS applies e^3 = -1 on top). Matches opticx's shg_sp output and a real-time propagation.
    pref = -4.0 * 0.25 * FEPS / (Auc * N)
    sig = pref * sig
    # Eq. (A3a) is not manifestly symmetric under the two (identical-frequency) field indices b<->c; the
    # paper states it must be symmetrised by permutation (text below Eq. A4). Same treatment as opticx's
    # excitonic SHG driver (get_sigma_shg_ex).
    return 0.5 * (sig + np.swapaxes(sig, 2, 3))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--tb', required=True); ap.add_argument('--nv', type=int, required=True)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--states'); g.add_argument('--omesp')
    ap.add_argument('--eta', type=float, default=0.15); ap.add_argument('--emin', type=float, default=3.0)
    ap.add_argument('--emax', type=float, default=6.0); ap.add_argument('--nw', type=int, default=300)
    ap.add_argument('--val-bands', help='comma-separated 0-based valence band indices (default: all below --nv)')
    ap.add_argument('--cond-bands', help='comma-separated 0-based conduction band indices (default: all from --nv)')
    ap.add_argument('--band-window', help='comma-separated 0-based band indices entering the n,m,l sums (default: all)')
    ap.add_argument('out')
    a = ap.parse_args()
    tb = read_tb(a.tb); K = k_from_states(a.states) if a.states else k_from_omesp(a.omesp)
    w = a.emin + (a.emax - a.emin) / a.nw * np.arange(a.nw)
    ints = lambda t: None if t is None else [int(x) for x in t.split(',')]
    sig = ipa_shg(tb, K, a.nv, w, a.eta, val_bands=ints(a.val_bands), cond_bands=ints(a.cond_bands),
                  band_window=ints(a.band_window))
    np.savez(a.out, w=w, sigma=sig)
    print(f"wrote {a.out}: {len(K)} k-points, {tb['norb']} bands, nv={a.nv}")


if __name__ == '__main__':
    main()
