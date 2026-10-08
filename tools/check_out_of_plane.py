#!/usr/bin/env python3
"""Out-of-plane second-order response and the absolute sign of the injection current, against a real-time propagation.

Model: non-interacting buckled hBN (buckled-hBN_tb.dat, point group C3v), 30x30, eta = 0.3 eV. Its positions are the
two Wannier centres only (z = -0.5, +0.5 Angstrom), so a hybrid gauge is exact: in-plane fields by minimal coupling,
k -> k + A(t) (electron charge -1), the z field by the length-gauge term +E_z(t) Z with Z = diag(z_i); currents
J_x = -<dH/dk_x>, J_z = -<i[H, Z]> per unit area. Out-of-plane quantities follow Quintela & Pedersen, PRB 107,
235416 (2023) and PRB 110, 085433 (2024): no k_z, the z connection is <n|Z|n>.

 1. Single-particle shift current along z (sigma^zxx, sigma^zzz) equals the excitonic shift-current route in the
    non-interacting limit (2e-2 of the component's maximum) and is not zero. The generalised derivative along the
    non-periodic direction is -i[xi^z, r]; before it was set to zero and every shift current along z vanished.
 2. Linear polarisations x, z and (x+z)/sqrt2 give the symmetric DC tensor D^abc = 2 Re sym sigma(0; w, -w) for
    xxx, zxx, xxz, zzz: the excitonic and the single-particle rectification equal the real-time value to 1e-2 at
    8 eV and 2e-2 at 4 eV (below the gap, where the signal is small). (At small eta the single-particle one-z components converge slowly with the k-mesh; see the docs.)
 3. Circular polarisation (x +- i z)/sqrt2: the difference of the two helicities gives the antisymmetric imaginary
    part A^{x,xz}, which holds the injection current. Its SIGN and size must match the real-time value (2e-2) for
    the excitonic rectification and for both single-particle methods. With term 3 of Taghizadeh & Pedersen
    Eq. (A9b) pairing its field factors with the wrong frequencies, the sign was reversed (covariant and excitonic).

Usage: python3 tools/check_out_of_plane.py --opticx bin/opticx --root . --workdir bin/check_out_of_plane
"""
import argparse, os, subprocess, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ipa_shift_numpy import read_tb, ANG, HA
from check_realtime_sign import mesh, C_AU
from ni_excitons import make_ni

IN = """# Periodic dimensions
2
# Wannier90_filename
{tb}
# Xatu_interface
true
{eig}
{sta}
# Exciton_cutoff
900
# Nfermi
1
# OME_sp
nonlinear
# OME_ex
nonlinear
# Response
{resp}
# Sp_method
{spm}
# Energy_variables
4.0 8.5 0.3 9
"""


class RealTime:
    def __init__(self, tbfile, n=30):
        tb = read_tb(tbfile); lat = tb['lat']; Rs = list(tb['H'])
        self.HR = np.array([tb['H'][R] for R in Rs]) / HA
        RP0 = tb['RP'][(0, 0, 0)]
        for R in Rs:
            for c in range(3):
                M = tb['RP'][R][c] if R in tb['RP'] else 0
                off = M - (np.diag(np.diag(M)) if R == (0, 0, 0) else 0)
                if abs(off).max() > 1e-12:
                    raise SystemExit('check_out_of_plane: positions are not centre-only; the hybrid gauge is not exact')
        tau = np.array([[np.diag(RP0[c]).real[i] for c in range(3)] for i in range(len(RP0[0]))]) / ANG
        Rv = np.array([(R[0] * lat[0] + R[1] * lat[1])[:2] / ANG for R in Rs])
        self.D = Rv[:, None, None, :] + tau[None, None, :, :2] - tau[None, :, None, :2]
        self.Z = np.diag(tau[:, 2]).astype(complex)
        self.K = mesh(lat, n); self.area = abs(np.linalg.det(lat[:2, :2])) / ANG ** 2

    def hk(self, Kp):
        ph = np.exp(1j * np.einsum('kc,rijc->krij', Kp, self.D)) * self.HR[None]
        return ph.sum(1), np.einsum('krij,rijc->ckij', 1j * ph, self.D)

    def J(self, w, eta, E0, phi, s, e_hat, t0f=6.0, dt=0.2):
        """E(t) = s E0 e^{eta t} Re[e_hat e^{-i(w t + phi)}], e_hat = (ex, ez) complex; returns (J_x, J_z) at t = 0."""
        ex, ez = e_hat
        t0 = -t0f / eta; n = int(round(-t0 / dt)); dt = -t0 / n
        _, U = np.linalg.eigh(self.hk(self.K)[0]); v = U[:, :, 0]
        rho = np.einsum('ki,kj->kij', v, v.conj())
        A = lambda t: -s * E0 * np.real(ex * np.exp((eta - 1j * w) * t - 1j * phi) / (eta - 1j * w))
        Ez = lambda t: s * E0 * np.real(ez * np.exp((eta - 1j * w) * t - 1j * phi))
        def H(t): return self.hk(self.K + np.array([A(t), 0.0])[None])[0] + Ez(t) * self.Z[None]
        def f(t, r):
            h = H(t); return -1j * (h @ r - r @ h)
        t = t0
        for _ in range(n):
            k1 = f(t, rho); k2 = f(t + dt / 2, rho + dt / 2 * k1); k3 = f(t + dt / 2, rho + dt / 2 * k2); k4 = f(t + dt, rho + dt * k3)
            rho = rho + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4); t += dt
        h, dH = H(0.0), self.hk(self.K + np.array([A(0.0), 0.0])[None])[1]
        jx = -np.einsum('kij,kji->', rho, dH[0]).real
        jz = -np.einsum('kij,kji->', rho, 1j * (h @ self.Z[None] - self.Z[None] @ h)).real
        return np.array([jx, jz]) / (len(self.K) * self.area)

    def dc(self, w_ev, eta_ev, e_hat, E0=1e-4):
        """DC part of the second-order current (even in s, phases 0 and pi/2 cancel the 2w part)."""
        w, eta = w_ev / HA, eta_ev / HA
        J2 = {p: 0.5 * (self.J(w, eta, E0, p, 1, e_hat) + self.J(w, eta, E0, p, -1, e_hat)) for p in (0.0, np.pi / 2)}
        return 0.5 * (J2[0.0] + J2[np.pi / 2])

    def tensor(self, w_ev, eta_ev, E0=1e-4):
        """With J(2) = 1/4 sum sigma E E: J_DC = 1/2 Re[sigma^abc e_b e_c*] E0^2. Linear: D = 2 Re sym sigma for
        (a; bc) = (x,z; xx, zz, xz); circular: A^{a,xz} = (Im sigma^axz - Im sigma^azx)/2 = (J(+) - J(-))/E0^2."""
        r = 2 ** -0.5
        jx = self.dc(w_ev, eta_ev, (1.0, 0.0), E0); jz = self.dc(w_ev, eta_ev, (0.0, 1.0), E0)
        jd = self.dc(w_ev, eta_ev, (r, r), E0)
        Dxx, Dzz, Dd = 4 * jx / E0 ** 2, 4 * jz / E0 ** 2, 4 * jd / E0 ** 2
        Dxz = Dd - 0.5 * (Dxx + Dzz)
        jp = self.dc(w_ev, eta_ev, (r, 1j * r), E0); jm = self.dc(w_ev, eta_ev, (r, -1j * r), E0)
        A = (jp - jm) / E0 ** 2
        return {'xxx': Dxx[0], 'zxx': Dxx[1], 'zzz': Dzz[1], 'xxz': Dxz[0], 'A_xxz': A[0]}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--opticx', required=True); ap.add_argument('--root', default='.')
    ap.add_argument('--workdir', default='check_out_of_plane')
    a = ap.parse_args()
    root, work, opticx = os.path.abspath(a.root), os.path.abspath(a.workdir), os.path.abspath(a.opticx)
    tb = os.path.join(root, 'buckled-hBN_tb.dat'); sta = os.path.join(root, 'hBN_N30.states')
    for f in (tb, sta):
        if not os.path.exists(f):
            print('SKIP: missing ' + f); sys.exit(0)
    os.makedirs(work, exist_ok=True)
    make_ni(tb, sta, os.path.join(work, 'buckNI30'))
    eig, stn = os.path.join(work, 'buckNI30.eigval'), os.path.join(work, 'buckNI30.states')
    nfail = 0

    def check(label, ok, msg=''):
        nonlocal nfail
        print(('  PASS   ' if ok else '  FAIL   ') + label + ('   ' + msg if msg else '')); nfail += 0 if ok else 1

    def run(resp, spm):
        d = os.path.join(work, f'{resp}_{spm}'); os.makedirs(d, exist_ok=True)
        name = f'buckNI30_{resp}_{spm}_eta0.3.txt'
        open(os.path.join(d, name), 'w').write(IN.format(tb=tb, eig=eig, sta=stn, resp=resp, spm=spm))
        r = subprocess.run([opticx, name], cwd=d, stdout=open(os.path.join(d, 'run.log'), 'w'), stderr=subprocess.STDOUT)
        if r.returncode != 0:
            print(f'FAIL: opticx {resp}/{spm} exited with {r.returncode}'); sys.exit(1)
        return d

    def rect(d, f):
        A = np.loadtxt(os.path.join(d, f), comments='#', ndmin=2)
        return A[:, 0], (A[:, 2::2] + 1j * A[:, 3::2]).reshape(-1, 3, 3, 3)

    def shift(d, f):
        A = np.loadtxt(os.path.join(d, f), comments='#', ndmin=2)
        return A[:, 0], A[:, 1:28].reshape(-1, 3, 3, 3)

    print('Out-of-plane components and the injection sign vs real time (NI buckled hBN 30x30, eta 0.3 eV)\n')
    ds = run('shift', 'covariant')
    _, ssp = shift(ds, 'shift_sp_lengthgauge_buckled-hBN.dat'); _, sex = shift(ds, 'shift_ex_lengthgauge_buckled-hBN.dat')
    for lab, (i, j, k) in (('zxx', (2, 0, 0)), ('zzz', (2, 2, 2))):
        e = abs(ssp[:, i, j, k] - sex[:, i, j, k]).max() / abs(sex[:, i, j, k]).max()
        check(f'sp shift {lab} = excitonic shift route, nonzero', e < 2e-2 and abs(ssp[:, i, j, k]).max() > 0.1 * abs(sex[:, i, j, k]).max(),
              f'max diff {e:.1e} of max (max 2e-2)')

    dc = run('rectification', 'covariant'); dp = run('rectification', 'per_band')
    w, sp = rect(dc, 'second_rectification_lengthgauge_buckled-hBN.dat'); _, ex = rect(dc, 'second_ex_rectification_lengthgauge_buckled-hBN.dat')
    _, pb = rect(dp, 'second_rectification_lengthgauge_buckled-hBN.dat')
    idx = {'x': 0, 'y': 1, 'z': 2}
    rt = RealTime(tb)
    for wv in (8.0, 4.0):
        R = rt.tensor(wv, 0.3); i = np.argmin(abs(w - wv))
        for lab in ('xxx', 'zxx', 'xxz', 'zzz'):
            a_, b_, c_ = (idx[ch] for ch in lab)
            for name, T in (('excitonic', ex), ('single-particle', sp)):
                v = 2 * 0.5 * (T[i, a_, b_, c_] + T[i, a_, c_, b_]).real / C_AU
                e = abs(v - R[lab]) / abs(R[lab]); tol = 1e-2 if wv > 7 else 2e-2       # 4 eV: small off-resonant signal
                check(f'{wv} eV {lab} {name} rectification = real time', e < tol, f'{v:+.4e} vs {R[lab]:+.4e} (rel {e:.1e}, max {tol:.0e})')
        if wv == 8.0:
            for name, T in (('excitonic', ex), ('single-particle covariant', sp), ('single-particle per_band', pb)):
                v = 0.5 * (T[i, 0, 0, 2].imag - T[i, 0, 2, 0].imag) / C_AU
                e = abs(v - R['A_xxz']) / abs(R['A_xxz'])
                check(f'8.0 eV injection channel A^(x,xz), {name}: sign and size = real time', e < 2e-2,
                      f'{v:+.4e} vs {R["A_xxz"]:+.4e} (rel {e:.1e})')
    if nfail:
        print(f'\nFAILED: {nfail} check(s)'); sys.exit(1)
    print('\nALL CHECKS PASSED')


if __name__ == '__main__':
    main()
