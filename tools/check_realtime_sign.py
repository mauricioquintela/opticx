#!/usr/bin/env python3
"""Absolute SIGN and NORMALISATION of opticx's second-order outputs against a real-time propagation.

Every other check compares opticx with itself or with another formula; this one compares with the PHYSICAL current.
Two-band hBN (hBN_tb.dat; position operator = Wannier centres only, so minimal coupling is exact):
    H(k) = sum_R H(R) exp(i k.(R + tau_j - tau_i)),  k -> k + A(t) (electron charge -1, a.u.),
    J = -(1/(N_k A)) sum_k Tr[rho dH/dk(k + A)],  spin factor dropped (as opticx),
E(t) = s E0 cos(w t + phi) e^{eta t} along x, propagated from t0 << 0 to t = 0 (each Fourier component then carries
the +i*eta of the perturbative formulas). With opticx's convention J(2)(t) = 1/4 sum sigma E E e^{-i(w_p+w_q)t} and
E(t) = 1/2 sum E(w_p) e^{-i w_p t}: J2(0) + J2(pi/2) = D E0^2/2 with D = sigma(0;w,-w) + sigma(0;-w,w), and
J2(0) - J2(pi/2) = Re sigma_SHG E0^2, 2 J2(pi/4) - [J2(0) + J2(pi/2)] = Im sigma_SHG E0^2 (J2 = even-in-s part).
Checked (tolerances 1e-3 relative):  SHG (Response = shg, default covariant) and rectification (Response =
rectification) at hw = 8 eV and 4 eV, eta = 0.3 eV, same 30x30 mesh; plus the SIGN of the shift current on resonance.
Validated by flipping the sign back (the pre-2026-10-06 code): every SHG/rectification check fails.
"""
import argparse, os, shutil, subprocess, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ipa_shift_numpy import read_tb, ANG, HA                               # noqa: E402

C_AU = 6.623618e-03 * 1e6 * (27.211386 ** -2) * 5.291772e-11 * 1e9      # |a.u. -> uA nm/V^2| (without the e^3 sign)
INP = """# Periodic dimensions
2
# Wannier90_filename
{tb}
# Xatu_interface
false
# Ncells
30
# Bandlist
0 1
# Nfermi
1
# OME_sp
nonlinear
# Response
{resp}
# Energy_variables
4.0 8.5 {eta} 9
"""


def mesh(lat, n):
    a = lat[:2, :2] / ANG; b = 2 * np.pi * np.linalg.inv(a).T
    u = np.arange(n) / n - 0.5 + (0.5 / n if n % 2 else 0.0)
    U1, U2 = np.meshgrid(u, u, indexing='ij')
    return U1.T.reshape(-1, 1) * b[0] + U2.T.reshape(-1, 1) * b[1]


class RT:
    def __init__(self, tbfile):
        tb = read_tb(tbfile); lat = tb['lat']; Rs = list(tb['H'])
        self.HR = np.array([tb['H'][R] for R in Rs]) / HA
        tau = np.array([[np.diag(tb['RP'][(0, 0, 0)][c]).real[i] for c in range(2)] for i in range(2)]) / ANG
        Rv = np.array([(R[0] * lat[0] + R[1] * lat[1])[:2] / ANG for R in Rs])
        self.D = Rv[:, None, None, :] + tau[None, None, :, :] - tau[None, :, None, :]
        self.K = mesh(lat, 30); self.area = abs(np.linalg.det(lat[:2, :2])) / ANG ** 2

    def hk(self, Kp):
        ph = np.exp(1j * np.einsum('kc,rijc->krij', Kp, self.D)) * self.HR[None]
        return ph.sum(1), np.einsum('krij,rijc->ckij', 1j * ph, self.D)

    def J(self, w, eta, E0, phi, s, t0f=6.0, dt=0.2):
        t0 = -t0f / eta; n = int(round(-t0 / dt)); dt = -t0 / n
        _, U = np.linalg.eigh(self.hk(self.K)[0]); v = U[:, :, 0]
        rho = np.einsum('ki,kj->kij', v, v.conj())
        A = lambda t: -s * E0 * np.real(np.exp((eta - 1j * w) * t - 1j * phi) / (eta - 1j * w))
        def f(t, r):
            H = self.hk(self.K + np.array([A(t), 0.0])[None])[0]; return -1j * (H @ r - r @ H)
        t = t0
        for _ in range(n):
            k1 = f(t, rho); k2 = f(t + dt / 2, rho + dt / 2 * k1); k3 = f(t + dt / 2, rho + dt / 2 * k2); k4 = f(t + dt, rho + dt * k3)
            rho = rho + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4); t += dt
        dH = self.hk(self.K + np.array([A(0.0), 0.0])[None])[1]
        return -np.einsum('kij,ckji->c', rho, dH).real[0] / (len(self.K) * self.area)

    def second(self, w_ev, eta_ev, E0=1e-4):
        w, eta = w_ev / HA, eta_ev / HA
        J2 = {p: 0.5 * (self.J(w, eta, E0, p, 1) + self.J(w, eta, E0, p, -1)) for p in (0.0, np.pi / 4, np.pi / 2)}
        dc = J2[0.0] + J2[np.pi / 2]
        return 2 * dc / E0 ** 2, (J2[0.0] - J2[np.pi / 2]) / E0 ** 2 + 1j * (2 * J2[np.pi / 4] - dc) / E0 ** 2


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--opticx', required=True); ap.add_argument('--root', default='.')
    ap.add_argument('--workdir', default='check_realtime_sign')
    a = ap.parse_args()
    root, work, opticx = os.path.abspath(a.root), os.path.abspath(a.workdir), os.path.abspath(a.opticx)
    tbf = os.path.join(root, 'hBN_tb.dat')
    shutil.rmtree(work, ignore_errors=True); os.makedirs(work)
    nfail = 0

    def check(label, ok, msg=''):
        nonlocal nfail
        print(('  PASS   ' if ok else '  FAIL   ') + label + ('   ' + msg if msg else '')); nfail += 0 if ok else 1

    eta = 0.3
    out = {}
    for resp in ('shg', 'rectification', 'shift'):
        d = os.path.join(work, f'hBN_N30_{resp}_eta{eta}'); os.makedirs(d)
        name = f'hBN_N30_{resp}_eta{eta}.txt'
        open(os.path.join(d, name), 'w').write(INP.format(tb=tbf, resp=resp, eta=eta))
        r = subprocess.run([opticx, name], cwd=d, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if r.returncode != 0:
            check(f'opticx {resp} ran', False, f'rc={r.returncode}')
        for f in os.listdir(d):
            if f.endswith('.omesp'):
                os.remove(os.path.join(d, f))
        out[resp] = d
    ld = lambda d, pre: np.loadtxt([os.path.join(d, f) for f in os.listdir(d) if f.startswith(pre)][0], comments='#')
    shg = ld(out['shg'], 'shg_sp_'); rect = ld(out['rectification'], 'second_rectification_'); sh = ld(out['shift'], 'shift_sp_')
    rt = RT(tbf)
    print(f'real-time reference vs opticx, hBN 30x30, eta = {eta} eV (a.u., physical charge)')
    for w in (8.0, 4.0):
        D, S = rt.second(w, eta)
        i = np.argmin(abs(shg[:, 0] - w)); s_ox = (shg[i, 1] + 1j * shg[i, 2]) / C_AU
        j = np.argmin(abs(rect[:, 0] - w)); d_ox = 2 * rect[j, 2] / C_AU
        e = abs(s_ox - S) / abs(S); check(f'SHG xxx at {w} eV = real time', e < 1e-3, f'{s_ox:.5e} vs {S:.5e} (rel {e:.1e})')
        e = abs(d_ox - D) / abs(D); check(f'rectification xxx at {w} eV = real time', e < 1e-2, f'{d_ox:.5e} vs {D:.5e} (rel {e:.1e})')
        if w == 8.0:
            k = np.argmin(abs(sh[:, 0] - w))
            check('shift current xxx has the sign of the real-time DC current (on resonance)', np.sign(sh[k, 1]) == np.sign(D),
                  f'shift {sh[k, 1] / C_AU:.4e} vs D/2 {D / 2:.4e}')
    if nfail:
        print(f'\nFAILED: {nfail} check(s)'); sys.exit(1)
    print('\nALL CHECKS PASSED')


if __name__ == '__main__':
    main()
