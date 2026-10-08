#!/usr/bin/env python3
"""The excitonic rectification sigma(0; w, -w): the whole causal response, injection current included.

Response = rectification evaluates Taghizadeh & Pedersen PRB 97, 205432 (2018), Eq. (B1a) (method A) with the
causal prescription (w + i eta, -w + i eta, w_2 = 2i eta), both field orderings explicit, and term 3 (the
exciton populations and coherences, Eq. A9b) on the bare exciton-exciton current. Re is the b<->c symmetric DC
response, Im the antisymmetric one, which holds the injection current. The shift current in the convention of
npj Comput. Mater. 11, 13 is Response = shift, a separate route not checked here.

 1. Non-interacting buckled hBN (C3v, injection allowed), 30x30: the in-plane block of Re equals the
    single-particle rectification (whose formula is independent) to < 5e-3 of its maximum. Components with a z
    index are NOT compared: there the single-particle rectification differs from every other route at the
    8.59 eV peak (xxz 1.829 against 1.598 from the excitonic rectification and from both shift-current routes;
    zxx 1.632 against 2.098 from both excitonic routes, the single-particle shift formula giving 0 for any
    current along z), an open single-particle question unrelated to this branch. The integrated weight of Im (the injection
    current) equals the single-particle one to 2% and eta x weight is constant to 2% (a rate times
    tau = hbar/2 eta). Line shapes differ at finite eta, so Im is compared by weight.
 2. Re is exactly b<->c symmetric and Im antisymmetric (the reality of the current; 1e-12).
 3. Non-interacting flat hBN (D3h, injection forbidden): max|Im| < 1e-4 of the buckled one.
 4. The formula re-evaluated densely in NumPy from the cached matrix elements reproduces the Fortran output
    (Re and Im, 24 frequencies) to 1e-7, on real hBN excitons (full off-diagonal P_nm) and on the
    non-interacting buckled set, with the SAME unit factor in both (same mesh and cell).
 5. Ex_rectification is obsolete: 'causal' is accepted with a note, 'shift' stops and names Response = shift.

Usage: python3 tools/check_ex_rectification.py --opticx bin/opticx --root . --workdir bin/check_ex_rectification
"""
import argparse, os, struct, subprocess, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ni_excitons import make_ni
import plot_test_outputs

HA = 27.211385
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
# Energy_variables
{grid}
# Cache_ome_ex
{cache}
{extra}"""


def read_cache(fn):
    """ome_second_ex_<mat>.omeex2, format v3: E (Ha), X_n (3,N), X_nm (3,N,N), P_nm (3,N,N)."""
    b = open(fn, 'rb').read(); p = 14
    ver, L = struct.unpack('<ii', b[p:p+8]); p += 8 + L
    nk, nv, nc, N = struct.unpack('<4i', b[p:p+16]); p += 16
    nb, = struct.unpack('<i', b[p:p+4]); p += 4 + 4*nb + 4
    E = np.frombuffer(b, '<f8', N, p).copy(); p += 8*N
    X = np.frombuffer(b, '<c16', 3*N, p).reshape(N, 3).T.copy(); p += 2*16*3*N
    Xi = np.frombuffer(b, '<c16', 3*N*N, p).reshape(N, N, 3).transpose(2, 1, 0).copy(); p += 16*3*N*N
    P = np.frombuffer(b, '<c16', 3*N*N, p).reshape(N, N, 3).transpose(2, 1, 0).copy()
    return E, X, Xi, P


def rect_numpy(E, X, Xi, P, w_ev, eta_ev):
    """Causal Eq. (B1a), term 3 on P, symmetrised over the two orderings; no prefactor/units, before the sign."""
    Pi1 = -1j*E[None, :]*X
    def S(zp, zq):
        z2 = zp + zq
        out = np.zeros((len(zp), 3, 3, 3), complex)
        A1 = Pi1[:, :, None]/(z2[None, None] - E[None, :, None]); B1 = X.conj()[:, :, None]/(zq[None, None] - E[None, :, None])
        A2 = Pi1.conj()[:, :, None]/(z2[None, None] + E[None, :, None]); B2 = X[:, :, None]/(zq[None, None] + E[None, :, None])
        A3 = X[:, :, None]/(zq[None, None] + E[None, :, None]); B3 = X.conj()[:, :, None]/(zp[None, None] - E[None, :, None])
        for b in range(3):
            for c in range(3):
                W1 = Xi[b] @ B1[c]; W2 = Xi[b].conj() @ B2[c]
                out[:, :, b, c] += (np.einsum('anw,nw->wa', A1, W1) + np.einsum('anw,nw->wa', A2, W2))
        for a in range(3):                          # term 3, Eq. (A9b) pairing: X_n (w_q + E_n) is the c field,
            for c in range(3):                      # X*_m (w_p - E_m) the b field
                W3 = P[a] @ B3[c]
                out[:, a, c, :] -= np.einsum('bnw,nw->wb', A3, W3)
        return out
    zp = (w_ev + 1j*eta_ev)/HA; zq = (-w_ev + 1j*eta_ev)/HA
    A, B = S(zp, zq), S(zq, zp)
    return 0.5*(A + B.transpose(0, 1, 3, 2))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--opticx', required=True); ap.add_argument('--root', default='.')
    ap.add_argument('--workdir', default='check_ex_rectification')
    a = ap.parse_args()
    root, work, opticx = os.path.abspath(a.root), os.path.abspath(a.workdir), os.path.abspath(a.opticx)
    for f in ('hBN_tb.dat', 'buckled-hBN_tb.dat', 'hBN_N30.states', 'hBN_N30.eigval'):
        if not os.path.exists(os.path.join(root, f)):
            print('SKIP: missing ' + os.path.join(root, f)); sys.exit(0)
    os.makedirs(work, exist_ok=True)
    nfail = 0

    def check(label, ok, msg=''):
        nonlocal nfail
        print(('  PASS   ' if ok else '  FAIL   ') + label + ('   ' + msg if msg else '')); nfail += 0 if ok else 1

    def cache_name(tb):
        return 'ome_second_ex_' + os.path.basename(tb).replace('_tb.dat', '') + '.omeex2'

    def run(tag, tb, eig, sta, grid, cache='off', link=None, extra='', resp='rectification', expect_ok=True):
        d = os.path.join(work, tag); os.makedirs(d, exist_ok=True)
        if link:
            src = os.path.join(work, link, cache_name(tb)); dst = os.path.join(d, cache_name(tb))
            if not os.path.exists(dst): os.symlink(src, dst)
        name = f'{tag}.txt'
        open(os.path.join(d, name), 'w').write(IN.format(tb=tb, eig=eig, sta=sta, resp=resp, grid=grid, cache=cache, extra=extra))
        r = subprocess.run([opticx, name], cwd=d, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log = r.stdout.decode('utf-8', 'replace'); open(os.path.join(d, 'run.log'), 'w').write(log)
        if expect_ok and r.returncode != 0:
            print(f'FAIL: opticx run {tag} exited with {r.returncode} (see {d}/run.log)'); sys.exit(1)
        return d, r.returncode, log

    def load(d, fn):
        A = np.loadtxt(os.path.join(d, fn), comments='#')
        return A[:, 0], (A[:, 2::2] + 1j*A[:, 3::2]).reshape(-1, 3, 3, 3)

    asym = lambda T: 0.5*(T - T.transpose(0, 1, 3, 2))
    sym = lambda T: 0.5*(T + T.transpose(0, 1, 3, 2))
    tb_b, tb_f = os.path.join(root, 'buckled-hBN_tb.dat'), os.path.join(root, 'hBN_tb.dat')
    sta = os.path.join(root, 'hBN_N30.states')
    make_ni(tb_b, sta, os.path.join(work, 'buckNI30')); make_ni(tb_f, sta, os.path.join(work, 'flatNI30'))
    eigb, stab = os.path.join(work, 'buckNI30.eigval'), os.path.join(work, 'buckNI30.states')
    eigf, staf = os.path.join(work, 'flatNI30.eigval'), os.path.join(work, 'flatNI30.states')
    grid = '3.0 19.0 {eta} 1601'

    print('Excitonic rectification: the whole causal sigma(0; w, -w), injection current in Im\n')
    d05, _, _ = run('buckNI30_rect_eta0.05', tb_b, eigb, stab, grid.format(eta=0.05), 'write')
    d10, _, _ = run('buckNI30_rect_eta0.1', tb_b, eigb, stab, grid.format(eta=0.1), 'read', link='buckNI30_rect_eta0.05')
    w, ex05 = load(d05, 'second_ex_rectification_lengthgauge_buckled-hBN.dat')
    _, sp05 = load(d05, 'second_rectification_lengthgauge_buckled-hBN.dat')
    _, ex10 = load(d10, 'second_ex_rectification_lengthgauge_buckled-hBN.dat')
    _, sp10 = load(d10, 'second_rectification_lengthgauge_buckled-hBN.dat')
    dw = w[1] - w[0]; c = (1, 1, 2)
    ip = (slice(None), slice(0, 2), slice(0, 2), slice(0, 2))      # in-plane block, see the docstring
    e = abs(ex05.real[ip] - sym(sp05.real)[ip]).max()/abs(sp05.real[ip]).max()
    check('NI buckled hBN: Re = single-particle rectification, in-plane (eta 0.05)', e < 5e-3, f'{e:.1e} (max 5e-3)')
    W = lambda T: T[(slice(None),) + c].sum()*dw
    r05 = W(ex05.imag)/W(asym(sp05.imag)); r10 = W(ex10.imag)/W(asym(sp10.imag))
    check('NI buckled hBN: injection weight (Im) = single-particle', abs(r05 - 1) < 0.02 and abs(r10 - 1) < 0.02,
          f'ratio {r05:.4f} (eta 0.05), {r10:.4f} (eta 0.1) (want 1 +- 0.02)')
    t = (0.05*W(ex05.imag))/(0.1*W(ex10.imag))
    check('NI buckled hBN: eta x injection weight constant (~ tau)', abs(t - 1) < 0.02, f'{t:.4f} (want 1 +- 0.02)')
    rd = max(abs(asym(ex05.real)).max(), abs(sym(ex05.imag)).max())/abs(ex05).max()
    check('Re b<->c symmetric, Im antisymmetric (reality of the current)', rd < 1e-12, f'{rd:.1e} (max 1e-12)')

    df, _, _ = run('flatNI30_rect_eta0.05', tb_f, eigf, staf, grid.format(eta=0.05))
    _, fl = load(df, 'second_ex_rectification_lengthgauge_hBN.dat')
    q = abs(fl.imag).max()/abs(ex05.imag).max()
    check('NI flat hBN (D3h): injection forbidden', q < 1e-4, f'max|Im| flat/buckled {q:.1e} (max 1e-4)')

    idx = np.arange(0, len(w), len(w)//24)[:24]
    E, X, Xi, P = read_cache(os.path.join(d05, cache_name(tb_b)))
    ref = -rect_numpy(E, X, Xi, P, w[idx], 0.05)                 # causal outputs carry the physical sign (-1)
    s_ni = np.vdot(ref.ravel(), ex05[idx].ravel()).real/np.vdot(ref.ravel(), ref.ravel()).real
    res_ni = np.linalg.norm(ex05[idx] - s_ni*ref)/np.linalg.norm(ex05[idx])
    check('NI buckled hBN: Fortran = dense NumPy evaluation', res_ni < 1e-7, f'residual {res_ni:.1e} (max 1e-7)')

    eig_r = os.path.join(root, 'hBN_N30.eigval')
    dr, _, _ = run('hBN_N30_rect_eta0.1', tb_f, eig_r, sta, '4.0 12.0 0.1 801', 'write')
    wr, re = load(dr, 'second_ex_rectification_lengthgauge_hBN.dat')
    E, X, Xi, P = read_cache(os.path.join(dr, cache_name(tb_f)))
    idr = np.arange(0, len(wr), len(wr)//24)[:24]
    ref = -rect_numpy(E, X, Xi, P, wr[idr], 0.1)
    s = np.vdot(ref.ravel(), re[idr].ravel()).real/np.vdot(ref.ravel(), ref.ravel()).real
    res = np.linalg.norm(re[idr] - s*ref)/np.linalg.norm(re[idr])
    check('real hBN 30x30: Fortran = dense NumPy evaluation (full P_nm, Re and Im)', res < 1e-7, f'residual {res:.1e} (max 1e-7)')
    check('same unit factor (-1/(Nk V) x a.u.->SI) in both comparisons', abs(s/s_ni - 1) < 1e-6, f'ratio {s/s_ni:.8f}')

    _, rc, log = run('keyword_causal', tb_f, eig_r, sta, '4.0 12.0 0.1 41', 'read', link='hBN_N30_rect_eta0.1',
                     extra='# Ex_rectification\ncausal\n')
    check("Ex_rectification = causal: accepted with a note", rc == 0 and 'obsolete' in log.lower(), f'rc={rc}')
    _, rc, log = run('keyword_shift', tb_f, eig_r, sta, '4.0 12.0 0.1 41', 'read', link='hBN_N30_rect_eta0.1',
                     extra='# Ex_rectification\nshift\n', expect_ok=False)
    check("Ex_rectification = shift: refused, names Response = shift", rc != 0 and 'response = shift' in log.lower(), f'rc={rc}')

    spec = os.path.join(work, 'check_ex_rectification_spectra.dat')
    np.savetxt(spec, np.c_[w, ex05.real[:, 0, 0, 0], sym(sp05.real)[:, 0, 0, 0], ex05.imag[:, 1, 1, 2], asym(sp05.imag)[:, 1, 1, 2]],
               header='kind: ex_rectification\ncolumns: E_eV ex_re_xxx sp_re_xxx ex_im_yyz sp_im_yyz', comments='# ')
    plot_test_outputs.plot(spec)
    if nfail:
        print(f'\nFAILED: {nfail} check(s)'); sys.exit(1)
    print('\nALL CHECKS PASSED')


if __name__ == '__main__':
    main()
