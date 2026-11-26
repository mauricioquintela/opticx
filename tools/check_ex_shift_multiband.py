#!/usr/bin/env python3
"""Excitonic shift current and rectification without interaction: more than one valence or conduction band, and a band
window that leaves out bands of the model.

TERM 3. The shift-current convention of npj Comput. Mater. 11, 13 (SI Eqs. 9-11, Response = shift) cancels term 3 of
Taghizadeh & Pedersen PRB 97, 205432 (2018), Eq. (B1a). Without interaction and with one valence and one conduction
band that term is absent from the shift current; with more bands it couples the non-interacting excitons (v,c,k) and
(v',c,k), (v,c',k) through v_vv', v_cc'. By partial fractions (tools/sumrule_shift.py, coherence_gen_der) the
same-type three-band pieces of the single-particle shift current split into a coherence part, which is term 3, and a
transition part, which terms 1-2 keep. So, on the same mesh and with a window that holds the whole model:
    excitonic shift route (B)   = single-particle shift (A) - coherence parts
    excitonic rectification (C) = single-particle rectification (D)

WINDOW. Restricted to a band window, the position components do not commute, and the excitonic routes lose
    comm_out^{ac}_{vc} = i sum_{l outside the window} (r^a_vl r^c_lc - r^c_vl r^a_lc),
the remote part of (r^c_vc);a - (r^a_vc);c = i sum_l (...), entering in the form of Eq. 9 (r^b_cv comm_out^{ac}_vc).
Window_correction (default true) adds the bands outside the window as non-interacting virtual states: the
rectification gets the out-of-window pieces of Eq. (B1a) (it is then D), the shift route the commutator itself (it is
then A minus the coherence parts inside the window only, i.e. A for a 1 x 1 window).

Models: the 2D SSH model of arXiv:2611.05990 (four bands, time-reversal symmetric, gap 1.09 eV, same-type bands
>= 2.4 eV apart) and the same model with the orbital centres buckled by +-0.5 Angstrom in z (H unchanged, out-of-plane
components); 21x21, synthetic non-interacting Xatu text files on opticx's mesh. Checks:
 1. SSH 2v x 2c (the whole model): A = the velocity sum rule; B = A - coherence parts, while A and B differ by O(1);
    C = D; Response = shift WARNs (term 3), the rectification does not; window correction none (log).
 2. Buckled SSH 1v x 1c (bands 2|3 of 4), Window_correction = false (the cause): B = A - comm_out and C = D - comm_out
    on every component, out-of-plane ones included, while comm_out is O(1) of A.
 3. Same, default (the fix): B = A; C = D in Re (symmetric) and Im (antisymmetric, the injection current); the log
    says the correction was applied; no term-3 WARNING.
    A 2D map with w_q of both signs (method A at every point): = the single-particle map, in both half-planes, and its
    anti-diagonal = Response = rectification on the same grid.
 4. The window-edge WARNING: silent on the SSH models (bands >= 2.4 eV apart); MoS2 (repo model, needs --root) with the
    window 0 1, which cuts the spin-degenerate pairs, warns that an outside band lies within 30 meV of a window band.
 5. The 2 eta report: logged on every run with outside bands (buckled SSH: 0% of the k-points, no WARNING); with an
    eta of 1.5 eV (2 eta beyond the 2.4 eV distance to the outside bands) it exceeds a quarter of the mesh and warns.

Usage: python3 tools/check_ex_shift_multiband.py --opticx bin/opticx --root . --workdir bin/check_ex_shift_multiband
"""
import argparse, os, subprocess, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sumrule_shift as S
from ipa_shift_numpy import ANG, HA
from check_realtime_sign import mesh, C_AU
import plot_test_outputs

NC = 21
EV = (1.0, 9.0, 0.1, 400)                     # opticx Energy_variables: E1 E2 eta nw (eV)
NF = 2                                        # filled bands of the SSH model
SP = """# Periodic dimensions
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
{resp}
# Energy_variables
{ev}
"""
EX = """# Periodic dimensions
2
# Wannier90_filename
{tb}
# Xatu_interface
true
{eig}
{sta}
# Nfermi
{nf}
# OME_sp
nonlinear
# OME_ex
nonlinear
# Response
{resp}
# Energy_variables
{ev}
"""
WARN = 'WARNING (shift)'
WC_ON = 'window correction applied'
EDGE = 'lies within 30 meV of a window band'
ETA_LINE = 'an outside band lies within 2 eta of a window band at'
ETA_WARN = 'within 2 eta of the window at more than a quarter'


def write_ni(M, val, con, prefix):
    """Non-interacting Xatu text files on opticx's mesh: every transition (v, c, k) of the window is an exciton with a
    unit envelope, energy eps_c - eps_v; basis valence fastest, then conduction, then k (Xatu's order)."""
    K = mesh(M['lat'], NC); nv, nc = len(val), len(con); dim = len(K)*nv*nc
    E = np.empty(dim); rows = []
    for ik, k in enumerate(K):
        e = np.linalg.eigvalsh(S.kpoint_ops(M, k)[0])*HA
        for ic, c in enumerate(con):
            for iv, v in enumerate(val):
                E[ik*nv*nc + ic*nv + iv] = e[c] - e[v]
                rows.append(f'{k[0]/ANG:.12f} {k[1]/ANG:.12f} 0.000000000000 {v} {c}')
    order = np.argsort(E, kind='stable')
    with open(prefix + '.eigval', 'w') as f:
        f.write(f'{NC}\n{dim}\n{dim}\n' + ''.join(f'{E[i]:.12f}\n' for i in order))
    zero = ['0 0']*dim
    with open(prefix + '.states', 'w') as f:
        f.write(f'{dim}\n' + '\n'.join(rows) + '\n')
        for i in order:
            zero[i] = '1 0'; f.write(' '.join(zero) + '\n'); zero[i] = '0 0'
    return prefix + '.eigval', prefix + '.states'


def window_comm(M, val, con, w, eta):
    """comm_out in the form of Eq. 9 (resonant delta, a.u.) for the window transitions (v in val, c in con): the
    commutator through the bands outside the window."""
    K = mesh(M['lat'], NC); out = np.zeros((3, 3, 3, len(w)), complex)
    for k in K:
        E, v, r, _ = S.gen_der(*S.kpoint_ops(M, k), 1e-6)
        for iv in val:
            for ic in con:
                C = np.array([[1j*sum(r[a, iv, l]*r[c, l, ic] - r[c, iv, l]*r[a, l, ic] for l in range(len(E))
                                      if l not in val + con) for c in range(3)] for a in range(3)])
                I = np.einsum('b,ac->abc', r[:, ic, iv], C); I = I + I.transpose(0, 2, 1)
                out += I[..., None]*(1j*(1/(E[ic] - E[iv] - w - 1j*eta)).imag)[None, None, None, :]
    return (out/(2*len(K)*M['area'])).real


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--opticx', required=True); ap.add_argument('--workdir', default='check_ex_shift_multiband')
    ap.add_argument('--root', default='.')
    a = ap.parse_args()
    work, opticx, root = os.path.abspath(a.workdir), os.path.abspath(a.opticx), os.path.abspath(a.root)
    os.makedirs(work, exist_ok=True)
    nfail = 0

    def check(label, ok, msg=''):
        nonlocal nfail
        print(('  PASS   ' if ok else '  FAIL   ') + label + ('   ' + msg if msg else '')); nfail += 0 if ok else 1

    def run(tag, text):
        d = os.path.join(work, tag); os.makedirs(d, exist_ok=True)
        open(os.path.join(d, f'{tag}.txt'), 'w').write(text)
        r = subprocess.run([opticx, f'{tag}.txt'], cwd=d, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log = r.stdout.decode('utf-8', 'replace'); open(os.path.join(d, 'run.log'), 'w').write(log)
        if r.returncode != 0:
            print(f'FAIL: opticx run {tag} exited with {r.returncode} (see {d}/run.log)'); sys.exit(1)
        return d, log

    def shift(d, ex, mat):
        A = np.loadtxt(os.path.join(d, f'shift_{"ex" if ex else "sp"}_lengthgauge_{mat}.dat'), comments='#')
        return A[:, 1:28].reshape(-1, 3, 3, 3).transpose(1, 2, 3, 0)/C_AU

    def rect(d, ex, mat):
        """b<->c symmetric Re and antisymmetric Im, a.u."""
        A = np.loadtxt(os.path.join(d, f'second_{"ex_" if ex else ""}rectification_lengthgauge_{mat}.dat'), comments='#')
        T = (A[:, 2::2] + 1j*A[:, 3::2]).reshape(-1, 3, 3, 3).transpose(1, 2, 3, 0)/C_AU
        return 0.5*(T + T.transpose(0, 2, 1, 3)).real, 0.5*(T - T.transpose(0, 2, 1, 3)).imag

    e1, e2, eta_ev, nw = EV; ev = f'{e1} {e2} {eta_ev} {nw}'
    w = (e1 + (e2 - e1)/nw*np.arange(nw))/HA; eta = eta_ev/HA
    nrm = np.linalg.norm
    rel = lambda X, Y: nrm(X - Y)/nrm(Y)
    print(f'SSH model of arXiv:2611.05990, {NC}x{NC}, eta {eta_ev} eV, non-interacting excitons\n')

    # 1. SSH, 2v x 2c: the whole model
    tb = os.path.join(work, 'ssh_domain1_tb.dat'); S.write_ssh_tb(tb); M = S.model(tb); mat = 'ssh_domain1'
    val, con = [0, 1], [2, 3]
    dA, _ = run('sp_shift_2v2c', SP.format(tb=tb, nc=NC, bl='-1 0 1 2', nf=NF, resp='shift', ev=ev))
    dD, _ = run('sp_rect_2v2c', SP.format(tb=tb, nc=NC, bl='-1 0 1 2', nf=NF, resp='rectification', ev=ev))
    eig, sta = write_ni(M, val, con, os.path.join(work, 'ssh_ni_2v2c'))
    dB, logB = run('ex_shift_2v2c', EX.format(tb=tb, eig=eig, sta=sta, nf=NF, resp='shift', ev=ev))
    dC, logC = run('ex_rect_2v2c', EX.format(tb=tb, eig=eig, sta=sta, nf=NF, resp='rectification', ev=ev))
    A, B, (C, _), (D, _) = shift(dA, False, mat), shift(dB, True, mat), rect(dC, True, mat), rect(dD, False, mat)
    full = S.spectra(M, NC, NF, val + con, w, eta)[0].real
    coh = S.spectra(M, NC, NF, val + con, w, eta, coherence=True)[0].real

    print('SSH, 2 valence x 2 conduction bands (the whole model)')
    e = rel(A, full); check('single-particle shift A = velocity sum rule', e < 1e-6, f'{e:.1e} (tol 1e-6)')
    g = rel(A, B); check('A and the excitonic shift route B differ by O(1)', g > 0.5, f'||A - B||/||B|| = {g:.3f} (> 0.5)')
    e = rel(B, A - coh); check('B = A - coherence parts', e < 1e-2, f'{e:.1e} (tol 1e-2)')
    e = rel(C, D); check('excitonic rectification C = single-particle rectification D', e < 1e-2, f'{e:.1e} (tol 1e-2)')
    check('Response = shift warns (term 3, 2v x 2c)', WARN in logB)
    check('Response = rectification does not warn', WARN not in logC)
    check('window correction: none (the window is the whole model)', 'Window correction: none' in logB and WC_ON not in logB)

    # 2.-3. buckled SSH, 1v x 1c: the window leaves out bands 1 and 4
    tbb = os.path.join(work, 'sshbuck_tb.dat'); S.write_ssh_tb(tbb, buckle=0.5); Mb = S.model(tbb); matb = 'sshbuck'
    val1, con1 = [1], [2]
    dA1, _ = run('sp_shift_buck_1v1c', SP.format(tb=tbb, nc=NC, bl='0 1', nf=NF, resp='shift', ev=ev))
    dD1, _ = run('sp_rect_buck_1v1c', SP.format(tb=tbb, nc=NC, bl='0 1', nf=NF, resp='rectification', ev=ev))
    eig1, sta1 = write_ni(Mb, val1, con1, os.path.join(work, 'sshbuck_ni_1v1c'))
    off = '# Window_correction\nfalse\n'
    dB0, _ = run('ex_shift_buck_1v1c_off', EX.format(tb=tbb, eig=eig1, sta=sta1, nf=NF, resp='shift', ev=ev) + off)
    dC0, _ = run('ex_rect_buck_1v1c_off', EX.format(tb=tbb, eig=eig1, sta=sta1, nf=NF, resp='rectification', ev=ev) + off)
    dB1, logB1 = run('ex_shift_buck_1v1c', EX.format(tb=tbb, eig=eig1, sta=sta1, nf=NF, resp='shift', ev=ev))
    dC1, logC1 = run('ex_rect_buck_1v1c', EX.format(tb=tbb, eig=eig1, sta=sta1, nf=NF, resp='rectification', ev=ev))
    # 2D map with w_q of both signs: it contains w_1 + w_2 = 0 points, so method A evaluates all of it
    evm, m2 = f'1 6 {eta_ev} 25', '# Energy_variables_2\n-5 5 50\n'
    dMs, _ = run('sp_map_buck_1v1c', SP.format(tb=tbb, nc=NC, bl='0 1', nf=NF, resp='general', ev=evm) + m2)
    dMe, logM = run('ex_map_buck_1v1c', EX.format(tb=tbb, eig=eig1, sta=sta1, nf=NF, resp='general', ev=evm) + m2)
    dMo, _ = run('ex_map_buck_1v1c_off', EX.format(tb=tbb, eig=eig1, sta=sta1, nf=NF, resp='general', ev=evm) + m2 + off)
    dMr, _ = run('ex_rect_buck_1v1c_mapgrid', EX.format(tb=tbb, eig=eig1, sta=sta1, nf=NF, resp='rectification', ev=evm))
    for f in (sta, sta1):
        os.remove(f)
    A1, D1, D1i = shift(dA1, False, matb), *rect(dD1, False, matb)
    B0, (C0, _) = shift(dB0, True, matb), rect(dC0, True, matb)
    B1, (C1, C1i) = shift(dB1, True, matb), rect(dC1, True, matb)
    comm = window_comm(Mb, val1, con1, w, eta)

    print('buckled SSH, 1 valence x 1 conduction band (window omits bands 1 and 4)')
    print('  Window_correction = false (the cause):')
    g = nrm(comm)/nrm(A1); check('comm_out is O(1) of the single-particle shift', g > 0.2, f'{g:.3f} (> 0.2)')
    e = rel(B0, A1 - comm); check('excitonic shift = single-particle shift - comm_out', e < 1e-2, f'{e:.1e} (tol 1e-2)')
    for c in ('xyy', 'zxx', 'xxz'):
        i = tuple('xyz'.index(x) for x in c)
        e = rel(B0[i], (A1 - comm)[i]); check(f'    {c}', e < 1e-2, f'{e:.1e} (tol 1e-2; uncorrected vs A: {rel(B0[i], A1[i]):.2f})')
    e = rel(C0, D1 - comm); check('excitonic rectification = single-particle rectification - comm_out', e < 1e-2, f'{e:.1e} (tol 1e-2)')
    print('  default (the fix):')
    e = rel(B1, A1); check('shift = single-particle shift (Eq. 9; no coherence parts in a 1 x 1 window)', e < 1e-2,
                           f'{e:.1e} (tol 1e-2; uncorrected {rel(B0, A1):.2f})')
    e = rel(C1, D1); check('rectification, Re symmetric = single-particle', e < 5e-3, f'{e:.1e} (tol 5e-3; uncorrected {rel(C0, D1):.2f})')
    e = rel(C1i, D1i); check('rectification, Im antisymmetric (injection) = single-particle', e < 1e-2, f'{e:.1e} (tol 1e-2)')
    check('log: window correction applied (shift and rectification)', WC_ON in logB1 and WC_ON in logC1)
    check('Response = shift does not warn (1v x 1c)', WARN not in logB1)

    print('  2D map (w_p 1-6 eV, w_q -5..5 eV, method A):')
    def gen(f):
        a = np.loadtxt(f, comments='#'); return a[:, 0], a[:, 1], (a[:, 2::2] + 1j*a[:, 3::2]).reshape(-1, 3, 3, 3)
    wpm, wqm, Ms = gen(os.path.join(dMs, f'second_general_lengthgauge_{matb}.dat'))
    _, _, Me = gen(os.path.join(dMe, f'second_ex_general_lengthgauge_{matb}.dat'))
    _, _, Mo = gen(os.path.join(dMo, f'second_ex_general_lengthgauge_{matb}.dat'))
    wr, _, Rr = gen(os.path.join(dMr, f'second_ex_rectification_lengthgauge_{matb}.dat'))
    pt = lambda Z: np.sqrt((np.abs(Z)**2).sum(axis=(1, 2, 3)))
    e, e0 = pt(Me - Ms).sum()/pt(Ms).sum(), pt(Mo - Ms).sum()/pt(Ms).sum()
    check('map = single-particle map (sum of per-point norms)', e < 1e-2, f'{e:.1e} (tol 1e-2; uncorrected {e0:.2f})')
    for lab, sel in (('w_q < 0', wqm < -1e-9), ('w_q > 0', wqm > 1e-9)):
        e = pt(Me[sel] - Ms[sel]).sum()/pt(Ms[sel]).sum(); check(f'    {lab}', e < 1e-2, f'{e:.1e} (tol 1e-2)')
    dc = np.abs(wpm + wqm) < 1e-8
    idx = [int(np.argmin(np.abs(wr - x))) for x in wpm[dc]]
    e = np.abs(Me[dc] - Rr[idx]).max()/np.abs(Rr).max()
    check('anti-diagonal = Response = rectification', dc.sum() > 10 and np.abs(wr[idx] - wpm[dc]).max() < 1e-9 and e < 1e-12,
          f'{dc.sum()} points, {e:.1e} (tol 1e-12)')
    check('log: window correction applied (map)', WC_ON in logM)

    print('window-edge warning')
    check('silent when the outside bands are far (buckled SSH)', EDGE not in logB1 and EDGE not in logC1)
    m = [l for l in logB1.splitlines() if ETA_LINE in l]
    check('2 eta report logged, 0% of the k-points, no WARNING (buckled SSH, eta 0.1)',
          len(m) == 1 and ' 0.0%' in m[0] and ETA_WARN not in logB1, m[0].strip() if m else 'no line')
    _, logw = run('sp_shift_buck_eta1p5', SP.format(tb=tbb, nc=NC, bl='0 1', nf=NF, resp='shift', ev=f'{e1} {e2} 1.5 {nw}'))
    check('eta 1.5 eV (2 eta beyond the outside bands): warns', ETA_WARN in logw)
    mos2 = os.path.join(root, 'MoS2_spin_wannier_07032024_tb.dat')
    if not os.path.exists(mos2):
        print(f'  SKIP   MoS2 window-edge check: {mos2} not found')
    else:
        _, logm = run('sp_shift_mos2_edge', SP.format(tb=mos2, nc=9, bl='0 1', nf=26, resp='shift', ev=ev))
        check('MoS2, window 0 1 (cuts the spin-degenerate pairs): warns', EDGE in logm and 'Widen Bandlist' in logm)

    out = os.path.join(work, 'ex_shift_multiband_spectra.dat')
    cols = {'sp_shift': A, 'ex_shift': B, 'sp_shift_minus_coh': A - coh, 'ex_rect': C, 'sp_rect': D}
    bcols = {'b_sp_shift': A1, 'b_ex_shift_off': B0, 'b_ex_shift': B1, 'b_sp_shift_minus_comm': A1 - comm,
             'b_sp_rect': D1, 'b_ex_rect_off': C0, 'b_ex_rect': C1}
    names, data = ['E_eV'], [w*HA]
    for c in ('xxx', 'xyy', 'zxx'):
        i = tuple('xyz'.index(x) for x in c)
        for n_, T in list(cols.items()) + list(bcols.items()):
            names.append(f'{c}_{n_}'); data.append(T[i]*C_AU)
    np.savetxt(out, np.array(data).T, fmt='%.8e', comments='',
               header='# kind: ex_shift_multiband\n# columns: ' + ' '.join(names))
    plot_test_outputs.plot(out)
    print(f'\n{"ALL PASS" if nfail == 0 else f"{nfail} FAIL"}')
    sys.exit(1 if nfail else 0)


if __name__ == '__main__':
    main()
