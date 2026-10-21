#!/usr/bin/env python3
"""Band phases of a non-orthonormal model follow Xatu's convention, so the exciton envelopes refer to the same states.

Xatu fixes the phase of each band state by making the sum of its eigenvector's coefficients real (ExcitonTB::fixGlobalPhase).
With an overlap it diagonalises the Lowdin Hamiltonian S^{-1/2} H S^{-1/2}, so the coefficients are S^{1/2} c, not the
c of the generalised problem H c = E S c. opticx now applies the same rule to S^{1/2} c (ome_sp fix_band_phases); before,
it applied it to c, which on MoSe2 SOC (2v x 4c, CRYSTAL model) put the excitonic absorbance 15% below Xatu's (1.3e-3 with
the convention).

Model: hBN_tb.dat rewritten with an overlap (on-site off-diagonal 0.35 plus a small random longer-range part) and positions
obeying r(R) = r(-R)^+ + R S(R); non-interacting Xatu-format excitons on hBN_N30's k mesh. Checks:
  * the run stores S(k) c consistent with S(k) rebuilt here (the k convention of the check is opticx's);
  * every stored window state has sum(S^{1/2} c) real (1e-10), while sum(c) is not (sensitivity: max |arg| > 1e-3; 0.06 rad here; the old rule fails
    the first check by that much);
  * the .omesp records the convention (bit 4 of the tagged tail); the orthonormal hBN run does not set it and keeps the
    old rule (sum(c) real);
  * an excitonic run with OME_sp = none and an .omesp without the record stops before computing; with it, it runs.

Usage: python3 tools/check_nonortho_phase.py --opticx bin/opticx --root . --workdir bin/check_nonortho_phase
"""
import argparse, os, shutil, struct, subprocess, sys
import numpy as np

ANG = 0.52917721067121        # bohr in Angstrom
OMESP_TAG = 1330464562


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--opticx', required=True); ap.add_argument('--root', default='.')
    ap.add_argument('--workdir', default='check_nonortho_phase')
    a = ap.parse_args()
    root, work, opticx = os.path.abspath(a.root), os.path.abspath(a.workdir), os.path.abspath(a.opticx)
    sys.path.insert(0, os.path.join(root, 'tools'))
    from check_tb_hermiticity import read_tb, write_tb
    from ni_excitons import make_ni
    tb0, states = os.path.join(root, 'hBN_tb.dat'), os.path.join(root, 'hBN_N30.states')
    for f in (tb0, states):
        if not os.path.exists(f):
            print('SKIP: missing ' + f); sys.exit(0)
    os.makedirs(work, exist_ok=True)
    nfail = 0

    def check(label, ok, msg=''):
        nonlocal nfail
        print(('  PASS   ' if ok else '  FAIL   ') + label + ('   ' + msg if msg else '')); nfail += 0 if ok else 1

    # non-orthonormal hBN: same H(R), an overlap S(R) and positions consistent with it
    m = read_tb(tb0); Rs, norb, lat = m['Rs'], m['norb'], m['latv']
    mR = lambda R: (-R[0], -R[1], -R[2])
    rng = np.random.default_rng(11); S = {}
    for R in Rs:
        if R == (0, 0, 0):
            S[R] = np.array([[1.0, 0.35], [0.35, 1.0]], complex)
        elif mR(R) in S:
            S[R] = S[mR(R)].conj().T
        else:
            S[R] = 0.02*(rng.standard_normal((norb, norb)) + 1j*rng.standard_normal((norb, norb)))*np.exp(-np.linalg.norm(np.array(R)))
    RP = {R: m['RP'][R] + 0.5*(np.array(R, float) @ lat)[:, None, None]*S[R][None] for R in Rs}
    tbn = os.path.join(work, 'hBN_nonortho_tb.dat'); write_tb(tbn, m, m['H'], RP, S)
    pre = os.path.join(work, 'hBN_NI30'); make_ni(tb0, states, pre)

    def run(tag, tb, ortho, omesp='nonlinear', extra=''):
        d = os.path.join(work, tag); os.makedirs(d, exist_ok=True)
        txt = (f'# Periodic dimensions\n2\n# Wannier90_filename\n{tb}\n# Xatu_interface\ntrue\n{pre}.eigval\n{pre}.states\n'
               f'# Nfermi\n1\n# OME_sp\n{omesp}\n# OME_ex\nnonlinear\n# Response\nshift\n# Energy_variables\n4 12 0.2 81\n'
               + ('' if ortho else '# Orthonormal\nfalse\n') + extra)
        open(os.path.join(d, f'{tag}.txt'), 'w').write(txt)
        r = subprocess.run([opticx, f'{tag}.txt'], cwd=d, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log = r.stdout.decode('utf-8', 'replace'); open(os.path.join(d, 'run.log'), 'w').write(log)
        return d, r.returncode, log

    def tail(path, nb):
        """k list (bohr^-1) and the tagged tail of a nonlinear .omesp: flags, norb, win_evec, win_sevec."""
        b = open(path, 'rb').read()
        nk = struct.unpack_from('<i', b, 4)[0]
        rk = np.frombuffer(b, '<f8', 3*nk, 12).reshape(3, nk).T
        off = 4 + 8 + 3*nk*8 + nk*(nb*8 + 2*3*nb*nb*16 + 9*nb*nb*8 + 9*nb*nb*16) + nk*9*nb*nb*8 + nk*9*nb*nb*16
        off0 = off
        tag, flags, no = struct.unpack_from('<3i', b, off); off += 12
        assert tag == OMESP_TAG, 'no tagged tail'
        if flags & 1: off += 16*nk*nb*nb
        ev = np.frombuffer(b, '<c16', no*nb*nk, off).reshape(nk, nb, no).transpose(2, 1, 0); off += 16*no*nb*nk
        sev = np.frombuffer(b, '<c16', no*nb*nk, off).reshape(nk, nb, no).transpose(2, 1, 0)
        return rk, flags, ev, sev, off0

    d, rc, log = run('nonortho', tbn, False)
    check('non-orthonormal hBN, excitonic shift from Xatu-format excitons: runs', rc == 0 and 'Opticx calculation ended' in log, f'rc {rc}')
    omesp = os.path.join(d, 'ome_nonlinear_sp_hBN_nonortho.omesp')
    rk, flags, ev, sev, flag_off = tail(omesp, 2)
    check('the .omesp records the convention (tail bit 4)', bool(flags & 16), f'flags {flags}')
    Sk = [sum(np.exp(1j*k[:2] @ ((np.array(R, float) @ lat)[:2]/ANG))*S[R] for R in Rs) for k in rk]
    dS = max(abs(Sk[i] @ ev[:, :, i] - sev[:, :, i]).max() for i in range(len(rk)))
    check('S(k) rebuilt here matches the run (S c stored by opticx)', dS < 1e-10, f'{dS:.1e}')
    ph_l, ph_c = 0.0, 0.0
    for i in range(len(rk)):
        w, U = np.linalg.eigh(Sk[i]); sh = (U*np.sqrt(w)) @ U.conj().T
        ph_l = max(ph_l, abs(np.angle((sh @ ev[:, :, i]).sum(0))).max())
        ph_c = max(ph_c, abs(np.angle(ev[:, :, i].sum(0))).max())
    check('every window state has sum(S^1/2 c) real (Xatu convention)', ph_l < 1e-10, f'max |arg| {ph_l:.1e}')
    check('sensitivity: sum(c) is not real (the old convention differs here)', ph_c > 1e-3, f'max |arg| {ph_c:.2f}')

    d, rc, log = run('ortho', tb0, True)
    rk, flags, ev, _, _ = tail(os.path.join(d, 'ome_nonlinear_sp_hBN.omesp'), 2)
    po = max(abs(np.angle(ev[:, :, i].sum(0))).max() for i in range(len(rk)))
    check('orthonormal hBN: convention bit not set, sum(c) real as before', rc == 0 and not flags & 16 and po < 1e-10,
          f'flags {flags}, max |arg| {po:.1e}')

    old = os.path.join(work, 'nonortho_oldomesp'); os.makedirs(old, exist_ok=True)
    b = bytearray(open(omesp, 'rb').read()); f = struct.unpack_from('<i', b, flag_off + 4)[0]
    struct.pack_into('<i', b, flag_off + 4, f & ~16)
    open(os.path.join(old, 'ome_nonlinear_sp_hBN_nonortho.omesp'), 'wb').write(bytes(b))
    d, rc, log = run('nonortho_oldomesp', tbn, False, omesp='none')
    check('OME_sp = none with an .omesp without the record: stops before computing', rc != 0 and
          'before opticx fixed the band' in log and 'Entering sigma' not in log, f'rc {rc}')
    shutil.copy(omesp, os.path.join(old, 'ome_nonlinear_sp_hBN_nonortho.omesp'))
    d, rc, log = run('nonortho_oldomesp', tbn, False, omesp='none')
    check('... with the record: runs', rc == 0 and 'Opticx calculation ended' in log, f'rc {rc}')

    if nfail:
        print(f'\nFAILED: {nfail} check(s)'); sys.exit(1)
    print('\nALL CHECKS PASSED')


if __name__ == '__main__':
    main()
