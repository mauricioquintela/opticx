#!/usr/bin/env python3
"""The band states opticx builds are the ones the exciton envelopes refer to, whatever phase rule Xatu used.

An exciton envelope A_vck is a set of coefficients on Xatu's band states, so opticx must rebuild exactly those states.
Xatu builds differ in how they fix a band's phase (ExcitonTB::fixGlobalPhase): the sum of the coefficients made real, or
the weighted sum sum_i w_i c_i with w_i = exp(2 pi i g (i-1)), g = (sqrt(5) - 1)/2. With a mismatched rule every
excitonic matrix element is scrambled (In2Se3 102x102: excitonic absorbance xx - yy 34% where C3 requires 0, 1.3% as in
Xatu once matched). opticx now takes the window states from the archive when Xatu stored them (/bands/eigenvectors), and
otherwise applies the rule named by Band_phase.

Model: hBN_tb.dat on hBN_N30's mesh, with INTERACTING excitons built here from a model BSE whose direct term
-V(k-k') <c_k|c_k'><v_k'|v_k> is written on the band states, so the same physics can be generated with either rule.
Checks (excitonic SHG, the response that uses every kind of excitonic element):
  * archives that store the band states: the weighted-rule and the plain-sum-rule excitons give the same SHG (the
    physics is the same); the log names the rule found and the largest angle to Xatu's states (rounding);
  * archives without them: Band_phase = weighted / sum on the matching excitons = the same SHG; the wrong rule
    does not (sensitivity: the defect this fixes);
  * the band states of a different model stop the run; an archive in the atomic Bloch gauge stops the run;
  * the .omesp records the band states it holds (tail bit 6 / '#BAND_PHASE 3' in the linear file): OME_sp = none with
    it gives the same SHG, with a file for other states stops before computing;
  * the excitonic cache records them (format 5) and a run with other states refuses it; plain-sum runs keep writing
    the old layout (format 3);
  * a non-orthonormal model: Xatu stores the Lowdin coefficients S^{1/2} c; aligned = Band_phase = weighted on text.

Usage: python3 tools/check_band_phase.py --opticx bin/opticx --root . --workdir bin/check_band_phase
"""
import argparse, glob, os, shutil, struct, subprocess, sys
import numpy as np

G = 0.5*(np.sqrt(5.0) - 1.0)
OMESP_TAG = 1330464562


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--opticx', required=True); ap.add_argument('--root', default='.')
    ap.add_argument('--workdir', default='check_band_phase')
    a = ap.parse_args()
    root, work, opticx = os.path.abspath(a.root), os.path.abspath(a.workdir), os.path.abspath(a.opticx)
    sys.path.insert(0, os.path.join(root, 'tools'))
    try:
        import h5py
    except ImportError:
        print('SKIP: h5py not available'); sys.exit(0)
    from check_tb_hermiticity import read_tb, write_tb
    tb0, states0 = os.path.join(root, 'hBN_tb.dat'), os.path.join(root, 'hBN_N30.states')
    for f in (tb0, states0):
        if not os.path.exists(f):
            print('SKIP: missing ' + f); sys.exit(0)
    shutil.rmtree(work, ignore_errors=True); os.makedirs(work)
    nfail = 0

    def check(label, ok, msg=''):
        nonlocal nfail
        print(('  PASS   ' if ok else '  FAIL   ') + label + ('   ' + msg if msg else '')); nfail += 0 if ok else 1

    m = read_tb(tb0); Rs, norb, lat = m['Rs'], m['norb'], m['latv']
    with open(states0) as f:
        nk = int(f.readline().split()[0]); kl = [f.readline().split() for _ in range(nk)]
    K = np.array([[float(x) for x in l[:3]] for l in kl]); tail = [' '.join(l[3:]) for l in kl]   # Angstrom^-1
    Rc = np.array([R[0]*lat[0] + R[1]*lat[1] + R[2]*lat[2] for R in Rs])
    deg = np.array(m['deg'], float)
    Bvec = 2*np.pi*np.linalg.inv(lat).T                                   # reciprocal vectors (rows)

    # non-orthonormal twin: same H(R), an overlap S(R) and positions consistent with it (as check_nonortho_phase)
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
    # a different model (stronger hoppings) for the mismatch check
    Hx = {R: (m['H'][R]*1.3 if R != (0, 0, 0) else m['H'][R]) for R in Rs}
    tbx = os.path.join(work, 'hBN_other_tb.dat'); write_tb(tbx, m, Hx, m['RP'])

    def bands(rule, overlap):
        """Energies (nk, 2) and the coefficient vectors Xatu would store (nk, 2, norb), phases by `rule`."""
        E = np.zeros((nk, norb)); V = np.zeros((nk, norb, norb), complex)
        w = np.exp(2j*np.pi*G*np.arange(norb))
        for j, k in enumerate(K):
            ph = np.exp(1j*(Rc @ k))/deg
            Hk = np.einsum('r,rij->ij', ph, np.array([m['H'][R] for R in Rs]))
            if overlap:
                Sk = np.einsum('r,rij->ij', ph, np.array([S[R] for R in Rs]))
                s, U = np.linalg.eigh(Sk); Sm = (U/np.sqrt(s)) @ U.conj().T
                Hk = Sm @ Hk @ Sm                                         # Lowdin Hamiltonian, as Xatu
            e, X = np.linalg.eigh(0.5*(Hk + Hk.conj().T))
            t = X.sum(0) if rule == 'sum' else w @ X
            E[j], V[j] = e, (X*np.exp(-1j*np.angle(t))).T
        return E, V

    def excitons(E, V):
        """Model BSE (TDA, direct term only) on these band states; energies (eV) ascending, envelopes (nk, nexc)."""
        dk = K[:, None, :] - K[None, :, :]
        frac = dk @ np.linalg.inv(Bvec); frac -= np.round(frac)
        q2 = ((frac @ Bvec)**2).sum(-1)
        Vq = 2.0*np.exp(-q2/0.5**2)*(30.0/np.sqrt(nk))**2/30.0
        oc = V[:, 1, :].conj() @ V[:, 1, :].T                            # <c_k|c_k'>
        ov = V[:, 0, :].conj() @ V[:, 0, :].T                            # <v_k|v_k'>
        Hb = np.diag(E[:, 1] - E[:, 0]) - Vq*oc*ov.T
        ee, A = np.linalg.eigh(0.5*(Hb + Hb.conj().T))
        return ee, A

    def write_text(prefix, ee, A):
        with open(prefix + '.eigval', 'w') as f:
            f.write(f'30\n{nk}\n{nk}\n' + ''.join(f'{x:.10f}\n' for x in ee))
        with open(prefix + '.states', 'w') as f:
            f.write(f'{nk}\n' + ''.join(f'{K[j,0]:12.7f} {K[j,1]:12.7f} {K[j,2]:12.7f} {tail[j]}\n' for j in range(nk)))
            for n in range(A.shape[1]):
                f.write(' '.join(f'{x.real:.15e} {x.imag:.15e}' for x in A[:, n]) + '\n')

    def write_archive(prefix, E, V, with_bands=True, gauge='lattice'):
        h5 = prefix + ('.h5' if with_bands else '_nobands.h5')
        subprocess.run([sys.executable, os.path.join(root, 'tools', 'xatu_text_to_h5.py'), prefix + '.eigval', prefix + '.states',
                        '--nfermi', '1', '--ndim', '2', '-o', h5], check=True, stdout=subprocess.DEVNULL)
        with h5py.File(h5, 'a') as f:
            if gauge is not None: f.attrs['gauge'] = gauge
            if with_bands:
                g = f.create_group('bands')
                g.create_dataset('indices', data=np.array([0, 1], np.int64))
                g.create_dataset('energies', data=E)
                g.create_dataset('eigenvectors', data=V)
        return h5

    sets = {}
    for rule in ('weighted', 'sum'):
        for ov in (False, True):
            E, V = bands(rule, ov); ee, A = excitons(E, V)
            pre = os.path.join(work, f'hBN_{rule}' + ('_nonortho' if ov else ''))
            write_text(pre, ee, A)
            sets[(rule, ov)] = dict(pre=pre, h5=write_archive(pre, E, V), h5n=write_archive(pre, E, V, with_bands=False))
    pre_at = sets[('weighted', False)]['pre']
    shutil.copy(pre_at + '.eigval', pre_at + 'A.eigval'); shutil.copy(pre_at + '.states', pre_at + 'A.states')
    atomic = write_archive(pre_at + 'A', *bands('weighted', False), gauge='atomic')

    def run(tag, src, tb=tb0, ortho=True, omesp='nonlinear', resp='shg', extra='', d=None):
        d = d or os.path.join(work, tag); os.makedirs(d, exist_ok=True)
        xin = src if isinstance(src, str) else f'{src[0]}\n{src[1]}'
        txt = (f'# Periodic dimensions\n2\n# Wannier90_filename\n{tb}\n# Xatu_interface\ntrue\n{xin}\n# Exciton_cutoff\n300\n'
               f'# Nfermi\n1\n# OME_sp\n{omesp}\n# OME_ex\n{"linear" if resp == "absorbance" else "nonlinear"}\n# Response\n{resp}\n'
               f'# Energy_variables\n2 8 0.2 61\n' + ('' if ortho else '# Orthonormal\nfalse\n') + extra)
        open(os.path.join(d, f'{tag}.txt'), 'w').write(txt)
        r = subprocess.run([opticx, f'{tag}.txt'], cwd=d, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log = r.stdout.decode('utf-8', 'replace'); open(os.path.join(d, f'{tag}.log'), 'w').write(log)
        return d, r.returncode, log

    def shg(d):
        x = np.loadtxt(glob.glob(os.path.join(d, 'shg_ex_lengthgauge_*.dat'))[0]); return x[:, 1::2] + 1j*x[:, 2::2]

    def rel(x, y):
        return np.abs(x - y).max()/np.abs(y).max()

    ended = lambda rc, log: rc == 0 and 'Opticx calculation ended' in log

    # 1. archives that store the band states
    dW, rc, log = run('aligned_weighted', sets[('weighted', False)]['h5'])
    ok = ended(rc, log)
    check('archive with band states (weighted rule): runs, states taken from it, rule recognised', ok and
          'taken from the Xatu archive (phases: weighted sum' in log, f'rc {rc}')
    dev = [l for l in log.splitlines() if 'largest angle' in l]
    check('... largest angle to Xatu''s states < 1e-6 rad (rounding)', bool(dev) and
          float(dev[0].split('angle')[1].split('rad')[0]) < 1e-6, dev[0].strip() if dev else 'no angle line')
    ref = shg(dW) if ok else None
    dS, rc, log = run('aligned_sum', sets[('sum', False)]['h5'])
    check('archive with band states (plain-sum rule) = weighted rule: same excitonic SHG (1e-8)', ended(rc, log) and ref is not None
          and rel(shg(dS), ref) < 1e-8 and 'phases: sum of the coefficients real' in log,
          f'{rel(shg(dS), ref):.1e}' if ended(rc, log) and ref is not None else f'rc {rc}')

    # 2. archives without them: Band_phase decides
    for rule, kw, same in (('weighted', 'weighted', True), ('sum', 'sum', True), ('weighted', 'sum', False)):
        d, rc, log = run(f'nobands_{rule}_bp{kw}', sets[(rule, False)]['h5n'], extra=f'# Band_phase\n{kw}\n')
        r_ = rel(shg(d), ref) if ended(rc, log) and ref is not None else np.inf
        if same:
            check(f'no band states, {rule}-rule excitons, Band_phase = {kw}: = aligned run (1e-8)', r_ < 1e-8, f'{r_:.1e}')
        else:
            check(f'sensitivity: {rule}-rule excitons read with Band_phase = {kw} differ from it (> 5%)',
                  ended(rc, log) and r_ > 5e-2, f'{r_:.2f}')
    d, rc, log = run('text_weighted', (sets[('weighted', False)]['pre'] + '.eigval', sets[('weighted', False)]['pre'] + '.states'),
                     extra='# Band_phase\nweighted\n')
    check('text files, Band_phase = weighted: = aligned run, log states the rule', ended(rc, log) and rel(shg(d), ref) < 1e-6 and
          'Band phases of the exciton basis: weighted sum real (Band_phase = weighted' in log,
          f'{rel(shg(d), ref):.1e}' if ended(rc, log) else f'rc {rc}')
    d, rc, log = run('aligned_bp_given', sets[('weighted', False)]['h5'], extra='# Band_phase\nsum\n')
    check('Band_phase given with an archive that stores the states: ignored with a NOTE', ended(rc, log) and
          'Band_phase = sum ignored' in log and rel(shg(d), ref) < 1e-12)

    # 3. guards
    d, rc, log = run('other_model', sets[('weighted', False)]['h5'], tb=tbx)
    check('band states of a different model: stops, nothing computed', rc != 0 and 'are not the eigenstates of this model' in log
          and 'Entering sigma' not in log, f'rc {rc}')
    d, rc, log = run('atomic_gauge', atomic)
    check('archive in the atomic Bloch gauge: stops before computing', rc != 0 and 'Bloch gauge' in log and 'Entering ome' not in log,
          f'rc {rc}')

    # 4. the .omesp record
    omesp = os.path.join(dW, 'ome_nonlinear_sp_hBN.omesp')
    b = open(omesp, 'rb').read(); nb = 2; nkk = struct.unpack_from('<i', b, 4)[0]
    off = 4 + 8 + 3*nkk*8 + nkk*(nb*8 + 2*3*nb*nb*16 + 9*nb*nb*8 + 9*nb*nb*16) + nkk*9*nb*nb*8 + nkk*9*nb*nb*16
    tag, flags = struct.unpack_from('<2i', b, off)
    check('.omesp records Xatu''s states (tail bit 6)', tag == OMESP_TAG and bool(flags & 64) and not flags & 32, f'flags {flags}')
    d, rc, log = run('aligned_weighted', sets[('weighted', False)]['h5'], omesp='none', d=dW)
    check('OME_sp = none with that .omesp: runs, = the computed run (1e-12)', ended(rc, log) and rel(shg(dW), ref) < 1e-12,
          f'rc {rc}')
    dO = os.path.join(work, 'omesp_other'); os.makedirs(dO)
    shutil.copy(omesp, dO)
    d, rc, log = run('omesp_other', sets[('weighted', False)]['h5n'], omesp='none', extra='# Band_phase\nweighted\n', d=dO)
    check('OME_sp = none, .omesp of Xatu''s states but a run with Band_phase = weighted: stops before computing', rc != 0 and
          'band states in the .omesp were made with the band states of the Xatu archive' in log and 'Entering sigma' not in log,
          f'rc {rc}')
    bb = bytearray(b); struct.pack_into('<i', bb, off + 4, flags & ~64)
    open(os.path.join(dO, 'ome_nonlinear_sp_hBN.omesp'), 'wb').write(bytes(bb))
    d, rc, log = run('omesp_other', sets[('weighted', False)]['h5'], omesp='none', d=dO)
    check('OME_sp = none, an older .omesp (plain-sum states) with an archive that stores the states: stops', rc != 0 and
          'were made with the plain-sum phase rule' in log and 'Entering sigma' not in log, f'rc {rc}')
    d, rc, log = run('linear_aligned', sets[('weighted', False)]['h5'], resp='absorbance', omesp='linear')
    lin = open(os.path.join(d, 'ome_linear_sp_hBN.omesp')).read() if ended(rc, log) else ''
    check('linear .omesp records it (#BAND_PHASE 3)', '#BAND_PHASE 3' in lin, f'rc {rc}')
    d, rc, log = run('linear_aligned', sets[('weighted', False)]['h5'], resp='absorbance', omesp='none', d=d)
    check('... and OME_sp = none reads it back', ended(rc, log), f'rc {rc}')

    # 5. the excitonic cache
    dC = os.path.join(work, 'cache'); os.makedirs(dC)
    d, rc, log = run('cache', sets[('weighted', False)]['h5'], extra='# Cache_ome_ex\nwrite\n', d=dC)
    cf = os.path.join(dC, 'ome_second_ex_hBN.omeex2')
    ver = struct.unpack_from('<i', open(cf, 'rb').read(), 14)[0] if os.path.exists(cf) else -1
    check('cache of a run on Xatu''s states: format 5', ended(rc, log) and ver == 5, f'version {ver}')
    d, rc, log = run('cache', sets[('weighted', False)]['h5n'], extra='# Cache_ome_ex\nread\n# Band_phase\nweighted\n', d=dC)
    check('... read by a run with other states (Band_phase = weighted): stops', rc != 0 and 'DIFFERENT band states' in log, f'rc {rc}')
    d, rc, log = run('cache', sets[('weighted', False)]['h5'], extra='# Cache_ome_ex\nread\n', d=dC)
    check('... read by a run on the same states: hit, = computed (1e-12)', ended(rc, log) and 'read from cache' in log.lower()
          and rel(shg(dC), ref) < 1e-12, f'rc {rc}')
    dC3 = os.path.join(work, 'cache_sum'); os.makedirs(dC3)
    d, rc, log = run('cache_sum', sets[('sum', False)]['h5n'], extra='# Cache_ome_ex\nwrite\n', d=dC3)
    cf = os.path.join(dC3, 'ome_second_ex_hBN.omeex2')
    ver = struct.unpack_from('<i', open(cf, 'rb').read(), 14)[0] if os.path.exists(cf) else -1
    check('plain-sum run without band states: cache layout unchanged (format 3)', ended(rc, log) and ver == 3, f'version {ver}')

    # 6. non-orthonormal model: Xatu stores the Lowdin coefficients
    d, rc, log = run('nonortho_aligned', sets[('weighted', True)]['h5'], tb=tbn, ortho=False)
    ok = ended(rc, log); refn = shg(d) if ok else None
    dev = [l for l in log.splitlines() if 'largest angle' in l]
    check('non-orthonormal, archive with (Lowdin) band states: runs, angle < 1e-6 rad', ok and bool(dev) and
          float(dev[0].split('angle')[1].split('rad')[0]) < 1e-6, dev[0].strip() if dev else f'rc {rc}')
    d, rc, log = run('nonortho_text', (sets[('weighted', True)]['pre'] + '.eigval', sets[('weighted', True)]['pre'] + '.states'),
                     tb=tbn, ortho=False, extra='# Band_phase\nweighted\n')
    check('non-orthonormal: aligned = text with Band_phase = weighted (1e-6)', ended(rc, log) and refn is not None and
          rel(shg(d), refn) < 1e-6, f'{rel(shg(d), refn):.1e}' if ended(rc, log) and refn is not None else f'rc {rc}')
    d, rc, log = run('nonortho_sum', sets[('sum', True)]['h5'], tb=tbn, ortho=False)
    check('non-orthonormal: plain-sum archive = weighted archive (1e-8)', ended(rc, log) and refn is not None and
          rel(shg(d), refn) < 1e-8, f'{rel(shg(d), refn):.1e}' if ended(rc, log) and refn is not None else f'rc {rc}')

    if nfail:
        print(f'\nFAILED: {nfail} check(s)'); sys.exit(1)
    print('\nALL CHECKS PASSED')


if __name__ == '__main__':
    main()
