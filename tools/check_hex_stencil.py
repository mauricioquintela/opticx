#!/usr/bin/env python3
"""The covariant envelope derivative of X_nm is C3-covariant on a hexagonal mesh.

X_nm is built from transported differences of the exciton envelopes between mesh neighbours. Along the two crystal axes
alone that difference operator is not mapped onto itself by C3, and for envelopes the mesh resolves poorly it makes X_nm
anisotropic (MoSe2 with spin-orbit coupling, 87x87, 2078 excitons: excitonic SHG C3 residual 38.6%). On a 2D hexagonal
mesh opticx uses the three star directions with least-squares weights instead, which C3 maps onto itself.

Test object: hBN on opticx's own 30x30 mesh with NON-INTERACTING excitons (one transition each, energies from the .omesp of a
single-particle run), the lowest ones up to a gap of the transition spectrum: maximally unresolved envelopes, a set closed
under C3. For any two multiplets A, B of degenerate excitons, T_AB = sum_{n in A, m in B, n != m} X^a_nm conj(X^b_nm) (a, b in
x, y) must lie in span{1, eps}; the weighted residual measures the anisotropy of X_nm (axis stencil 0.45, star 0.002). The
diagonal X_nn is left out: for a single-transition envelope it is a band Berry-connection difference at one k, which depends on
the band-phase convention (the same reason a truncated non-interacting SHG is not a C3 test).
Checks: log line and pair residual < 1e-2 on the hexagonal mesh; the excitonic cache records the stencil (X_nm method code 3);
a cache with code 2 (two-axis stencil) stops the run and says why; a stretched (non-hexagonal) lattice keeps the axis stencil
(no log line, code 2).

Usage: python3 tools/check_hex_stencil.py --opticx bin/opticx --root . --workdir bin/check_hex_stencil
"""
import argparse, glob, os, shutil, struct, subprocess, sys
import numpy as np

NCELLS = 30


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--opticx', required=True); ap.add_argument('--root', default='.')
    ap.add_argument('--workdir', default='check_hex_stencil')
    a = ap.parse_args()
    root, work, opticx = os.path.abspath(a.root), os.path.abspath(a.workdir), os.path.abspath(a.opticx)
    sys.path.insert(0, os.path.join(root, 'tools'))
    try:
        import h5py
    except ImportError:
        print('SKIP: h5py not available'); sys.exit(0)
    from check_tb_hermiticity import read_tb, write_tb
    tb0 = os.path.join(root, 'hBN_tb.dat')
    if not os.path.exists(tb0):
        print('SKIP: missing ' + tb0); sys.exit(0)
    shutil.rmtree(work, ignore_errors=True); os.makedirs(work)
    nfail = 0

    def check(label, ok, msg=''):
        nonlocal nfail
        print(('  PASS   ' if ok else '  FAIL   ') + label + ('   ' + msg if msg else '')); nfail += 0 if ok else 1

    def ni_archive(omesp, out):
        """Non-interacting excitons, the lowest transitions up to the largest gap among the first ~200."""
        b = open(omesp, 'rb').read(); nk, nb = struct.unpack_from('<2i', b, 4)
        kb = np.frombuffer(b, '<f8', 3*nk, 12).reshape(3, nk).T
        off = 12 + 3*nk*8; rec = nb*8 + 2*3*nb*nb*16 + 9*nb*nb*8 + 9*nb*nb*16
        ek = np.array([np.frombuffer(b, '<f8', nb, off + i*rec) for i in range(nk)])*27.211385
        E = ek[:, 1] - ek[:, 0]; order = np.argsort(E, kind='stable'); Es = E[order]
        ncut = 150 + 1 + int(np.argmax(np.diff(Es[150:220])))
        with h5py.File(out, 'w') as f:
            A = f.attrs
            A['format'] = 'xatu-excitons'; A['format_version'] = 1; A['n_kpoints'] = nk; A['ncells'] = NCELLS
            A['submesh_factor'] = 1; A['fermi_level'] = 0; A['n_excitons'] = ncut; A['n_excitons_computed'] = ncut
            A['exciton_basis_dim'] = nk; A['tamm_dancoff'] = np.bool_(True)
            A['valence_bands'] = np.array([0], np.int64); A['conduction_bands'] = np.array([1], np.int64)
            f.create_dataset('basis', data=np.array([(0, 1, ik) for ik in range(nk)], np.int64))
            f.create_dataset('kpoints', data=kb/0.52917721067121)
            f.create_dataset('summary/energies', data=Es[:ncut])
            g = f.create_group('excitons')
            for n in range(ncut):
                st = np.zeros(nk, complex); st[order[n]] = 1.0
                s_ = g.create_group(f'{n+1:04d}'); s_.create_dataset('eigval', data=Es[n]); s_.create_dataset('state', data=st)
        return ncut, 1e3*(Es[ncut] - Es[ncut - 1])

    def run_pair(tag, tb):
        d = os.path.join(work, tag); os.makedirs(os.path.join(d, 'ex'))
        common = f'# Periodic dimensions\n2\n# Wannier90_filename\n{tb}\n'
        tail = '# Nfermi\n1\n# Response\nshg\n# Energy_variables\n2 8 0.2 61\n'
        sp = common + f'# Xatu_interface\nfalse\n# Ncells\n{NCELLS}\n# Bandlist\n0 1\n# OME_sp\nnonlinear\n' + tail
        open(os.path.join(d, f'{tag}_sp.txt'), 'w').write(sp)
        r = subprocess.run([opticx, f'{tag}_sp.txt'], cwd=d, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        omesp = glob.glob(os.path.join(d, 'ome_nonlinear_sp_*.omesp'))
        if r.returncode != 0 or not omesp:
            return None
        h5 = os.path.join(d, 'ni.h5'); ncut, gap = ni_archive(omesp[0], h5)
        shutil.copy(omesp[0], os.path.join(d, 'ex'))
        ex = common + f'# Xatu_interface\ntrue\n{h5}\n# OME_sp\nnone\n# OME_ex\nnonlinear\n# Cache_ome_ex\nwrite\n' + tail
        open(os.path.join(d, 'ex', f'{tag}_ex.txt'), 'w').write(ex)
        r = subprocess.run([opticx, f'{tag}_ex.txt'], cwd=os.path.join(d, 'ex'), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log = r.stdout.decode('utf-8', 'replace'); open(os.path.join(d, 'ex', 'run.log'), 'w').write(log)
        cache = glob.glob(os.path.join(d, 'ex', '*.omeex2'))
        return dict(rc=r.returncode, log=log, cache=cache[0] if cache else None, d=os.path.join(d, 'ex'), ex=ex, ncut=ncut, gap=gap)

    def read_cache(f):
        """Header offset of the X_nm method code, the code, energies and X[m, n, a] = X^a_nm."""
        b = open(f, 'rb').read(); o = 14
        ver, nlen = struct.unpack_from('<2i', b, o); o += 8 + nlen
        nk, nv, nc, N = struct.unpack_from('<4i', b, o); o += 16
        nb = struct.unpack_from('<i', b, o)[0]; o += 4 + 4*nb
        mpos = o; code = struct.unpack_from('<i', b, o)[0]; o += 4
        if ver >= 4: o += 12
        if ver == 5: o += 4
        E = np.frombuffer(b, '<f8', N, o); o += 8*N + 2*16*3*N
        X = np.frombuffer(b, '<c16', 3*N*N, o).reshape(N, N, 3)
        return mpos, code, E, X

    def pair_residual(E, X, tol=1e-4/27.211385):
        N = len(E); grp = np.concatenate([[0], np.cumsum(np.diff(E) > tol)]); ng = grp[-1] + 1
        Gm = np.zeros((ng, N)); Gm[grp, np.arange(N)] = 1
        Xa = X[:, :, :2].copy(); Xa[np.arange(N), np.arange(N), :] = 0
        T = np.zeros((ng, ng, 2, 2), complex)
        for p in range(2):
            for q in range(2):
                T[:, :, p, q] = Gm @ (Xa[:, :, p]*Xa[:, :, q].conj()).T @ Gm.T
        tr = 0.5*(T[:, :, 0, 0] + T[:, :, 1, 1]); an = 0.5*(T[:, :, 0, 1] - T[:, :, 1, 0])
        PT = np.zeros_like(T); PT[:, :, 0, 0] = PT[:, :, 1, 1] = tr; PT[:, :, 0, 1] = an; PT[:, :, 1, 0] = -an
        return np.linalg.norm((T - PT).reshape(ng, ng, 4), axis=2).sum()/np.linalg.norm(T.reshape(ng, ng, 4), axis=2).sum()

    # hexagonal mesh
    h = run_pair('hex', tb0)
    ok = h is not None and h['rc'] == 0 and h['cache']
    check('hBN, non-interacting excitons: runs and writes the cache', bool(ok),
          f"{h['ncut']} excitons, cut at a {h['gap']:.0f} meV gap" if h else 'single-particle run failed')
    if ok:
        check('log: X_nm on the hexagonal star', 'hexagonal star of mesh directions' in h['log'])
        mpos, code, E, X = read_cache(h['cache'])
        r = pair_residual(E, X)
        check('off-diagonal X_nm C3-covariant per multiplet pair (< 1e-2; two-axis stencil 0.45)', r < 1e-2, f'{r:.1e}')
        check('cache records the hexagonal stencil (X_nm method code 3)', code == 3, f'code {code}')
        b = bytearray(open(h['cache'], 'rb').read()); struct.pack_into('<i', b, mpos, 2)
        open(h['cache'], 'wb').write(bytes(b))
        txt = h['ex'].replace('# OME_ex\nnonlinear\n', '# OME_ex\nnone\n').replace('# Cache_ome_ex\nwrite\n', '# Cache_ome_ex\nread\n')
        open(os.path.join(h['d'], 'read_old.txt'), 'w').write(txt)
        rr = subprocess.run([opticx, 'read_old.txt'], cwd=h['d'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        lg = rr.stdout.decode('utf-8', 'replace')
        check('a cache made with the two-axis stencil (code 2) stops the run and says why', rr.returncode != 0 and
              'used the two-axis stencil' in lg and 'Entering sigma' not in lg, f'rc {rr.returncode}')

    # stretched lattice: not hexagonal
    m = read_tb(tb0); lat = m['lat'][:]
    v = [float(x) for x in lat[1].split()]; lat[1] = ' '.join(f'{1.1*x:.10f}' for x in v)
    m2 = dict(m); m2['lat'] = lat
    tbs = os.path.join(work, 'hBN_stretched_tb.dat'); write_tb(tbs, m2, m['H'], m['RP'])
    s = run_pair('stretched', tbs)
    ok = s is not None and s['rc'] == 0 and s['cache']
    check('stretched lattice: runs', bool(ok))
    if ok:
        check('stretched lattice: two-axis stencil (no log line, cache code 2)', 'hexagonal star' not in s['log'] and
              read_cache(s['cache'])[1] == 2, f"code {read_cache(s['cache'])[1]}")

    if nfail:
        print(f'\nFAILED: {nfail} check(s)'); sys.exit(1)
    print('\nALL CHECKS PASSED')


if __name__ == '__main__':
    main()
