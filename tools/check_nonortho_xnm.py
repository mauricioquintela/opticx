#!/usr/bin/env python3
"""Excitonic inter-exciton elements X_nm on a non-orthonormal model whose overlap depends on k.

X_nm is built with a covariant derivative of the exciton envelopes: the neighbour's envelope is carried into the basis at
k with the overlaps W = (S(k) c(k))^H c(k') of the window states, plus the Wannier-position connection c^H A c. When S
depends on k, W is not unitary, and the valence part of the connection must then be built from A^H (ome_ex build_xcov).
With A there the excitonic responses of such a model were wrong (this test: excitonic SHG 56% off at 30x30 and at 60x60);
orthonormal models and overlaps without intercell terms were unaffected.

Model: hBN_tb.dat and a REWRITE of it in a non-orthogonal basis with intercell mixing, phi'_i(0) = sum_{R,j} T_ji(R) phi_j(R)
(T on R = 0 and the four in-plane neighbours): H' = T^H H T, S' = T^H T, r' = T^H (r + R) T cell by cell. The two describe
the same crystal, so every response must agree. Excitons: the COMPLETE set of non-interacting transitions on opticx's own
30x30 mesh (energies from the .omesp of a single-particle run), for which the excitonic route must reproduce the
single-particle one. Checks: the rewrite's overlap reaches between cells; single particle rewrite = plain (1e-6);
excitonic = single particle for the plain model and for the rewrite (5e-3, the mesh error of X_nm at 30x30: 2.6e-3);
excitonic rewrite = excitonic plain (5e-3).

Usage: python3 tools/check_nonortho_xnm.py --opticx bin/opticx --root . --workdir bin/check_nonortho_xnm
"""
import argparse, glob, os, shutil, struct, subprocess, sys
import numpy as np

NCELLS = 30


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--opticx', required=True); ap.add_argument('--root', default='.')
    ap.add_argument('--workdir', default='check_nonortho_xnm')
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

    # the rewrite: <k R1|r|l R3> = <k 0|r|l R3-R1> + R1 delta_kl delta_{R1,R3} in the orthonormal basis
    m = read_tb(tb0); lat = m['latv']
    H0, r0 = m['H'], m['RP']
    rng = np.random.default_rng(5)
    T = {(0, 0, 0): np.array([[1.0, 0.3 + 0.1j], [-0.2j, 1.0]])}
    for R in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)):
        T[R] = 0.12*(rng.standard_normal((2, 2)) + 1j*rng.standard_normal((2, 2)))
    add = lambda x, y: tuple(p + q for p, q in zip(x, y)); sub = lambda x, y: tuple(p - q for p, q in zip(x, y))
    span = sorted(set(sub(add(R, R2), R1) for R in m['Rs'] for R1 in T for R2 in T) | set(m['Rs']))
    H, S, RP = {}, {}, {}
    for R in span:
        h = np.zeros((2, 2), complex); s = np.zeros((2, 2), complex); r = np.zeros((3, 2, 2), complex)
        for R1, T1 in T.items():
            for R2, T2 in T.items():
                R3 = sub(add(R, R2), R1); A = T1.conj().T
                if R3 in H0:
                    h += A @ H0[R3] @ T2; r += np.einsum('ij,cjk,kl->cil', A, r0[R3], T2)
                if R3 == (0, 0, 0):
                    s += A @ T2; r += (np.array(R1, float) @ lat)[:, None, None]*(A @ T2)[None]
        H[R], S[R], RP[R] = h, s, r
    m2 = dict(m); m2['Rs'] = span; m2['deg'] = [1]*len(span)
    tbr = os.path.join(work, 'hBN_rewrite_tb.dat'); write_tb(tbr, m2, H, RP, S)
    smax = max(abs(S[R]).max() for R in span if R != (0, 0, 0))
    check('the rewrite has an intercell overlap (S depends on k)', smax > 0.1, f'max |S(R != 0)| {smax:.2f}')

    def ni_archive(omesp, out):
        """Complete non-interacting excitons: every transition (k, v, c), unit envelope, E = E_c - E_v."""
        b = open(omesp, 'rb').read(); nk, nb = struct.unpack_from('<2i', b, 4)
        kb = np.frombuffer(b, '<f8', 3*nk, 12).reshape(3, nk).T
        off = 12 + 3*nk*8; rec = nb*8 + 2*3*nb*nb*16 + 9*nb*nb*8 + 9*nb*nb*16
        ek = np.array([np.frombuffer(b, '<f8', nb, off + i*rec) for i in range(nk)])*27.211385
        E = ek[:, 1] - ek[:, 0]; order = np.argsort(E, kind='stable')
        with h5py.File(out, 'w') as f:
            A = f.attrs
            A['format'] = 'xatu-excitons'; A['format_version'] = 1; A['n_kpoints'] = nk; A['ncells'] = NCELLS
            A['submesh_factor'] = 1; A['fermi_level'] = 0; A['n_excitons'] = nk; A['n_excitons_computed'] = nk
            A['exciton_basis_dim'] = nk; A['tamm_dancoff'] = np.bool_(True)
            A['valence_bands'] = np.array([0], np.int64); A['conduction_bands'] = np.array([1], np.int64)
            f.create_dataset('basis', data=np.array([(0, 1, ik) for ik in range(nk)], np.int64))
            f.create_dataset('kpoints', data=kb/0.52917721067121)
            f.create_dataset('summary/energies', data=E[order])
            g = f.create_group('excitons')
            for n in range(nk):
                st = np.zeros(nk, complex); st[order[n]] = 1.0
                s_ = g.create_group(f'{n+1:04d}'); s_.create_dataset('eigval', data=E[order[n]]); s_.create_dataset('state', data=st)

    def run_pair(tag, tb, ortho):
        d = os.path.join(work, tag); os.makedirs(os.path.join(d, 'ex'))
        common = f'# Periodic dimensions\n2\n# Wannier90_filename\n{tb}\n'
        tail = f'# Nfermi\n1\n# Response\nshg\n# Energy_variables\n2 8 0.2 121\n# Orthonormal\n{ortho}\n'
        sp = common + f'# Xatu_interface\nfalse\n# Ncells\n{NCELLS}\n# Bandlist\n0 1\n# OME_sp\nnonlinear\n' + tail
        open(os.path.join(d, f'{tag}_sp_shg.txt'), 'w').write(sp)
        r1 = subprocess.run([opticx, f'{tag}_sp_shg.txt'], cwd=d, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        omesp = glob.glob(os.path.join(d, 'ome_nonlinear_sp_*.omesp'))
        if r1.returncode != 0 or not omesp:
            return None, None
        h5 = os.path.join(d, 'ni.h5'); ni_archive(omesp[0], h5)
        shutil.copy(omesp[0], os.path.join(d, 'ex'))
        ex = common + f'# Xatu_interface\ntrue\n{h5}\n# OME_sp\nnone\n# OME_ex\nnonlinear\n' + tail
        open(os.path.join(d, 'ex', f'{tag}_ex_shg.txt'), 'w').write(ex)
        r2 = subprocess.run([opticx, f'{tag}_ex_shg.txt'], cwd=os.path.join(d, 'ex'), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        open(os.path.join(d, 'ex', 'run.log'), 'wb').write(r2.stdout)
        if r2.returncode != 0:
            return None, None
        ld = lambda f: (lambda x: x[:, 1::2] + 1j*x[:, 2::2])(np.loadtxt(glob.glob(f)[0]))
        return ld(os.path.join(d, 'shg_sp_lengthgauge_*.dat')), ld(os.path.join(d, 'ex', 'shg_ex_lengthgauge_*.dat'))

    rel = lambda x, y: np.abs(x - y).max()/np.abs(y).max()
    sp0, ex0 = run_pair('plain', tb0, 'true')
    sp1, ex1 = run_pair('rewrite', tbr, 'false')
    ok = all(x is not None for x in (sp0, ex0, sp1, ex1))
    check('all four runs completed', ok)
    if ok:
        check('single particle: rewrite = plain (1e-6)', rel(sp1, sp0) < 1e-6, f'{rel(sp1, sp0):.1e}')
        check('plain hBN: excitonic (complete NI) = single particle (5e-3)', rel(ex0, sp0) < 5e-3, f'{rel(ex0, sp0):.1e}')
        check('rewrite: excitonic (complete NI) = single particle (5e-3)', rel(ex1, sp1) < 5e-3, f'{rel(ex1, sp1):.1e}')
        check('excitonic: rewrite = plain (5e-3)', rel(ex1, ex0) < 5e-3, f'{rel(ex1, ex0):.1e}')

    if nfail:
        print(f'\nFAILED: {nfail} check(s)'); sys.exit(1)
    print('\nALL CHECKS PASSED')


if __name__ == '__main__':
    main()
