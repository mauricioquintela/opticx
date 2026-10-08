#!/usr/bin/env python3
"""Band structure along a k-path (bands_<material>.dat), the Kpath input and Response = bands.

Kpath: one vertex per line, "f1 f2 f3 n", reduced coordinates along the reciprocal lattice vectors and the number of
points from this vertex to the next; n = 1 gives the vertex alone and a jump to the next one; the last vertex carries
1 (included) or 0. Without Kpath a default for the lattice is used (hexagonal 2D: Gamma-M-K-Gamma, with K found as
the zone corner, which is (1/3,1/3) or (2/3,1/3) depending on the angle between the reciprocal vectors).

Checks (a few seconds):
 1. Eigenvalues at every path point equal an independent NumPy diagonalisation of the same _tb.dat (hBN,
    buckled hBN, the 34-orbital MoS2 model), 1e-8 eV.
 2. hBN default path: labels G M K G, the K vertex is a zone corner, the gap along the path is the K gap
    (smallest direct gap at the K vertex).
 3. Kpath with a jump: the row count is the sum of the counts, the vertex rows in the header are right, the path
    length does not advance across the jump, Kpath_labels appear in the header.
 4. Response = bands stops after the band structure (no matrix-element or spectrum files).
 5. A malformed Kpath line stops with a message.

Usage: python3 tools/check_bands.py --opticx bin/opticx --root . --workdir bin/check_bands
"""
import argparse, os, subprocess, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ipa_shift_numpy import read_tb, ANG, HA

IN = """# Periodic dimensions
2
# Wannier90_filename
{tb}
# Xatu_interface
false
# Ncells
6
# Bandlist
{bands}
# Nfermi
{nf}
# Response
bands
{extra}"""


def numpy_bands(tbfile, K):
    tb = read_tb(tbfile); lat = tb['lat']
    E = []
    for k in K:
        Hk = sum(np.exp(1j * (k @ ((R[0] * lat[0] + R[1] * lat[1] + R[2] * lat[2]) / ANG))) * M for R, M in tb['H'].items())
        E.append(np.linalg.eigvalsh(Hk))
    return np.array(E)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--opticx', required=True); ap.add_argument('--root', default='.')
    ap.add_argument('--workdir', default='check_bands')
    a = ap.parse_args()
    root, work, opticx = os.path.abspath(a.root), os.path.abspath(a.workdir), os.path.abspath(a.opticx)
    models = {'hBN': (os.path.join(root, 'hBN_tb.dat'), '0 1', 1),
              'buckled-hBN': (os.path.join(root, 'buckled-hBN_tb.dat'), '0 1', 1),
              'MoS2_spin_wannier_07032024': (os.path.join(root, 'MoS2_spin_wannier_07032024_tb.dat'), '-1 0 1 2', 26)}
    for tb, _, _ in models.values():
        if not os.path.exists(tb):
            print('SKIP: missing ' + tb); sys.exit(0)
    os.makedirs(work, exist_ok=True)
    nfail = 0

    def check(label, ok, msg=''):
        nonlocal nfail
        print(('  PASS   ' if ok else '  FAIL   ') + label + ('   ' + msg if msg else '')); nfail += 0 if ok else 1

    def run(tag, mat, extra=''):
        tb, bands, nf = models[mat]
        d = os.path.join(work, tag); os.makedirs(d, exist_ok=True)
        for f in os.listdir(d):
            os.remove(os.path.join(d, f))
        name = f'{tag}.txt'
        open(os.path.join(d, name), 'w').write(IN.format(tb=tb, bands=bands, nf=nf, extra=extra))
        r = subprocess.run([opticx, name], cwd=d, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log = r.stdout.decode('utf-8', 'replace'); open(os.path.join(d, 'run.log'), 'w').write(log)
        f = os.path.join(d, f'bands_{mat}.dat')
        if not os.path.exists(f):
            return r.returncode, log, d, None, []
        head = [l for l in open(f) if l.startswith('#')]
        return r.returncode, log, d, np.loadtxt(f, comments='#'), head

    def vertices(head):
        out = []
        for l in head:
            if l.startswith('# vertex ') and 'label' not in l:
                p = l.split(); out.append((p[3], int(p[4]), float(p[5]), [float(x) for x in p[6:9]]))
        return out

    print('Band structure along a k-path\n')
    for mat in models:
        rc, log, d, B, head = run(f'{mat}_default', mat)
        if B is None:
            check(f'{mat}: default path written', False, f'rc={rc}'); continue
        E = numpy_bands(models[mat][0], B[:, :3]) * 1.0
        e = abs(np.sort(B[:, 4:], axis=1) - E).max()
        check(f'{mat}: eigenvalues on the path = independent NumPy diagonalisation', rc == 0 and e < 1e-8,
              f'max diff {e:.1e} eV over {len(B)} points (max 1e-8)')
        if mat == 'hBN':
            v = vertices(head)
            check('hBN default path: G M K G', [x[0] for x in v] == ['G', 'M', 'K', 'G'], str([x[0] for x in v]))
            kK = B[v[2][1] - 1, :2]; tbv = read_tb(models[mat][0])
            corner = abs(np.linalg.norm(kK) - np.linalg.norm(np.linalg.solve((tbv['lat'][:2, :2] / ANG).T,
                         np.eye(2))[0, :]) * 2 * np.pi / np.sqrt(3))
            gaps = B[:, 5] - B[:, 4]
            check('hBN default path: K vertex is the zone corner and the direct-gap minimum',
                  corner < 1e-6 and np.argmin(gaps) == v[2][1] - 1, f'|K| - |b|/sqrt3 = {corner:.1e}, gap at K {gaps[v[2][1]-1]:.4f} eV')

    path = """# Kpath
 0.0 0.0 0.0 40
 0.5 0.0 0.0 1
 0.6666667 0.3333333 0.0 30
 0.0 0.0 0.0 1
# Kpath_labels
G M K G
"""
    rc, log, d, B, head = run('hBN_kpath_jump', 'hBN', path)
    v = vertices(head)
    ok = rc == 0 and B is not None and len(B) == 40 + 1 + 30 + 1
    check('Kpath: row count = sum of the counts (40 + 1 + 30 + 1)', ok, f'{0 if B is None else len(B)} rows')
    if ok:
        check('Kpath: vertex rows and labels in the header', [(x[0], x[1]) for x in v] == [('G', 1), ('M', 41), ('K', 42), ('G', 72)],
              str([(x[0], x[1]) for x in v]))
        check('Kpath: the path length does not advance across the jump (count 1)', abs(B[41, 3] - B[40, 3]) < 1e-12,
              f's(M) {B[40, 3]:.6f}, s(K) {B[41, 3]:.6f}')
        check('Response = bands: no matrix-element or spectrum files',
              sorted(f for f in os.listdir(d) if not f.endswith('.txt') and f != 'run.log') == ['bands_hBN.dat'],
              str(sorted(os.listdir(d))))
    rc, log, d, B, head = run('hBN_kpath_bad', 'hBN', '# Kpath\n 0.0 0.0 40\n 0.5 0.0 0.0 1\n')
    check('malformed Kpath line: stops with a message', rc != 0 and 'Kpath line' in log, f'rc={rc}')
    if nfail:
        print(f'\nFAILED: {nfail} check(s)'); sys.exit(1)
    print('\nALL CHECKS PASSED')


if __name__ == '__main__':
    main()
