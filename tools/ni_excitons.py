#!/usr/bin/env python3
"""Synthetic NON-INTERACTING exciton files (Xatu format) for a two-band tight-binding model.

Exciton n is a single transition (v, c, k_n): its energy is E_n = e_c(k_n) - e_v(k_n) (eV, sorted ascending, as Xatu
does) and its wavefunction is a unit vector on the (v,c,k) basis. Fed to opticx, the excitonic shift/SHG code then
evaluates the independent-particle limit of the excitonic formulas on exactly the same k-grid as the single-particle
code, which makes the two directly comparable (sign, normalisation, index roles).

The k list (and the band-index columns) are taken from a real Xatu .states file, so the grid ordering is guaranteed to
be the one opticx expects. The 2-band restriction is only for the energies: E_n = w[1]-w[0].

Usage:  python3 ni_excitons.py hBN_tb.dat hBN_N30.states hBN_NI30      # writes hBN_NI30.eigval and hBN_NI30.states
"""
import sys
import numpy as np


def make_ni(tb_file, states_in, out_prefix):
    lines = open(tb_file).read().split('\n')
    lat = np.array([list(map(float, lines[i].split())) for i in (1, 2, 3)])
    norb, nR = int(lines[4].split()[0]), int(lines[5].split()[0])
    if norb != 2:
        raise SystemExit(f"ni_excitons: only two-band models are supported (norb = {norb})")
    tok = iter(' '.join(lines[6:]).split())
    deg = [int(next(tok)) for _ in range(nR)]
    H = {}
    for iR in range(nR):
        R = (int(next(tok)), int(next(tok)), int(next(tok))); M = np.zeros((2, 2), complex)
        for _ in range(4):
            i, j = int(next(tok)) - 1, int(next(tok)) - 1
            M[i, j] = float(next(tok)) + 1j * float(next(tok))
        H[R] = M / deg[iR]
    with open(states_in) as f:
        n = int(f.readline().split()[0]); kl = [f.readline().split() for _ in range(n)]
    K = np.array([[float(x) for x in l[:3]] for l in kl]); tail = [l[3:] for l in kl]      # Angstrom^-1
    E = np.zeros(n)
    for j, k in enumerate(K):
        Hk = sum(np.exp(1j * (k @ (R[0] * lat[0] + R[1] * lat[1] + R[2] * lat[2]))) * M for R, M in H.items())
        w = np.linalg.eigvalsh(Hk); E[j] = w[1] - w[0]
    order = np.argsort(E, kind='stable')          # exciton m lives on basis element (k index) order[m]
    with open(out_prefix + '.eigval', 'w') as f:
        f.write(f"{int(round(n ** 0.5))}\n{n}\n{n}\n")
        for m in range(n):
            f.write(f"{E[order[m]]:14.8f}\n")
    with open(out_prefix + '.states', 'w') as f:
        f.write(f"{n}\n")
        for j in range(n):
            f.write(f"{K[j, 0]:12.7f} {K[j, 1]:12.7f} {K[j, 2]:12.7f} {' '.join(tail[j])}\n")
        row = np.zeros(2 * n)
        for m in range(n):
            row[:] = 0; row[2 * order[m]] = 1.0
            f.write(' '.join(f"{x:.1f}" for x in row) + '\n')
    return n, E.min(), E.max()


if __name__ == '__main__':
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    n, emin, emax = make_ni(*sys.argv[1:])
    print(f"wrote {sys.argv[3]}.eigval/.states: {n} excitons, transition energies {emin:.4f}..{emax:.4f} eV")
