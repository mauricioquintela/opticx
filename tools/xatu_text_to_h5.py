#!/usr/bin/env python3
"""Convert a Xatu .eigval/.states text pair into a Xatu exciton archive (the layout 'xatu -H' writes).

The result holds what opticx reads from an archive: the band window, the k mesh, the exciton energies and
the resonant envelopes, with every number exactly as parsed from the text files (so opticx gives the same
answer from either). The text files do not record the Xatu filling or the mesh, so they are arguments:
--nfermi is the number of filled bands of the Xatu run (opticx's Nfermi), --ndim the periodic dimension
(ncells = nk^(1/ndim); submesh runs cannot be converted). Calculation parameters that the text files do not
carry (potential, dielectric, ...) are absent; the root attribute 'converted_from' names the source files.

Usage: python3 tools/xatu_text_to_h5.py hBN_N30.eigval hBN_N30.states --nfermi 1 --ndim 2 -o hBN_N30.h5
"""
import argparse
import os
import sys
from datetime import datetime

import h5py
import numpy as np


def read_eigval(path):
    tokens = open(path).read().split()
    ncell, dim, nprint = int(tokens[0]), int(tokens[1]), int(tokens[2])
    energies = np.array([float(x) for x in tokens[3:3 + nprint]])
    if energies.size != nprint:
        sys.exit(f'ERROR: {path} lists {nprint} energies but holds {energies.size}')
    return ncell, dim, energies


def read_states(path, nstates):
    with open(path) as f:
        dim = int(f.readline())
        basis = [f.readline().split() for _ in range(dim)]
        k = np.array([[float(x) for x in row[:3]] for row in basis])
        vc = np.array([[int(row[3]), int(row[4])] for row in basis], dtype=np.int64)
        states = np.empty((nstates, dim), dtype=complex)
        for i in range(nstates):
            vals = np.array([float(x) for x in f.readline().split()])
            if vals.size != 2 * dim:
                sys.exit(f'ERROR: {path}: state {i + 1} has {vals.size // 2} coefficients, expected {dim}')
            states[i] = vals[0::2] + 1j * vals[1::2]
    return dim, k, vc, states


def band_window(vc):
    """Valence and conduction bands in basis order: valence fastest, then conduction, then k."""
    vbands = []
    for v in vc[:, 0]:
        if v in vbands:
            break
        vbands.append(int(v))
    nv = len(vbands)
    cbands = []
    for c in vc[::nv, 1]:
        if c in cbands:
            break
        cbands.append(int(c))
    return vbands, cbands


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('eigval')
    ap.add_argument('states')
    ap.add_argument('--nfermi', type=int, required=True, help='filled bands of the Xatu run (opticx Nfermi)')
    ap.add_argument('--ndim', type=int, required=True, help='periodic dimension')
    ap.add_argument('-o', '--output', required=True)
    a = ap.parse_args()

    ncell_file, dim_e, energies = read_eigval(a.eigval)
    dim, k, vc, states = read_states(a.states, energies.size)
    if dim != dim_e:
        sys.exit(f'ERROR: basis size {dim} in {a.states} but {dim_e} in {a.eigval}')
    vbands, cbands = band_window(vc)
    npairs = len(vbands) * len(cbands)
    if dim % npairs:
        sys.exit(f'ERROR: basis size {dim} is not a multiple of nv*nc = {npairs}')
    nk = dim // npairs
    kidx = np.arange(dim) // npairs
    expect = np.array([[vbands[r % len(vbands)], cbands[(r % npairs) // len(vbands)]] for r in range(dim)])
    if not np.array_equal(expect, vc):
        sys.exit('ERROR: the basis is not ordered valence fastest, then conduction, then k')
    kpoints = k[::npairs]
    if not np.array_equal(k, kpoints[kidx]):
        sys.exit('ERROR: k points within a k block of the basis differ')
    ncells = int(round(nk ** (1.0 / a.ndim)))
    if ncells ** a.ndim != nk:
        sys.exit(f'ERROR: {nk} k points is not a full ncells^{a.ndim} mesh (Xatu submesh runs cannot be converted)')
    if ncells != ncell_file:
        print(f'note: .eigval header says ncells = {ncell_file}, the mesh has {ncells} per axis; using {ncells}')

    with h5py.File(a.output, 'w', track_order=True) as f:
        f.attrs['format'] = 'xatu-excitons'
        f.attrs['format_version'] = np.int64(1)
        f.attrs['created'] = datetime.now().astimezone().isoformat(timespec='seconds')
        f.attrs['converted_from'] = f'{os.path.abspath(a.eigval)}, {os.path.abspath(a.states)}'
        f.attrs['units'] = 'energies in eV; lengths in the units of the system file; k in inverse length'
        f.attrs['ncells'] = np.int64(ncells)
        f.attrs['submesh_factor'] = np.int64(1)
        f.attrs['n_kpoints'] = np.int64(nk)
        f.attrs['valence_bands'] = np.array(vbands, dtype=np.int64)
        f.attrs['conduction_bands'] = np.array(cbands, dtype=np.int64)
        f.attrs['fermi_level'] = np.int64(a.nfermi - 1)
        f.attrs['tamm_dancoff'] = np.bool_(True)  # the text files hold the resonant block only
        f.attrs['exciton_basis_dim'] = np.int64(dim)
        f.attrs['n_excitons_computed'] = np.int64(energies.size)
        f.attrs['n_excitons'] = np.int64(energies.size)
        f.create_dataset('kpoints', data=kpoints)
        f.create_dataset('basis', data=np.column_stack([vc, kidx]).astype(np.int64))
        s = f.create_group('summary', track_order=True)
        s.create_dataset('energies', data=energies)
        ex = f.create_group('excitons', track_order=True)
        width = max(4, len(str(energies.size)))
        for i in range(energies.size):
            g = ex.create_group(str(i + 1).zfill(width))
            g.attrs['index'] = np.int64(i + 1)
            g.create_dataset('eigval', data=energies[i])
            g.create_dataset('state', data=states[i])
    print(f'{a.output}: {energies.size} excitons, {nk} k points, valence {vbands}, conduction {cbands}')


if __name__ == '__main__':
    main()
