#!/usr/bin/env python3
"""Convert a Xatu .eigval/.states text pair into a Xatu exciton archive (the layout 'xatu -H' writes).

The result holds what opticx reads from an archive: the band window, the k mesh, the exciton energies and
the resonant envelopes, with every number exactly as parsed from the text files (so opticx gives the same
answer from either). The text files do not record the Xatu filling or the mesh, so they are arguments:
--nfermi is the number of filled bands of the Xatu run (opticx's Nfermi), --ndim the periodic dimension
(ncells = nk^(1/ndim); submesh runs cannot be converted). Calculation parameters that the text files do not
carry (potential, dielectric, ...) are absent; the root attribute 'converted_from' names the source files.

The states are streamed: the main process reads the .states lines in order and writes each dataset, while
--jobs worker processes parse the text (and gzip it, with --compression) on a bounded window of lines, so
memory stays at a few times the basis size per worker whatever the number of excitons (a complete 102x102
basis of 41616 states is 28 GB of coefficients). Parsing the text is the whole cost of a conversion; the
archive does not depend on --jobs. --compression N (gzip level 1-9) compresses the state datasets; opticx
reads either. If the conversion stops on an error, the partial archive is deleted.

Usage: python3 tools/xatu_text_to_h5.py hBN_N30.eigval hBN_N30.states --nfermi 1 --ndim 2 -o hBN_N30.h5 [-j 16]
"""
import argparse
import os
import sys
import zlib
from collections import deque
from concurrent.futures import ProcessPoolExecutor
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


def read_basis(f):
    """Header of an open .states file: the basis size and, per basis row, k and the (v, c) band pair."""
    dim = int(f.readline())
    basis = [f.readline().split() for _ in range(dim)]
    k = np.array([[float(x) for x in row[:3]] for row in basis])
    vc = np.array([[int(row[3]), int(row[4])] for row in basis], dtype=np.int64)
    return dim, k, vc


# numpy >= 1.23 parses with a C reader; older versions get the vectorised str -> float conversion
_C_LOADTXT = tuple(int(x) for x in np.__version__.split('.')[:2]) >= (1, 23)


def parse_state(line, path, i, dim):
    """One exciton envelope from its .states line (re, im pairs), as complex numbers."""
    text = line.decode('ascii')
    vals = np.loadtxt([text], ndmin=1) if _C_LOADTXT and text.strip() else np.array(text.split(), dtype=float)
    if vals.size != 2 * dim:
        raise ValueError(f'{path}: state {i + 1} has {vals.size // 2} coefficients, expected {dim}')
    return vals[0::2] + 1j * vals[1::2]


def parse_batch(block, path, first, dim, level):
    """Worker task: parse a block of whole .states lines, states first, first+1, ...; with level > 0 return
    them as zlib streams (the format of HDF5's gzip filter)."""
    out = []
    for j, line in enumerate(block.split(b'\n')[:-1]):
        state = parse_state(line, path, first + j, dim)
        out.append(zlib.compress(state.tobytes(), level) if level else state)
    return out


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
    ap.add_argument('--compression', type=int, default=0, metavar='N',
                    help='gzip level 1-9 for the state datasets (default 0: uncompressed)')
    ap.add_argument('-j', '--jobs', type=int, default=0, metavar='N',
                    help='worker processes that parse the states (default: the CPUs available to the job, at most 32)')
    a = ap.parse_args()
    if not 0 <= a.compression <= 9:
        sys.exit('ERROR: --compression must be between 0 and 9')
    if a.jobs < 0:
        sys.exit('ERROR: --jobs must be positive')
    if a.jobs == 0:
        a.jobs = min(32, len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else os.cpu_count() or 1)

    ncell_file, dim_e, energies = read_eigval(a.eigval)
    fs = open(a.states, 'rb')
    dim, k, vc = read_basis(fs)
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

    try:
        write_archive(a, fs, dim, k, vc, kpoints, kidx, nk, ncells, vbands, cbands, energies)
    except BaseException as e:      # incl. KeyboardInterrupt: never leave a truncated archive behind
        fs.close()
        if os.path.exists(a.output):
            os.remove(a.output)
        if isinstance(e, ValueError):
            sys.exit(f'ERROR: {e}; no archive written')
        raise
    fs.close()
    print(f'{a.output}: {energies.size} excitons, {nk} k points, valence {vbands}, conduction {cbands}')


def write_archive(a, fs, dim, k, vc, kpoints, kidx, nk, ncells, vbands, cbands, energies):
    gz = dict(compression='gzip', compression_opts=a.compression) if a.compression else {}
    step = max(1, energies.size // 10)
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
        for i, state in enumerate(parsed_states(a, fs, dim, energies.size)):
            g = ex.create_group(str(i + 1).zfill(width))
            g.attrs['index'] = np.int64(i + 1)
            g.create_dataset('eigval', data=energies[i])
            if a.compression:       # one chunk per state, deflated by a worker: the gzip filter's own format
                d = g.create_dataset('state', shape=(dim,), dtype=np.complex128, chunks=(dim,), **gz)
                d.id.write_direct_chunk((0,), state)
            else:
                g.create_dataset('state', data=state)
            if energies.size >= 1000 and (i + 1) % step == 0:
                print(f'  {i + 1} / {energies.size} states written', file=sys.stderr, flush=True)


def parsed_states(a, fs, dim, n):
    """The n states of the open .states file in order, parsed by a.jobs workers over a bounded window.

    The file is read in blocks cut at a line end (much faster than reading line by line: one state line is
    several MB on large meshes); each block is one task, a few states at least."""
    block_bytes = max(8 << 20, 64 * dim)          # >= ~2 states per block (~32 characters per coefficient)
    window = 3 * a.jobs                            # blocks in flight: bounds the memory
    with ProcessPoolExecutor(max_workers=a.jobs) as pool:
        pending = deque()
        first, rest = 0, b''
        while first < n or pending:
            while first < n and len(pending) < window:
                block = rest + fs.read(block_bytes)
                cut = block.rfind(b'\n') + 1
                if cut == 0 and len(block) == len(rest):          # end of file inside a state
                    cut, block = len(block) + 1, block + b'\n'
                if cut == 0:                                      # no line end yet: read more
                    rest = block
                    continue
                count = min(block.count(b'\n', 0, cut), n - first)
                if count == 0:                                    # end of file, states missing
                    raise ValueError(f'{a.states}: state {first + 1} is missing (end of file)')
                end = 0
                for _ in range(count):                            # never parse past state n
                    end = block.index(b'\n', end) + 1
                pending.append(pool.submit(parse_batch, block[:end], a.states, first, dim, a.compression))
                rest = block[end:] if first + count < n else b''
                first += count
            try:
                yield from pending.popleft().result()
            except BaseException:
                for fut in pending:
                    fut.cancel()
                raise

if __name__ == '__main__':
    main()
