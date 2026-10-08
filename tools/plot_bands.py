#!/usr/bin/env python3
"""Plot bands_<material>.dat (written by every OptiX run, or by Response = bands).

Uses the header: vertex labels and rows for the x ticks, the Bandlist window bands (drawn in colour), and the
VBM along the path to put the Fermi level at zero (--absolute keeps the Wannier-model energies).

Usage: python3 tools/plot_bands.py bands_MoS2.dat [out.png] [--emin -3 --emax 3] [--absolute]
"""
import argparse
import numpy as np


def read(fn):
    verts, window, vbm = [], [], None
    for line in open(fn):
        if not line.startswith('#'):
            break
        p = line.split()
        if line.startswith('# vertex ') and 'label' not in line:
            verts.append((p[3], int(p[4])))
        elif 'Bandlist window bands:' in line:
            window = [int(x) for x in line.split(':')[1].split()]
        elif 'VBM (band Nfermi)' in line:
            vbm = float(line.split("Nfermi)")[1].split()[0])
    return np.loadtxt(fn, comments='#'), verts, window, vbm


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('file'); ap.add_argument('out', nargs='?')
    ap.add_argument('--emin', type=float); ap.add_argument('--emax', type=float)
    ap.add_argument('--absolute', action='store_true', help='do not shift the VBM to zero')
    a = ap.parse_args()
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    B, verts, window, vbm = read(a.file)
    s, E = B[:, 3], B[:, 4:]
    shift = 0.0 if (a.absolute or vbm is None) else vbm
    fig, ax = plt.subplots(figsize=(6, 4.5))
    for n in range(E.shape[1]):
        inwin = (n + 1) in window
        ax.plot(s, E[:, n] - shift, color='#2a78d6' if inwin else '#6b6a64', lw=1.6 if inwin else 0.8)
    ticks = [s[r - 1] for _, r in verts]
    labels = ['Γ' if l.upper() in ('G', 'GAMMA') else l for l, _ in verts]
    for t in ticks:
        ax.axvline(t, color='#d6dbe3', lw=0.8)
    if shift:
        ax.axhline(0, color='#c4521f', lw=0.8, ls='--')
    ax.set_xticks(ticks, labels); ax.set_xlim(s[0], s[-1])
    if a.emin is not None or a.emax is not None:
        ax.set_ylim(a.emin, a.emax)
    ax.set_ylabel('E - E_VBM (eV)' if shift else 'E (eV)')
    ax.set_title(a.file.split('/')[-1] + ('   (Bandlist window in colour)' if window else ''), fontsize=9, loc='left')
    fig.tight_layout()
    out = a.out or a.file.rsplit('.', 1)[0] + '.png'
    fig.savefig(out, dpi=130)
    print('plot written:', out)


if __name__ == '__main__':
    main()
