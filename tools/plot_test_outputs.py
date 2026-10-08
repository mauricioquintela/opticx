#!/usr/bin/env python3
"""Simple plots of the spectra written by the OPTICX tests.

Every test that has a spectrum writes a plain-text file whose first lines are
    # kind: <name>
    # columns: <space separated column names>
followed by numbers. This script reads such a file and saves a PNG next to it (or to the given path):

    python3 tools/plot_test_outputs.py bin/test_run/shg_real_data_spectra.dat [out.png]

kinds: shg_real_data, shift_real_data, shg_consistency, sp_shift, ex_rectification. Plotting is best-effort: a missing matplotlib only
prints a message (the Makefile calls this with a leading '-', so it never fails a test).
"""
import sys
import numpy as np


def read(fn):
    kind, cols = None, None
    with open(fn) as f:
        for line in f:
            if line.startswith('# kind:'):
                kind = line.split(':', 1)[1].strip()
            elif line.startswith('# columns:'):
                cols = line.split(':', 1)[1].split()
            elif not line.startswith('#'):
                break
    d = np.loadtxt(fn, comments='#', ndmin=2)
    return kind, {c: d[:, i] for i, c in enumerate(cols)}


def complex_abs(D, name):
    return np.hypot(D[name + '_re'], D[name + '_im'])


def plot(fn, out=None):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('plot_test_outputs: matplotlib not available, no plot written'); return None
    kind, D = read(fn)
    out = out or fn.rsplit('.', 1)[0] + '.png'
    if kind == 'shg_real_data':
        fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
        x = D['Ew_eV'] if 'Ew_eV' in D else D['E2w_eV'] / 2   # files written before 2026-09-24 stored 2*hbar*omega
        ax[0].plot(x, complex_abs(D, 'xxx_B'), 'k-', lw=2, label='method B (production)')
        ax[0].plot(x, complex_abs(D, 'xxx_A_Pi'), 'C1--', lw=1.5, label=r'method A, $\Pi=-iEX$')
        ax[0].plot(x, complex_abs(D, 'xxx_A_P'), 'C3:', lw=2, label='method A on the code\'s bare P (wrong)')
        ax[0].set_title(r'$|\sigma_{xxx}|$'); ax[1].set_title(r'$|\sigma_{xyy}|$ (b,c-symmetrised)')
        ax[1].plot(x, complex_abs(D, 'xyy_B'), 'k-', lw=2, label='method B')
        ax[1].plot(x, complex_abs(D, 'xyy_A_Pi'), 'C1--', lw=1.5, label=r'method A, $\Pi=-iEX$')
        for a in ax: a.set_xlabel(r'$\hbar\omega$ (eV, fundamental)'); a.set_ylabel(r'$\mu$A nm / V$^2$'); a.legend(fontsize=8)
        fig.suptitle('Excitonic SHG, hBN 30x30, real data (test_shg_real_data)')
    elif kind == 'shift_real_data':
        fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
        for a, c in zip(ax, ('xxx', 'xyy', 'yxy')):
            a.plot(D['E_eV'], D[c + '_code'], 'k-', lw=2, label='opticx (Lorentzian)')
            a.plot(D['E_eV'], D[c + '_ref'], 'C1--', lw=1.5, label='paper Eq. (11), independent')
            a.set_title(r'$\sigma_{%s}$' % c); a.set_xlabel(r'$\hbar\omega$ (eV)'); a.legend(fontsize=8)
        ax[0].set_ylabel(r'$\mu$A nm / V$^2$')
        fig.suptitle('Excitonic shift conductivity, hBN 30x30, 100 excitons (test_shift_real_data)')
    elif kind == 'shg_consistency':
        fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
        x = D['omega']
        ax[0].plot(x, complex_abs(D, 'xxx_B'), 'k-', lw=2, label='method B'); ax[0].plot(x, complex_abs(D, 'xxx_A_scalar'), 'C1--', lw=1.5, label='method A scalar')
        ax[0].plot(x, complex_abs(D, 'xxx_A_matrix'), 'C2:', lw=2, label='method A matrix')
        ax[1].plot(x, D['xxx_B_re'], 'k-', label='Re B'); ax[1].plot(x, D['xxx_B_im'], 'C0-', label='Im B')
        ax[0].set_title(r'$|\sigma_{xxx}|$'); ax[1].set_title('method B, real / imaginary part')
        for a in ax: a.set_xlabel(r'$\omega$ (arb. units, synthetic model)'); a.legend(fontsize=8)
        fig.suptitle('SHG consistency test (synthetic complex-Hermitian data, eta = 0.05)')
    elif kind == 'sp_shift':
        fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
        for a, c in zip(ax, ('xxx', 'xyy', 'yxy')):
            a.plot(D['E_eV'], D[c + '_sp'], 'k-', lw=2, label='opticx sp (shift_shiftvector)')
            a.plot(D['E_eV'], D[c + '_exact'], 'C1--', lw=1.5, label='exact IPA, Eq. (9) (NumPy)')
            a.plot(D['E_eV'], D[c + '_ex'], 'C2:', lw=2, label='excitonic IPA limit, Eq. (10)')
            a.set_title(r'$\sigma_{%s}$' % c); a.set_xlabel(r'$\hbar\omega$ (eV)'); a.legend(fontsize=7)
        ax[0].set_ylabel(r'$\mu$A nm / V$^2$')
        fig.suptitle('Single-particle shift, non-interacting hBN 30x30 (check_sp_shift)')
    elif kind == 'ex_rectification':
        fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
        ax[0].plot(D['E_eV'], D['sp_re_xxx'], 'k-', lw=2, label='single-particle rectification')
        ax[0].plot(D['E_eV'], D['ex_re_xxx'], 'C1--', lw=1.5, label='excitonic rectification, NI limit')
        ax[0].set_title(r'Re $\sigma^{xxx}(0;\omega,-\omega)$')
        ax[1].plot(D['E_eV'], D['sp_im_yyz'], 'k-', lw=2, label='single-particle, antisym Im')
        ax[1].plot(D['E_eV'], D['ex_im_yyz'], 'C1--', lw=1.5, label='excitonic, antisym Im')
        ax[1].set_title(r'Im $\sigma^{yyz}$ (injection channel)')
        for a in ax: a.set_xlabel(r'$\hbar\omega$ (eV)'); a.legend(fontsize=7)
        ax[0].set_ylabel(r'$\mu$A nm / V$^2$')
        fig.suptitle('Excitonic rectification, non-interacting buckled hBN 30x30, eta 0.05 eV (check_ex_rectification)')
    else:
        print(f'plot_test_outputs: unknown kind {kind!r} in {fn}'); return None
    fig.tight_layout(); fig.savefig(out, dpi=110); plt.close(fig)
    print('plot written:', out)
    return out


if __name__ == '__main__':
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    plot(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
