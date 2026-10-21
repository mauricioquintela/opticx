#!/usr/bin/env python3
"""Scissor correction (keywords Scissor, eV, and Scissor_mode = r_unchanged | v_unchanged).

A scissor shifts every conduction band rigidly, E_c -> E_c + Delta. r_unchanged keeps the interband position elements and
rescales v_cv by (w_cv + Delta)/w_cv (Levine & Allan, PRL 63, 1719 (1989); Nastos et al., PRB 72, 045223 (2005));
v_unchanged keeps v_cv and rescales r_cv by w_cv/(w_cv + Delta) (the convention of Xatu's Kubo conductivity).

Checks (hBN, 30x30, a few seconds each):
  * Scissor = 0 gives the output of no Scissor keyword (to round-off) (shift and absorbance, single particle and excitonic).
  * r_unchanged, single-particle shift current: r and its generalised derivative are unchanged, so the spectrum is the
    unshifted one moved rigidly by Delta (Delta = 100 grid steps; 1e-10).
  * Non-interacting excitons (synthetic, energies + Delta) against the single-particle code with the same Scissor, in both
    modes: absorbance, shift, rectification and SHG agree as well as they do at Delta = 0 (the two routes build their
    matrix elements independently, so a scaling applied on one side only would show).
  * Scissor_mode changes the answer (both modes are really applied).
  * The Xatu archive: a recorded scissor different from Scissor stops the run before anything is computed; a matching
    one runs and says so; an archive without the record runs with a NOTE.
  * The excitonic cache: Scissor = 0 writes the version-3 layout; Scissor /= 0 writes version 4 with the scissor; a
    cache built with a different scissor or mode, or a version-3 cache read with Scissor /= 0, stops the run.
  * Guards: Sp_method = per_band with a Scissor stops; an unknown Scissor_mode stops.

Usage: python3 tools/check_scissor.py --opticx bin/opticx --root . --workdir bin/check_scissor
"""
import argparse, os, shutil, struct, subprocess, sys
import numpy as np

DELTA = 1.0          # eV


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--opticx', required=True); ap.add_argument('--root', default='.')
    ap.add_argument('--workdir', default='check_scissor')
    a = ap.parse_args()
    root, work, opticx = os.path.abspath(a.root), os.path.abspath(a.workdir), os.path.abspath(a.opticx)
    sys.path.insert(0, os.path.join(root, 'tools'))
    from ni_excitons import make_ni
    tb, states = os.path.join(root, 'hBN_tb.dat'), os.path.join(root, 'hBN_N30.states')
    for f in (tb, states):
        if not os.path.exists(f):
            print('SKIP: missing ' + f); sys.exit(0)
    os.makedirs(work, exist_ok=True)
    nfail = 0

    def check(label, ok, msg=''):
        nonlocal nfail
        print(('  PASS   ' if ok else '  FAIL   ') + label + ('   ' + msg if msg else '')); nfail += 0 if ok else 1

    def ni(delta):                       # synthetic non-interacting excitons, energies + delta (eV)
        pre = os.path.join(work, f'hBN_NI30_d{delta}')
        if not os.path.exists(pre + '.eigval'):
            make_ni(tb, states, pre)
            L = open(pre + '.eigval').read().split('\n'); E = [float(x) for x in L[3:] if x.strip()]
            open(pre + '.eigval', 'w').write('\n'.join(L[:3] + [f'{e + delta:14.8f}' for e in E]) + '\n')
        return pre

    def run(tag, resp, delta=None, mode=None, ev='4.0 13.0 0.1 900', exc=None, extra='', keep=False):
        d = os.path.join(work, tag); os.makedirs(d, exist_ok=True)
        for f in os.listdir(d):
            if not (keep and f.endswith('.omeex2')):
                os.remove(os.path.join(d, f))
        order = 'linear' if resp == 'absorbance' else 'nonlinear'
        xatu = exc if exc else f'{ni(delta or 0.0)}.eigval\n{ni(delta or 0.0)}.states'
        txt = (f'# Periodic dimensions\n2\n# Wannier90_filename\n{tb}\n# Xatu_interface\ntrue\n{xatu}\n'
               f'# Nfermi\n1\n# OME_sp\n{order}\n# OME_ex\n{order}\n# Response\n{resp}\n# Energy_variables\n{ev}\n')
        if delta is not None: txt += f'# Scissor\n{delta}\n'
        if mode is not None: txt += f'# Scissor_mode\n{mode}\n'
        name = f'{tag}.txt'; open(os.path.join(d, name), 'w').write(txt + extra)
        r = subprocess.run([opticx, name], cwd=d, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log = r.stdout.decode('utf-8', 'replace'); open(os.path.join(d, 'run.log'), 'w').write(log)
        return d, r.returncode, log

    def load(d, stem):
        f = [x for x in os.listdir(d) if x.startswith(stem)]
        if not f: return None, None
        A = np.loadtxt(os.path.join(d, f[0]), comments='#')
        if stem.startswith(('second', 'shg')):
            B = A[:, -54:]; return A[:, 0], B[:, 0::2] + 1j*B[:, 1::2]
        return A[:, 0], A[:, 1:28] if stem.startswith('shift') else A[:, 1:10]

    def files(d):
        return {f: np.loadtxt(os.path.join(d, f), comments='#') for f in os.listdir(d)
                if f.endswith('.dat') and not f.startswith('bands')}

    def same(A, B, tol=1e-12):     # same spectra to round-off (threaded sums are not bit-reproducible)
        return sorted(A) == sorted(B) and all(abs(A[f] - B[f]).max() <= tol*max(abs(B[f]).max(), 1e-300) for f in A)

    print('\nScissor = 0 is the identity\n')
    for resp in ('shift', 'absorbance'):
        d0, r0, _ = run(f'{resp}_nokey', resp)
        d1, r1, _ = run(f'{resp}_zero', resp, delta=0.0, mode='v_unchanged')
        f0, f1 = files(d0), files(d1)
        check(f'{resp}: Scissor = 0 gives the output of no Scissor keyword ({len(f0)} files, 1e-12)',
              r0 == 0 and r1 == 0 and same(f1, f0) and len(f0) >= 2, f'rc {r0}/{r1}')

    print('\nr_unchanged: the single-particle shift current moves rigidly by Delta\n')
    d0, r0, _ = run('rigid_0', 'shift', ev='3.0 13.0 0.1 1000')
    d1, r1, log = run('rigid_d', 'shift', delta=DELTA, mode='r_unchanged', ev='3.0 13.0 0.1 1000')
    w, s0 = load(d0, 'shift_sp'); _, s1 = load(d1, 'shift_sp')
    n = int(round(DELTA/(w[1] - w[0])))
    err = abs(s1[n:] - s0[:-n]).max()/abs(s0).max() if r0 == 0 and r1 == 0 else np.inf
    check(f'sp shift with Scissor = {DELTA} eV = unshifted spectrum moved by {n} grid steps', err < 1e-10, f'{err:.1e}')
    check('the log names the shift and the mode', 'Scissor: conduction bands shifted by' in log and 'r_unchanged' in log)

    print('\nNon-interacting excitons with energies + Delta against the single-particle code, both modes\n')
    pairs = {'absorbance': ('sigma_first_sp_real_', 'sigma_first_ex_real_'), 'shift': ('shift_sp_', 'shift_ex_'),
             'rectification': ('second_rectification_', 'second_ex_rectification_'), 'shg': ('shg_sp_', 'shg_ex_')}
    for resp, (ssp, sex) in pairs.items():
        base = None
        for delta, mode in ((None, None), (DELTA, 'r_unchanged'), (DELTA, 'v_unchanged')):
            d, rc, log = run(f'ni_{resp}_{mode}', resp, delta=delta, mode=mode)
            w, sp = load(d, ssp); _, ex = load(d, sex)
            if rc != 0 or sp is None or ex is None:
                check(f'{resp}, {mode}: run', False, f'rc {rc}'); continue
            diff = abs(ex - sp).max()/abs(sp).max()
            if mode is None:
                base = diff; print(f'         {resp:13s} baseline (no scissor): {diff:.2e}'); continue
            check(f'{resp:13s} {mode}: excitonic = single particle as at Delta = 0', diff < 1.5*base + 1e-4,
                  f'{diff:.2e} vs baseline {base:.2e}')
    _, spr = load(os.path.join(work, 'ni_shift_r_unchanged'), 'shift_sp_')
    _, spv = load(os.path.join(work, 'ni_shift_v_unchanged'), 'shift_sp_')
    rel = abs(spr - spv).max()/abs(spr).max()
    check('Scissor_mode changes the answer (shift, r_unchanged vs v_unchanged)', rel > 0.05, f'{rel:.2f}')

    print('\nThe Xatu archive records the scissor\n')
    conv = os.path.join(root, 'tools', 'xatu_text_to_h5.py')
    try:
        import h5py
    except ImportError:
        h5py = None
    if h5py is None:
        print('  SKIP   archive checks (no h5py)')
    else:
        pre = ni(DELTA); arch = os.path.join(work, 'hBN_NI30_scissor.h5')
        r = subprocess.run([sys.executable, conv, pre + '.eigval', pre + '.states', '--nfermi', '1', '--ndim', '2',
                            '-o', arch], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        d, rc, log = run('h5_norecord', 'absorbance', delta=DELTA, exc=arch)
        if 'was built without HDF5' in log:
            print('  SKIP   archive checks (opticx built without HDF5)')
        else:
            check('archive without a scissor record: runs, with a NOTE', rc == 0 and 'NOTE: the archive does not record a scissor' in log,
                  f'rc {rc}')
            for val, ok_expected in ((DELTA, True), (DELTA + 0.25, False)):
                a2 = os.path.join(work, f'hBN_NI30_s{val}.h5'); shutil.copy(arch, a2)
                with h5py.File(a2, 'r+') as f:
                    f.attrs['scissor'] = np.float64(val)
                d, rc, log = run(f'h5_s{val}', 'absorbance', delta=DELTA, exc=a2)
                if ok_expected:
                    check('archive scissor = Scissor: runs and reports it', rc == 0 and 'as in the Xatu archive' in log, f'rc {rc}')
                else:
                    check('archive scissor /= Scissor: stops before computing, naming both values',
                          rc != 0 and 'used scissor =' in log and 'Entering ome' not in log, f'rc {rc}')

    print('\nThe excitonic cache records the scissor\n')
    cache = 'ome_second_ex_hBN.omeex2'
    def version(d):
        p = os.path.join(d, cache)
        return struct.unpack_from('<i', open(p, 'rb').read(), 14)[0] if os.path.exists(p) else None
    d, rc, _ = run('cache', 'shift', extra='# Cache_ome_ex\nwrite\n')
    check('Scissor = 0 writes the version-3 cache layout', rc == 0 and version(d) == 3, f'version {version(d)}')
    d, rc, log = run('cache', 'shift', delta=DELTA, extra='# Cache_ome_ex\nread\n', keep=True)
    check('version-3 cache read with Scissor /= 0: stops', rc != 0 and 'DIFFERENT scissor' in log, f'rc {rc}')
    d, rc, _ = run('cache', 'shift', delta=DELTA, mode='r_unchanged', extra='# Cache_ome_ex\nwrite\n')
    ref = files(d)
    check('Scissor /= 0 writes version 4', rc == 0 and version(d) == 4, f'version {version(d)}')
    d, rc, log = run('cache', 'shift', delta=DELTA, mode='r_unchanged', extra='# Cache_ome_ex\nread\n', keep=True)
    check('same scissor and mode: cache hit, output identical (1e-12)', rc == 0 and 'from cache' in log and same(files(d), ref),
          f'rc {rc}')
    d, rc, log = run('cache', 'shift', delta=DELTA, mode='v_unchanged', extra='# Cache_ome_ex\nread\n', keep=True)
    check('different mode: stops', rc != 0 and 'DIFFERENT scissor' in log, f'rc {rc}')
    d, rc, log = run('cache', 'shift', delta=DELTA + 0.5, mode='r_unchanged', extra='# Cache_ome_ex\nread\n', keep=True)
    check('different scissor: stops', rc != 0 and 'DIFFERENT scissor' in log, f'rc {rc}')

    print('\nGuards\n')
    d, rc, log = run('perband', 'shift', delta=DELTA, extra='# Sp_method\nper_band\n')
    check('Sp_method = per_band with a Scissor: stops with a message', rc != 0 and 'covariant single-particle routes only' in log, f'rc {rc}')
    d, rc, log = run('badmode', 'shift', delta=DELTA, mode='both')
    check('unknown Scissor_mode: stops', rc != 0 and 'Scissor_mode' in log and 'not recognised' in log, f'rc {rc}')

    if nfail:
        print(f'\nFAILED: {nfail} check(s)'); sys.exit(1)
    print('\nALL CHECKS PASSED')


if __name__ == '__main__':
    main()
