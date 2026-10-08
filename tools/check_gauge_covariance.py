#!/usr/bin/env python3
"""Gauge covariance of the second-order matrix elements and responses.

Physical results must not depend on the phases of the Bloch eigenvectors. Every gauge defect this project
has found survived for a while because NO test varied the gauge: a test that does not vary the gauge cannot
detect one. This one does: with the environment variable OPTICX_GAUGE_TEST set, opticx
multiplies every eigenvector by a smooth but rapidly varying phase exp(i theta_n(k)) -- about 1 rad between
neighbours of the 30x30 hBN mesh -- after its own phase convention, and carries the change into the exciton
envelopes. Each run is repeated with and without it, on hBN_N30, and the outputs compared:

  * covariant X_nm (Xnm_derivative = covariant, the default): X_nm, the excitonic shift and SHG, and the
    single-particle shift and SHG must be unchanged to round-off.
  * finite_difference X_nm: X_nm MUST change by O(1). This is the sensitivity check -- it proves the test
    can see a gauge defect, so a pass above is not a test that cannot fail.
  * the single-particle shift used to change by 9.1% at ANY scramble amplitude, down to 0.015 rad: its phase
    derivative was a difference of two principal values, which jumps by 2 pi across the negative real axis
    (fixed 2026-10-04). The amplitude-independence is checked too, because it is what identified that defect.

Degenerate-block basis, on a model WITH exactly degenerate bands -- hBN has none -- namely
monolayer MoS2 on a 15x15 mesh (MoS2_N15.eigval/.states, from MoS2_spin_wannier_07032024_tb.dat): with
OPTICX_DEGEN_SCRAMBLE set, the envelopes inside every degenerate block are rotated by a k-dependent unitary,
simulating Xatu having picked a different basis there. With Exciton_basis_repair = true the excitonic SHG
must come back unchanged; with it off it must change by O(1) (sensitivity). Skipped, loudly, if the MoS2_N15
files are not present.
"""
import argparse, os, shutil, subprocess, sys
import numpy as np

IN = """# Periodic dimensions
2
# Wannier90_filename
{tb}
# Xatu_interface
true
{eig}
{sta}
# Exciton_cutoff
100
# Nfermi
1
# OME_sp
nonlinear
# OME_ex
nonlinear
# Response
{resp}
# Energy_variables
4.6 7.6 0.05 120
# Xnm_derivative
{method}
# Cache_ome_ex
write
"""
CACHE = 'ome_second_ex_hBN.omeex2'


def read_xnm(path):
    """X_nm from a format-3 .omeex2 cache: (3, N, N) complex."""
    f = open(path, 'rb')
    r = lambda dt, n: np.fromfile(f, dtype=dt, count=n)
    f.read(14); ver = r('<i4', 1)[0]
    assert ver == 3, f'cache format {ver}, expected 3'
    nlen = r('<i4', 1)[0]; f.read(nlen)
    npt, nv, nc, N = r('<i4', 4); nb = r('<i4', 1)[0]; r('<i4', nb); r('<i4', 1); r('<f8', N)
    r('<c16', 3 * N); r('<c16', 3 * N)
    X = np.empty((3, N, N), complex)
    for m in range(N):
        X[:, :, m] = r('<c16', 3 * N).reshape(N, 3).T
    return X


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--opticx', required=True)
    ap.add_argument('--root', default='.')
    ap.add_argument('--workdir', default='check_gauge_covariance')
    ap.add_argument('--degen-root', default=None, help='directory holding MoS2_N15.eigval/.states (default: --root)')
    a = ap.parse_args()
    root, work, opticx = os.path.abspath(a.root), os.path.abspath(a.workdir), os.path.abspath(a.opticx)
    tb, eig, sta = (os.path.join(root, f) for f in ('hBN_tb.dat', 'hBN_N30.eigval', 'hBN_N30.states'))
    for f in (tb, eig, sta):
        if not os.path.exists(f):
            sys.exit(f'missing {f}')
    shutil.rmtree(work, ignore_errors=True); os.makedirs(work)

    nfail = 0

    def check(label, ok, msg=''):
        nonlocal nfail
        print(('  PASS   ' if ok else '  FAIL   ') + label + ('   ' + msg if msg else ''))
        nfail += 0 if ok else 1

    def run(method, resp, gauge):
        d = os.path.join(work, f'{method}_{resp}_g{gauge}')
        os.makedirs(d)
        open(os.path.join(d, 'in.txt'), 'w').write(IN.format(tb=tb, eig=eig, sta=sta, resp=resp, method=method))
        env = dict(os.environ)
        env.pop('OPTICX_GAUGE_TEST', None)
        if gauge != '0':
            env['OPTICX_GAUGE_TEST'] = gauge
        r = subprocess.run([opticx, 'in.txt'], cwd=d, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log = r.stdout.decode('utf-8', 'replace')
        open(os.path.join(d, 'run.log'), 'w').write(log)
        for f in os.listdir(d):
            if f.endswith('.omesp'):
                os.remove(os.path.join(d, f))
        return r.returncode, log, d

    rel = lambda x, y: abs(x - y).max() / max(abs(x).max(), 1e-300)
    load = lambda d, f: np.loadtxt(os.path.join(d, f))[:, 1:]

    print('Gauge covariance under scrambled eigenvector phases (hBN_N30, 100 excitons)\n')
    out = {}
    for method in ('covariant', 'finite_difference'):
        for resp in ('shift_shiftvector', 'shg'):
            for g in ('0', '1'):
                rc, log, d = run(method, resp, g)
                out[method, resp, g] = d
                if rc != 0:
                    check(f'{method} {resp} gauge={g} ran', False, f'rc={rc}')
                if g == '1':
                    check(f'{method:17s} {resp:17s} log announces the gauge test', 'gauge test' in log.lower())

    c0, c1 = out['covariant', 'shift_shiftvector', '0'], out['covariant', 'shift_shiftvector', '1']
    f0, f1 = out['finite_difference', 'shift_shiftvector', '0'], out['finite_difference', 'shift_shiftvector', '1']
    e = rel(read_xnm(os.path.join(c0, CACHE)), read_xnm(os.path.join(c1, CACHE)))
    check('covariant          X_nm unchanged', e < 1e-12, f'{e:.2e} (tol 1e-12)')
    e = rel(load(c0, 'shift_ex_lengthgauge_hBN.dat'), load(c1, 'shift_ex_lengthgauge_hBN.dat'))
    check('covariant          excitonic shift unchanged', e < 1e-12, f'{e:.2e} (tol 1e-12)')
    s0, s1 = out['covariant', 'shg', '0'], out['covariant', 'shg', '1']
    e = rel(load(s0, 'shg_ex_lengthgauge_hBN.dat'), load(s1, 'shg_ex_lengthgauge_hBN.dat'))
    check('covariant          excitonic SHG unchanged', e < 1e-10, f'{e:.2e} (tol 1e-10)')
    e = rel(load(s0, 'shg_sp_lengthgauge_hBN.dat'), load(s1, 'shg_sp_lengthgauge_hBN.dat'))
    check('single-particle    SHG unchanged', e < 1e-8, f'{e:.2e} (tol 1e-8)')
    e = rel(load(c0, 'shift_sp_lengthgauge_hBN.dat'), load(c1, 'shift_sp_lengthgauge_hBN.dat'))
    check('single-particle    shift unchanged (branch-cut fix)', e < 1e-8,
          f'{e:.2e} (tol 1e-8; 9.1e-02 before the 2026-10-04 fix)')

    # sensitivity: the finite_difference form is NOT covariant, so the scramble must move it
    e = rel(read_xnm(os.path.join(f0, CACHE)), read_xnm(os.path.join(f1, CACHE)))
    check('SENSITIVITY        finite_difference X_nm DOES change', e > 1e-2,
          f'{e:.2e} (must exceed 1e-2, else the test cannot see a gauge defect)')

    # amplitude independence of the single-particle shift: a branch-cut defect shows up at ANY amplitude
    rc, log, d = run('covariant', 'shift_shiftvector', '0.01')
    e = rel(load(c0, 'shift_sp_lengthgauge_hBN.dat'), load(d, 'shift_sp_lengthgauge_hBN.dat'))
    check('single-particle    shift unchanged at 0.015 rad amplitude', e < 1e-8, f'{e:.2e} (tol 1e-8)')

    # -- degenerate-block basis repair, needs a model with exact degeneracies ----------------
    droot = os.path.abspath(a.degen_root or root)
    mtb = os.path.join(root, 'MoS2_spin_wannier_07032024_tb.dat')
    meig, msta = os.path.join(droot, 'MoS2_N15.eigval'), os.path.join(droot, 'MoS2_N15.states')
    if not all(os.path.exists(f) for f in (mtb, meig, msta)):
        print('\n  SKIP   degenerate-block checks: MoS2_N15.eigval/.states or the MoS2 tb file not found'
              f' (looked in {droot}). Generate with xatu --w90 26 <tb> exciton.config -ec -n 899 at ncells 15.')
    else:
        def mrun(tag, repair, scramble):
            d = os.path.join(work, f'mos2_{tag}'); os.makedirs(d)
            txt = IN.format(tb=mtb, eig=meig, sta=msta, resp='shg', method='covariant')
            txt = txt.replace('# Exciton_cutoff\n100', '# Exciton_cutoff\n400').replace('# Nfermi\n1', '# Nfermi\n26')
            txt = txt.replace('4.6 7.6 0.05 120', '0.8 3.2 0.025 120').replace('# Cache_ome_ex\nwrite', '# Cache_ome_ex\noff')
            open(os.path.join(d, 'in.txt'), 'w').write(txt + f'# Exciton_basis_repair\n{repair}\n')
            env = dict(os.environ); env.pop('OPTICX_GAUGE_TEST', None); env.pop('OPTICX_DEGEN_SCRAMBLE', None)
            if scramble:
                env['OPTICX_DEGEN_SCRAMBLE'] = '1'
            r = subprocess.run([opticx, 'in.txt'], cwd=d, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            log = r.stdout.decode('utf-8', 'replace'); open(os.path.join(d, 'run.log'), 'w').write(log)
            for f in os.listdir(d):
                if f.endswith('.omesp'):
                    os.remove(os.path.join(d, f))
            return r.returncode, log, d
        out2 = {}
        for repair in ('true', 'false'):
            for scr in (False, True):
                rc, log, d = mrun(f'{repair}_{int(scr)}', repair, scr)
                out2[repair, scr] = d
                if rc != 0:
                    check(f'MoS2_N15 repair={repair} scramble={scr} ran', False, f'rc={rc}')
                if scr:
                    check(f'MoS2_N15 repair={repair:5s} log announces the scramble', 'degenerate-block scramble' in log.lower())
                if repair == 'true' and not scr:
                    check('MoS2_N15 repair=true  log reports the repaired k-points', 'basis repair:' in log.lower()
                          and 'exactly degenerate' in log.lower())
        f = 'shg_ex_lengthgauge_MoS2_spin_wannier_07032024.dat'
        e = rel(load(out2['true', False], f), load(out2['true', True], f))
        check('MoS2_N15 repair=true  excitonic SHG unchanged by a different Xatu basis', e < 1e-9, f'{e:.2e} (tol 1e-9)')
        e = rel(load(out2['false', False], f), load(out2['false', True], f))
        check('SENSITIVITY  repair=false excitonic SHG DOES change', e > 1e-2,
              f'{e:.2e} (must exceed 1e-2, else the scramble does not reach the degenerate blocks)')

    if nfail:
        print(f'\nFAILED: {nfail} check(s)'); sys.exit(1)
    print('\nALL CHECKS PASSED')


if __name__ == '__main__':
    main()
