#!/usr/bin/env python3
"""Regression check of the single-particle SHG conductivity (get_sigma_shg_sp / get_shg_intens_sp), with
tolerances, on GeS's real 27-band model (opticx's own k-grid, Xatu_interface = false).

What it does (~1 minute): runs opticx twice with Response = shg, OME_sp = nonlinear, on two band windows
around the gap (bands 11-27 and 1-27, 1-based, Nfermi = 20) and checks that sigma^{abc} agrees between the
two windows within a tolerance -- i.e. that adding the 10 Ge 3d semicore bands (1-10) barely changes the result, the
basic sanity check for a length-gauge SHG calculation (paper Fig. 7: method A should already be close to
converged with enough bands; in a 6/13/27-band window study, 6->13 bands changed the result by ~100x because of
near-degenerate bands, 13->27 by <2%). The small window was 15-27 until 2026-10-05; with
the correct filling (Nfermi 20) that window drops four s/p valence bands (11-14) and xyy moved by 5.9%,
so it is now 11-27 (all s/p valence + conduction; agreement 0.1%).

Why GeS and not hBN: this check is about band-window convergence, which needs real remote bands, and hBN's
two-orbital model has none. (hBN itself is a valid SHG target: with the gauge fixes of 2026-09-24 its two-band
SHG is gauge invariant to machine precision, satisfies D3h to 6e-6 on 30x30, 60x60 and 120x120, and equals the
independent NumPy evaluation of Eq. A3a in tools/ipa_shg_numpy.py to 5e-7.)

Usage: python3 tools/check_sp_shg.py --opticx bin/opticx --root . --workdir bin/check_sp_shg
"""
import argparse, os, subprocess, sys
import numpy as np


def bandlist(lo, hi, nf):
    return ' '.join(str(b - nf) for b in range(lo, hi + 1))


def run_shg(opticx, root, work, tag, lo, hi, nf, ncells, w_ev, eta_ev):
    tb = os.path.join(root, 'wannier90_files_input', 'GeS_wannier_04062024_tb.dat')
    inp = f"""# Periodic dimensions
2
# Wannier90_filename
{tb}
# Xatu_interface
false
# Ncells
{ncells}
# Bandlist
{bandlist(lo, hi, nf)}
# Nfermi
{nf}
# OME_sp
nonlinear
# OME_ex
false
# Response
shg
# Energy_variables
{w_ev} {w_ev} {eta_ev} 1
"""
    wdir = os.path.join(work, tag)
    os.makedirs(wdir, exist_ok=True)
    open(os.path.join(wdir, f'GeS_shg_{tag}.txt'), 'w').write(inp)
    r = subprocess.run([os.path.abspath(opticx), f'GeS_shg_{tag}.txt'], cwd=wdir,
                        stdout=open(os.path.join(wdir, 'out.log'), 'w'), stderr=subprocess.STDOUT)
    if r.returncode != 0:
        print(f'FAIL: opticx exited with status {r.returncode} for {tag} (see {wdir}/out.log)'); sys.exit(1)
    d = np.loadtxt(os.path.join(wdir, 'shg_sp_lengthgauge_GeS_wannier_04062024.dat'), comments='#')
    vals = d[1:].reshape(3, 3, 3, 2)
    return vals[..., 0] + 1j * vals[..., 1]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--opticx', required=True); ap.add_argument('--root', default='.'); ap.add_argument('--workdir', default='check_sp_shg')
    a = ap.parse_args()
    root = os.path.abspath(a.root); work = os.path.abspath(a.workdir)
    nf, ncells, w_ev, eta_ev = 20, 10, 3.0, 0.15   # 20: gap 2.76 eV between bands 20 and 21 (21 was metallic)

    print(f'sp SHG on GeS (27 bands, {ncells}x{ncells} grid, hbar*omega = {w_ev} eV, eta = {eta_ev} eV)')
    sig17 = run_shg(a.opticx, root, work, 'w17', 11, 27, nf, ncells, w_ev, eta_ev)
    sig27 = run_shg(a.opticx, root, work, 'w27', 1, 27, nf, ncells, w_ev, eta_ev)

    nfail = 0
    def check(label, ok, msg):
        nonlocal nfail
        print(('  PASS   ' if ok else '  FAIL   ') + label + '   ' + msg); nfail += 0 if ok else 1

    for name, (i, j, k) in {'xxx': (0, 0, 0), 'xyy': (0, 1, 1), 'yxy': (1, 0, 1)}.items():
        v17, v27 = sig17[i, j, k], sig27[i, j, k]
        rel = abs(v27 - v17) / abs(v27)
        check(f'convergence (17->27 bands) {name}', rel < 0.05,
              f'17 bands: {v17:.4f}, 27 bands: {v27:.4f}, rel diff {rel:.3f} (max 0.05)')

    check('finite (no NaN/Inf)', np.isfinite(sig27).all(), '')

    # Filling guard: Nfermi must put the Fermi level in a gap everywhere. 20 does (gap 2.76 eV); 21 is metallic
    # in this model (band 21 reaches 1.14 eV, band 22 comes down to -0.43 eV), and OptiX must say so.
    fill_warn = lambda log: 'does not put the Fermi level in a gap' in log
    log20 = open(os.path.join(work, 'w27', 'out.log')).read()
    check('filling guard: Nfermi = 20 (insulating) gives no warning', not fill_warn(log20), '')
    wdir = os.path.join(work, 'nf21'); os.makedirs(wdir, exist_ok=True)
    inp21 = open(os.path.join(work, 'w27', 'GeS_shg_w27.txt')).read()
    inp21 = inp21.replace('# Nfermi\n20', '# Nfermi\n21').replace('# OME_sp\nnonlinear', '# OME_sp\nlinear') \
                 .replace('# Response\nshg', '# Response\nnone').replace(bandlist(1, 27, 20), bandlist(1, 27, 21))
    open(os.path.join(wdir, 'GeS_filling_nf21.txt'), 'w').write(inp21)
    r = subprocess.run([os.path.abspath(a.opticx), 'GeS_filling_nf21.txt'], cwd=wdir, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT)
    log21 = r.stdout.decode('utf-8', 'replace'); open(os.path.join(wdir, 'out.log'), 'w').write(log21)
    check('filling guard: Nfermi = 21 (metallic) is warned about', fill_warn(log21), f'rc={r.returncode}')
    check('intrinsic permutation symmetry (sigma^abc = sigma^acb, by construction)',
          np.allclose(sig27, np.swapaxes(sig27, 1, 2), atol=1e-10), '')

    if nfail:
        print(f'\nFAILED: {nfail} check(s)'); sys.exit(1)
    print('\nALL CHECKS PASSED')


if __name__ == '__main__':
    main()
