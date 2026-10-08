#!/usr/bin/env python3
"""Behaviour AND observability of the second-order excitonic OME cache (Cache_ome_ex).

Driven through Response = shift_shiftvector, which exercises the same excitonic OME path and exists in
every version of the code, so the test is independent of which second-order responses are built.

Exercises the five modes on hBN_N30 and asserts, for each, both what the run DID (was the exciton
k-loop run or skipped? was the cache written?) and what it SAID about it.

Also the OME_ex = none route (OME_sp = none's excitonic counterpart): with Cache_ome_ex = read|readwrite a
second-order run takes the excitonic matrix elements from the cache, and since every OME_sp = nonlinear run stores
all the block-covariant data, the .omesp written by the shift_shiftvector runs here serves Response = shg too.
That run must hit the cache, skip the k-loop and the envelope read, and equal a direct shg run; OME_ex = none
without a cache, or with a cache file that is absent, must stop before computing anything.

The log assertions are the point. An earlier ad-hoc version of these checks asserted only behaviour --
k-point line counts, file existence, mtimes -- and all five passed for days while the cache's "read from
cache" message sat behind an unreachable `return`. A cache that cannot be seen in a log
cannot be debugged from one, so silence is a failure here, not a cosmetic issue.
"""
import argparse
import struct, os, re, shutil, subprocess, sys, time
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
{ncut}
# Nfermi
1
# OME_sp
nonlinear
# OME_ex
nonlinear
# Response
shift_shiftvector
# Energy_variables
2.0 8.0 0.05 60
# Cache_ome_ex
{mode}
"""
CACHE = 'ome_second_ex_hBN.omeex2'
OMESP = 'ome_nonlinear_sp_hBN.omesp'
OMESP_TAG = 1330464562
# a second-order run with every matrix-element switch free (OME_ex = none route)
IN2 = """# Periodic dimensions
2
# Wannier90_filename
{tb}
# Xatu_interface
true
{eig}
{sta}
# Exciton_cutoff
120
# Nfermi
1
# OME_sp
{osp}
# OME_ex
{oex}
# Response
{resp}
# Energy_variables
2.0 8.0 0.05 60
# Cache_ome_ex
{mode}
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--opticx', required=True)
    ap.add_argument('--root', default='.')
    ap.add_argument('--workdir', default='check_ome_cache')
    a = ap.parse_args()
    root, work = os.path.abspath(a.root), os.path.abspath(a.workdir)
    opticx = os.path.abspath(a.opticx)
    tb = os.path.join(root, 'hBN_tb.dat')
    eig = os.path.join(root, 'hBN_N30.eigval')
    sta = os.path.join(root, 'hBN_N30.states')
    for f in (tb, eig, sta):
        if not os.path.exists(f):
            print('SKIP: missing ' + f); sys.exit(0)
    shutil.rmtree(work, ignore_errors=True); os.makedirs(work)

    nfail = 0

    def check(label, ok, msg=''):
        nonlocal nfail
        print(('  PASS   ' if ok else '  FAIL   ') + label + ('   ' + msg if msg else ''))
        nfail += 0 if ok else 1

    def hdr_offsets(blob):
        """Locate the header fields. Layout: MAGIC(14) VERSION(4) | nlen(4) name | np nv nc nex(16)
        | nband_ex(4) | nband_index(4*nband_ex) [v2+] | xnm_method(4) [v3+].
        Returns (version, off_nband, nband, method); fields absent from that version are None."""
        off = 14
        ver = struct.unpack_from('<i', blob, off)[0]; off += 4
        nlen = struct.unpack_from('<i', blob, off)[0]; off += 4 + nlen
        off += 16
        if ver < 2:
            return ver, off, None, None
        nband = struct.unpack_from('<i', blob, off)[0]
        method = struct.unpack_from('<i', blob, off + 4 + 4 * nband)[0] if ver >= 3 else None
        return ver, off, nband, method

    def mangle_cache(path, what):
        """'swap' reverses the cached band list (nv/nc unchanged: two bands stored in the wrong order);
        'v1' rewrites the file in the old format that carried no band list at all;
        'v2' rewrites it in the format that carried the band list but not the X_nm method;
        'method' flips the recorded X_nm method (2 = covariant <-> 1 = finite_difference)."""
        b = bytearray(open(path, 'rb').read())
        ver, off, nband, method = hdr_offsets(b)
        assert ver == 3, f'expected cache format version 3, got {ver}'
        idx = list(struct.unpack_from(f'<{nband}i', b, off + 4))
        moff = off + 4 + 4 * nband                      # offset of the v3 method field
        if what == 'swap':
            struct.pack_into(f'<{nband}i', b, off + 4, *idx[::-1])
            out = bytes(b)
        elif what == 'method':
            struct.pack_into('<i', b, moff, 1 if method == 2 else 2)
            out = bytes(b)
        elif what == 'v2':
            out = bytes(b[:14] + struct.pack('<i', 2) + b[18:moff] + b[moff + 4:])
        else:
            out = bytes(b[:14] + struct.pack('<i', 1) + b[18:off] + b[moff + 4:])
        open(path, 'wb').write(out)
        return idx

    def run(tag, mode, ncut=120, cache_from=None, truncate=None, mangle=None):
        """One opticx run in its own directory. Returns (returncode, log, dirpath)."""
        d = os.path.join(work, tag); os.makedirs(d, exist_ok=True)
        open(os.path.join(d, 'in.txt'), 'w').write(
            IN.format(tb=tb, eig=eig, sta=sta, ncut=ncut, mode=mode))
        if cache_from:
            shutil.copy2(os.path.join(work, cache_from, CACHE), os.path.join(d, CACHE))
            if mangle is not None:
                mangle_cache(os.path.join(d, CACHE), mangle)
            if truncate is not None:
                p = os.path.join(d, CACHE)
                with open(p, 'r+b') as fh:
                    fh.truncate(max(1, os.path.getsize(p) // truncate))
        r = subprocess.run([opticx, 'in.txt'], cwd=d, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT)
        log = r.stdout.decode('utf-8', 'replace')
        open(os.path.join(d, 'run.log'), 'w').write(log)
        return r.returncode, log, d

    kloop = lambda log: len(re.findall(r'OME \(ex\): k-point', log))
    said_read = lambda log: 'read from cache' in log.lower()
    said_wrote = lambda log: 'cached to' in log.lower()
    # The exciton envelopes are the second most expensive part of start-up and a cache HIT never uses
    # them, so they are not read until the cache has been consulted. Asserted through
    # the log for the same reason the read/write messages are: a saving that cannot be seen in a log
    # cannot be shown to still be happening.
    read_env = lambda log: 'reading exciton wavefunctions' in log.lower()
    deferred = lambda log: 'not read yet' in log.lower()

    print('Second-order excitonic OME cache: behaviour and log observability\n')

    # -- write: must compute AND write AND say it wrote -------------------------------------------
    rc, log, d = run('write', 'write')
    check('write      computed (k-loop ran)', rc == 0 and kloop(log) > 0, f'{kloop(log)} k-points')
    check('write      cache file created', os.path.exists(os.path.join(d, CACHE)))
    check('write      log announces the write', said_wrote(log))
    check('write      log does NOT claim a read', not said_read(log))
    check('write      envelopes read (no hit was possible)', read_env(log))
    ref = np.loadtxt(os.path.join(d, 'shift_ex_lengthgauge_hBN.dat'))

    # -- read with a cache present: must skip, not rewrite, and SAY it read ------------------------
    t0 = os.path.getmtime(os.path.join(work, 'write', CACHE))
    rc, log, d = run('read', 'read', cache_from='write')
    check('read       k-loop skipped', rc == 0 and kloop(log) == 0, f'{kloop(log)} k-points')
    check('read       log announces the read', said_read(log),
          '' if said_read(log) else 'cache used but silent in the log')
    # the copy carries the source mtime, so 'read' must leave it untouched
    check('read       cache not rewritten',
          abs(os.path.getmtime(os.path.join(d, CACHE)) - t0) < 1e-6,
          f'mtime delta {os.path.getmtime(os.path.join(d, CACHE)) - t0:.3f} s')
    got = np.loadtxt(os.path.join(d, 'shift_ex_lengthgauge_hBN.dat'))
    rel = abs(got - ref).max() / max(abs(ref).max(), 1e-300)
    check('read       result matches the computed one', rel < 1e-12, f'max rel diff {rel:.2e}')
    check('read       envelopes NOT read on a hit', not read_env(log),
          '' if not read_env(log) else 'the .states wavefunctions were read for nothing')
    check('read       log says the read was deferred', deferred(log))

    # -- read with NO cache: clean miss, silent about reading, and it recomputes --------------------
    rc, log, d = run('read_nocache', 'read')
    check('read       absent cache -> recomputes', rc == 0 and kloop(log) > 0, f'{kloop(log)} k-points')
    check('read       absent cache -> no read claimed', not said_read(log))
    check('read       absent cache -> envelopes loaded lazily', read_env(log),
          '' if read_env(log) else 'the k-loop would run on unread envelopes')

    # -- readwrite from cold: computes then writes --------------------------------------------------
    rc, log, d = run('readwrite', 'readwrite')
    check('readwrite  cold: computed then wrote',
          rc == 0 and kloop(log) > 0 and os.path.exists(os.path.join(d, CACHE)) and said_wrote(log))

    # -- off: cache present but must be ignored ------------------------------------------------------
    rc, log, d = run('off', 'off', cache_from='write')
    check('off        cache present but ignored', rc == 0 and kloop(log) > 0, f'{kloop(log)} k-points')
    check('off        log does not mention a read', not said_read(log))

    # -- a truncated cache is a clean miss, and says so ----------------------------------------------
    rc, log, d = run('truncated', 'read', cache_from='write', truncate=3)
    check('truncated  falls back to computing', rc == 0 and kloop(log) > 0)
    check('truncated  log explains the fallback',
          'truncated' in log.lower() or not said_read(log))

    # -- asking for more excitons than were cached is a miss, not a silent wrong answer --------------
    rc, log, d = run('bigger', 'read', ncut=160, cache_from='write')
    check('bigger     larger cutoff -> recomputes', rc == 0 and kloop(log) > 0)
    check('bigger     larger cutoff -> no read claimed', not said_read(log))

    # -- an unrecognised value must stop with the valid list -------------------------------------------
    rc, log, d = run('bogus', 'sometimes')
    check('bogus      refused', rc != 0 and 'not recognised' in log.lower())

    # -- asking for more excitons than the .eigval holds must say so, not die on EOF -------------------
    # hBN_N30 has 900 states. This used to end in a bare "Fortran runtime error: End of file" plus a
    # backtrace, naming neither the keyword nor the limit.
    rc, log, d = run('overrun', 'off', ncut=1000)
    named = 'exciton_cutoff' in log.lower() and '900' in log
    check('overrun    Exciton_cutoff past the end of .eigval is refused', rc != 0 and named,
          '' if named else 'stopped without naming the keyword and the limit available')
    check('overrun    no Fortran runtime backtrace', 'backtrace' not in log.lower())

    # -- the band list is in the fingerprint ------------------------------------------
    # Version 1 recorded only the COUNTS nv_ex and nc_ex. Both [60,61] and [61,60] give nv=4, nc=2, so
    # a count-based header could not see two conduction bands stored in swapped order, which once silently
    # invalidated a whole ReS2 result set (excitonic shift 2.23x too large). These checks are the reason the format was bumped.
    blob = open(os.path.join(work, 'write', CACHE), 'rb').read()
    ver, _, nband, method = hdr_offsets(blob)
    check('bandlist   cache is format version 3', ver == 3, f'version {ver}')
    check('bandlist   header records the band count', nband == 2, f'nband_ex = {nband}')
    check('method     header records the X_nm method (2 = covariant, the default)', method == 2,
          f'method = {method}')

    # a REVERSED band list with identical nv/nc must STOP -- the swapped-band-order case exactly
    rc, log, d = run('bandswap', 'read', cache_from='write', mangle='swap')
    low = log.lower()
    check('bandswap   reversed band list is REFUSED', rc != 0,
          '' if rc != 0 else 'ran with the bands in the wrong order')
    check('bandswap   message names the band set', 'band set' in low or 'band list' in low)
    check('bandswap   message flags the ORDER specifically', 'order' in low)
    check('bandswap   nothing was computed', kloop(log) == 0, f'{kloop(log)} k-points')

    # a genuine v1 file must be a CLEAN MISS -- it cannot be upgraded, its order is unknowable
    rc, log, d = run('v1cache', 'read', cache_from='write', mangle='v1')
    low = log.lower()
    check('v1cache    old format -> recomputes', rc == 0 and kloop(log) > 0, f'{kloop(log)} k-points')
    check('v1cache    old format -> no read claimed', not said_read(low))
    check('v1cache    log names the version and the hazard',
          'version 1' in low and ('band' in low))

    # -- the X_nm method is in the fingerprint ----------------------------------------
    # The covariant and finite_difference forms of X_nm differ by O(1) on multi-band spinful models
    # (MoS2: excitonic SHG C3 residual 62% vs 5%), so elements built by one must not be read by the other.
    rc, log, d = run('v2cache', 'read', cache_from='write', mangle='v2')
    low = log.lower()
    check('v2cache    format-2 file (no X_nm method) -> recomputes', rc == 0 and kloop(log) > 0, f'{kloop(log)} k-points')
    check('v2cache    format-2 file (no X_nm method) -> no read claimed', not said_read(low))
    check('v2cache    log names the version and X_nm', 'version 2' in low and 'x_nm' in low)
    rc, log, d = run('methodswap', 'read', cache_from='write', mangle='method')
    low = log.lower()
    check('method     cache built by the OTHER X_nm method is REFUSED', rc != 0,
          '' if rc != 0 else 'mixed X_nm methods silently')
    check('method     message names the X_nm method', 'x_nm method' in low)
    check('method     nothing was computed', kloop(log) == 0, f'{kloop(log)} k-points')

    # -- OME_ex = none: the excitonic elements from the cache, the single-particle ones from any .omesp ----
    def run2(tag, osp, oex, resp, mode, files=()):
        d = os.path.join(work, tag); os.makedirs(d, exist_ok=True)
        for f in files:
            shutil.copy2(os.path.join(work, 'write', f), os.path.join(d, f))
        open(os.path.join(d, 'in.txt'), 'w').write(
            IN2.format(tb=tb, eig=eig, sta=sta, osp=osp, oex=oex, resp=resp, mode=mode))
        r = subprocess.run([opticx, 'in.txt'], cwd=d, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log = r.stdout.decode('utf-8', 'replace')
        open(os.path.join(d, 'run.log'), 'w').write(log)
        return r.returncode, log, d

    blob = open(os.path.join(work, 'write', OMESP), 'rb').read()
    i = blob.find(struct.pack('<i', OMESP_TAG))
    flags = struct.unpack_from('<i', blob, i + 4)[0] if i >= 0 else -1
    check('omesp      a shift run stores every block-covariant section (flags 15)', flags == 15, f'flags {flags}')
    rc, log, dn = run2('none_shg', 'none', 'none', 'shg', 'read', files=(OMESP, CACHE))
    check('none_shg   OME_sp = none + OME_ex = none + read: runs, says it read the cache',
          rc == 0 and said_read(log), f'rc={rc}')
    check('none_shg   nothing computed and the envelopes not read',
          kloop(log) == 0 and not read_env(log), f'{kloop(log)} k-points, envelopes read: {read_env(log)}')
    rc2, _, dd = run2('direct_shg', 'nonlinear', 'nonlinear', 'shg', 'off')
    for f in ('shg_sp_lengthgauge_hBN.dat', 'shg_ex_lengthgauge_hBN.dat'):
        ok = rc == 0 and rc2 == 0 and all(os.path.exists(os.path.join(x, f)) for x in (dn, dd))
        e = (abs(np.loadtxt(os.path.join(dn, f)) - np.loadtxt(os.path.join(dd, f))).max()
             / abs(np.loadtxt(os.path.join(dd, f))).max()) if ok else 1.0
        check(f'none_shg   = direct run: {f}', ok and e < 1e-12, f'{e:.1e} (tol 1e-12)')
    rc, log, _ = run2('none_nocache', 'none', 'none', 'shg', 'off', files=(OMESP,))
    check('none_nocache OME_ex = none without Cache_ome_ex: stops before the response, says why',
          rc != 0 and 'set cache_ome_ex = read' in log.lower() and 'entering optical_response' not in log.lower(), f'rc={rc}')
    rc, log, _ = run2('none_miss', 'none', 'none', 'shg', 'read', files=(OMESP,))
    check('none_miss  OME_ex = none + read without a cache file: stops, says why',
          rc != 0 and 'ome_ex = none reads the excitonic matrix elements' in log.lower()
          and 'entering optical_response' not in log.lower(), f'rc={rc}')

    if nfail:
        print(f'\nFAILED: {nfail} check(s)'); sys.exit(1)
    print('\nALL CHECKS PASSED')


if __name__ == '__main__':
    main()
