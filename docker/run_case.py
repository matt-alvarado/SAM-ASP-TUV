#!/usr/bin/env python3
"""Container entry point: build one SAM-ASP-TUV fire case and run it.

Every setting can be given as an option or as an environment variable
(docker run -e NAME=value ...); options win.

  --lat / LAT              fire latitude (deg N)                     required
  --lon / LON              fire longitude (deg E)                    required
  --date / DATE            start date, UTC (YYYY-MM-DD)              required
  --start / START_TIME     start time, UTC (HH:MM)                   required
  --hours / RUN_HOURS      run length (h)                            required
  --decks / ASP_DECK_DIR   directory with the ASP input decks (mount it,
                           e.g. -v /host/decks:/decks -e ASP_DECK_DIR=/decks);
                           default: the decks built into the image
  --massemis / MASSEMIS    fire mass emission flux (prm massemis), default 4.86e-6
  --fire-start-step / FIRE_START_STEP
                           step after which fire emissions start (spin-up);
                           default: the template's (1080 steps)
  --dt / DT                SAM time step (s), default 8
  --case / CASE            case name, default "fire"
  --output / OUTPUT_DIR    where results go, default /output (mount it)
  --set NAME=VALUE / SAM_SET="NAME=VALUE;NAME=VALUE"
                           extra SAM &PARAMETERS settings for the case's prm
                           (Fortran syntax, e.g. perturb_type=2); --set may repeat
  --check-decks            only check the decks against the model's tracer
                           count (ntracers) and exit
  --era5-pl / ERA5_PL_FILE, --era5-sl / ERA5_SL_FILE
                           use these ERA5 netCDF files (pressure levels,
                           single levels) instead of downloading
  ERA5 download: set CDSAPI_KEY (your Copernicus CDS key) or mount a
  .cdsapirc at /root/.cdsapirc.

Before any work the decks' SAM tracer count is checked against the model's
ntracers; a mismatch stops the run with the NTRACERS value to build with.

Outputs in OUTPUT_DIR: OUT_STAT, OUT_3D, OUT_2D, OUT_MOMENTS, OUT_MOVIES,
RESTART (SAM), ASP/ (ASP output files), case/ (the case's input files,
including met_source.txt), ERA5_cache/, run.out (model log).
"""
import argparse
import os
import resource
import shutil
import subprocess
import sys

MODEL = os.environ.get('SAMASP_HOME', '/opt/sam-asp-tuv')
EXE = 'SAM_RAD_CAM_MICRO_M2005_w_aero'
DECK_FILES = ('ASPInputDeck.in', 'GasPhaseChems.in', 'GasChemicalMechanism.in',
              'PhotoRates.in', 'PhotoMap.in', 'AerosolModes.in', 'EnvAerosolModes.in')
SAM_OUT = ('OUT_2D', 'OUT_3D', 'OUT_MOMENTS', 'OUT_STAT', 'OUT_MOVIES', 'RESTART')


def env(name, default=None):
    v = os.environ.get(name)
    return v if v not in (None, '') else default


def fail(msg):
    sys.stderr.write('ERROR: %s\n' % msg)
    sys.exit(2)


def _deck_rows(path, skip_flag_line=False):
    """data rows (';'-separated fields, comments removed) of an ASP deck file"""
    rows = []
    for ln in open(path, errors='replace'):
        s = ln.replace('\r', '').split('!')[0].strip()
        if s and ';' in s and not s.startswith('END-OF-FILE'):
            rows.append([f.strip() for f in s.split(';')])
    return rows[1:] if skip_flag_line else rows


def deck_tracers(decks):
    """SAM tracers the decks need, counted as samasp_layout_init does with ASP's species counts:
    evolving gases + H2O, then per size bin the aqueous species and ions + water, the organic
    species + 1, and the aqueous organics, then one number tracer per bin.
    -> (needed, description), or (None, reason) if it cannot be counted here"""
    j = lambda f: os.path.join(decks, f)   # noqa: E731
    gas = sum(1 for r in _deck_rows(j('GasPhaseChems.in'), True) if len(r) > 1 and r[1] == '1') + 1
    aq = len(_deck_rows(j('AqPhaseChems.in'))) + len(_deck_rows(j('AqPhaseIons.in'))) + 1
    org = len(_deck_rows(j('OrgPhaseChems.in'))) + 1
    aqorg = len(_deck_rows(j('HydrophilicOrgChems.in')))
    vals = [ln.replace('\r', '').split('!')[0].strip() for ln in open(j('AerosolModes.in'), errors='replace')]
    vals = [v for v in vals if v]
    if vals[0].upper() not in ('F', '.FALSE.'):
        return None, 'AerosolModes.in allows Lagrangian particles; the model checks the count at start-up'
    nbins = int(float(vals[3]))          # flag, mass flag, upper cutoff, number of bins
    spec = aq + org + aqorg
    return gas + nbins * (spec + 1), ('%d gases + %d size bins x (%d species + 1 number)'
                                      % (gas, nbins, spec))


def compiled_ntracers():
    import re
    for ln in open(os.path.join(MODEL, 'SRC_TUV', 'domain.f90')):
        m = re.match(r'^\s*integer, parameter :: ntracers\s*=\s*(\d+)', ln)
        if m:
            return int(m.group(1))
    return None


def check_tracers(decks, quiet=False):
    """stop before any work if the decks need a different tracer count than the model was built with"""
    try:
        need, how = deck_tracers(decks)
    except (OSError, ValueError, IndexError) as e:
        print('WARNING: could not count the tracers of %s (%s); the model checks at start-up' % (decks, e))
        return
    have = compiled_ntracers()
    if need is None or have is None:
        print('WARNING: tracer count not checked here (%s); the model checks at start-up'
              % (how if need is None else 'ntracers not found in SRC_TUV/domain.f90'))
        return
    if need != have:
        fail('the ASP decks in %s need %d SAM tracers (%s), but this model was built with ntracers = %d.\n'
             '       Rebuild the image with --build-arg NTRACERS=%d (or use decks with the same species and '
             'size bins). Nothing was run.' % (decks, need, how, have, need))
    if not quiet:
        print('tracers: %d (%s) = ntracers of this build' % (need, how))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--lat', type=float, default=env('LAT'))
    ap.add_argument('--lon', type=float, default=env('LON'))
    ap.add_argument('--date', default=env('DATE'))
    ap.add_argument('--start', default=env('START_TIME'))
    ap.add_argument('--hours', type=float, default=env('RUN_HOURS'))
    ap.add_argument('--decks', default=env('ASP_DECK_DIR', os.path.join(MODEL, 'SRC_TUV/ASP/InorganicInputDecks')))
    ap.add_argument('--massemis', type=float, default=float(env('MASSEMIS', '4.86e-6')))
    ap.add_argument('--fire-start-step', type=int, default=env('FIRE_START_STEP'))
    ap.add_argument('--dt', type=float, default=float(env('DT', '8')))
    ap.add_argument('--case', default=env('CASE', 'fire'))
    ap.add_argument('--output', default=env('OUTPUT_DIR', '/output'))
    ap.add_argument('--era5-pl', default=env('ERA5_PL_FILE'))
    ap.add_argument('--era5-sl', default=env('ERA5_SL_FILE'))
    ap.add_argument('--set', action='append', default=[], metavar='NAME=VALUE')
    ap.add_argument('--check-decks', action='store_true')
    a = ap.parse_args(argv)

    if a.check_decks:
        decks = os.path.abspath(a.decks)
        absent = [f for f in DECK_FILES if not os.path.isfile(os.path.join(decks, f))]
        if absent:
            fail('ASP deck directory %s lacks %s (is it mounted?)' % (decks, ', '.join(absent)))
        check_tracers(decks)
        return 0

    missing = [n for n, v in (('lat/LAT', a.lat), ('lon/LON', a.lon), ('date/DATE', a.date),
                              ('start/START_TIME', a.start), ('hours/RUN_HOURS', a.hours)) if v is None]
    if missing:
        fail('set ' + ', '.join(missing) + ' (see --help)')
    decks = os.path.abspath(a.decks)
    absent = [f for f in DECK_FILES if not os.path.isfile(os.path.join(decks, f))]
    if absent:
        fail('ASP deck directory %s lacks %s (is it mounted?)' % (decks, ', '.join(absent)))
    check_tracers(decks)
    if bool(a.era5_pl) != bool(a.era5_sl):
        fail('give both ERA5_PL_FILE and ERA5_SL_FILE, or neither')
    if not a.era5_pl and not (env('CDSAPI_KEY') or os.path.exists(os.path.expanduser('~/.cdsapirc'))):
        fail('no ERA5 source: set CDSAPI_KEY (Copernicus CDS key), mount a .cdsapirc at '
             '/root/.cdsapirc, or give ERA5_PL_FILE and ERA5_SL_FILE')

    out = os.path.abspath(a.output)
    os.makedirs(os.path.join(out, 'ASP'), exist_ok=True)
    nstop = int(round(a.hours * 3600.0 / a.dt))
    with open(os.path.join(MODEL, 'NRANKS')) as f:
        nranks = int(f.read())

    # 1. the case (meteorology from ERA5, prm from the template)
    sets = ["asp_deck_dir='%s/'" % decks.rstrip('/'),
            "asp_output_dir='%s/'" % os.path.join(out, 'ASP')]
    if nstop < 100:                                    # SAM: nstat must not exceed nstop
        sets.append('nstat=%d' % nstop)
    if a.fire_start_step is not None:
        sets.append('fire_nstep_start=%d' % a.fire_start_step)
    sets += [s.strip() for s in (env('SAM_SET') or '').split(';') if s.strip()] + a.set
    cmd = [sys.executable, os.path.join(MODEL, 'create_run.py'), '--case', a.case,
           '--lat', str(a.lat), '--lon', str(a.lon), '--time', '%sT%s' % (a.date, a.start),
           '--massemis', repr(a.massemis), '--dt', repr(a.dt), '--hours', repr(a.hours),
           '--rundir', MODEL, '--cache', os.path.join(out, 'ERA5_cache')]
    if a.era5_pl:
        cmd += ['--pl-file', a.era5_pl, '--sl-file', a.era5_sl]
    for s in sets:
        cmd += ['--set', s]
    subprocess.run(cmd, check=True)
    casedir = os.path.join(MODEL, 'RUNS', a.case)
    shutil.copytree(casedir, os.path.join(out, 'case'), dirs_exist_ok=True)

    # fire timing (SAM: emissions on steps fire_nstep_start+1 .. + 2*(plume_y-dy)/vg/dt)
    with open(os.path.join(casedir, 'prm')) as f:
        prm = f.read()

    def prm_value(name, default):
        for l in prm.splitlines():
            l = l.split('!')[0]
            if '=' in l and l.split('=')[0].strip().lower() == name:
                return float(l.split('=')[1].strip(' ,').replace('d', 'e'))
        return default
    fstart = int(prm_value('fire_nstep_start', 1080))     # params.f90 default
    vg, plume_y, dy = prm_value('vg', 0.0), prm_value('plume_y', 0.0), prm_value('dy', 0.0)
    fend = fstart + (2 * (plume_y - dy) / vg / a.dt if vg > 0 else float('inf'))
    print('run: %d steps of %g s; fire emissions from step %d to ~%.0f; %d MPI ranks'
          % (nstop, a.dt, fstart, fend, nranks))
    if nstop <= fstart:
        print('WARNING: the run ends before the fire starts (no fire emissions); '
              'increase RUN_HOURS or set FIRE_START_STEP')
    elif nstop <= fend:
        print('WARNING: the run ends before the fire emissions end; ASP chemistry starts after them')

    # 2. run SAM in the model directory, outputs to OUTPUT_DIR
    os.chdir(MODEL)
    with open('CaseName', 'w') as f:
        f.write(a.case + '\n')
    for d in SAM_OUT:
        os.makedirs(os.path.join(out, d), exist_ok=True)
        if os.path.islink(d) or os.path.isfile(d):
            os.remove(d)
        elif os.path.isdir(d):
            shutil.rmtree(d)
        os.symlink(os.path.join(out, d), d)
    try:
        resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))
    except (ValueError, OSError):
        print('WARNING: could not lift the stack limit; run the container with --ulimit stack=-1')
    # Ubuntu names it mpirun.openmpi; elsewhere (e.g. the RHEL 9 image) it is mpirun
    mpirun = shutil.which('mpirun.openmpi') or shutil.which('mpirun') or 'mpirun'
    with open(os.path.join(out, 'run.out'), 'w') as log:
        rc = subprocess.run([mpirun, '--allow-run-as-root', '--oversubscribe',
                             '-np', str(nranks), './' + EXE], stdout=log, stderr=subprocess.STDOUT).returncode
    last = [l for l in open(os.path.join(out, 'run.out'), errors='replace') if 'NSTEP =' in l]
    step = int(last[-1].split()[2]) if last else 0
    print('SAM finished at step %d of %d (exit code %d); log: %s' % (step, nstop, rc, os.path.join(out, 'run.out')))
    return 0 if step >= nstop else (rc or 1)


if __name__ == '__main__':
    sys.exit(main())
