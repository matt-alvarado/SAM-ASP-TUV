#!/usr/bin/env python3
"""Container entry point: build one SAM-ASP-TUV fire case and run it.

Every setting is an option or an environment variable (docker run -e
NAME=value); options win. The ASP decks' tracer count is checked against the
model's ntracers before any work. ERA5 is downloaded with CDSAPI_KEY or a
mounted /root/.cdsapirc, unless ERA5_PL_FILE and ERA5_SL_FILE are given.

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


def require_decks(decks):
    """absolute path of the ASP deck directory; stops if a deck file is missing"""
    decks = os.path.abspath(decks)
    absent = [f for f in DECK_FILES if not os.path.isfile(os.path.join(decks, f))]
    if absent:
        fail('ASP deck directory %s lacks %s (is it mounted?)' % (decks, ', '.join(absent)))
    return decks


def prm_value(text, name):
    """value of the uncommented entry `name = value` in prm text, or None"""
    for line in text.splitlines():
        line = line.split('!')[0]
        if '=' in line and line.split('=')[0].strip().lower() == name:
            return float(line.split('=')[1].strip(' ,').replace('d', 'e'))
    return None


def require_prm_value(text, name, where):
    v = prm_value(text, name)
    if v is None:
        fail('%s does not set %s' % (where, name))
    return v


def parse_args(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    def opt(flag, var, help, **kw):
        kw.setdefault('default', env(var))
        ap.add_argument(flag, help='%s [env %s]' % (help, var), **kw)
    opt('--lat', 'LAT', 'fire latitude (deg N); required', type=float)
    opt('--lon', 'LON', 'fire longitude (deg E); required', type=float)
    opt('--date', 'DATE', 'start date, UTC (YYYY-MM-DD); required')
    opt('--start', 'START_TIME', 'start time, UTC (HH:MM); required')
    opt('--hours', 'RUN_HOURS', 'run length (h); required', type=float)
    opt('--decks', 'ASP_DECK_DIR', 'directory with the ASP input decks (mount it, e.g. -v /host/decks:/decks '
        '-e ASP_DECK_DIR=/decks); default: the decks built into the image',
        default=env('ASP_DECK_DIR', os.path.join(MODEL, 'SRC_TUV/ASP/InorganicInputDecks')))
    opt('--massemis', 'MASSEMIS', 'fire mass emission flux (prm massemis); default %(default)s',
        type=float, default=float(env('MASSEMIS', '4.86e-6')))
    opt('--fire-start-step', 'FIRE_START_STEP', 'step after which fire emissions start (spin-up); '
        'default: the prm template value', type=int)
    opt('--dt', 'DT', 'SAM time step (s); default %(default)s', type=float, default=float(env('DT', '8')))
    opt('--case', 'CASE', 'case name; default %(default)s', default=env('CASE', 'fire'))
    opt('--output', 'OUTPUT_DIR', 'where results go (mount it); default %(default)s',
        default=env('OUTPUT_DIR', '/output'))
    opt('--era5-pl', 'ERA5_PL_FILE', 'ERA5 pressure-level netCDF file to use instead of downloading')
    opt('--era5-sl', 'ERA5_SL_FILE', 'ERA5 single-level netCDF file to use instead of downloading')
    ap.add_argument('--set', action='append', default=[], metavar='NAME=VALUE',
                    help='extra SAM &PARAMETERS setting for the case\'s prm (Fortran syntax, e.g. '
                    'perturb_type=2); repeatable [env SAM_SET="NAME=VALUE;NAME=VALUE"]')
    ap.add_argument('--check-decks', action='store_true',
                    help='only check the decks against the model\'s tracer count (ntracers) and exit')
    a = ap.parse_args(argv)
    a.set = [s.strip() for s in (env('SAM_SET') or '').split(';') if s.strip()] + a.set
    return a


def validate_inputs(a):
    """stop before any work if a setting, the decks or the ERA5 source is missing,
    or if the decks need another tracer count; -> absolute deck directory"""
    missing = [n for n, v in (('lat/LAT', a.lat), ('lon/LON', a.lon), ('date/DATE', a.date),
                              ('start/START_TIME', a.start), ('hours/RUN_HOURS', a.hours)) if v is None]
    if missing:
        fail('set ' + ', '.join(missing) + ' (see --help)')
    decks = require_decks(a.decks)
    check_tracers(decks)
    if bool(a.era5_pl) != bool(a.era5_sl):
        fail('give both ERA5_PL_FILE and ERA5_SL_FILE, or neither')
    if not a.era5_pl and not (env('CDSAPI_KEY') or os.path.exists(os.path.expanduser('~/.cdsapirc'))):
        fail('no ERA5 source: set CDSAPI_KEY (Copernicus CDS key), mount a .cdsapirc at '
             '/root/.cdsapirc, or give ERA5_PL_FILE and ERA5_SL_FILE')
    return decks


def case_settings(a, decks, out, nstop, template):
    """&PARAMETERS settings passed to create_run.py; the user's --set/SAM_SET come last and win"""
    sets = ["asp_deck_dir='%s/'" % decks.rstrip('/'),
            "asp_output_dir='%s/'" % os.path.join(out, 'ASP')]
    with open(template) as f:
        nstat = require_prm_value(f.read(), 'nstat', template)
    if nstop < nstat:                                  # SAM: nstat must not exceed nstop
        sets.append('nstat=%d' % nstop)
    if a.fire_start_step is not None:
        sets.append('fire_nstep_start=%d' % a.fire_start_step)
    return sets + a.set


def build_case(a, decks, out, nstop):
    """write RUNS/<case>/ with create_run.py (ERA5 meteorology, prm from the template),
    copy it to OUTPUT_DIR/case; -> the case's prm text"""
    template = os.path.join(MODEL, 'RUNDATA', 'prm_template')
    cmd = [sys.executable, os.path.join(MODEL, 'create_run.py'), '--case', a.case,
           '--lat', str(a.lat), '--lon', str(a.lon), '--time', '%sT%s' % (a.date, a.start),
           '--massemis', repr(a.massemis), '--dt', repr(a.dt), '--hours', repr(a.hours),
           '--rundir', MODEL, '--cache', os.path.join(out, 'ERA5_cache')]
    if a.era5_pl:
        cmd += ['--pl-file', a.era5_pl, '--sl-file', a.era5_sl]
    for s in case_settings(a, decks, out, nstop, template):
        cmd += ['--set', s]
    subprocess.run(cmd, check=True)
    casedir = os.path.join(MODEL, 'RUNS', a.case)
    shutil.copytree(casedir, os.path.join(out, 'case'), dirs_exist_ok=True)
    with open(os.path.join(casedir, 'prm')) as f:
        return f.read()


def report_fire_window(prm, nstop, dt, nranks):
    """print the steps with fire emissions (SAM: fire_nstep_start+1 .. + 2*(plume_y-dy)/vg/dt)
    and warn if the run ends before they start or end; -> (first step, last step)"""
    fstart = int(require_prm_value(prm, 'fire_nstep_start', 'the case prm'))
    vg, plume_y, dy = (require_prm_value(prm, n, 'the case prm') for n in ('vg', 'plume_y', 'dy'))
    fend = fstart + (2 * (plume_y - dy) / vg / dt if vg > 0 else float('inf'))
    print('run: %d steps of %g s; fire emissions from step %d to ~%.0f; %d MPI ranks'
          % (nstop, dt, fstart, fend, nranks))
    if nstop <= fstart:
        print('WARNING: the run ends before the fire starts (no fire emissions); '
              'increase RUN_HOURS or set FIRE_START_STEP')
    elif nstop <= fend:
        print('WARNING: the run ends before the fire emissions end; ASP chemistry starts after them')
    return fstart, fend


def link_outputs(case, out):
    """in the model directory: set CaseName and point SAM's output directories at OUTPUT_DIR"""
    with open(os.path.join(MODEL, 'CaseName'), 'w') as f:
        f.write(case + '\n')
    for d in SAM_OUT:
        os.makedirs(os.path.join(out, d), exist_ok=True)
        link = os.path.join(MODEL, d)
        if os.path.islink(link) or os.path.isfile(link):
            os.remove(link)
        elif os.path.isdir(link):
            shutil.rmtree(link)
        os.symlink(os.path.join(out, d), link)


def run_sam(out, nranks, nstop):
    """run SAM with MPI in the model directory, log to OUTPUT_DIR/run.out; -> exit code"""
    try:
        resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))
    except (ValueError, OSError):
        print('WARNING: could not lift the stack limit; run the container with --ulimit stack=-1')
    # Ubuntu names it mpirun.openmpi; elsewhere (e.g. the RHEL 9 image) it is mpirun
    mpirun = shutil.which('mpirun.openmpi') or shutil.which('mpirun') or 'mpirun'
    runlog = os.path.join(out, 'run.out')
    with open(runlog, 'w') as log:
        rc = subprocess.run([mpirun, '--allow-run-as-root', '--oversubscribe', '-np', str(nranks), './' + EXE],
                            cwd=MODEL, stdout=log, stderr=subprocess.STDOUT).returncode
    last = [l for l in open(runlog, errors='replace') if 'NSTEP =' in l]
    step = int(last[-1].split()[2]) if last else 0
    print('SAM finished at step %d of %d (exit code %d); log: %s' % (step, nstop, rc, runlog))
    return 0 if step >= nstop else (rc or 1)


def main(argv=None):
    a = parse_args(argv)
    if a.check_decks:
        check_tracers(require_decks(a.decks))
        return 0
    decks = validate_inputs(a)
    out = os.path.abspath(a.output)
    os.makedirs(os.path.join(out, 'ASP'), exist_ok=True)
    nstop = int(round(a.hours * 3600.0 / a.dt))
    with open(os.path.join(MODEL, 'NRANKS')) as f:
        nranks = int(f.read())
    prm = build_case(a, decks, out, nstop)
    report_fire_window(prm, nstop, a.dt, nranks)
    link_outputs(a.case, out)
    return run_sam(out, nranks, nstop)


if __name__ == '__main__':
    sys.exit(main())
