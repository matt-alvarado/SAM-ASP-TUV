#!/usr/bin/env python3
"""Create a SAM-ASP-TUV case directory (RUNS/<case>/) with the non-fire
meteorological inputs from ERA5.

For the fire's location and time this writes
  snd  sounding: potential temperature, water vapour mixing ratio and wind on
       pressure levels (ERA5, grid point nearest the fire, nearest hour)
  lsf  large-scale forcing: the same winds (no large-scale tendencies)
  sfc  surface sensible and latent heat flux and surface stress (u*^2)
  prm  namelist, from a template (default RUNDATA/prm_template)
  grd, lst, and met_source.txt (where the meteorology came from)
The meteorology is held constant over the run (the same profile at 3 times,
one day apart). SAM's y axis is aligned with the boundary-layer wind
(--wind-axis Ay, as before), and vg (the plume transport speed in prm) is the
mean boundary-layer wind along y.

Examples
  python3 create_run.py --case W1101 --lat 39.69583 --lon -122.6 \\
                        --time 2009-12-17T10:00 --massemis 4.86e-6
  python3 create_run.py --case W1101 --params params.txt     # old input file

ERA5 is downloaded from the Copernicus Climate Data Store; this needs the
cdsapi package and a CDS key in ~/.cdsapirc (see samasp_case/era5.py). Use
--pl-file/--sl-file to use ERA5 netCDF files you already have instead.
Replaces the NARR/NCEP version of this script (in git history before
this change).
"""
import argparse
import datetime as dt
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from samasp_case import era5, case_files          # noqa: E402


def parse_time(s):
    for fmt in ('%Y-%m-%dT%H:%M', '%Y-%m-%d %H:%M', '%Y-%m-%dT%H', '%Y-%m-%d'):
        try:
            return dt.datetime.strptime(s, fmt)
        except ValueError:
            pass
    raise argparse.ArgumentTypeError('time must be YYYY-MM-DDTHH:MM (UTC), got %r' % s)


def read_params(path):
    """Cases from the original params.txt: a header line, then
    lat,lon,year,month,day,hour,minute,massemis per line"""
    cases = []
    with open(path) as f:
        lines = [l for l in f.read().splitlines() if l.strip()]
    for l in lines[1:]:
        v = [x.strip() for x in l.split(',') if x.strip()]
        cases.append(dict(lat=float(v[0]), lon=float(v[1]),
                          time=dt.datetime(int(v[2]), int(v[3]), int(v[4]), int(v[5]), int(v[6])),
                          massemis=float(v[7])))
    return cases


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--case', required=True, help='case name (RUNS/<case>/, caseid in prm)')
    ap.add_argument('--lat', type=float, help='fire latitude (deg N)')
    ap.add_argument('--lon', type=float, help='fire longitude (deg E, -180..360)')
    ap.add_argument('--time', type=parse_time, help='start time, UTC (YYYY-MM-DDTHH:MM)')
    ap.add_argument('--massemis', type=float, help='mass emission flux of the fire (prm massemis)')
    ap.add_argument('--params', help='take lat, lon, time and massemis from an old-style params.txt '
                    '(several lines give cases <case>_01, <case>_02, ...)')
    ap.add_argument('--rundir', default=HERE, help='SAM-ASP-TUV checkout (default: this script\'s directory)')
    ap.add_argument('--prm-template', help='prm template (default <rundir>/RUNDATA/prm_template)')
    ap.add_argument('--dt', type=float, default=8.0, help='SAM time step (s), default 8')
    ap.add_argument('--hours', type=float, default=2.0, help='run length (h), default 2')
    ap.add_argument('--wind-axis', default='Ay', choices=case_files.WIND_AXES,
                    help='orientation of SAM\'s axes (default Ay: y along the boundary-layer wind)')
    ap.add_argument('--cache', help='ERA5 download cache (default <rundir>/ERA5_cache)')
    ap.add_argument('--pl-file', help='use this ERA5 pressure-level netCDF file (no download)')
    ap.add_argument('--sl-file', help='use this ERA5 single-level netCDF file (no download)')
    ap.add_argument('--set', action='append', default=[], metavar='NAME=VALUE',
                    help='set a &PARAMETERS entry in prm (Fortran syntax, e.g. '
                    "--set fire_nstep_start=0 --set \"asp_deck_dir='decks/'\"); repeatable")
    a = ap.parse_args(argv)

    if a.params:
        cases = read_params(a.params)
        names = [a.case] if len(cases) == 1 else ['%s_%02d' % (a.case, i + 1) for i in range(len(cases))]
    else:
        missing = [k for k in ('lat', 'lon', 'time', 'massemis') if getattr(a, k) is None]
        if missing:
            ap.error('give --%s (or --params)' % ', --'.join(missing))
        cases = [dict(lat=a.lat, lon=a.lon, time=a.time, massemis=a.massemis)]
        names = [a.case]
    if bool(a.pl_file) != bool(a.sl_file):
        ap.error('give both --pl-file and --sl-file, or neither')
    if a.pl_file and len(cases) > 1:
        ap.error('--pl-file/--sl-file can only be used for one case')

    template = a.prm_template or os.path.join(a.rundir, 'RUNDATA', 'prm_template')
    lst = os.path.join(a.rundir, 'RUNDATA', 'lst')
    cache = a.cache or os.path.join(a.rundir, 'ERA5_cache')
    nstop = int(round(a.hours * 3600.0 / a.dt))
    settings = []
    for kv in a.set:
        if '=' not in kv:
            ap.error('--set needs NAME=VALUE, got %r' % kv)
        k, v = kv.split('=', 1)
        settings.append((k.strip(), v.strip()))

    for name, c in zip(names, cases):
        lon = c['lon'] - 360.0 if c['lon'] > 180.0 else c['lon']
        if a.pl_file:
            pl, sl = a.pl_file, a.sl_file
        else:
            pl, sl = era5.download(c['lat'], lon, c['time'], cache)
        m = era5.read_column(pl, sl, c['lat'], lon, c['time'])
        d = case_files.derived(m, a.wind_axis)
        casedir = os.path.join(a.rundir, 'RUNS', name)
        case_files.write_case(casedir, m, d, template,
                              dict(caseid=name, latitude0=c['lat'], longitude0=lon, dt=a.dt,
                                   nstop=nstop, massemis=c['massemis'],
                                   requested_time=c['time'].strftime('%Y-%m-%d %H:%M'),
                                   wind_axis=a.wind_axis), lst, settings)
        print('%s: ERA5 %.2fN %.2fE %s UTC; psfc %.1f hPa, PBL %.0f m, H %.1f, LE %.1f W/m2, '
              'vg %.2f m/s -> %s' % (name, m.lat, m.lon, m.time.strftime('%Y-%m-%d %H:%M'),
                                     m.psfc, m.hpbl, m.shf, m.lhf, d['vg'], casedir))
    return 0


if __name__ == '__main__':
    sys.exit(main())
