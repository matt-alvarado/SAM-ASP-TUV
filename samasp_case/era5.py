"""ERA5 meteorology for a SAM-ASP case: download (Copernicus CDS) and read
the column at the grid point nearest the fire.

Data: ERA5 hourly reanalysis (ECMWF), 0.25 deg, global, from the Copernicus
Climate Data Store (https://cds.climate.copernicus.eu):
  reanalysis-era5-pressure-levels: temperature, specific humidity, u and v
      wind on PRESSURE_LEVELS (hPa)
  reanalysis-era5-single-levels: surface pressure, boundary layer height,
      surface sensible and latent heat flux (accumulated over the hour ending
      at the valid time, J m-2, positive downward), instantaneous eastward
      and northward turbulent surface stress (N m-2) and 2 m temperature

Only a 3 x 3 box of grid points around the fire and one hour are requested
(a few kB). Retrieval needs the cdsapi package and a CDS account: put
    url: https://cds.climate.copernicus.eu/api
    key: <your personal access token>
in ~/.cdsapirc (see the CDS "How to API" page) and accept the ERA5 licence
on the dataset pages once. Downloaded files are cached and reused.
"""
import datetime as dt
import os

import numpy as np

PRESSURE_LEVELS = [1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 700,
                   650, 600, 550, 500, 450, 400, 350, 300, 250, 225, 200, 175,
                   150, 125, 100]                       # hPa, surface first
PL_VARIABLES = ['temperature', 'specific_humidity', 'u_component_of_wind',
                'v_component_of_wind']
SL_VARIABLES = ['surface_pressure', 'boundary_layer_height',
                'surface_sensible_heat_flux', 'surface_latent_heat_flux',
                'instantaneous_eastward_turbulent_surface_stress',
                'instantaneous_northward_turbulent_surface_stress',
                '2m_temperature']
GRID = 0.25                                             # deg
RD = 287.05                                             # J kg-1 K-1, dry air


class MetColumn:
    """ERA5 values at one grid point and hour.

    level  pressure levels (hPa), surface first
    temp   temperature (K); shum specific humidity (kg/kg); u, v wind (m/s)
    psfc   surface pressure (hPa); hpbl boundary layer height (m)
    shf, lhf   sensible and latent heat flux (W m-2), positive UPWARD (into
               the atmosphere), the hour ending at `time`
    taux, tauy surface stress (N m-2); t2m 2 m temperature (K)
    lat, lon, time  the grid point and valid time used
    """

    def __init__(self, **kw):
        self.__dict__.update(kw)


def nearest_hour(t):
    """t rounded to the nearest whole hour"""
    t0 = t.replace(minute=0, second=0, microsecond=0)
    return t0 + dt.timedelta(hours=1) if t - t0 >= dt.timedelta(minutes=30) else t0


def _area(lat, lon):
    """[N, W, S, E] of a 3 x 3 box of ERA5 grid points around (lat, lon)"""
    la = round(lat / GRID) * GRID
    lo = round(lon / GRID) * GRID
    return [min(90.0, la + GRID), lo - GRID, max(-90.0, la - GRID), lo + GRID]


def requests(lat, lon, time):
    """The two CDS requests (dataset, request dict) for one case"""
    t = nearest_hour(time)
    common = {'product_type': ['reanalysis'], 'year': ['%04d' % t.year],
              'month': ['%02d' % t.month], 'day': ['%02d' % t.day],
              'time': ['%02d:00' % t.hour], 'area': _area(lat, lon),
              'data_format': 'netcdf', 'download_format': 'unarchived'}
    pl = dict(common, variable=PL_VARIABLES,
              pressure_level=[str(p) for p in PRESSURE_LEVELS])
    sl = dict(common, variable=SL_VARIABLES)
    return [('reanalysis-era5-pressure-levels', pl), ('reanalysis-era5-single-levels', sl)]


def cache_names(cache_dir, lat, lon, time):
    t = nearest_hour(time)
    a = _area(lat, lon)
    tag = '%s_%.2f_%.2f' % (t.strftime('%Y%m%d%H'), (a[0] + a[2]) / 2, (a[1] + a[3]) / 2)
    return (os.path.join(cache_dir, 'era5_pl_%s.nc' % tag),
            os.path.join(cache_dir, 'era5_sl_%s.nc' % tag))


def download(lat, lon, time, cache_dir, client=None):
    """Fetch the pressure-level and single-level files (if not cached)"""
    os.makedirs(cache_dir, exist_ok=True)
    files = cache_names(cache_dir, lat, lon, time)
    todo = [(r, f) for r, f in zip(requests(lat, lon, time), files) if not os.path.exists(f)]
    if todo:
        if client is None:
            try:
                import cdsapi
            except ImportError:
                raise SystemExit('ERA5 download needs the cdsapi package (pip install cdsapi) '
                                 'and a CDS key in ~/.cdsapirc; see samasp_case/era5.py')
            client = cdsapi.Client()
        for (dataset, req), f in todo:
            print('downloading %s -> %s' % (dataset, f))
            client.retrieve(dataset, req, f + '.part')
            os.replace(f + '.part', f)
    return files


# ---------------------------------------------------------------- reading
def _dim(ds, names):
    for n in names:
        if n in ds.dimensions or n in ds.variables:
            return n
    raise KeyError('none of %s in %s' % (names, ds.filepath()))


def _times(ds):
    import netCDF4
    name = _dim(ds, ('valid_time', 'time'))
    v = ds.variables[name]
    out = netCDF4.num2date(v[:], v.units, getattr(v, 'calendar', 'standard'),
                           only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    return name, [dt.datetime(x.year, x.month, x.day, x.hour, x.minute) for x in np.atleast_1d(out)]


def _point(ds, lat, lon):
    la = ds.variables[_dim(ds, ('latitude', 'lat'))][:]
    lo = ds.variables[_dim(ds, ('longitude', 'lon'))][:]
    j = int(np.argmin(np.abs(la - lat)))
    i = int(np.argmin(np.abs(((lo - lon + 180.0) % 360.0) - 180.0)))
    return j, i, float(la[j]), float(lo[i])


def _get(ds, var, it, j, i, tname):
    """Values of var at time index it and point (j, i), any other dims kept
    (pressure level) and singleton dims (expver, number) squeezed"""
    v = ds.variables[var]
    idx = []
    for d in v.dimensions:
        if d == tname:
            idx.append(it)
        elif d in ('latitude', 'lat'):
            idx.append(j)
        elif d in ('longitude', 'lon'):
            idx.append(i)
        else:
            idx.append(slice(None))
    a = np.ma.filled(np.ma.asarray(v[tuple(idx)]).astype(float), np.nan)
    return np.squeeze(a)


def _check_point(path, lat, lon, glat, glon):
    dlon = abs(((glon - lon + 180.0) % 360.0) - 180.0)
    if abs(glat - lat) > GRID or dlon > GRID:
        raise SystemExit('%s does not cover %.3f N %.3f E (nearest grid point %.2f N %.2f E)'
                         % (path, lat, lon, glat, glon))


def _time_index(path, times, t):
    if t not in times:
        raise SystemExit('%s has no data for %s UTC (has %s .. %s)'
                         % (path, t, min(times), max(times)))
    return times.index(t)


def read_column(pl_file, sl_file, lat, lon, time):
    """MetColumn at the grid point nearest (lat, lon) and the hour nearest time"""
    import netCDF4
    t = nearest_hour(time)
    with netCDF4.Dataset(pl_file) as ds:
        tname, times = _times(ds)
        it = _time_index(pl_file, times, t)
        j, i, glat, glon = _point(ds, lat, lon)
        _check_point(pl_file, lat, lon, glat, glon)
        lname = _dim(ds, ('pressure_level', 'level'))
        level = np.asarray(ds.variables[lname][:], dtype=float)
        order = np.argsort(-level)                     # surface (high p) first
        col = {k: _get(ds, k, it, j, i, tname)[order] for k in ('t', 'q', 'u', 'v')}
        level = level[order]
    with netCDF4.Dataset(sl_file) as ds:
        tname, times = _times(ds)
        it = _time_index(sl_file, times, t)
        j2, i2, slat, slon = _point(ds, lat, lon)
        _check_point(sl_file, lat, lon, slat, slon)
        s = {k: float(_get(ds, k, it, j2, i2, tname))
             for k in ('sp', 'blh', 'sshf', 'slhf', 'iews', 'inss', 't2m')}
    return MetColumn(level=level, temp=col['t'], shum=col['q'], u=col['u'], v=col['v'],
                     psfc=s['sp'] / 100.0, hpbl=s['blh'],
                     shf=-s['sshf'] / 3600.0, lhf=-s['slhf'] / 3600.0,   # J m-2 per hour, down -> W m-2, up
                     taux=s['iews'], tauy=s['inss'], t2m=s['t2m'],
                     lat=glat, lon=glon, time=t)
