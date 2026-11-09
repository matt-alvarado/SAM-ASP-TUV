"""Offline tests for create_run.py / samasp_case (no network needed).

Synthetic ERA5 files are written in the current CDS netCDF layout
(valid_time, pressure_level) and in the older layout (time, level, packed
shorts). Run:  python3 -m unittest discover -s tests -v
"""
import datetime as dt
import os
import re
import shutil
import sys
import tempfile
import unittest

import numpy as np
import netCDF4

HERE = os.path.dirname(os.path.abspath(__file__))
TOP = os.path.dirname(HERE)
sys.path.insert(0, TOP)
import create_run                                   # noqa: E402
from samasp_case import era5, case_files            # noqa: E402

LEVELS = np.array(era5.PRESSURE_LEVELS, dtype=float)
T_VALID = dt.datetime(2009, 11, 17, 20, 0)


def profile(level):
    """A smooth sounding: T (K), q (kg/kg), u, v (m/s) on pressure levels"""
    t = 288.0 - 6.5e-3 * 7.4e3 * np.log(1000.0 / level)
    q = 0.006 * (level / 1000.0) ** 3
    u = 2.0 + 10.0 * (1000.0 - level) / 900.0
    v = -1.0 + 3.0 * (1000.0 - level) / 900.0
    return t, q, u, v


def write_era5(path_pl, path_sl, layout='new', lat0=39.75, lon0=-122.5):
    """3 x 3 points, 3 hours (19, 20, 21 UTC); the centre point at 20 UTC holds
    profile(); other points/hours are offset so a wrong pick is detected"""
    lats = np.array([lat0 + 0.25, lat0, lat0 - 0.25])     # ERA5: descending
    lons = np.array([lon0 - 0.25, lon0, lon0 + 0.25])
    times = [T_VALID - dt.timedelta(hours=1), T_VALID, T_VALID + dt.timedelta(hours=1)]
    tname, lname = ('valid_time', 'pressure_level') if layout == 'new' else ('time', 'level')
    tunits = 'seconds since 1970-01-01' if layout == 'new' else 'hours since 1900-01-01 00:00:00.0'
    lev = LEVELS if layout == 'new' else LEVELS[::-1]      # old files: ascending order
    t, q, u, v = profile(lev)

    def axes(ds):
        ds.createDimension(tname, 3)
        ds.createDimension('latitude', 3)
        ds.createDimension('longitude', 3)
        vt = ds.createVariable(tname, 'i8' if layout == 'new' else 'i4', (tname,))
        vt.units = tunits
        vt.calendar = 'proleptic_gregorian'
        vt[:] = netCDF4.date2num(times, tunits, 'proleptic_gregorian')
        ds.createVariable('latitude', 'f4', ('latitude',))[:] = lats
        ds.createVariable('longitude', 'f4', ('longitude',))[:] = lons

    def var(ds, name, dims, data):
        if layout == 'new':
            x = ds.createVariable(name, 'f4', dims)
        else:                                               # packed short
            x = ds.createVariable(name, 'i2', dims)
            lo, hi = float(np.min(data)), float(np.max(data))
            x.scale_factor = (hi - lo) / 60000.0 if hi > lo else 1.0
            x.add_offset = (hi + lo) / 2.0
        x[:] = data

    with netCDF4.Dataset(path_pl, 'w') as ds:
        axes(ds)
        ds.createDimension(lname, len(lev))
        ds.createVariable(lname, 'f8', (lname,))[:] = lev
        shape = (3, len(lev), 3, 3)
        for name, prof in (('t', t), ('q', q), ('u', u), ('v', v)):
            a = np.empty(shape)
            for it in range(3):
                for j in range(3):
                    for i in range(3):
                        a[it, :, j, i] = prof * (1 + 0.01 * (abs(it - 1) + abs(j - 1) + abs(i - 1)))
            var(ds, name, (tname, lname, 'latitude', 'longitude'), a)
    sl = dict(sp=98877.0, blh=850.0, sshf=-244.57 * 3600, slhf=-43.339 * 3600,
              iews=0.15, inss=-0.12, t2m=290.0)
    with netCDF4.Dataset(path_sl, 'w') as ds:
        axes(ds)
        for name, val in sl.items():
            a = np.full((3, 3, 3), val)
            a[0] *= 1.5
            a[2] *= 0.5                                    # other hours differ
            a[1, 0, 0] *= 3.0                              # another point differs
            var(ds, name, (tname, 'latitude', 'longitude'), a)
    return sl


def read(path):
    with open(path) as f:
        return f.read()


class Era5Read(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.d)

    def check_layout(self, layout):
        pl, sl = os.path.join(self.d, 'pl.nc'), os.path.join(self.d, 'sl.nc')
        s = write_era5(pl, sl, layout)
        # a fire 0.08 deg from the centre point, 12 min after the hour
        m = era5.read_column(pl, sl, 39.70, -122.58, dt.datetime(2009, 11, 17, 20, 12))
        t, q, u, v = profile(LEVELS)
        tol = 1e-5 if layout == 'new' else 2e-3              # packed shorts
        self.assertEqual(m.time, T_VALID)
        self.assertAlmostEqual(m.lat, 39.75, 5)
        self.assertAlmostEqual(m.lon, -122.5, 5)
        np.testing.assert_allclose(m.level, LEVELS)          # surface first
        np.testing.assert_allclose(m.temp, t, rtol=tol)
        np.testing.assert_allclose(m.u, u, rtol=tol, atol=tol)
        self.assertAlmostEqual(m.psfc, s['sp'] / 100.0, delta=0.02)
        self.assertAlmostEqual(m.hpbl, s['blh'], delta=0.1)
        self.assertAlmostEqual(m.shf, 244.57, delta=0.02)    # upward positive
        self.assertAlmostEqual(m.lhf, 43.339, delta=0.02)
        self.assertAlmostEqual(m.taux, 0.15, delta=1e-4)

    def test_new_cds_layout(self):
        self.check_layout('new')

    def test_old_cds_layout_packed(self):
        self.check_layout('old')

    def test_wrong_file_is_an_error(self):
        pl, sl = os.path.join(self.d, 'pl.nc'), os.path.join(self.d, 'sl.nc')
        write_era5(pl, sl)
        with self.assertRaises(SystemExit):                  # far from the file's box
            era5.read_column(pl, sl, 34.7, -120.0, T_VALID)
        with self.assertRaises(SystemExit):                  # hour not in the file
            era5.read_column(pl, sl, 39.75, -122.5, T_VALID + dt.timedelta(hours=5))

    def test_hour_rounding(self):
        self.assertEqual(era5.nearest_hour(dt.datetime(2009, 1, 1, 5, 29)), dt.datetime(2009, 1, 1, 5))
        self.assertEqual(era5.nearest_hour(dt.datetime(2009, 1, 1, 23, 30)), dt.datetime(2009, 1, 2, 0))


class Derived(unittest.TestCase):
    def column(self):
        t, q, u, v = profile(LEVELS)
        return era5.MetColumn(level=LEVELS.copy(), temp=t, shum=q, u=u, v=v, psfc=988.77,
                              hpbl=850.0, shf=244.57, lhf=43.339, taux=0.15, tauy=-0.12,
                              t2m=290.0, lat=39.75, lon=-122.5, time=T_VALID)

    def test_thermo_and_rotation(self):
        m = self.column()
        d = case_files.derived(m, 'Ay')
        np.testing.assert_allclose(d['ptemp'], m.temp * (1000.0 / m.level) ** case_files.RC)
        np.testing.assert_allclose(d['wmix'], m.shum / (1 - m.shum) * 1000.0)
        # boundary-layer levels: above ground (1000 hPa is below the 988.77 hPa
        # surface here) and below the boundary-layer top
        ptop = min(d['hpblp'], m.psfc - 50.0)
        bl = (m.level <= m.psfc) & (m.level > ptop)
        np.testing.assert_array_equal(d['inbl'], bl)
        self.assertFalse(d['inbl'][0])
        self.assertTrue(bl.sum() >= 2)
        # Ay: SAM y along the mean boundary-layer wind: after rotation its mean
        # x component is 0 and its mean y component is the mean wind's speed
        mu, mv = np.mean(m.u[bl]), np.mean(m.v[bl])
        self.assertAlmostEqual(float(np.mean(d['uls'][bl])), 0.0, places=10)
        self.assertAlmostEqual(float(np.mean(d['vls'][bl])), float(np.hypot(mu, mv)), places=10)
        # rotation keeps wind speed at every level
        np.testing.assert_allclose(np.hypot(d['uls'], d['vls']), np.hypot(m.u, m.v))
        # vg: mean of vls over the boundary-layer levels
        self.assertAlmostEqual(d['vg'], float(np.mean(d['vls'][bl])), places=12)
        # tau = |stress| / rho with rho from surface pressure and 2 m temperature
        rho = m.psfc * 100 / (case_files.RD * m.t2m)
        self.assertAlmostEqual(d['tau'], np.hypot(0.15, -0.12) / rho, places=12)

    def test_no_level_in_boundary_layer(self):
        # the layer is at least 50 hPa deep (ptop <= psfc - 50), so with ERA5's
        # 25-50 hPa spacing it always holds a level; with sparser levels it can
        # be empty, and then the lowest level above the ground is used
        lev = np.array([1000.0, 850.0, 700.0, 500.0])
        t, q, u, v = profile(lev)
        m = era5.MetColumn(level=lev, temp=t, shum=q, u=u, v=v, psfc=990.0, hpbl=300.0,
                           shf=0.0, lhf=0.0, taux=0.1, tauy=0.0, t2m=290.0,
                           lat=0.0, lon=0.0, time=T_VALID)
        d = case_files.derived(m, 'Ay')
        np.testing.assert_array_equal(d['inbl'], [False, True, False, False])
        self.assertAlmostEqual(d['uls'][1], 0.0, places=10)
        self.assertAlmostEqual(d['vg'], float(np.hypot(u[1], v[1])), places=10)

    def test_other_axes(self):
        m = self.column()
        for ax, (uu, vv) in {'E': (m.u, m.v), 'W': (-m.u, -m.v), 'N': (m.v, -m.u),
                             'S': (-m.v, m.u)}.items():
            d = case_files.derived(m, ax)
            np.testing.assert_allclose(d['uls'], uu)
            np.testing.assert_allclose(d['vls'], vv)
        d = case_files.derived(m, 'U')
        np.testing.assert_allclose(d['uls'], np.hypot(m.u, m.v))

    def test_calendar_day(self):
        self.assertAlmostEqual(case_files.calendar_day(dt.datetime(2009, 1, 1)), 1.0)
        self.assertAlmostEqual(case_files.calendar_day(T_VALID), 321.0 + 20 / 24.0)


class FakeClient:
    def __init__(self, pl, sl):
        self.calls, self.src = [], {'reanalysis-era5-pressure-levels': pl,
                                    'reanalysis-era5-single-levels': sl}

    def retrieve(self, dataset, request, target):
        self.calls.append((dataset, request))
        shutil.copy(self.src[dataset], target)


class DownloadAndCase(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.pl, self.sl = os.path.join(self.d, 'pl.nc'), os.path.join(self.d, 'sl.nc')
        write_era5(self.pl, self.sl)
        self.run = os.path.join(self.d, 'run')
        os.makedirs(os.path.join(self.run, 'RUNDATA'))
        shutil.copy(os.path.join(TOP, 'RUNDATA', 'prm_template'), os.path.join(self.run, 'RUNDATA'))
        shutil.copy(os.path.join(TOP, 'RUNDATA', 'lst'), os.path.join(self.run, 'RUNDATA'))

    def tearDown(self):
        shutil.rmtree(self.d)

    def test_requests_and_cache(self):
        c = FakeClient(self.pl, self.sl)
        cache = os.path.join(self.d, 'cache')
        t = dt.datetime(2009, 11, 17, 19, 50)
        f1 = era5.download(39.70, -122.58, t, cache, client=c)
        self.assertEqual(len(c.calls), 2)
        (ds1, r1), (ds2, r2) = c.calls
        self.assertEqual(ds1, 'reanalysis-era5-pressure-levels')
        self.assertEqual(r1['time'], ['20:00'])
        self.assertEqual(r1['day'], ['17'])
        self.assertEqual(r1['pressure_level'][0], '1000')
        self.assertEqual(r1['area'], [40.0, -122.75, 39.5, -122.25])
        self.assertIn('boundary_layer_height', r2['variable'])
        era5.download(39.70, -122.58, t, cache, client=c)   # cached: no new calls
        self.assertEqual(len(c.calls), 2)
        self.assertTrue(all(os.path.exists(f) for f in f1))

    def test_case_files(self):
        create_run.main(['--case', 'TESTCASE', '--lat', '39.70', '--lon', '-122.58',
                         '--time', '2009-11-17T20:12', '--massemis', '4.86e-6',
                         '--rundir', self.run, '--pl-file', self.pl, '--sl-file', self.sl])
        cd = os.path.join(self.run, 'RUNS', 'TESTCASE')
        for f in ('snd', 'lsf', 'sfc', 'prm', 'grd', 'lst', 'met_source.txt'):
            self.assertTrue(os.path.exists(os.path.join(cd, f)), f)
        snd = read(os.path.join(cd, 'snd')).splitlines()
        nlev = len(LEVELS)
        self.assertEqual(len(snd), 1 + 3 * (nlev + 1))
        head = snd[1].replace(',', ' ').split()
        self.assertAlmostEqual(float(head[0]), 321.0 + 20 / 24.0, places=5)
        self.assertEqual(int(head[1]), nlev)
        self.assertAlmostEqual(float(head[2]), 988.77, places=1)
        row = snd[2].split()
        self.assertEqual(row[0], '-999.')
        self.assertEqual(float(row[1]), 1000.0)
        sfc = read(os.path.join(cd, 'sfc')).splitlines()
        self.assertEqual(len(sfc), 4)
        self.assertAlmostEqual(float(sfc[1].split()[2]), 244.57, places=1)
        prm = read(os.path.join(cd, 'prm'))
        self.assertNotIn('{', prm)
        self.assertIn("caseid = 'TESTCASE'", prm)
        self.assertRegex(prm, r'day0 = 321\.833333')
        self.assertIn('year0 = 2009', prm)
        self.assertIn('nstop = 900', prm)                    # 2 h at dt = 8 s
        self.assertIn('massemis = 4.86e-06', prm)
        self.assertIn('longitude0 = -122.58', prm)

    def test_params_file(self):
        p = os.path.join(self.d, 'params.txt')
        with open(p, 'w') as f:
            f.write('lat,lon,year,month,day,hour,minute,Massemis,\n39.70,-122.58,2009,11,17,20,0,4.86e-06,\n')
        create_run.main(['--case', 'P', '--params', p, '--rundir', self.run,
                         '--pl-file', self.pl, '--sl-file', self.sl])
        self.assertTrue(os.path.exists(os.path.join(self.run, 'RUNS', 'P', 'prm')))


if __name__ == '__main__':
    unittest.main()
