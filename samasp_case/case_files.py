"""Write a SAM-ASP case (snd, lsf, sfc, prm, grd, lst) from a MetColumn.

The meteorology is held constant: the same profile is written at NPTS times
one day apart (SAM needs input times spanning the run). The conversions are
those of the original create_run.py (NARR version); the wind axes ('Ax', 'Ay')
and vg use the mean wind over the boundary-layer levels above the ground.
"""
import datetime as dt
import os
import shutil

import numpy as np

NPTS = 3                                  # times written to snd/lsf/sfc
P0 = 1000.0                               # hPa, reference for potential temperature
RC = 8.314472 / 29.19                     # R / cp (as in the original script)
HOVERT = 8.31432 / (0.0289644 * 9.80665)  # scale height / temperature (m/K)
RD = 287.05                               # J kg-1 K-1
WIND_AXES = ('U', 'N', 'S', 'E', 'W', 'T', 'Ax', 'Ay')


def calendar_day(t):
    """Calendar day of t's year as in SAM's CAM radiation (1.0 = 1 Jan 00 UTC)"""
    d = t - dt.datetime(t.year, 1, 1)
    return d.days + 1 + d.seconds / 86400.0


def derived(m, wind_axis='Ay', theta=np.pi * 0.40):
    """Quantities written to the case files.

    Returns dict with ptemp (K), wmix (g/kg), uls, vls (m/s, rotated to the
    SAM axes), vg (mean boundary-layer wind along SAM y, m/s), hpblp (hPa),
    tau (u*^2, m2 s-2), theta (rotation angle used, rad), inbl (the
    boundary-layer levels: above ground, below the boundary-layer top).
    """
    if wind_axis not in WIND_AXES:
        raise ValueError('wind_axis must be one of %s' % (WIND_AXES,))
    level, temp, u, v = m.level, m.temp, m.u, m.v
    ptemp = temp * (P0 / level) ** RC
    wmix = m.shum / (1.0 - m.shum) * 1000.0
    # boundary-layer top in pressure (lowest-level temperature, < 2% error)
    hpblp = m.psfc * np.exp(-m.hpbl / (HOVERT * temp[0]))
    ptop = min(hpblp, m.psfc - 50.0)
    # boundary-layer levels: above the ground (ERA5 also gives values below
    # it, extrapolated) and below the boundary-layer top; if no level lies
    # in that layer, the lowest level above the ground
    above = level <= m.psfc
    inbl = above & (level > ptop)
    if not inbl.any():
        inbl = np.zeros_like(above)
        inbl[np.where(above)[0][0]] = True
    if wind_axis in ('Ax', 'Ay'):
        # direction of the mean boundary-layer wind (vector mean of u and v)
        theta = float(np.arctan2(np.mean(v[inbl]), np.mean(u[inbl])))
    if wind_axis == 'N':
        uls, vls = v, -u
    elif wind_axis == 'S':
        uls, vls = -v, u
    elif wind_axis == 'E':
        uls, vls = u, v
    elif wind_axis == 'W':
        uls, vls = -u, -v
    elif wind_axis in ('T', 'Ax'):
        uls = u * np.cos(theta) + v * np.sin(theta)
        vls = -u * np.sin(theta) + v * np.cos(theta)
    elif wind_axis == 'Ay':
        uls = u * np.sin(theta) - v * np.cos(theta)
        vls = u * np.cos(theta) + v * np.sin(theta)
    else:                                  # 'U': all wind along SAM x
        uls, vls = np.sqrt(u ** 2 + v ** 2), v * 0.0
    # mean boundary-layer wind along SAM y (plume transport speed)
    vg = float(np.mean(vls[inbl]))
    rho = m.psfc * 100.0 / (RD * m.t2m)
    tau = float(np.hypot(m.taux, m.tauy) / rho)
    return dict(ptemp=ptemp, wmix=wmix, uls=uls, vls=vls, vg=vg, hpblp=hpblp,
                tau=tau, theta=theta, rho=rho, inbl=inbl)


def write_case(casedir, m, d, prm_template, prm_values, lst_file):
    """Write snd, lsf, sfc, grd, lst, prm and met_source.txt into casedir"""
    os.makedirs(casedir, exist_ok=True)
    day0 = calendar_day(m.time)
    days = [day0 + i for i in range(NPTS)]
    nlev = len(m.level)
    with open(os.path.join(casedir, 'snd'), 'w') as f:
        f.write(' z[m] p[mb] tp[K] q[g/kg] u[m/s] v[m/s]\n')
        for day in days:
            f.write(' %.6f,  %d, %.2f  day,levels,pres0\n' % (day, nlev, m.psfc))
            for k in range(nlev):
                f.write('-999.\t%.1f\t%.3f\t%.5f\t%.4f\t%.4f\n'
                        % (m.level[k], d['ptemp'][k], d['wmix'][k], d['uls'][k], d['vls'][k]))
    with open(os.path.join(casedir, 'lsf'), 'w') as f:
        f.write(' z[m] p[mb]  tls[K/s] qls[kg/kg/s] uls  vls  wls[m/s]\n')
        for day in days:
            f.write(' %.6f,  %d, %.2f  day,levels,pres0\n' % (day, nlev, m.psfc))
            for k in range(nlev):
                f.write('-999.\t%.1f\t0.\t0.\t%.4f\t%.4f\t0.\n' % (m.level[k], d['uls'][k], d['vls'][k]))
    with open(os.path.join(casedir, 'sfc'), 'w') as f:
        f.write('day sst(K) H(W/m2) LE(W/m2) TAU(m2/s2)\n')
        for day in days:
            f.write(' %.6f\t0.\t%.3f\t%.3f\t%.5f\n' % (day, m.shf, m.lhf, d['tau']))
    with open(os.path.join(casedir, 'grd'), 'w') as f:
        f.write('20\n60\n100')
    shutil.copy(lst_file, os.path.join(casedir, 'lst'))
    values = dict(prm_values, day0='%.6f' % day0, year0=m.time.year,
                  latitude0=prm_values['latitude0'], longitude0=prm_values['longitude0'],
                  vg='%.4f' % d['vg'])
    with open(prm_template) as f:
        text = f.read()
    try:
        prm = text.format(**values)
    except KeyError as e:
        raise SystemExit('prm template %s uses {%s}, which create_run.py does not provide'
                         % (prm_template, e.args[0]))
    with open(os.path.join(casedir, 'prm'), 'w') as f:
        f.write(prm)
    with open(os.path.join(casedir, 'met_source.txt'), 'w') as f:
        f.write('ERA5 hourly reanalysis (ECMWF, Copernicus Climate Data Store)\n')
        f.write('grid point: %.2f N, %.2f E; valid time %s UTC (requested %s)\n'
                % (m.lat, m.lon, m.time.strftime('%Y-%m-%d %H:%M'),
                   prm_values.get('requested_time', '')))
        f.write('surface pressure %.2f hPa; boundary layer height %.1f m (%.1f hPa)\n'
                % (m.psfc, m.hpbl, d['hpblp']))
        f.write('sensible heat flux %.2f W/m2, latent heat flux %.2f W/m2 (upward, hour ending '
                'at the valid time)\n' % (m.shf, m.lhf))
        f.write('surface stress %.4f N/m2, air density %.3f kg/m3, tau = u*^2 = %.5f m2/s2\n'
                % (float(np.hypot(m.taux, m.tauy)), d['rho'], d['tau']))
        f.write('wind axis %s, rotation %.4f rad; mean boundary-layer wind vg = %.4f m/s\n'
                % (prm_values.get('wind_axis', ''), d['theta'], d['vg']))
        f.write('meteorology held constant at %d times one day apart from day0 = %.6f (%d)\n'
                % (NPTS, day0, m.time.year))
