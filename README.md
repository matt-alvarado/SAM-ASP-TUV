SAM-ASP-TUV © 2024 by Matthew Alvarado, Chantelle Lonsdale, Archana Dayalu, and Marikate Mountain of Atmospheric and Environmental Research, Inc. and Nicole June and Jeffrey Pierce of Colorado State University is licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International. To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/


## Creating a case (non-fire meteorology from ERA5)

`create_run.py` writes `RUNS/<case>/` (`snd`, `lsf`, `sfc`, `prm`, `grd`, `lst`
and `met_source.txt`) for a fire's location and start time, with the
meteorology from the ERA5 reanalysis column nearest the fire:

    pip install -r requirements.txt          # numpy, netCDF4, cdsapi
    python3 create_run.py --case MYFIRE --lat 39.70 --lon -122.58 \
        --time 2009-11-17T20:00 --massemis 4.86e-6

ERA5 is downloaded from the Copernicus Climate Data Store (a 3 x 3 box of
grid points and one hour, a few kB) and cached in `ERA5_cache/`. This needs a
free CDS account: put your API key in `~/.cdsapirc` (`url:
https://cds.climate.copernicus.eu/api` and `key: <token>`) and accept the ERA5
licence on the dataset pages once. `--pl-file`/`--sl-file` use ERA5 netCDF
files you already have instead. The `prm` comes from `RUNDATA/prm_template`
(edit it, or pass `--prm-template`); `--dt`, `--hours` and `--wind-axis` set
the time step, run length and axis orientation. `--params params.txt` reads
the old input file format. Tests (offline): `python3 -m unittest discover -s tests`.
