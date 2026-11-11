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

## Running in a container

Prepared fire cases (Williams Flats 2019, Williams 2009, BBOP 2013) are in
[cases/](cases/README.md): `docker run ... --env-file cases/<case>/case.env sam-asp-tuv`.

See [docs/Container_User_Guide.md](docs/Container_User_Guide.md) for the full
guide (settings, ASP decks, ERA5, outputs, troubleshooting, and running on AWS),
and [docs/OnPrem_RHEL_CMMC_Guide.md](docs/OnPrem_RHEL_CMMC_Guide.md) for an on-prem
RHEL VM in a CMMC Level 2 (CUI) enclave, with rootless Podman. `Dockerfile.ubi9`
builds the same image on Red Hat's Universal Base Image 9 (RHEL 9), and
`Dockerfile.ubi9-gnu` builds it with gfortran for ARM64 (Apple silicon; see
[docs/Mac_AppleSilicon_CMMC_Guide.md](docs/Mac_AppleSilicon_CMMC_Guide.md)).
`build_gnu.sh` builds the model with the GNU compilers outside a container.

The `Dockerfile` builds the model (Intel ifx, Open MPI, NetCDF) and runs one
fire case per container. Build once (the common_ASP submodule must be checked
out):

    git submodule update --init
    docker build -t sam-asp-tuv .            # --build-arg NSUBDOMAINS_X=8 for 8 MPI ranks

Run a case, with the ASP input decks and an output directory mounted:

    docker run --rm --ulimit stack=-1 \
      -v /path/to/my_decks:/decks  -v /path/to/results:/output \
      -e ASP_DECK_DIR=/decks \
      -e LAT=39.70 -e LON=-122.58 -e DATE=2009-11-17 -e START_TIME=20:00 \
      -e RUN_HOURS=4 -e CDSAPI_KEY=<your CDS key> \
      sam-asp-tuv

| Setting | Variable (or option) | |
|---|---|---|
| fire latitude, longitude | `LAT`, `LON` (`--lat`, `--lon`) | required |
| start date and time, UTC | `DATE` (YYYY-MM-DD), `START_TIME` (HH:MM) | required |
| run length (hours) | `RUN_HOURS` (`--hours`) | required |
| ASP input decks | `ASP_DECK_DIR` (`--decks`), a mounted directory | default: decks in the image |
| fire emission flux | `MASSEMIS` | default 4.86e-6 |
| fire start (spin-up) step | `FIRE_START_STEP` | default 1080 |
| time step (s), case name | `DT` (8), `CASE` (fire) | |
| ERA5 | `CDSAPI_KEY`, or mount `/root/.cdsapirc`, or `ERA5_PL_FILE` + `ERA5_SL_FILE` (mounted files) | |
| extra SAM namelist settings | `SAM_SET="NAME=VALUE;..."` (`--set NAME=VALUE`) | |

The deck directory must hold a full SAM-ASP deck set (`ASPInputDeck.in`,
`GasPhaseChems.in`, `GasChemicalMechanism.in`, `PhotoRates.in`, `PhotoMap.in`,
the aerosol mode files, ...), like `SRC_TUV/ASP/InorganicInputDecks/`.
Results in the output directory: SAM's `OUT_STAT`, `OUT_3D`, `OUT_2D`, ...;
`ASP/` (ASP output files); `case/` (the case's input files and
`met_source.txt`); `ERA5_cache/`; `run.out` (model log). The run uses the MPI
rank count the image was built with (16 by default), oversubscribing the cores
if there are fewer; fire emissions start after `FIRE_START_STEP` steps
(spin-up), so `RUN_HOURS` must cover that, and the run warns if it does not.
