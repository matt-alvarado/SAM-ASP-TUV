# SAM-ASP-TUV container: user's guide

This guide covers building the SAM-ASP-TUV container image, running fire cases
with it on any Linux machine with Docker, and running it on Amazon Web
Services (AWS). For an on-prem RHEL virtual machine in a CMMC Level 2 (CUI)
enclave, with rootless Podman, see also the
[on-prem guide](OnPrem_RHEL_CMMC_Guide.md).

SAM-ASP-TUV couples three models:

- **SAM** (System for Atmospheric Modeling): large-eddy dynamics of the smoke plume
- **ASP** (Aerosol Simulation Program): gas chemistry and aerosol in the plume
- **TUV** (Tropospheric Ultraviolet and Visible model): photolysis rates, including the plume's own aerosol

One container runs one fire case. You give the fire's location, date, start
time and run length, and (optionally) your own ASP input decks. The container:

1. gets the background meteorology from the ERA5 reanalysis;
2. writes the SAM case files;
3. runs the model with MPI;
4. leaves all results in a directory you mount.

Contents:
1. [What is in the image](#1-what-is-in-the-image)
2. [Before you start](#2-before-you-start)
3. [Building the image](#3-building-the-image)
4. [Running a case](#4-running-a-case)
5. [Settings reference](#5-settings-reference)
6. [ASP input decks](#6-asp-input-decks)
7. [Meteorology (ERA5)](#7-meteorology-era5)
8. [Choosing the run length: spin-up and fire timing](#8-choosing-the-run-length-spin-up-and-fire-timing)
9. [Outputs](#9-outputs)
10. [Performance and machine size](#10-performance-and-machine-size)
11. [Running on AWS](#11-running-on-aws)
12. [Troubleshooting](#12-troubleshooting)
13. [Reproducibility](#13-reproducibility)

---

## 1. What is in the image

| Path in the image | What |
|---|---|
| `/opt/sam-asp-tuv` | The SAM-ASP-TUV source tree, the compiled model `SAM_RAD_CAM_MICRO_M2005_w_aero`, TUV data, the default ASP input decks (`SRC_TUV/ASP/InorganicInputDecks/`) and the post-processing tools (`UTIL/`). |
| `/opt/sam-asp-tuv/docker/run_case.py` | The entry point: it builds the case and runs it. |
| `/opt/sam-asp-tuv/create_run.py` | The case builder: ERA5 meteorology, then `snd`, `lsf`, `sfc` and `prm`. |
| `/opt/intel-ifx/lib`, `/opt/netcdf-ifx/lib` | Only the Intel Fortran runtime and NetCDF-Fortran libraries the model binary needs. |
| System packages | Ubuntu 24.04, Open MPI, NetCDF-C/HDF5, Python 3 with numpy, netCDF4 and cdsapi. |

The image is about 1.1 GB on disk. The model is compiled for a fixed number
of MPI ranks: 16 by default. See [3.2](#32-build-options).

The model grid is 240 × 1 × 80 points (a 2-D slice through the plume),
with 500 m horizontal spacing (`dx`, `dy` in the `prm` template). The domain
moves with the boundary-layer wind, following the plume.

---

## 2. Before you start

You need:

- **A Linux machine with Docker** (Docker Engine 23 or newer, for BuildKit
  and `docker build --secret`). Any x86-64 machine works. For CPU and memory,
  see [section 10](#10-performance-and-machine-size).
- **About 10 GB of free disk** to build the image, plus space for outputs.
  A few hours of model time writes a few GB.
- **The SAM-ASP-TUV repository** with the common_ASP submodule (the ASP
  source), to build the image.
- **A Copernicus Climate Data Store (CDS) account**, for ERA5 meteorology.
  This is free; set it up as below. Alternatively, supply ERA5 files you
  downloaded yourself (see [7.3](#73-using-era5-files-you-already-have)).

### 2.1 Setting up ERA5 access (once)

1. Register at <https://cds.climate.copernicus.eu> and log in.
2. Accept the ERA5 licence on both dataset pages. Open each one, go to its
   *Download* tab, and accept the licence near the bottom:
   - "ERA5 hourly data on pressure levels from 1940 to present"
     (`reanalysis-era5-pressure-levels`)
   - "ERA5 hourly data on single levels from 1940 to present"
     (`reanalysis-era5-single-levels`)

   Without this, downloads fail with a "required licences not accepted" error.
3. Copy your **personal access token (API key)** from your CDS user profile
   page.
4. Either pass the key to the container as `CDSAPI_KEY=<token>`, or write a
   file `~/.cdsapirc` and mount it into the container (see
   [7.1](#71-how-the-download-works)):

   ```
   url: https://cds.climate.copernicus.eu/api
   key: <your personal access token>
   ```
   Protect it with `chmod 600 ~/.cdsapirc`.

Treat the key like a password. Don't put it in scripts you share or commit.

---

## 3. Building the image

### 3.1 Basic build

```bash
git clone https://github.com/matt-alvarado/SAM-ASP-TUV.git
cd SAM-ASP-TUV
git checkout asp_v2.2.2
git submodule update --init          # fetches the ASP source (common_ASP)
docker build -t sam-asp-tuv .
```

The build takes roughly 5–10 minutes on a 4-core machine:

- downloading the Intel Fortran compiler from conda-forge, about 1 minute;
- building NetCDF-Fortran with it, about 30 s;
- compiling ASP, TUV and SAM, about 2–4 minutes.

The compiler (about 4 GB) is used only in the first build stage. The final
image keeps only what the model needs at run time.

If the build stops with
`ASP source missing: run "git submodule update --init"`, the submodule wasn't
checked out.

### 3.2 Build options

| Option | Meaning |
|---|---|
| `--build-arg NSUBDOMAINS_X=8` | Number of MPI ranks (subdomains in x). It must divide 240, so 4, 6, 8, 10, 12, 16, 20, 24, … work; the default is 16. Choose it to match the physical cores of the machine you'll run on ([section 10](#10-performance-and-machine-size)). |
| `--secret id=ca,src=/path/to/ca.crt` | An extra CA certificate, for building behind a corporate proxy that intercepts TLS. It's used only during the build and not kept in the image. |
| `--build-arg HTTPS_PROXY=http://proxy:port` | Standard Docker proxy settings, if your network needs them. |

It's useful to tag the image with the rank count and the code version, for example:

```bash
docker build --build-arg NSUBDOMAINS_X=16 \
  -t sam-asp-tuv:16r-$(git rev-parse --short HEAD) .
```

### 3.3 Checking the build

```bash
docker run --rm sam-asp-tuv --help           # lists all settings
docker run --rm --entrypoint cat sam-asp-tuv /opt/sam-asp-tuv/NRANKS   # MPI ranks
```

---

## 4. Running a case

### 4.1 Quick start

```bash
mkdir -p ~/samasp/results

docker run --rm --ulimit stack=-1 --shm-size=1g \
  -v ~/samasp/results:/output \
  -e LAT=39.70 -e LON=-122.58 \
  -e DATE=2009-11-17 -e START_TIME=20:00 \
  -e RUN_HOURS=4 \
  -e CDSAPI_KEY=<your CDS key> \
  sam-asp-tuv
```

This uses the ASP input decks built into the image. With your own decks,
mount them and set `ASP_DECK_DIR`:

```bash
docker run --rm --ulimit stack=-1 --shm-size=1g \
  -v ~/samasp/my_decks:/decks:ro \
  -v ~/samasp/results:/output \
  -e ASP_DECK_DIR=/decks \
  -e LAT=39.70 -e LON=-122.58 -e DATE=2009-11-17 -e START_TIME=20:00 \
  -e RUN_HOURS=4 -e CDSAPI_KEY=<your CDS key> \
  sam-asp-tuv
```

The Docker flags:

- `--ulimit stack=-1`: TUV and ASP use large local arrays and need an
  unlimited stack. The entry point tries to lift the limit itself and warns
  if it can't.
- `--shm-size=1g`: Open MPI ranks exchange data through shared memory. 1 GB
  avoids running out with 16 ranks.
- `-v host_dir:/output`: where the results go. Without this mount the
  results stay inside the container and are lost when it's removed (`--rm`).
- `:ro` on the deck mount: the model only reads the decks.

### 4.2 What you'll see

```
fire: ERA5 39.75N -122.50E 2009-11-17 20:00 UTC; psfc 988.8 hPa, PBL 850 m, H 244.6, LE 43.3 W/m2, vg 2.81 m/s -> /opt/sam-asp-tuv/RUNS/fire
run: 1800 steps of 8 s; fire emissions from step 1080 to ~1125; 16 MPI ranks
SAM finished at step 1800 of 1800 (exit code 9); log: /output/run.out
```

- **Line 1:** the ERA5 grid point and hour used, and the key meteorology:
  surface pressure, boundary-layer depth, sensible and latent heat flux, and
  `vg`, the boundary-layer wind that carries the plume.
- **Line 2:** the run's step count, when the fire emits, and the MPI rank
  count. Warnings appear here if the run is too short; see
  [section 8](#8-choosing-the-run-length-spin-up-and-fire-timing).
- **Line 3:** the result. SAM's exit code 9 means a normal finish. The
  container exits with 0 if the run reached its last step.

Follow progress from another terminal with
`tail -f ~/samasp/results/run.out` (look for `NSTEP =` lines).

### 4.3 Running several cases

Use a separate output directory and case name for each run:

```bash
while IFS=, read -r name lat lon date start hours; do
  mkdir -p ~/samasp/results/$name
  docker run --rm --ulimit stack=-1 --shm-size=1g \
    -v ~/samasp/results/$name:/output \
    -e CASE=$name -e LAT=$lat -e LON=$lon -e DATE=$date -e START_TIME=$start \
    -e RUN_HOURS=$hours --env-file ~/samasp/cds.env \
    sam-asp-tuv > ~/samasp/results/$name.log 2>&1
done < cases.csv
```

Here `cases.csv` holds lines like `williams,39.70,-122.58,2009-11-17,20:00,4`,
and `~/samasp/cds.env` holds the line `CDSAPI_KEY=<key>` (mode 600). Run
cases one after another on one machine: each run already uses all its cores.

---

## 5. Settings reference

Every setting can be given as an environment variable (`-e NAME=value`) or as
a command-line option after the image name (`sam-asp-tuv --lat 39.7 ...`).
Options win over variables.

| Variable | Option | Meaning | Default |
|---|---|---|---|
| `LAT` | `--lat` | Fire latitude, degrees north | **required** |
| `LON` | `--lon` | Fire longitude, degrees east (−180…180 or 0…360) | **required** |
| `DATE` | `--date` | Start date, UTC, `YYYY-MM-DD` | **required** |
| `START_TIME` | `--start` | Start time, UTC, `HH:MM` (ERA5 is hourly; the nearest hour is used) | **required** |
| `RUN_HOURS` | `--hours` | Run length in hours, including spin-up (see [section 8](#8-choosing-the-run-length-spin-up-and-fire-timing)) | **required** |
| `ASP_DECK_DIR` | `--decks` | Directory, inside the container, with the ASP input decks; mount it with `-v` | the image's `SRC_TUV/ASP/InorganicInputDecks` |
| `MASSEMIS` | `--massemis` | Fire mass emission flux (namelist `massemis`) | `4.86e-6` |
| `FIRE_START_STEP` | `--fire-start-step` | Fire emissions start after this time step (the spin-up) | `1080` (from the `prm` template) |
| `DT` | `--dt` | SAM time step, seconds | `8` |
| `CASE` | `--case` | Case name (used in output file names) | `fire` |
| `OUTPUT_DIR` | `--output` | Output directory inside the container | `/output` |
| `CDSAPI_KEY` | | Your CDS API key, for downloading ERA5 | |
| `ERA5_PL_FILE`, `ERA5_SL_FILE` | `--era5-pl`, `--era5-sl` | ERA5 pressure-level and single-level netCDF files to use instead of downloading | |

Other model settings, such as the fire box, the CO threshold for running
ASP, output intervals, and SAM physics switches, come from the namelist
template `RUNDATA/prm_template`. To change them for all runs, edit the
template and rebuild the image. For a single run, mount your own template
and use the case builder directly; see [12](#12-troubleshooting), "changing
other namelist settings".

---

## 6. ASP input decks

`ASP_DECK_DIR` must point to a complete SAM-ASP deck set: a single directory
with all the ASP input files. The entry point checks for these before starting:

| File | Contents |
|---|---|
| `ASPInputDeck.in` | ASP run settings |
| `GasPhaseChems.in` | Gas species: molecular mass, **fire emission factors (g/kg burned)**, background mixing ratios, deposition |
| `GasChemicalMechanism.in` | Gas-phase reactions |
| `PhotoRates.in` | Photolysis rate table layout (number of rate columns) |
| `PhotoMap.in` | How TUV's photolysis rates map to ASP's rate columns (SAM-ASP-TUV specific) |
| `AerosolModes.in` | Fire (plume) aerosol: size distribution and composition |
| `EnvAerosolModes.in` | Background aerosol |

ASP also reads the other files of a deck set: `AqPhaseChems.in`,
`AqPhaseIons.in`, `AqEquilibriumReactions.in`, `Dissolution.in`,
`OrgPhaseChems.in`, `HydrophilicOrgChems.in`, `OrganicDissolution.in`,
`HydrophilicOrganicDissolution.in`, `AerosolThermoKusikMeissner.in`, and
others the deck refers to. These must be in the same directory.

The simplest way to make your own decks is to start from the image's set and
change what you need. For a new fire, that's usually the emission factors in
`GasPhaseChems.in` and the composition and size in `AerosolModes.in`:

```bash
docker run --rm --entrypoint tar sam-asp-tuv \
  -C /opt/sam-asp-tuv/SRC_TUV/ASP -cf - InorganicInputDecks | tar -xf - -C ~/samasp
mv ~/samasp/InorganicInputDecks ~/samasp/my_decks
# edit ~/samasp/my_decks/*.in, then run with -v ~/samasp/my_decks:/decks:ro -e ASP_DECK_DIR=/decks
```

Notes:

- **Tracer count.** SAM allocates one tracer per gas species and per aerosol
  species and size bin. The image is compiled for the default deck's tracer
  count, and the model stops at start-up with the required value if your
  decks need a different number. To use decks with a different species list,
  set `ntracers` in `SRC_TUV/domain.f90` to that value and rebuild the image.
- **Mechanisms.** The SAM-ASP-TUV interface reads species, bins and
  photolysis mapping from the decks at run time. Decks for other mechanisms
  (the ASP repository's CB6r5, SAPRC07, CRACMM, GEOS-Chem, MCM decks) also need
  a matching `PhotoMap.in` and the right tracer count.
- `NitrateGroups.in` (used by the standalone ASP parcel model) isn't needed
  for SAM-ASP runs.

---

## 7. Meteorology (ERA5)

### 7.1 How the download works

For each run, `create_run.py` asks the CDS for a **3 × 3 box of ERA5 grid
points (0.25°) around the fire, for one hour**: the hour nearest the start
time. That's two small requests, a few kB each:

- `reanalysis-era5-pressure-levels`: temperature, specific humidity, u and v
  wind on 27 levels from 1000 to 100 hPa;
- `reanalysis-era5-single-levels`: surface pressure, boundary-layer height,
  surface sensible and latent heat flux, eastward and northward surface
  stress, and 2 m temperature.

CDS requests are queued, so a download takes from a few seconds to several
minutes depending on CDS load. Downloaded files are cached in
`/output/ERA5_cache/`, so rerunning a case with the same output directory
doesn't download again.

To use a `.cdsapirc` file instead of `CDSAPI_KEY`, mount it:
`-v ~/.cdsapirc:/root/.cdsapirc:ro`.

### 7.2 What is written

From the grid point nearest the fire, the case builder writes:

- **`snd`:** the sounding: potential temperature, water vapour mixing ratio,
  and winds rotated so SAM's y axis follows the mean boundary-layer wind.
- **`lsf`:** the same winds as large-scale forcing (no large-scale tendencies).
- **`sfc`:** sensible and latent heat flux (positive upward, the ERA5 hour
  ending at the start time) and the surface stress u*².
- **`vg`** in `prm`: the mean boundary-layer wind along y. It sets how fast
  the domain moves with the plume and how long the fire emits.

The meteorology is held **constant** during the run. The boundary layer is
taken as the ERA5 levels above the ground and below the boundary-layer top
(or the lowest 50 hPa, whichever is deeper).

`case/met_source.txt` in the output records the grid point, valid time and
values used.

### 7.3 Using ERA5 files you already have

If you have ERA5 netCDF files for the fire, or the machine can't reach the
CDS, mount them and set both file variables. They must contain the variables
above, for an area and hour that includes the fire:

```bash
docker run --rm --ulimit stack=-1 --shm-size=1g \
  -v ~/samasp/era5:/era5:ro -v ~/samasp/results:/output \
  -e ERA5_PL_FILE=/era5/pl.nc -e ERA5_SL_FILE=/era5/sl.nc \
  -e LAT=39.70 -e LON=-122.58 -e DATE=2009-11-17 -e START_TIME=20:00 -e RUN_HOURS=4 \
  sam-asp-tuv
```

The run stops with a clear message if the files don't cover the fire's
location or hour. To download such files on another machine, the requests
are in `samasp_case/era5.py` (`requests()`), or use the CDS web form with
the same variables, levels, hour and area, in netCDF format.

---

## 8. Choosing the run length: spin-up and fire timing

A run has three phases:

| Phase | Steps | What happens |
|---|---|---|
| Spin-up | 1 … `FIRE_START_STEP` | Turbulence develops in the boundary layer; no fire. |
| Emission | `FIRE_START_STEP` + 1 … ≈ `FIRE_START_STEP` + 2·(plume_y − dy)/vg/dt | The fire emits into the fire box. With the template (`plume_y` = 1000 m, `dy` = 500 m) this is 1000/(vg·dt) steps: about 45 steps at vg = 2.8 m/s and dt = 8 s. |
| Chemistry | after the emission | ASP runs in every grid cell where CO exceeds 130 ppb. |

So `RUN_HOURS` must cover the spin-up, the emission and the plume age you
want to study:

```
RUN_HOURS ≈ (FIRE_START_STEP + emission steps) × DT / 3600 + plume age (h)
```

With the defaults (1080 steps × 8 s = 2.4 h of spin-up), a 2-hour-old plume
needs about `RUN_HOURS=4.5`. The entry point prints the computed fire window
and warns if the run ends before the fire starts or before its emissions
end.

For quick tests you can shorten the spin-up, for example
`-e FIRE_START_STEP=0 -e RUN_HOURS=0.2`. This gives a few minutes of plume
in unspun-up turbulence: fine for checking a setup, not for science.

---

## 9. Outputs

Everything goes to the mounted output directory:

| Path | Contents |
|---|---|
| `run.out` | The model log (SAM, TUV, ASP). Check here first if something goes wrong. |
| `OUT_STAT/<case>_<case>.stat` | SAM statistics: horizontal-mean profiles and time series of every field and tracer, every `nstat` steps (binary). |
| `OUT_3D/<case>_<case>_<n>.bin2D` | Snapshots of the 2-D (x–z) fields, including every gas tracer and the J-values, every `nsave3D` steps (binary). |
| `OUT_2D`, `OUT_MOMENTS`, `OUT_MOVIES`, `RESTART` | Other SAM output, as configured in `prm` (mostly empty with the template). |
| `ASP/` | ASP's own output: `Transcript.TXT` (ASP log), `GasPhaseODESet.txt` and `GasPhaseJacobian.txt` (the mechanism as ASP read it). |
| `case/` | The case's input files: `prm`, `snd`, `lsf`, `sfc`, `grd`, `lst`, and `met_source.txt`. |
| `ERA5_cache/` | The downloaded ERA5 files. |

### 9.1 Converting outputs to netCDF

The image includes SAM's conversion tools in `/opt/sam-asp-tuv/UTIL`
(`bin2D2nc`, `stat2nc`, `bin3D2nc`, `com3D2nc`, `2Dbin2nc`, and the rest).
They are built from source (`UTIL/SRC`) with the model; the repository keeps
no prebuilt binaries. Outside the container, build them with
`UTIL/build_util.sh` (all tools), or for example `UTIL/build_util.sh
stat2nc`; it uses the same environment file as `build_ifx.sh`.

To convert the 2-D field snapshots:

```bash
docker run --rm -v ~/samasp/results/OUT_3D:/data -w /data \
  --entrypoint /opt/sam-asp-tuv/UTIL/bin2D2nc sam-asp-tuv fire_fire_16.bin2D
# writes fire_fire_16.bin2D_1.nc: every field on (time, z, x)
```

To convert the statistics file (horizontal-mean profiles and time series):

```bash
docker run --rm -v ~/samasp/results/OUT_STAT:/data -w /data \
  --entrypoint /opt/sam-asp-tuv/UTIL/stat2nc sam-asp-tuv fire_fire.stat
# writes fire_fire.nc: every profile on (time, z), plus the time series
```

`stat2nc` handles SAM-ASP-TUV files with thousands of profiles. Notes on
its output:

- netCDF names can't start with `(`. Species such as `(NH4)2SO401` are
  written as `_(NH4)2SO401`, and their original name is kept in the
  attribute `sam_name`.
- The file is netCDF-3 with 64-bit offsets, so it can be larger than 2 GB.

---

## 10. Performance and machine size

- **CPU.** The model runs `NRANKS` MPI processes (16 by default). Give each
  rank its own **physical core**. Cloud "vCPUs" are usually hyper-threads,
  two per core, so 16 ranks want 16 cores (32 vCPUs, or 16 vCPUs with
  hyper-threading turned off). With fewer cores the ranks share them, which
  works but is slower. Match the build option `NSUBDOMAINS_X` to the cores
  you'll have.
- **Memory.** About 0.7 GB per rank, so 11–12 GB for 16 ranks. 32 GB gives
  comfortable headroom.
- **Disk.** Outputs grow with run length and output frequency: roughly 1–2
  GB per model hour with the template's settings.
- **Run time.** Most of the cost is ASP chemistry in the plume, and it grows
  as the plume spreads. As a guide: a 72-step test, with ASP active for the
  last ~25 steps, took 2.7 minutes with 16 ranks sharing only 4 cores. On 16
  dedicated cores expect roughly 3–5× faster, but measure a short run on
  your machine before planning long ones.

---

## 11. Running on AWS

The simplest way is a single EC2 virtual machine running Docker, as in
section 4. It's described step by step below, followed by options for
storing images, keys and data on AWS.

The commands use the AWS CLI (`aws`), version 2, with your account
configured. Replace the region, names and IDs with yours. AWS products,
instance types and prices change, so check the current AWS documentation
for anything that doesn't match.

### 11.1 Choose an instance

| Use | Suggested instance | Notes |
|---|---|---|
| 16-rank runs (default image) | **c7i.8xlarge**: 32 vCPUs = 16 physical cores, 64 GiB | Launch with hyper-threading off (1 thread per core): 16 vCPUs on 16 cores, same price. |
| Smaller or cheaper runs, or testing | c7i.4xlarge: 16 vCPUs = 8 cores, 32 GiB | Build with `NSUBDOMAINS_X=8`, or run 16 ranks oversubscribed. |
| Bigger domains or more ranks | c7i.16xlarge (32 cores) or larger | Build with a matching `NSUBDOMAINS_X`; it must divide 240. |

Compute-optimised AMD instances (c7a) work too: the Intel Fortran runtime
runs on AMD processors. Spot instances cost much less but can be stopped by
AWS at short notice. Use them for runs you can repeat.

**AWS WorkSpaces** (virtual desktops) aren't a good fit for model runs: they
offer few cores, and Docker isn't set up by default. Use one, if you like,
as the desktop you work from, and run the model on EC2.

### 11.2 One-time AWS setup

1. **Key pair**, to log in with SSH:
   ```bash
   aws ec2 create-key-pair --key-name samasp --query KeyMaterial --output text > samasp.pem
   chmod 400 samasp.pem
   ```
2. **Security group** allowing SSH only from your address:
   ```bash
   aws ec2 create-security-group --group-name samasp-ssh --description "SSH for SAM-ASP"
   aws ec2 authorize-security-group-ingress --group-name samasp-ssh \
     --protocol tcp --port 22 --cidr $(curl -s https://checkip.amazonaws.com)/32
   ```
3. **S3 bucket** for decks and results (bucket names are global, so pick a unique one):
   ```bash
   aws s3 mb s3://my-samasp-bucket
   ```
4. **The CDS key in AWS Systems Manager Parameter Store** (encrypted), so it
   never sits in scripts:
   ```bash
   aws ssm put-parameter --name /samasp/cdsapi_key --type SecureString --value '<your CDS key>'
   ```
5. **An IAM role for the instance** (an "instance profile"), so the machine
   can read the key and use the bucket without stored credentials. In the IAM
   console, create a role for EC2 with:
   - `ssm:GetParameter` on `arn:aws:ssm:<region>:<account>:parameter/samasp/*`
     (plus `kms:Decrypt` for the default key);
   - read and write on `arn:aws:s3:::my-samasp-bucket` and
     `arn:aws:s3:::my-samasp-bucket/*`;
   - if you use ECR ([11.5](#115-storing-the-image-in-amazon-ecr-optional)), the
     managed policy `AmazonEC2ContainerRegistryReadOnly`.

   Call the role `samasp-ec2`.

### 11.3 Launch the instance

Use an Ubuntu 24.04 LTS image (AMI) for x86-64. Look up the current AMI ID
for your region in the EC2 console, or with:

```bash
aws ssm get-parameter --query Parameter.Value --output text \
  --name /aws/service/canonical/ubuntu/server/24.04/stable/current/amd64/hvm/ebs-gp3/ami-id
```

Then launch, here a c7i.8xlarge with hyper-threading off, a 100 GB disk and
the role attached:

```bash
aws ec2 run-instances \
  --image-id <ami-id> --instance-type c7i.8xlarge \
  --cpu-options CoreCount=16,ThreadsPerCore=1 \
  --key-name samasp --security-groups samasp-ssh \
  --iam-instance-profile Name=samasp-ec2 \
  --block-device-mappings 'DeviceName=/dev/sda1,Ebs={VolumeSize=100,VolumeType=gp3}' \
  --tag-specifications 'ResourceType=instance,Tags=[{Key=Name,Value=samasp}]'
```

Find its public address, and log in once it's running:

```bash
aws ec2 describe-instances --filters Name=tag:Name,Values=samasp \
  --query 'Reservations[].Instances[].PublicDnsName' --output text
ssh -i samasp.pem ubuntu@<public-dns-name>
```

### 11.4 Set up the instance

On the instance:

```bash
sudo apt-get update && sudo apt-get install -y docker.io git awscli
sudo usermod -aG docker ubuntu     # log out and back in for this to take effect
exit
```

Log in again, then build the image (section 3), choosing the rank count for
this instance:

```bash
git clone https://github.com/matt-alvarado/SAM-ASP-TUV.git
cd SAM-ASP-TUV && git checkout asp_v2.2.2 && git submodule update --init
docker build --build-arg NSUBDOMAINS_X=16 -t sam-asp-tuv .
```

The repositories may be private. In that case, clone over HTTPS with a
GitHub personal access token, or build once and use ECR
([11.5](#115-storing-the-image-in-amazon-ecr-optional)) so later instances
just pull the image.

### 11.5 Storing the image in Amazon ECR (optional)

Push the built image once, then pull it on any instance without rebuilding:

```bash
REGION=us-east-1; ACCOUNT=$(aws sts get-caller-identity --query Account --output text)
REPO=$ACCOUNT.dkr.ecr.$REGION.amazonaws.com/sam-asp-tuv
aws ecr create-repository --repository-name sam-asp-tuv --region $REGION
aws ecr get-login-password --region $REGION | docker login --username AWS --password-stdin $ACCOUNT.dkr.ecr.$REGION.amazonaws.com
docker tag sam-asp-tuv $REPO:16r
docker push $REPO:16r
# on another instance:
docker pull $REPO:16r
```

### 11.6 Run a case on the instance

Copy your decks from S3, run, and copy the results back:

```bash
aws s3 sync s3://my-samasp-bucket/decks/williams ~/samasp/decks
mkdir -p ~/samasp/results/williams

CDSAPI_KEY=$(aws ssm get-parameter --name /samasp/cdsapi_key --with-decryption \
             --query Parameter.Value --output text)

docker run --rm --ulimit stack=-1 --shm-size=1g \
  -v ~/samasp/decks:/decks:ro -v ~/samasp/results/williams:/output \
  -e ASP_DECK_DIR=/decks -e CASE=williams \
  -e LAT=39.70 -e LON=-122.58 -e DATE=2009-11-17 -e START_TIME=20:00 \
  -e RUN_HOURS=4.5 -e CDSAPI_KEY="$CDSAPI_KEY" \
  sam-asp-tuv 2>&1 | tee ~/samasp/results/williams.log

aws s3 sync ~/samasp/results/williams s3://my-samasp-bucket/results/williams
```

For long runs, start the run inside `tmux` (or with `nohup ... &`), so it
keeps going if your SSH connection drops:

```bash
tmux new -s run      # run the docker command above; detach with Ctrl-b d, reattach with: tmux attach -t run
```

### 11.7 Running without logging in (optional)

To run a case and shut the machine down automatically, put the steps in a
script and pass it as **user data** when launching. It runs once, as root,
at first boot. A minimal example, using an image in ECR:

```bash
#!/bin/bash
set -e
apt-get update && apt-get install -y docker.io awscli
REGION=us-east-1; ACCOUNT=<account>; REPO=$ACCOUNT.dkr.ecr.$REGION.amazonaws.com/sam-asp-tuv:16r
aws ecr get-login-password --region $REGION | docker login --username AWS --password-stdin $ACCOUNT.dkr.ecr.$REGION.amazonaws.com
docker pull $REPO
mkdir -p /data/decks /data/out
aws s3 sync s3://my-samasp-bucket/decks/williams /data/decks
KEY=$(aws ssm get-parameter --region $REGION --name /samasp/cdsapi_key --with-decryption --query Parameter.Value --output text)
docker run --rm --ulimit stack=-1 --shm-size=1g -v /data/decks:/decks:ro -v /data/out:/output \
  -e ASP_DECK_DIR=/decks -e CASE=williams -e LAT=39.70 -e LON=-122.58 \
  -e DATE=2009-11-17 -e START_TIME=20:00 -e RUN_HOURS=4.5 -e CDSAPI_KEY="$KEY" \
  $REPO > /data/out/container.log 2>&1 || true
aws s3 sync /data/out s3://my-samasp-bucket/results/williams
shutdown -h now
```

Launch with `--user-data file://run_williams.sh` and
`--instance-initiated-shutdown-behavior terminate`, so the instance is
deleted, and stops costing money, when the run ends.

For many runs, **AWS Batch** can manage the machines for you:

- a compute environment with your instance type;
- a job queue;
- a job definition using the ECR image, with the environment variables as
  job parameters, `ulimits` stack −1 and `sharedMemorySize`.

Batch jobs don't get host directories, so the decks must be in the image (a
small Dockerfile `FROM sam-asp-tuv` with `COPY my_decks /decks`), and the
results must be copied to S3 inside the job (a wrapper entry point that
calls `run_case.py` and then `aws s3 sync`). This is more setup than the
EC2 route above; consider it once the EC2 workflow works.

### 11.8 Costs and cleanup

- A running instance costs money every hour, whether or not the model is
  running. **Stop** it between work sessions (`aws ec2 stop-instances
  --instance-ids <id>`): the disk and your built image are kept, at a small
  storage cost. **Terminate** it when you're done.
- S3 and ECR storage are charged per GB-month. Delete old results you no
  longer need.
- Set a billing alert (AWS Budgets) when you start.

---

## 12. Troubleshooting

| Symptom | Cause and fix |
|---|---|
| `ERROR: set lon/LON, ...` | A required setting is missing. |
| `ASP deck directory /decks lacks ...` | The deck directory isn't mounted, or it's the wrong directory. Check the `-v host:/decks` path and `ASP_DECK_DIR`. |
| `no ERA5 source` | Give `CDSAPI_KEY`, mount `/root/.cdsapirc`, or give both ERA5 files. |
| ERA5 download fails with 401 / "invalid key" | The CDS key is wrong or expired; copy it again from your CDS profile. |
| ERA5 download fails with "required licences not accepted" | Accept the ERA5 licence on both dataset pages ([2.1](#21-setting-up-era5-access-once)). |
| ERA5 download hangs | The CDS queue is busy; requests can wait several minutes. Check the queue status on the CDS website. |
| `... does not cover ...` or `has no data for ...` | The mounted ERA5 files don't include the fire's location or hour. |
| `WARNING: the run ends before the fire starts` | `RUN_HOURS` is shorter than the spin-up; see [section 8](#8-choosing-the-run-length-spin-up-and-fire-timing). |
| The run stops early; `run.out` shows `Segmentation fault` near the start | The stack limit is too small: add `--ulimit stack=-1`. |
| MPI errors about shared memory or `vader` | Add or increase `--shm-size` (e.g. `--shm-size=2g`). |
| `forrtl: error (78): process killed (SIGTERM)` at the very end of `run.out` | Harmless: Open MPI stops the other ranks after SAM's normal exit (code 9). The run is complete if it reached its last step. |
| `ntracers ... does not match` at start-up | Your decks need a different tracer count; see [section 6](#6-asp-input-decks). |
| "No space left on device" | Free disk, or mount `/output` on a larger volume. |
| Very slow run | More ranks than physical cores (check `NRANKS` and the machine), or ASP chemistry in a large plume. See [section 10](#10-performance-and-machine-size). |

**Changing other namelist settings for one run.** Mount an edited copy of
the template and call the case builder and model yourself in an interactive
container:

```bash
docker run --rm -it --ulimit stack=-1 --shm-size=1g \
  -v ~/samasp/results:/output --entrypoint bash sam-asp-tuv
# inside: python3 /opt/sam-asp-tuv/docker/run_case.py --help
# create_run.py accepts --set NAME=VALUE to set any &PARAMETERS entry, e.g.
#   --set asp_co_threshold_ppb=200.
```

**Inspecting a finished case.** `case/prm` shows the exact namelist used.
`case/met_source.txt` shows the meteorology. `run.out` and `ASP/Transcript.TXT`
hold the logs.

---

## 13. Reproducibility

To be able to repeat a run exactly, keep:

- **the image:** tag it with the code version (`git rev-parse --short HEAD`
  in the checkout you built from) and keep it, for example in ECR;
- **the ASP decks you mounted;**
- **`case/`** from the output, holding the namelist and the meteorology
  actually used, and **`ERA5_cache/`**, holding the ERA5 data;
- **the `docker run` command**, including every `-e` setting.

The model build uses strict floating-point settings (`-fp-model precise`)
for SAM, ASP and TUV, so the same image, inputs and rank count give the same
results run after run.
