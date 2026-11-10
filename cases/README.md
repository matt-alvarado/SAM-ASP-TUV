# Fire cases

Each folder holds what is specific to one fire: the run settings
(`case.env`) and, where available, the ASP input decks (`decks/`). The
meteorology comes from ERA5 when the case is run (see
`docs/Container_User_Guide.md`), so no SAM sounding or forcing files are kept.

| Folder | Fire | Decks | SAM tracers (`NTRACERS`) |
|---|---|---|---|
| `williams_flats_2019/` | Williams Flats, WA, 3 Aug 2019 (FIREX-AQ) | SAM-ASP-TUV v1.0.0 decks (the image's default decks) | 1488 (default image) |
| `williams_2009/` | Williams fire, CA, 17 Nov 2009 | ASP repository `Williams/` baseline deck, 10 size bins | 1486 |
| `bbop_colockum_2013/` | Colockum Tarps, WA, 9 Aug 2013 (BBOP) | none yet | depends on the decks |

The cases come from the SAM-ASP-TUV v1.0.0 release (Zenodo,
[10.5281/zenodo.13385240](https://doi.org/10.5281/zenodo.13385240)), whose
`RUNS/` folders are also in this repository: `SND1` (Williams Flats),
`W1101`, `WILLIAMS` and `tuvtest` (Williams 2009), and `ADTEST` and
`SAMASP0001` (BBOP). Each case's README lists the original settings.

## Running a case

The image contains `cases/`, so a case runs with its env file:

```bash
docker run --rm --ulimit stack=-1 --shm-size=1g \
  --env-file cases/williams_flats_2019/case.env \
  -v ~/.cdsapirc:/root/.cdsapirc:ro -v ~/samasp/runs/wf2019:/output \
  sam-asp-tuv
```

(`podman run` takes the same options.) Any setting can be overridden with a
command-line option, for example `--hours 1`.

**Defaults for every case.** All case decks use **10 aerosol size bins**
(sectional cutoffs 0.025 and 1.01 um, as in the v1.0.0 decks), and
**photolysis is computed online by TUV** at every step: SAM-ASP-TUV overrides
the decks' `PHOTO FLAG` and uses `PhotoRates.in` only for its layout. ASP
repository decks have 40 bins; with 40 bins the Williams deck fails in SAM
(sectional coagulation, "Total Volume concentration is lost"), so convert such
decks to 10 bins (the three sectional lines of `AerosolModes.in`) before use.

**Tracer count.** SAM is compiled for one deck layout: `ntracers` in
`SRC_TUV/domain.f90` must equal the number the decks need. The entry point
checks this first and stops with the required value before any work (also
`--check-decks` to check only); SAM checks again at start-up. The default image fits
`williams_flats_2019`. For `williams_2009`, build an image with
`--build-arg NTRACERS=1486`, for example tagged `sam-asp-tuv:nt1486`.

`SAM_SET` in an env file (or `--set NAME=VALUE`) passes extra settings to the
case's SAM namelist, such as the warm-bubble start of Williams Flats.
