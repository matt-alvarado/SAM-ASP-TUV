# Changes since SAM-ASP-TUV v1.0.0

Changes that alter results or how the model is built and run, newest group
last. Details are in the git history; the current design is described in
[SAM_ASP_Interface.md](SAM_ASP_Interface.md).

## Coupling (SAM ↔ ASP ↔ TUV)

- Tracer layout, names and gas indices come from ASP (`samasp_layout.f90`);
  start-up stops with the required `ntracers`. Bugs this fixed: number
  concentrations were stored at 1476+bin and overlapped the last two
  aqueous-organic species of bin 10 (the mechanism has 618 gases, not 616);
  number-concentration names and 3-D outputs were offset by 2; aqueous tracer
  names carried the wrong bin number, so the 3-D aqueous totals used the
  wrong bins. 3-D gas position 167 (meant to be PAN) had become H13CO2CHO.
- Photolysis by TUV label (`PhotoMap.in`, `samasp_photo.f90`), replacing the
  `sammap` table and its out-of-bounds writes; output byte-identical.
- TUV site, date and vertical grid from SAM (before: fixed 2019-08-03 at
  47.9 N, and a hard-wired 40 m / 1 km grid that mapped SAM's top 5 layers
  to 3–8 km). SAM's solar zenith angle uses longitude. Results change.
- common_ASP (`asp_v2.2.2`) as the submodule `SRC_TUV/ASP/src`; bundled ASP
  copies and prebuilt ifort ASP/TUV libraries removed. Not carried over
  from SAM's old ASP copy: its electrolyte-equilibrium edits and the
  downgraded coagulation volume-loss check.
- ASP time step = SAM's substep (was 10 s); CO threshold, fire box and start
  step, deck and output directories are namelist settings. Fire emissions
  on rank 7 wrote outside the tracer array; the box is now in global indices.
- Initial plume CO for the particle-number conversion comes from
  `GasPhaseChems.in` (10000 ppb; was a hard-coded 3400 ppb): emitted POA and
  particle numbers scale by 0.34.
- `day0`/`day` follow SAM's CAM radiation convention (1.0 = 1 January 00 UTC).
- Photolysis is always online (deck `PHOTO FLAG` overridden).
- 3-D field count follows the number of size bins; IVOC4Gas–IVOC8Gas were
  declared (nz,ny,nzm) instead of (nx,ny,nzm), so POA4Gas and
  SOA1Gas–SOA4Gas 3-D output was garbage (model fields unaffected).
- SAM statistic `MCR` (core mass flux) renamed `MCRCORE`: it clashed with the
  ASP tracer `MCR` (methacrolein), and both profiles held the sum. Duplicate
  statistics names now stop the run.
- End-of-run move of the 2-D output used a fixed rank count of 12; it now
  uses `nsubdomains` and `caseid`, and only moves into an existing
  `COMPLETE_RUNS_AEROSOL/`.

## Performance (asptest 579 s → 52 s, same physics)

- TUV fixed set-up done once (`tuv_setup`); TUV compiled with precise FP.
- common_ASP: 8-slot DPREP ordering cache; flat-array gas-chemistry
  derivative and Jacobian evaluators (byte-identical).
- ASP load balancing across MPI ranks (`samasp_balance.f90`) with a per-cell
  state reset (`SetASPCellReset`), which changes MPI results once relative to
  warm-started cells.
- SAM compiled with precise FP by default (`SAM_FP`): results change once at
  rounding level; ~3% slower.
- Hash index for statistics names; TUV work arrays sized for 4 streams;
  per-call log messages only with `samasp_verbose`; TUV column reuse.

## Cases, containers and builds

- `create_run.py` takes the meteorology from ERA5 (Copernicus CDS) instead of
  NARR/NCEP R1; surface stress uses the actual air density (was 1.2 kg/m³);
  SAM's axes and `vg` use the vector-mean boundary-layer wind. The CDS
  returns single-level data as a zip of netCDF files; `era5.py` merges it.
- Containers: `Dockerfile` (Ubuntu, ifx), `Dockerfile.ubi9` (RHEL 9, ifx),
  `Dockerfile.ubi9-gnu` (RHEL 9, gfortran, x86-64 and ARM64); entry point
  `docker/run_case.py` with a tracer pre-check; fire cases in `cases/`.
- GNU build (`build.sh gnu`): plume means agree with ifx within ~0.5%, not
  bit-identical.
- `UTIL/` tools built from source only (`UTIL/build_util.sh`). `stat2nc`
  reads up to 100,000 profiles with 20-character names (was 1,000 and 8);
  `bin2D2nc`/`com2D2nc` overflowed `varid(99)` with SAM-ASP's 143 fields.
- Python scripts converted to Python 3.
