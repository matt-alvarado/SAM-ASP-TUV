# The SAM ↔ ASP ↔ TUV interface

How SAM (the host model, `SRC_TUV/*.f90`), ASP (chemistry and aerosol, the
common_ASP submodule in `SRC_TUV/ASP/src`) and TUV (photolysis,
`SRC_TUV/TUV/`) are coupled. The design goal: a change to ASP's mechanism,
species lists, number of size bins or photolysis set (the input decks) needs
no edits to SAM code. SAM takes every count, offset and index from ASP at run
time and stops at start-up with a clear message when something does not match.
The only compile-time value that must follow the decks is `ntracers`.

The dated development history is in [CHANGELOG.md](CHANGELOG.md).

## 1. Start-up (`tracers_init`, `SRC_TUV/tracers.f90`)

1. `InitializeASP(DeckDir=asp_deck_dir, OutputDir=asp_output_dir)` reads the
   decks. `SetPhotoFlag(1)` then overrides the deck's `PHOTO FLAG`:
   photolysis rates always come from TUV online (`PhotoRates.in` is used only
   for its column layout).
2. `samasp_layout_init` (`SRC_TUV/samasp_layout.f90`) computes the tracer
   layout (§2) and stops if `ntracers` does not match.
3. The initial plume CO (`GasPhaseCOInit`, from `GasPhaseChems.in`) converts
   ASP's initial particle numbers to particles per kg burned.
4. `samasp_photo_init` (`SRC_TUV/samasp_photo.f90`) reads `PhotoMap.in` (§3).
5. `tuv4samasp_dims` gives TUV's array sizes (`kz`, `kw`, `kj`), which size
   SAM's TUV work arrays.
6. With MPI (`dompi`), `SetASPCellReset(.true.)` makes each ASP cell call
   independent of the cells computed before it (needed for §5).

## 2. Tracer layout and the `ntracers` check

`SRC_TUV/samasp_layout.f90` holds all counts and offsets, taken from ASP
after `InitializeASP`:

| Quantity | From ASP |
|---|---|
| `nGas` | `HowManyEvolveGasChems` |
| `nAq` | `HowManyAqChems + HowManyAqCations + HowManyAqAnions` |
| `nOrg`, `nAqOrg` | `HowManyOrgChems`, `HowManyAqOrgChems` |
| `nSpec` | `nAq + nOrg + nAqOrg` (aerosol species per bin) |
| `nBins` | `ModelParameters::HowManyBins` |

Tracer index functions (use these; no numeric offsets anywhere else):

```
itr_gas(j)       = j                             gas j = 1..nGas
itr_mass(bin,j)  = nGas + (bin-1)*nSpec + j      aerosol species j in bin:
                                                 aqueous | organic | aqueous organic
itr_num(bin)     = nGas + nBins*nSpec + bin      particle number in bin
nTracersNeeded   = nGas + nBins*(nSpec + 1)
```

Tracer names come from the same functions (`aero_species_name`: species
name + 2-digit bin, e.g. `K+01`; `num_conc_name`), so names and values cannot
drift apart. `samasp_gas_index` finds H2O and CO (and any other gas) by name
and stops if the name is missing; so does `Find_Chem`.

`ntracers` in `SRC_TUV/domain.f90` is a compile-time constant (SAM's tracer
arrays are static). If it differs from `nTracersNeeded`, SAM prints the layout
and stops with "Set ntracers to N in domain.f90 and rebuild". The container
build sets it with `--build-arg NTRACERS=N`, and `docker/run_case.py` counts
the tracers the decks need before any work (`deck_tracers`, the same
formula), so a mismatch is caught before ERA5 is downloaded.

3-D output (`SRC_TUV/write_fields3D.f90`) chooses gases by name
(`gas3d_names`) and J-values by TUV label (`jval3d_labels`); names not in the
current mechanism are skipped with a message. The field count follows nBins,
the gases found and the J-values found.

## 3. Photolysis: TUV → ASP (`PhotoMap.in`)

ASP photolysis reactions in `GasChemicalMechanism.in` (`; 0.; N.`) use
column N of `PhotoRates.in`; common_ASP allocates `PhotoRates` with
`HowManyPhotoRxns` columns read from that file. TUV computes J-values for
the reactions switched on in `SRC_TUV/TUV/INPUTS/samasp`, identified by
their labels. `PhotoMap.in`, in the deck directory, connects the two:

```
ASP column ; scale factor ; TUV reaction label
```

Lines with the same column are summed. Columns 49, 50 and 51 have no TUV
reaction of their own and are scaled from others: 49 = 0.33·J(C2H5CHO),
50 = 0.16·J(C2H5CHO) (butanal's two channels), 51 = 100·J(MACR) (HPALD).

Checks: at start-up (`samasp_photo_init`) a column outside
`1..HowManyPhotoRxns`, or a column used by the mechanism with no
`PhotoMap.in` line, stops the run; at the first TUV call
(`samasp_photo_set_labels`) a label that is not among TUV's selected
reactions stops it. `samasp_photo_map` applies the map every call.

To add a photolysis reaction: add its column to the mechanism and to
`PhotoRates.in` (header and data), switch the TUV reaction on in
`INPUTS/samasp`, and add a `PhotoMap.in` line. No code changes.

## 4. TUV site, date and grids (`compute_column_photolysis`, `TUV4SAMASP.f`)

- `tuv4samasp` takes SAM's `latitude0`, `longitude0`, the date (namelist
  `year0`; `day` is the calendar day of `year0` as in SAM's CAM radiation,
  1.0 = 1 January 00 UTC, converted by `day_to_date`) and the UTC time. The
  site and date in `INPUTS/samasp` are not used. SAM's own solar zenith angle
  (passed to ASP, in `tracers_physics`) also uses longitude.
- Vertical grid (`tuv_setup`, once): TUV's lowest levels are SAM's layer
  interfaces `zi`, then levels about 1 km apart up to 120 km, so TUV layer k
  is SAM layer k. J-values are interpolated to SAM's layer mid-heights;
  aerosol optical depth uses each layer's thickness.
- Wavelengths: ASP's Fast-TUV optics use `NFastTUVWvl` bins
  (`FastTUVWavelengths`); `check_tuv_wavelengths` stops at the first TUV call
  if TUV's bin centres do not match.
- The fixed TUV set-up (input file, grids, spectra, cross sections) is done
  once in `tuv_setup`; each call does only the time, aerosol and radiative
  transfer. The discrete-ordinates work arrays are sized for 4 streams
  (`nstr` in `INPUTS/samasp`; `psndo` stops if it is raised).
- Columns whose aerosol optics are identical to an already computed column of
  the same step (every column outside the plume) reuse its J-values
  (`tuv_key`).
- TUV's data paths (`DATAE1/`, `DATAJ1/`, `DATAS1/`, `INPUTS/samasp`) are
  stock TUV and resolve through the symlinks at the repository root, so SAM
  must run in the repository directory.

## 5. ASP calls and load balancing (`tracers_physics`)

Per (sub)step:

1. `add_fire_emissions`: on steps `fire_nstep_start+1 .. fire_nstep_end`
   (`fire_nstep_end = fire_nstep_start + 2*(plume_y-dy)/vg/dt`), in the box of
   global x indices `fire_i1..fire_i2`, levels `fire_k1..fire_k2`.
2. `collect_plume_cells`: after the emission period, cells with CO above
   `asp_co_threshold_ppb` are packed with their ASP inputs.
3. `balanced_compute` (`SRC_TUV/samasp_balance.f90`): the plume cells of all
   ranks, in global order, are split into equal contiguous blocks, computed
   (`asp_cells`, ASP time step = SAM's substep `dtn`) and returned with
   `MPI_Alltoallv`. Because of the cell reset (§1, step 6), the results do
   not depend on the rank a cell is computed on.
4. `store_asp_results` writes the results back to the tracers and the
   aerosol optics arrays.
5. `compute_column_photolysis` runs TUV per column; the J-values are used by
   ASP on the next step.

SAM and TUV are compiled with precise floating-point semantics (`SAM_FP`,
`TUV_FP`; gfortran `-ffp-contract=off`). Under the fast FP model, harmless
code rearrangements change rounding of `waterconc`, which ASP's inorganic
equilibrium amplifies to large differences in a few cells. Keep it precise
unless results need not be reproducible.

## 6. Namelist settings and deck directories

SAM-ASP-TUV settings in `&PARAMETERS` (`SRC_TUV/params.f90`, read in
`setparm.f90`; `RUNDATA/prm_template` lists them):

| Setting | Meaning |
|---|---|
| `asp_deck_dir`, `asp_output_dir` | ASP deck and output directories, relative to the run directory, ending in `/` (default `SRC_TUV/ASP/InorganicInputDecks/`, `SRC_TUV/ASP/OutputFiles/`) |
| `fire_nstep_start` | fire emissions start after this step (spin-up) |
| `fire_i1`, `fire_i2`, `fire_k1`, `fire_k2` | fire box: global x indices and levels |
| `asp_co_threshold_ppb` | ASP runs where CO exceeds this |
| `year0` | calendar year of `day0`/`day` (TUV's date) |
| `samasp_verbose` | per-call progress messages of SAM, TUV and ASP |

## 7. Build

ASP (common_ASP submodule: `git submodule update --init`) and TUV are built
out of tree in `_build/` by `build.sh`, then SAM by `./Build`; the Makefiles
take compilers and library paths from the environment (see README, "Build
variables"). TUV's copies of `ErrMsg`/`WrtBad`/`WrtDim`/`TstBad`/`D1MACH` in
`rtrans.f` are renamed (`ErrMsT`, …, `D1MACT`) so they do not clash with
common_ASP's at link time. The rotor job scripts (`run_tuvtest*.sh`,
`sub_*.sh`) and the `*-fromrotor`, `.no_tuv` and `.tuv` Makefile variants
still contain rotor paths.

## 8. Status

| Phase | Work | Status |
|---|---|---|
| 0 | Baseline build and short run (`RUNS/asptest`); `ParticleAttributes.h` restored | done |
| 1 | Bug fixes: number-concentration name offset, aqueous name suffix, `Find_Chem` not-found stop, `sammap` out-of-bounds writes | done |
| 2 | Tracer layout from ASP (`samasp_layout.f90`), `ntracers` check, gas outputs by name | done |
| 3 | Photolysis by TUV label (`PhotoMap.in`, `samasp_photo.f90`), allocatable `PhotoRates` | done |
| 4 | TUV site, date and grid from SAM; wavelength check; layer thickness | done |
| 5 | Deck/output dirs, ASP time step, CO threshold and fire box in the namelist; build paths from the environment | done |
| 6 | Interface in common_ASP `StepASP4SAM.f90`; common_ASP as a submodule; bundled ASP copies removed | done |

Acceptance test for the goal: add a dummy gas, an organic aerosol species
and a photolysis reaction (with its `PhotoMap.in` line) to the decks, change
only `ntracers` as the start-up message says, rebuild and run; the new
tracers appear with correct names and existing species are unchanged. Then
reorder two species in `GasPhaseChems.in`: outputs are still labelled
correctly.
