# Plan: remove hard-coded SAM ↔ ASP (↔ TUV) links

Goal: a change to ASP's chemical mechanism, species lists, bin count or photolysis set
(the files in `SRC_TUV/ASP/InorganicInputDecks/`) should not need hand edits to SAM
code. SAM should take every count, offset and index from ASP at run time, and it
should stop with a clear message when something cannot be matched.

Status (see §5 for details): phases 0, 1 (partly) and 2 are done on branch
`asp_v2.2.2` (developed as `claude/sam-asp-hardcoded-links-pstj1y`); the code
builds with Intel ifx
(`build_ifx.sh`) and the short `RUNS/asptest` case runs (`run_asptest_ifx.sh`).

---

## 1. What is hard-coded today

### 1.1 Tracer array layout (`SRC_TUV/tracers.f90`, `SRC_TUV/domain.f90`)

The tracer vector is `[gases | bin1 species | … | bin10 species | NumConc bin1..10]`,
with 616 gases (in the comments; the mechanism actually has 618), 86 aerosol species per bin and 10 bins. These numbers appear as literals:

| Where | Literal | Should come from |
|---|---|---|
| `domain.f90:24` | `ntracers = 1489` | `nGas + nBins*(nSpec+1)` (= 1488 today) |
| `tracers.f90:63` | `ngas = 617` | `HowManyEvolveGasChems` |
| `tracers.f90:69` | module `HowManyBins = 10` (hides ASP's run-time `HowManyBins`) | `ModelParameters::HowManyBins` |
| `tracers.f90` `tracers_init` | `GasConc(616)`, `AeroMassConc(10,86)`, `EnvMass(10,86)`, `NumConcWeight(10)` … | ASP counts |
| `tracers.f90` `tracers_physics` | `NumConc(10)`, `TermVel(10)`, … | `HowManyBins` |
| `tracers.f90` init, fire injection, physics | number-conc. tracer = `1476+bin` | computed offset |
| `tracers.f90:483` | number-conc. **name** = `1478+bin` | same computed offset |
| `tracers.f90:179` | `CO_index = 2` | `FindChem("CO",0)` |
| `tracers.f90:792` | `GasConc(1) = waterconc` (H2O assumed first) | `FindChem("H2O",0)` |
| `tracers.f90:815` | ASP is called only where `GasConc(2) > 130` (CO in ppb) | CO index plus a namelist threshold |
| `tracers.f90` debug blocks | species `56..64`, `ntrtrack2 = 673` | name lists |

**Existing bugs this exposes** (to fix first, in their own commits):

- **Number-concentration names are off by 2.** Values go to `1476+bin` but names go
  to `1478+bin`. Tracers 1477 and 1478 have no names, and the stat/3D labels `NoCo01…`
  sit on the wrong tracers.
- **Aqueous tracer names have the wrong bin suffix.** `write(tmp_name,…)` is commented
  out for the aqueous loop (`tracers.f90:422,425`). Bin 1 gets an uninitialised suffix,
  and each later bin gets the previous bin's number. `write_fields3D.f90` then calls
  `Find_Chem("K+01")` and similar names, which may not match.
- **`Find_Chem` returns an uninitialised `Index`** when a name is not found
  (`tracers.f90:105`). Wrong names fail silently.

### 1.2 Output (`SRC_TUV/write_fields3D.f90`)

- Gas outputs are selected by position: `do n=1,92` and tracer `167`. A reordered
  `GasPhaseChems.in` silently changes what is written.
- Per-bin stride `cur_index = 86*(binnum-1)` and `do binnum=1,10`.
- `nfields = 144` is a fixed total that has to be kept in step by hand.
- J-value output uses `PhotoArr(…,2|5|12|14)` by ASP column number.

### 1.3 Photolysis TUV → ASP (`SRC_TUV/TUV/TUV4SAMASP.f`, `tracers.f90`, ASP `StepASP.f90`, `Chemistry.f90`)

- The size 51 is fixed in `PhotoArr(…,51)`, `PhotolysisRates(51)` (tracers.f90,
  `SAM_wrapper`, `ASPInterface`, `InputConcentrations`), `PhotoRates(96,51)`
  (Chemistry.f90) and `samvalj(150,200)`.
- `TUV4SAMASP.f:993`: `sammap(1:59)` is a fixed table from TUV v5.3 reaction numbers to
  ASP `PhotoRates.in` columns. It breaks when the TUV reaction list, the selection in
  `TUV/INPUTS/samasp`, or the column order of `PhotoRates.in` changes. I spot-checked
  it against the current labels and it is correct today.
- **Bug:** the `-1` entries in `sammap` make the loop write `samvalj(-1,iz)`, which is
  out of bounds.
- **Bug:** `samvalj(49:51,iz)` is assigned after the `do iz` loop has ended, so `iz = nz+1`
  (`TUV4SAMASP.f:1023-1025`). `tracers.f90:997-999` recomputes 49–51 from ASP columns
  22 and 25. The `PhotoRates.in` header calls 49 "Butanal" and 50 "HPALD", and the
  mechanism uses 49 and 50 for butanal-type species. The derived-rate definitions live
  in two places and disagree with the header.
- The ASP reaction file already names its photolysis column by number (`; 0.; 22.`), and the
  `PhotoRates.in` header already carries the TUV label text for each column. So the
  mapping can be built by matching labels.

### 1.4 Aerosol optics and the TUV grids

- There are 7 wavelength bins: `ExtCoeff/SingScat/Assym(7)` (tracers.f90, StepASP.f90),
  `ExtCoeffArr(…,7)`, and the FastTUV centre wavelengths set in ASP
  (`LagrangianParticles/ParticleAttributes*.h`, `wvln_array`). These must match TUV
  `params kw=8` and `INPUTS/samasp nwint=-7`. The source says "SAM-ASP alt and wv grid
  and TUV alt and wv grid must match exactly", but nothing checks it.
- `samAOD(200,8)` is sized by hand to TUV's `kz`/`kw`.
- Vertical grid: TUV's `gridz` (grids.f) is hard-wired to a user grid of 193
  levels, every 40 m to 3 km and then every 1 km to 120 km (the `nz=197` in
  `INPUTS/samasp` is ignored). It was made to match tuvtest's SAM grid (80
  layers of 40 m, 0–3.2 km), but only by hand: SAM layer `kk` is copied into TUV
  layer `kk`, so SAM's top 5 layers (3.0–3.2 km) went into TUV's 1 km layers at
  3–8 km, J-values were taken at each layer's bottom rather than its middle, and
  any other SAM grid would be mapped wrongly. `samAOD` uses `dz` and not the
  layer thickness. The top layer index `193` is hard-coded.
- TUV's site and date are hard-coded in `TUV4SAMASP.f` (`lat = 47.9392`,
  `lon = -118.6183`, 2019-08-03, `tmzone = 0`), overriding `INPUTS/samasp`
  (34.7 N, 120.2 E, 2009-11-17); neither is taken from SAM
  `latitude0/longitude0/day`. Only the time comes from SAM. SAM's own solar
  zenith angle (passed to ASP) ignores longitude. So TUV and SAM can disagree
  about day and night: tuvtest (11 UTC) is night for TUV and all J-values
  are zero. `rdinp.f:132` hard-codes `inpfil='samasp'`.

### 1.5 Deck and build paths

- `ASP/src/ReadInputDeck.f90:78` opens `SRC_TUV/ASP/InorganicInputDecks/ASPInputDeck.in`,
  and the transcript goes to `SRC_TUV/ASP/OutputFiles/`. Both paths are relative to the
  run directory. common_ASP uses `InorganicInputDecks/…` instead.
- `Makefile` and `SRC_TUV/ASP/Makefile` use absolute paths `/nas/data_cellar/rotor/...`.

### 1.6 Code divergence

- SAM builds its own copy of ASP (`SRC_TUV/ASP/src`, plus older copies in `SRC/ASP` and
  `…/UNIFAC/`). That copy contains `SAM_wrapper`. common_ASP's `StepASP4SAM.f90` has no
  `SAM_wrapper`, and the two differ by about 300 lines ignoring whitespace. Mechanism
  and ASP fixes made in common_ASP never reach SAM.
- `SRC_TUV/ASP/src/LagrangianParticles/ParticleAttributes.h` is committed as a
  symlink (mode 120000) whose "target" is 64 KB of text. It does not check out as a file.

---

## 2. Design

### 2.1 One layout object, computed once after `InitializeASP`

Add a module `SAMASPLayout`. It belongs in common_ASP's `StepASP4SAM.f90` if §2.5 is
adopted, and in SAM otherwise. It holds:

```fortran
integer :: nGas, nAq, nOrg, nAqOrg, nSpec, nBins, nPhoto, nWvl
integer :: iGas0          ! = 0
integer :: iAero0         ! = nGas
integer :: iNum0          ! = nGas + nBins*nSpec
integer :: nTracersNeeded ! = iNum0 + nBins
integer :: iH2Og, iCOg    ! from FindChem("H2O",0), FindChem("CO",0)
contains
  integer function itr_gas(j)          ! = j
  integer function itr_mass(bin, j)    ! = iAero0 + (bin-1)*nSpec + j   (j over Aq|Org|AqOrg)
  integer function itr_num(bin)        ! = iNum0 + bin
```

- All counts come from ASP (`HowManyEvolveGasChems`, `HowManyAqChems+Cations+Anions`,
  `HowManyOrgChems`, `HowManyAqOrgChems`, `ModelParameters::HowManyBins`).
- Every `1476+…`, `86*…`, `617`, `10`, `1`, `2` literal in `tracers.f90` and
  `write_fields3D.f90` goes through these functions.
- Tracer names are generated from the same functions, so names and values cannot drift apart.
- **Init check:** if `nTracersNeeded /= ntracers`, stop with
  `"set ntracers = <N> in domain.f90"`. SAM's tracer arrays are static, so this is the
  smallest safe change. Making them allocatable is possible later (see Q2).
- `Find_Chem` and `FindChem` lookups used by the interface must abort when the name is
  missing. They must not return garbage.

### 2.2 Interface arrays sized by ASP

- `SAM_wrapper`, `ASPInterface`, `InputConcentrations` and `OutputConcentrations`
  declare arrays from ASP counts (most already do) or as assumed-shape `(:)`/`(:,:)`.
  This includes `PhotolysisRates(nPhoto)` and `ExtCoeff(nWvl)`.
- In `tracers.f90`, the local work arrays become allocatable and are sized from
  `SAMASPLayout`.
- In ASP `Chemistry.f90`, `PhotoRates(96,51)` becomes allocatable. Its size is read from
  `PhotoRates.in`: the number of data rows and the number of columns minus 2.

### 2.3 Photolysis mapping by label, not by number

- Add a small deck, `InorganicInputDecks/PhotoMap.in`, with one line per ASP photolysis
  column, for example:

  ```
  ! ASPcol ; scale ; TUV label (exact text of jlabel), may repeat to sum channels
  23 ; 1.00 ; HOCH2CHO -> CH2OH + HCO
  23 ; 1.00 ; HOCH2CHO -> CH3OH + CO
  23 ; 1.00 ; HOCH2CHO -> CH2CHO + OH
  49 ; 0.33 ; C2H5CHO -> C2H5 + HCO
  50 ; 0.16 ; C2H5CHO -> C2H5 + HCO
  51 ; 100. ; CH2=C(CH3)CHO -> Products
  ```

- At the first TUV call, each label is matched against `jlabel(imj(:))`. The result is an
  integer map plus scale factors. A missing label, or an ASP column referenced in
  `GasChemicalMechanism.in` with no source, stops the run with the column number and name.
- This replaces `sammap`, the `49/50/51` lines in both TUV4SAMASP.f and tracers.f90, and
  the `-1` out-of-bounds writes.
- `tuv4samasp` returns `samvalj(nPhoto, nzm)` already mapped, and `PhotoArr` becomes
  `(…, nPhoto)`.
- The J-value outputs in `write_fields3D.f90` are chosen by label (a list in the SAM
  namelist) and not by column number.

### 2.4 Optics and TUV grid consistency

- ASP gets a function that exposes `nWvl` and `wvln_array`. TUV gets one that exposes
  `nw`, `wc(:)`, `nz`, `z(:)`, `kz`, `kw`, `kj`. At init, check that the ASP wavelength
  centres match the TUV bin centres within a tolerance, and stop if they don't.
- Interpolate the SAM column onto the TUV altitude grid instead of copying index by
  index. Compute AOD from the real layer thickness `zi(k)-zi(k-1)`. Remove the
  hard-coded index `193`.
- Pass `latitude0`, `longitude0`, the date and the time from SAM into `tuv4samasp`, replacing
  the values hard-coded in `TUV4SAMASP.f`; use longitude in SAM's own solar zenith angle. Make the TUV input file name an argument.

### 2.5 Deck paths, build paths and code source

- Add SAM namelist entries `asp_deck_dir` and `asp_output_dir`. Pass them to
  `ReadMainInputDeck` and `SetInputDeckSubDir`, and remove the literal
  `SRC_TUV/ASP/...` paths.
- Makefiles: derive `ASP_SRC_DIR`/`TUV_SRC_DIR` from `$(SAM_DIR)` or the environment,
  with no `/nas/...` paths.
- **Recommended:** move `SAM_wrapper` and `SAMASPLayout` into common_ASP
  `StepASP4SAM.f90`. Then use common_ASP from SAM-ASP-TUV as a git submodule, the way
  the ASP and STILT-ASP repos already do, and delete the bundled `SRC_TUV/ASP/src`,
  `SRC/ASP` and `UNIFAC/` copies. After that, a mechanism change is made once in
  common_ASP and SAM picks it up.

---

## 3. Phases

Each phase ends with a regression run against the Phase 0 baseline.

| Phase | Repo(s) | Work | Expected output change |
|---|---|---|---|
| 0 | SAM | Build the current code and save a short-run baseline (`OUT_STAT`, a few 3D steps) from an existing case (`RUNS/tuvtest` or `W1101`). Fix the broken `ParticleAttributes.h` symlink. | none |
| 1 | SAM, TUV | Bug fixes only: NoCo name offset, aqueous name suffix, `Find_Chem` not-found abort, `sammap` `-1` writes, 49–51 after-loop assignment. | labels change; values change only where the out-of-bounds write mattered |
| 2 | SAM (or common_ASP) | `SAMASPLayout`, the `ntracers` check, and replacing every numeric offset in `tracers.f90`/`write_fields3D.f90`. Gas outputs selected by name list. | bit-identical to Phase 1 |
| 3 | ASP, TUV, SAM | Allocatable `PhotoRates`, `nPhoto`, `PhotoMap.in`, label-based mapping, remove `sammap`. | bit-identical (same map) |
| 4 | ASP, TUV, SAM | Wavelength and grid checks, vertical interpolation, real layer thickness, lat/lon/date from SAM. | changes (physics fix); review separately |
| 5 | SAM, Make | Deck and output dirs in the namelist, relative build paths. | none |
| 6 | common_ASP, SAM | Move the interface into `StepASP4SAM.f90` and use common_ASP as a submodule. | should match Phase 5 if both have the same ASP version |

**Acceptance test for the goal:** add one dummy gas species, one organic aerosol species
and one photolysis reaction to the decks (new `PhotoMap.in` line). Change only
`ntracers` as the init message instructs, rebuild and run. The new tracers must appear
with correct names, and existing species values must be unchanged. Then reorder two
species in `GasPhaseChems.in`: outputs must still be labelled correctly.

---

## 4. Open questions

1. **Source of ASP for SAM:** should SAM-ASP-TUV switch to common_ASP as a submodule
   (Phase 6)? If yes, which branch: `master` or `asp_v_2.2.1`, which the common_ASP
   patch notes refer to? Today the two ASP codebases differ.
2. **`ntracers`:** is a run-time check plus an "edit domain.f90" message enough, or should
   SAM's tracer arrays become allocatable so that no rebuild is needed?
3. **Derived photolysis columns 49/50/51:** the `PhotoRates.in` header says 49 = Butanal
   and 50 = HPALD, and 51 is not documented. The code sets 49 = 0.33·J(C2H5CHO),
   50 = 0.16·J(C2H5CHO) and 51 = 100·J(MACR). Which definitions are intended?
4. **Fire-case settings in `tracers_physics`** (`rank.eq.7`, `iden 100..115`,
   `zden 32..55`, `nstep>1080`, the `GasConc(2)>130` trigger, `tstep=10`): these are
   case set-up, not the ASP link. Should they move to the namelist as part of this
   work, or be handled separately?

---

## 5. Progress

**Build and test set-up**
- Makefiles take paths and compilers from the environment (defaults = rotor
  ifort build); `build_ifx.sh` builds with Intel ifx + Open MPI + NetCDF.
- `RUNS/asptest` (fire from step 0, 30 steps) exercises SAM, TUV and ASP in
  about 12 minutes; repeated runs give bit-identical output.

**Phase 1 (bug fixes) — done so far**
- TUV4SAMASP.f: `sammap` −1 entries no longer written (out of bounds); J-value
  loop bound `nj` → `nmj`; "fix top layer" applied for every wavelength
  instead of `iw = kw+1`. savout.f: `outfil` length mismatch.
- Fire emissions: rank 7 wrote `tracer(100..115,…)` out of bounds; the box is
  now in global x indices and set in the namelist (`fire_*`).
- `ParticleAttributes.h` restored (was a truncated symlink).
- Still open: the 49–51 derived J-values after the TUV loop (question 3).

**Phase 2 (tracer layout) — done**
- New `SRC_TUV/samasp_layout.f90`: counts, offsets, tracer names and the
  H2O/CO gas indices come from ASP after `InitializeASP`. SAM stops at start-up
  with the required `ntracers` if `domain.f90` does not match (tested by
  changing the bin count in `AerosolModes.in`).
- `tracers.f90` and `write_fields3D.f90` use the layout everywhere (no 616,
  617, 86, 10-bin, 1476/1478/1479 literals left); 3D gas outputs are chosen by
  name (`gas3d_names`); `Find_Chem` stops on a missing name.
- Bugs this fixed: the mechanism has 618 gases, but number concentrations were
  stored at 1476+bin, overlapping the last two aqueous-organic species of
  bin 10; number-concentration names and 3D outputs were offset by 2; aqueous
  tracer names carried the wrong bin number, so the 3D aqueous totals (TotK,
  TotNa, …) used the wrong bins; `GasConc(616)` was passed where 618 were
  expected.
- 3D gas output position 167 was meant to be PAN (per the code comment) but
  had become H13CO2CHO; the name list now writes PAN.


**Phase 6 (common_ASP as a submodule) — done (decisions: 1a, 2a, 3a)**
- `SRC_TUV/ASP/src` is now the common_ASP submodule, on branch `asp_v2.2.2`
  of common_ASP. The SAM-ASP-TUV changes were developed on common_ASP branch
  `claude/sam-asp-hardcoded-links-pstj1y` and fast-forwarded into
  `asp_v2.2.2` (at `bf6da0b`); the ASP model repository's `asp_v2.2.2`
  pins the same common_ASP commit.
  After cloning SAM-ASP-TUV: `git submodule update --init`.
- Moved into common_ASP: the SAM-ASP-TUV interface (`StepASP4SAM.f90`:
  photolysis rates from TUV, 7-bin Fast-TUV aerosol optics, M3_TO_CM3 = 1e6),
  `ShellRefIndAndRadFASTTUV`/`AerosolOptPropFASTTUV` (alongside the
  full-spectrum routines), the `PhotoRates` export, and `SetMainDeckDirs`.
- Not carried over (decisions 2a and 3a): SAM's edits to the electrolyte
  equilibrium in Thermodynamics.h (asp_v2.2.2's rewritten solver is used),
  and SAM's downgrade of the coagulation volume-loss check to a warning
  (asp_v2.2.2 keeps it fatal and regrids out-of-order sections first).
- SAM side: deck and output directories are namelist settings
  (`asp_deck_dir`, `asp_output_dir`); `ASPInputDeck.in` has the three
  lights on/off lines that asp_v2.2.2 requires; TUV's copies of
  ErrMsg/WrtBad/WrtDim/TstBad/D1MACH in rtrans.f are renamed (ErrMsT, …,
  D1MACT; same length, so fixed-form columns are unchanged) so they do not
  clash with common_ASP's at link time.
- The stale prebuilt ifort ASP files (SRC_TUV/ASP libAspLibrary.a,
  ObjectFiles/, ModuleFiles/, asp_v2.1.1_linux_ifort_dbl, and the root
  libAspLibrary.a link) and the old ASP copy under SRC/ASP have been removed.
  The ASP Makefile creates ObjectFiles/ and ModuleFiles/ when it builds.
- The stale ifort TUV objects (SRC_TUV/TUV/*.o, libTuvLibrary.a and the root
  libTuvLibrary.a link) have been removed as well.

**Question 3 (derived J-values 49–51) — answered**
- Intended: 49 = 0.33·J(C2H5CHO) (butanal's first channel and 22 other
  reactions in GasChemicalMechanism.in), 50 = 0.16·J(C2H5CHO) (butanal's
  second channel), 51 = 100·J(MACR) (HPALD1, HPALD2); TUV has no rates for
  these species.
  They are currently set in both TUV4SAMASP.f and tracers.f90; phase 3
  should keep these definitions and set them in one place.

**Phase 3 (photolysis by reaction label) — done**
- `SRC_TUV/ASP/InorganicInputDecks/PhotoMap.in`: one line per source,
  `ASP column ; scale ; TUV reaction label`; it reproduces the old `sammap`
  table and the derived columns 49–51.
- `SRC_TUV/samasp_photo.f90` reads it; TUV4SAMASP.f returns the J-values and
  labels of its selected reactions (no `sammap`, no 49–51 assignments, no
  per-call debug prints); `tracers.f90` sizes `PhotoArr`/`PhotolysisRates`
  by the number of ASP columns; J-value outputs are chosen by TUV label.
- common_ASP (d038ff2): `PhotoRates` is allocated with `HowManyPhotoRxns`
  columns read from `PhotoRates.in` (rows stay 96 = 24 h at 15 min).
- Tested: daylight `asptest` output (3D and statistics files) is
  byte-identical to the run before phase 3. Start-up stops with a clear
  message for an unknown TUV label, a column outside `PhotoRates.in`, and a
  mechanism photolysis column with no source.
- To add a photolysis reaction now: add its column to the mechanism and to
  `PhotoRates.in` (header and data), select the TUV reaction in
  `TUV/INPUTS/samasp`, and add a `PhotoMap.in` line — no code changes.

**Phase 4 (TUV site, date and grids from SAM) — done**
- `tuv4samasp` takes `latitude0`, `longitude0`, the calendar date (new namelist
  setting `year0`; `day` = days since 1 January 00 UTC) and SAM's layer
  interfaces. TUV's lowest levels are SAM's interfaces, then about 1 km apart to
  120 km, so TUV layer k is SAM layer k; J-values are interpolated to SAM's layer
  mid-heights; aerosol optical depth uses each layer's thickness; the top-layer
  fix uses the real top layer. (Correction: before this, TUV used the 40 m / 1 km
  grid hard-wired in `gridz`, which matched SAM's 40 m layers in tuvtest/asptest
  up to 3 km; only the top 5 layers were mapped to 3–8 km. The phase 4 commit
  message (~49 km) and an earlier version of this note (2x stretch, which
  assumed 20 m layers) were wrong.)
- SAM's own solar zenith angle (passed to ASP) now uses longitude.
- The tuv4samasp arrays are sized from TUV's params; ASP exports its optics
  wavelengths (`NFastTUVWvl`, `FastTUVWavelengths`), and SAM checks TUV's bins
  against them at the first TUV call.
- Test: asptest at 20 UTC: TUV's solar zenith angle 53.96° (expected ~54° for
  34.7 N, 18 November, local noon), 198 TUV levels, wavelength check passed,
  no runtime errors.
- Check against standalone TUV (TUV.f main program, the same TUV library,
  default atmosphere from `INPUTS/samasp`, asptest's site and time): in
  clear-sky columns the coupled J(O1D) and J(HONO) are within 1.2% of
  standalone TUV without aerosol at every SAM level from 20 m to 3.2 km
  (the coupled column has no aerosol in SAM's layers); both rise ~58%
  (O1D) and ~33% (HONO) from the surface to 3.2 km. Standalone TUV gives the
  same values on SAM's grid as on its own, and repeated `tuv4samasp` calls
  give identical results.

**Phase 5 (remaining hard-coded settings) — done**
- ASP's time step is SAM's substep `dtn` (was a literal 10 s).
- CO threshold for calling ASP: namelist `asp_co_threshold_ppb` (default 130).
- `create_run.py`: `rundir` defaults to the script's directory; writes `year0`.
- Left as they are: TUV's data paths (`DATAE1/...`, `INPUTS/samasp`) are stock
  TUV and resolve through the symlinks at the repository root; the rotor job
  scripts (`run_tuvtest*.sh`, `sub_*.sh`) and the `*-fromrotor` / `.no_tuv` /
  `.tuv` Makefile variants still contain rotor paths.
- Open question: `tracers.f90` converts ASP's initial plume particle number to
  particles per kg burned with `CO_init = 3400` ppb, while ASP's
  `GasPhaseChems.in` sets the initial plume CO to 10000 ppb.


**Python scripts — converted to Python 3**
- `create_run.py`, and both copies (SRC/, SRC_TUV/) of
  `InorganicInputDecks/defaults/MappingPrelim/process_photo_rate_table.py`
  and `_getListOfCompounds.py`. Conversion with fissix (lib2to3 fork), then
  reviewed: whole-hour integer divisions in `create_run.py` keep `//`; mixed
  tab/space indentation expanded as Python 2 read it (tab = 8 columns).
- Checked: `_getListOfCompounds.py` reproduces the committed
  `_GaseousSpeciesBlankList.in` exactly; `process_photo_rate_table.py`
  reproduces the committed `PHOTO.DAT` line for line, with the reactions in
  PhotoRates.in order (Python 2 used arbitrary dict order). `create_run.py`
  compiles but was not run (it downloads NCEP data).

**Initial plume CO from the deck and CAM day convention — done**
- `tracers_init` uses ASP's initial plume CO from `GasPhaseChems.in`
  (`GasPhaseCOInit`, 10000 ppb) instead of the hard-coded 3400 ppb when
  converting ASP's initial plume particle numbers to particles per kg burned.
  Test (asptest, before ASP's first call): emitted POA and particle numbers
  scale by 0.340 (= 3400/10000); gases and dynamics are bit-identical.
- `day0`/`day` follow SAM's CAM radiation (1.0 = 1 January 00 UTC);
  `create_run.py` writes them that way. asptest is 17 November 20 UTC again
  (TUV solar zenith angle 53.96°).

**Performance (same physics and chemistry) — items 1, 3, 2**

Timings are asptest (30 steps, 16 MPI ranks oversubscribed on this machine).

- *Item 1, TUV fixed set-up once* (`TUV4SAMASP.f`): input file, grids,
  spectra, cross sections and fixed optical depths were re-read and
  recomputed on every call (every column, every step), ~90% of a non-plume
  rank's time. Now done on the first call; per call only time, aerosol and the
  radiative transfer. 560 s -> 233 s.
- *TUV compiled with `-fp-model precise`* (`Makefile_SAM`): under the fast FP
  model the item-1 change altered TUV's single-precision rounding (J-values
  up to ~2e-5 relative; old and new code are identical under precise). With
  precise semantics TUV's results no longer depend on how unrelated code is
  arranged. J-values change once, relative to the fast-FP build.
- *Item 3, DPREP ordering cache with 8 slots* (common_ASP `opkda1.f`): ASP
  alternates between ODE systems of 5 sparsity structures; the one-slot cache
  hit on ~81% of calls, 8 slots on 99.99%. Output byte-identical;
  233 s -> 197 s.
- *Item 2, ASP load balancing across MPI ranks*:
  - Order independence (common_ASP `StepASP4SAM.f90`, option A): ASP calls
    warm-start from the particle state of the previous call, and the
    tabulated Gibbs-Duhem integral's range grows with the ionic strengths
    seen, so a cell's result depended on the cells computed before it.
    `SetASPCellReset(.TRUE.)` saves the particle state; every later
    `ASPInterface` call restores it and resets the Gibbs-Duhem tables.
    SAM turns this on only when `dompi` is true; it is off by default, so
    ASP's Lagrangian parcel runs are unchanged. (This changes SAM-ASP-TUV
    MPI results once relative to the warm-started runs.)
  - Tests: every cell recomputed in reverse order within its rank (0
    differences in 42 rank-steps), and every rank recomputing the next
    rank's cells (0 differences in 10,558 cells).
  - `SRC_TUV/samasp_balance.f90` + `tracers_physics`: the plume cells of all
    ranks, in global order, are split into equal contiguous blocks; inputs
    and outputs are exchanged with `MPI_Alltoallv`. `tracers_physics` now
    (1) collects each rank's plume-cell inputs, (2) runs ASP on the balanced
    blocks, (3) writes results back, (4) runs TUV per column (as before,
    after that column's ASP results are stored).
  - Result: with SAM compiled with `-fp-model=precise`, the balanced run is
    bit-identical to the unbalanced run (all 478 bin2D fields and the stat
    file); 170 s -> 114 s. With SAM's default (fast) FP model the two differ
    because the restructured loop changes the compiler's rounding of the
    water concentration (`waterconc`) by one ulp in some cells; ASP's
    inorganic equilibrium amplifies that to up to ~60% in a few cells of
    TotNO3/NH3/TotSO4 (most fields ~3e-6 relative). Without precise FP,
    any edit of `tracers.f90` can do the same.
- *SAM compiled with `-fp-model precise` by default* (`Makefile`, `SAM_FP`;
  override with `SAM_FP=` on the make command line or in the environment):
  for the reason above. Changes SAM's results once at the rounding level;
  ~3% slower. asptest from the repository build: 117 s, output identical to
  the verified precise balanced run.
- *Statistics lookup and tracer initialisation* (`hbuffer.f90`,
  `tracers_init`): `hbuf_put` scanned all 9168 statistics names on each
  call (1.9 s per statistics call, ~25% of asptest); a hash index makes it
  0.02 s. The start-up duplicate-name check (all pairs, ~2.5 s) and the
  per-cell recomputation of per-bin values in `tracers_init` (~1 s) are
  gone. Output bit-identical; asptest 117 s -> 78 s; start-up to the end of
  initialisation 6.4 s -> 2.9 s.
- *TUV radiative-transfer arrays and log output*: a profile showed ~12% of
  CPU in `memset`, nearly all of it TUV's discrete-ordinates solver
  (`psndo`, `rtrans.f`) clearing work arrays sized for 32 streams on every
  call, while TUV uses 4 (`nstr` in `INPUTS/samasp`). They are now sized
  for 4 streams (psndo stops with an error if `nstr` is raised without
  raising them). Per-call progress messages of SAM, TUV and ASP (~42 MB
  of log per asptest run) are printed only with the new namelist option
  `samasp_verbose`; ASP's switch (`StepPrints`, common_ASP) stays on for
  ASP's own drivers. UNIFAC opens `unifac.err` once instead of on every
  call. Output bit-identical; asptest 79 s -> 66 s; log 42 MB -> 1.8 MB.
- *Gas-phase chemistry evaluators* (common_ASP patch 38): the derivative
  and Jacobian routines called by the gas-chemistry ODE solver read flat
  arrays copied once from the mechanism's linked lists, with the same
  operations in the same order, and the Jacobian skips its empty entries
  (616^2 checked per Jacobian before, 5,929 non-empty). The solver itself is
  unchanged. Tested byte-identical on ASP's Lagrangian parcel decks
  Williams, Otavi, Alaska, Timbavati (baseline mechanism) and
  Williams_CRACMM3 (all 108 numerical output files each) and on asptest
  (bit-identical); asptest 66 s -> 54 s.
- *TUV column reuse and once-per-step optics checks* (`tracers_physics`):
  every column has the same site, time and grid, so a column whose aerosol
  optics are bit-identical to an already computed column of the same step
  (all columns outside the plume) reuses its photolysis rates instead of
  calling TUV; the whole-array bad-value checks of the optics arrays run once
  per step instead of after every plume cell. Output bit-identical; asptest
  54 s -> 52 s (TUV was ~9% of CPU, ~3 s of wall time here).
- *Meteorology from ERA5* (`create_run.py`, `samasp_case/`,
  `RUNDATA/prm_template`): the case files are made from the ERA5 column
  nearest the fire (Copernicus CDS; only a 3 x 3 point box and one hour are
  downloaded and cached) instead of whole-month NARR and whole-year NCEP
  files from NOAA's old FTP server. One global 0.25 deg hourly source for the
  sounding, boundary-layer height, heat fluxes and surface stress (NARR is
  32 km, North America only, with NCEP R1 at ~1.9 deg for the stress).
  Command-line interface; the old params.txt still works; the prm comes from
  a template holding the previous settings plus the newer SAM-ASP options as
  comments. Same conversions as before; surface stress divided by the
  actual air density (was 1.2 kg/m3). The SAM axes ('Ay'/'Ax') follow the
  mean boundary-layer wind (vector mean over the levels above the ground and
  below the boundary-layer top), and vg is the mean along y over the same
  levels; before, the axis followed the wind at a single level (often
  1000 hPa, which can be below the ground) and vg also averaged levels below
  the ground. Offline tests in
  `tests/`; a case written from synthetic ERA5 files ran in SAM (100 steps).
- *Container* (`Dockerfile`, `docker/run_case.py`, README "Running in a
  container"): a two-stage image (build: Ubuntu 24.04, Open MPI, NetCDF-C,
  Intel ifx/icx from conda-forge, NetCDF-Fortran rebuilt with ifx, then
  `build_ifx.sh`; runtime: 1.1 GB, only the Intel/NetCDF-Fortran libraries
  the binary links). One container runs one fire case: latitude, longitude,
  date, start time, run length and the ASP deck directory (mounted) are set
  at `docker run` (environment variables or options), plus optional emission
  flux, fire start step, time step and the ERA5 source (CDS key or mounted
  files); results go to a mounted /output. `create_run.py --set NAME=VALUE`
  sets any &PARAMETERS entry. Tested: image built here; a 72-step case with
  mounted decks and synthetic ERA5 files ran (fire, then ASP chemistry: O3
  titration, PAN and HONO in the plume); errors for a missing deck or ERA5
  source.
- Overall, asptest against the pre-submodule baseline (`2ecba69`, same
  machine and compilers): 579 s -> 52 s, about 11 times faster.
- Statistics name clash fixed: the ASP gas tracer `MCR` (methacrolein) had
  the same name as SAM's statistic `MCR` (core mass flux), so both
  stat-file profiles named `MCR` held the sum of the two quantities (since
  the first commit). SAM's statistic is now `MCRCORE` (`RUNDATA/lst`, the
  `RUNS/*/lst` copies, `statistics.f90`), and a duplicate statistics name now
  stops the run at start-up with a message (`hbuffer.f90`) instead of only
  printing a warning. asptest: 3D output identical; in the stat file only
  that entry's name changes (no cloud cores there, so its values were zero).
- `UTIL/stat2nc` fixed for SAM-ASP-TUV statistics files. It crashed because
  the readers (`UTIL/SRC/hbuf.inc`, `hbuf_lib.f`) allowed 1,000 profiles and
  8-character names, while SAM-ASP writes about 9,000 profiles with
  20-character names (truncation would also have merged, e.g., `1-ButeneFLX`
  into `1-Butene`). Now: 100,000 profiles and 20-character names; profiles are
  defined at once and read in one pass (0.5 s for asptest, instead of one full
  file scan per profile); names netCDF rejects (`(NH4)2SO401`) get a leading
  `_` with the original in `sam_name`; 64-bit-offset netCDF. asptest: all
  9,144 profiles at all 3 times equal the raw file exactly. Built from source
  in the container. The old prebuilt `UTIL/stat2nc` binary was removed from
  the repository.
- All UTIL post-processing tools are now built from source only
  (`UTIL/build_util.sh`, run by the Dockerfile; binaries are git-ignored via
  `UTIL/.gitignore`). The 22 old prebuilt binaries (and a stray `a.out`) were
  removed; `2Dcom2nc_sep` had a source but no binary and is now built too.
  Building from source exposed a bug: `bin2D2nc` and `com2D2nc` kept the
  netCDF variable ids in `varid(99)`, and SAM-ASP files have 143 fields, so
  the array overflowed. With the old binary this silently corrupted the first
  time value (1.9e-43 instead of 321.83 days) and `TotCa` at the 2nd and 3rd
  output times; with the new build it corrupted `x`. Now `varid(10000)`.
  asptest: the new `bin2D2nc` output equals the raw `.bin2D` data for all 143
  fields at all 3 times. The image no longer needs the `libz.so` link the old
  binaries wanted. Only `bin2D2nc` and `stat2nc` were tested on model output
  (asptest writes no 3D, com, 2D-bin or satellite-simulator files); all 23
  tools compile and load their libraries in the image.
- On-prem guide for a RHEL 9 VM in a CMMC Level 2 (CUI) enclave
  (`docs/OnPrem_RHEL_CMMC_Guide.md`): items to agree with the ISSO (with the
  NIST SP 800-171 requirements they touch), administrator setup for rootless
  Podman, building on the VM or transferring the image, running with `:Z`
  mounts, the CDS key as a mounted file, long runs as systemd user services,
  data handling, records and patching. Tested here: the image under rootless
  Podman 4.9 (72-step case, outputs byte-identical to Docker, owned by the
  user) and peak memory (about 1.1 GB for 16 ranks). The Dockerfile's base
  image is now fully qualified (`docker.io/library/ubuntu:24.04`) so Podman
  needs no short-name lookup.
- RHEL 9 container image (`Dockerfile.ubi9`, UBI 9.8; `BASE` build-arg for
  a hardened UBI 9): OS, Python 3.9, NumPy, requests and OpenSSL from RHEL;
  Intel ifx/icx 2026.2.0, Open MPI 4.1.6, HDF5 1.14.6 and NetCDF-C 4.9.3 from
  conda-forge (pinned; not in the UBI repositories); NetCDF-Fortran 4.6.2
  from a SHA-256-checked source tarball; cdsapi/netCDF4/cftime pinned from
  PyPI. Open MPI uses `plm=isolated` (no ssh in the image). `run_case.py`
  finds `mpirun.openmpi` or `mpirun`. The 72-step test case gives
  byte-identical stat and 2-D output to the Ubuntu image, as root and as a
  rootless Podman user. Fixed in both Dockerfiles: the UTIL tools' runtime
  libraries were looked up without the `UTIL/` prefix. The on-prem guide now
  assumes a 16-vCPU AMD VM with 64 GB and an NVIDIA L4 (not used: the model
  is CPU/MPI only).
- Fire cases (`cases/`) from the SAM-ASP-TUV v1.0.0 release (Zenodo
  10.5281/zenodo.13385240): `williams_flats_2019` (SND1; the v1.0.0 decks =
  the image's default decks, old 2x4-bin organic scheme; warm-bubble start via
  `SAM_SET`), `williams_2009` (W1101/tuvtest; ASP repository `Williams/`
  baseline deck + `PhotoMap.in`, converted to 10 size bins; start 19:00 UTC), `bbop_colockum_2013`
  (ADTEST/SAMASP0001; 2013-08-09 18:49 UTC; no decks). Each has `case.env`
  for `--env-file`. New: `run_case.py` `--set` / `SAM_SET` (extra namelist
  settings) and the `NTRACERS` build-arg in both Dockerfiles (Williams 2009
  needs 1486, the default decks 1488). All cases use 10 size bins: with 40
  bins (4096 tracers) the Williams deck failed in SAM's sectional coagulation.
  SAM-ASP-TUV now always sets ASP's PhotoFlag after reading the decks
  (`tracers.f90`): photolysis is TUV's online rates, and ASP parcel decks
  with `PHOTO FLAG = 0` no longer stop the run.
  Fixed `write_fields3D.f90`: the 3D field count assumed 10 size bins (one
  particle-number field per bin); now 37 + nBins + gases + J-values.
- The ASP repository (asp_v2.2.2) gained `WilliamsFlats` decks for the baseline
  and all eight other mechanisms, built from the v1.0.0 decks with the 9-bin
  volatility scheme (`WilliamsFlatsTranslation/`).
- Tracer pre-check (`docker/run_case.py`): before building the case or
  downloading ERA5, the entry point counts the tracers the ASP decks need
  (evolving gases + H2O, then per size bin the aqueous species and ions +
  water, the organic species + 1 and the aqueous organics, plus one number
  tracer per bin, as `samasp_layout_init`) and compares them with `ntracers`
  in `SRC_TUV/domain.f90`; a mismatch stops with the `NTRACERS` value to
  build with. `--check-decks` checks only. Counts verified against SAM's own
  layout output: 1488 (Williams Flats), 1486 (Williams 10 bins), 4096
  (Williams 40 bins); `tests/test_run_case.py`.
