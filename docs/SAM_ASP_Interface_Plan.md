# Plan: remove hard-coded SAM ↔ ASP (↔ TUV) links

Goal: a change to ASP's chemical mechanism, species lists, bin count or photolysis set
(the files in `SRC_TUV/ASP/InorganicInputDecks/`) should not need hand edits to SAM
code. SAM should take every count, offset and index from ASP at run time, and it
should stop with a clear message when something cannot be matched.

Status: plan only. No code has been changed. Nothing here has been built or run.

---

## 1. What is hard-coded today

### 1.1 Tracer array layout (`SRC_TUV/tracers.f90`, `SRC_TUV/domain.f90`)

The tracer vector is `[gases | bin1 species | … | bin10 species | NumConc bin1..10]`,
with 616 gases, 86 aerosol species per bin and 10 bins. These numbers appear as literals:

| Where | Literal | Should come from |
|---|---|---|
| `domain.f90:24` | `ntracers = 1489` | `nGas + nBins*(nSpec+1)` (= 1486 today) |
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
- Vertical grid: `INPUTS/samasp` has `nz=197, zstop=120 km`, but SAM has `nzm=80`
  stretched levels. SAM layer `kk` is copied into TUV layer `kk` without regridding.
  `samAOD` uses `dz` and not the layer thickness `zi(k)-zi(k-1)`. The top layer index
  `193` is hard-coded.
- Lat/lon/date in `INPUTS/samasp` are fixed (34.7 N, 2009-11-17) and are not taken from
  SAM `latitude0/longitude0/day`. `rdinp.f:132` hard-codes `inpfil='samasp'`.

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
- Pass `latitude0`, `longitude0`, the date and the time from SAM into `tuv4samasp`, so they
  override the values in `INPUTS/samasp`. Make the TUV input file name an argument.

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
