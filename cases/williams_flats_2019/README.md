# Williams Flats fire, Washington, 3 August 2019

Sampled in FIREX-AQ. SAM-ASP-TUV v1.0.0 case `SND1` (run `SND1E1S1` and its
variants `_v1` ... `_v10`, whose saved SAM settings are identical).

## Settings (`case.env`)

| Setting | Value | From the original `RUNS/SND1/prm` |
|---|---|---|
| Location | 47.9392 N, 118.6183 W | `latitude0`, `longitude0` |
| Start | 2019-08-03 19:36 UTC | `day0 = 215.817` |
| Run length, time step | 8 h, 10 s | `nstop = 2880`, `dt = 10` |
| Emission flux | `MASSEMIS=8.3e-07` | `massemis` |
| Warm-bubble start | `SAM_SET`: `perturb_type=2`, bubble at x = 60 km, z = 1.5 km, 1 km x 100 m, +10 K | `perturb_type`, `bubble_*` |
| Fire start | default (step 1080, 3 h) | not set (default) |

The plume deposition velocity `vg` now comes from the ERA5 boundary-layer
wind (it was 0.34 m/s with the original sounding), and the meteorology from
ERA5 instead of the `SND1` `snd`, `lsf` and `sfc` files.

## ASP decks (`decks/`)

The SAM-ASP-TUV v1.0.0 deck set, with its original organic-aerosol scheme:
two 4-bin volatility sets (IVOC1-4 emitted, IVOC5-8 for SOA) without OH aging,
and the lumped SOA precursors OXYARO and OTHERO. Gas emissions are emission
factors (g/kg), with the Williams 2009 values for most species and NOx emitted
as NO; the smoke mode is "Williams Flats Mod" (634,000 cm-3, 0.09 um,
sigma 1.5, 10 size bins). These are the image's default decks
(`SRC_TUV/ASP/InorganicInputDecks`), copied here with `PhotoMap.in`.

A version for standalone ASP with the current 9-bin scheme, for the baseline
mechanism and eight others, is in the ASP repository
(`InorganicInputDecks/WilliamsFlats*`, built by
`WilliamsFlatsTranslation/build_williamsflats_decks.py`).

Tracers: 1488 (the default image).
