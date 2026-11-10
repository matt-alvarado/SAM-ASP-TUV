# Colockum Tarps fire, Washington, 9 August 2013

Sampled in the BBOP campaign. SAM-ASP-TUV v1.0.0 cases `ADTEST` and
`SAMASP0001` (identical), made by the old `create_run.py` from `params.txt`
(47.20317 N, 120.25683 W, 2013-08-09 18:49 UTC, `massemis = 4.86e-6`).
Their `prm` had `day0 = 220.75` (8 Aug 18:00 UTC); this case uses the
`params.txt` date and time.

| Setting | Value |
|---|---|
| Location | 47.20317 N, 120.25683 W |
| Start | 2013-08-09 18:49 UTC |
| Run length, time step | 2 h, 8 s (`nstop = 900`) |
| Emission flux | `MASSEMIS=4.86e-06` |

**No ASP decks.** The v1.0.0 release kept only one deck set (the Williams Flats
one), so this case's emissions and smoke aerosol must be supplied: put a deck
set in `decks/`, set `ASP_DECK_DIR`, and build an image with the `NTRACERS`
value the model reports.

The original run (900 steps) ends before the default fire start (step 1080),
so set `RUN_HOURS` longer or `FIRE_START_STEP` (for example
`FIRE_START_STEP=0` with a separate spin-up).
