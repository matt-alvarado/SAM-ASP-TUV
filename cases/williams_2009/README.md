# Williams fire, California, 17 November 2009

The Williams fire studied by Akagi et al. (2012) and modeled with ASP by
Alvarado et al. (2015, ACP). SAM-ASP-TUV v1.0.0 cases `W1101` and `tuvtest`
(same settings) and `WILLIAMS` (shorter, `nstop = 380`, `massemis = 4e-5`).

## Settings (`case.env`)

| Setting | Value | From the original `RUNS/W1101/prm` |
|---|---|---|
| Location | 34.69583 N, 120.2064 W | `latitude0`, `longitude0` |
| Start | 2009-11-17 19:00 UTC | `day0 = 321.458` was 11:00; 19:00 UTC (about 11:00 local) is the observation start of the ASP Williams deck |
| Run length, time step | 5 h, 10 s | `nstop = 1800`, `dt = 10` |
| Emission flux | `MASSEMIS=3e-05` | `massemis` |
| Fire start | default (step 1080, 3 h) | not set (default) |

## ASP decks (`decks/`)

The ASP repository's `InorganicInputDecks/Williams/` baseline deck (asp_v2.2.2,
commit a83e398; emissions from the ACP 2015 supplement, 9-bin volatility
scheme with OH aging, TUV 5.4 `PhotoRates.in`), plus SAM-ASP-TUV's
`PhotoMap.in`, with one change: **10 aerosol size bins** instead of 40
(`AerosolModes.in`: cutoff 1.01 um, 10 bins, lower cutoff 0.025 um, as in the
other SAM-ASP-TUV decks). With 40 bins the deck needs 4096 SAM tracers and
fails in SAM at the start of ASP chemistry (sectional coagulation, "Total
Volume concentration is lost"); it runs in the standalone ASP parcel model.

Photolysis is computed online by TUV (SAM-ASP-TUV ignores the deck's
`PHOTO FLAG = 0`, which is for the parcel model).

**Tracers: 1486.** Build an image for it:

```bash
docker build --build-arg NTRACERS=1486 -t sam-asp-tuv:nt1486 .
```
