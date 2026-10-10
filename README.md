# famescan

Automated FAME area percent reports from Agilent GC/MSD `.D` data.

This does in code what was done by hand in ChemStation (G1701EA):

1. find the peaks of the total ion chromatogram (TIC),
2. check each peak's mass spectrum for a fatty acid methyl ester (FAME)
   molecular ion and main ions,
3. integrate the FAME peaks,
4. export an "Area Percent Report" per run, in ChemStation's layout.

## Setup

Tools come from Nix, Python packages from uv:

```sh
nix develop   # or direnv: cp envrc.example .envrc && direnv allow
uv sync       # creates .venv from pyproject.toml / uv.lock
```

## Usage

```sh
uv run famescan report PATH... [-o output] [--fames LABELS] [-v]
uv run famescan validate --reference REPORTS.txt --data-root DIR [-o output] [--rt-tolerance 0.03]
```

For example:

```sh
uv run famescan report DATA_DIR              # all .D runs below DATA_DIR
uv run famescan report DATA_DIR/SUBDIR X.D   # some folders and/or runs
```

All paths are given explicitly. The runs are only read, never written; the
outputs go to `output/`.

- A `PATH` is a `.D` folder or a folder searched recursively for `.D` runs.
  Blanks and alkane standards are processed like any other run.
- `--fames C14:0,C16:0,C16:1,C18:0,C18:1,C18:2,C19:0` reports only these FAMEs.
  % max and % total are then relative to these. By default every FAME of the
  table is reported.
- `-v` also logs why FAME candidates were rejected, e.g.
  `X.D 14.146 min: M+ 240 C14:1? main ions low: 69:18% 83:14% 97:14%`.

### Output (`output/`)

| File                  | Content                                                                                             |
| --------------------- | --------------------------------------------------------------------------------------------------- |
| `<run>.txt`           | Percent report of one run. When two runs share a name, the name is prefixed with the parent folder. |
| `percent_reports.txt` | All reports concatenated, the format `famescan validate --reference` reads                          |
| `fame_peaks.csv`      | One row per FAME peak (see below)                                                                   |
| `validation.csv`      | `famescan validate`: every reference peak with our match, plus our extra FAME peaks                 |

The text reports have the same layout as ChemStation's:

- CRLF line endings, the same peak table columns and widths, and a blank line
  after every 5th peak.
- Scan numbers are 1-based; R.T. is the apex scan.
- Areas are abundance · s, heights are above the baseline.
- The parameter block lists famescan's own parameters instead of `rteint.p`.
- Peak types are `BB`/`BV`/`VB`/`VV`, without ChemStation's `r` prefix, because
  they don't come from the RTE integrator.

`fame_peaks.csv` has these columns:

- run: `data_file`, `sample`, `acq_on`
- peak: `peak`, `rt_min`, `first_scan`, `max_scan`, `last_scan`, `type`,
  `height`, `area`
- percentages: `pct_max`, `pct_total`
- identification: `fame` (e.g. `C18:1`), `name`, `m_plus`, `m_plus_rel_pct` (M+
  in % of the base peak), `ecl` (equivalent chain length), `is_standard`

## How it works

All thresholds are in `src/famescan/config.py`.

### Peaks and integration (`peaks.py`)

- **Detection** works on the TIC (sum over m/z 40–550), smoothed with a 5-point
  Savitzky–Golay filter ("Filtering: 5").
  - A maximum is a peak if its prominence is ≥ 5 × the TIC noise and its width
    is ≥ 1.3 scans. Narrower maxima are spikes, such as the single scan a
    saturated peak can end with.
  - Apexes are merged if they are less than 5 scans apart ("Peak separation: 5")
    or separated by a dip of less than 5 %. They are also merged if the signal
    between them stays above 80 % of the higher apex: the flat, dented top of a
    peak that saturated the detector.
- **Boundaries** are local minima of the raw TIC, reached by walking downhill
  from the apex. One-scan noise bumps are walked over, and so are small bumps on
  steep flanks. A minimum shared with the neighbouring peak is a valley (`V`);
  otherwise the boundary is a baseline point (`B`).
- **Baseline**: peaks joined by valleys share one straight baseline between the
  outer `B` points and are separated by vertical drop lines. If a valley lies
  below that line, the cluster is split there into two clusters, each with its
  own baseline.
- **Area** is the trapezoidal integral of TIC − baseline over time in seconds,
  as ChemStation's RTE integrator computes it.

### FAME identification (`fames.py`)

The spectrum of a peak is the mean of apex ± 1 scan minus the background,
interpolated between the peak boundaries, in % of the base peak. A peak is a
FAME if all three rules hold:

1. **Molecular ion.** The highest significant ion that co-elutes with the peak
   must be the M⁺ of a FAME in `src/famescan/data/fame_molecular_ions.csv`, or
   one of its isotopes M+1/M+2 (at most half as abundant as M⁺).
   - The table covers C10–C24 with 0, 1 or 2 double bonds, M⁺ = 14n + 46 − 2u,
     and C19:0 as the internal standard.
   - Significant means m/z ≥ 150 and ≥ 2 % of the base peak. M⁺ itself must
     reach ≥ 2 %.
   - Co-eluting means the ion's chromatogram peaks within ± 2 scans of the apex,
     or in the upper half of the TIC peak.
   - Column bleed is ignored: polysiloxane ions, i.e. the series 73, 75, 89,
     207, 327, 341 plus k · 74 (e.g. 207, 281, 355, 429), and their isotopes.
     They are ignored only while they are not more abundant than m/z 73
     (Me₃Si⁺).
2. **Main ions** of the class implied by M⁺, in % of the base peak:

   | double bonds | required ions                                                                        |
   | ------------ | ------------------------------------------------------------------------------------ |
   | 0            | 74 ≥ 50, 87 ≥ 20                                                                     |
   | 1            | 55, 69, 83, 97 ≥ 25 and 74 ≥ 15 (the methyl ester ion; separates FAMEs from alkenes) |
   | 2            | 67, 81, 95 ≥ 40 and 55 ≥ 15                                                          |

   The ester must also be a methyl ester, not an ethyl ester. An ethyl ester has
   the M⁺ of the methyl ester with one more carbon: ethyl oleate (C₂₀H₃₈O₂,
   M⁺ 310) looks like "C19:1". Its ester ions are 88 and 101 instead of 74
   and 87. A peak is rejected as an ethyl ester if 88 and 101 each reach ≥ 10 %
   and at least half of 74 and 87. In methyl esters, 88 is only the isotope of
   87; across all FAME peaks, 88/74 stays below 0.25 (99th percentile).

3. **Retention.** The peak's equivalent chain length (ECL), interpolated between
   the saturated FAMEs of the same run, must be close to its chain length: ECL −
   n within −1 … +0.5, or −1 … +1 for two double bonds. The anchors are the
   saturated FAMEs of the largest total area whose RT increases with chain
   length.

FAME peaks below 0.1 % of the largest FAME peak are dropped. As in the manual
reports:

- % max is relative to the largest FAME peak;
- % total is relative to the sum of all FAME peaks, including the C19:0
  standard;
- only the largest C19:0 peak is marked `is_standard`.

## Validation

`famescan validate` compares the pipeline with manual ChemStation reports of the
same runs, given as one text file of concatenated "Area Percent Report" exports.
The numbers below are for the 210 manual reports of the thesis runs.

The manual reports contain only the peaks the analyst kept as FAMEs: 4–19 per
report, 2221 in total. They were made in two steps:

1. ChemStation's RTE integrator integrated the TIC automatically, and the
   analyst deleted the peaks that were no FAMEs.
2. The analyst integrated by hand the FAMEs the integrator had missed. 89 % of
   these are below its minimum area of 3 % of the largest peak. A few large
   peaks were also re-integrated by hand, e.g. detector-saturated ones.

The reports name no FAMEs. Whether a peak is one, and which, is the analyst's
judgement.

famescan integrates every peak automatically. The reference peaks are matched to
its peaks by apex RT (± 0.03 min).

| peaks of the manual reports                      | integrated automatically (`rBV`, `rVB`, `rVV`) | integrated by hand (`rM`) |
| ------------------------------------------------ | ---------------------------------------------: | ------------------------: |
| total                                            |                                           1249 |                       972 |
| also FAME peaks in famescan's report             |                                  1148 (91.9 %) |              677 (69.7 %) |
| detected, but not a FAME by the rules above      |                                            101 |                       285 |
| not detected                                     |                                              0 |                        10 |
| start / end scan identical                       |                                    78 % / 83 % |               34 % / 23 % |
| area ratio famescan / reference, median (p5–p95) |                            0.997 (0.975–1.016) |       1.007 (0.908–1.180) |
| areas within ± 5 %                               |                                           97 % |                      63 % |

- The retention times agree: median |ΔRT| is 0.000 min, the maximum 0.011 min.
- Sample, Acq On and ALS Vial match in all 210 reports.
- For the automatically integrated peaks, famescan reproduces the RTE
  integrator's boundaries and areas.
- The peak type letters often differ even where the boundaries agree. For
  example, 361 `rBV` peaks are `BB` in famescan.
- For the hand-integrated peaks, famescan's boundaries take the place of the
  analyst's:
  - famescan's peaks are a median of 3 scans wider, mostly on the tail.
  - Where famescan's area is more than 5 % larger (26 % of these peaks), its
    peak is wider in 92 % of the cases.
  - The agreement depends on peak size. Within ± 5 % are 96 % of the peaks of ≥
    5 % of the largest peak, 62 % of those of 1–5 %, and 47 % of those below 1
    %.
- 7 of the 10 undetected peaks are shoulders (1–2 % of the largest peak) on the
  tail of the large C18:1 peak (24.73–24.79 min) or C18:0 peak (25.33 min).
  - The analyst split them off by hand. famescan finds no maximum of their own,
    so their area stays in the large peak, which becomes 2–3 % larger.
  - Of the other three, one lies within a peak famescan rejects (11.07 min). Two
    are too small to detect: 0.33 % of the largest peak at 19.43 min, and 0.16 %
    at 35.50 min.

Mean |Δ % total|, in percentage points, for the main FAMEs:

|                                                                    | C16:0 | C18:0 | C18:1 | C18:2 | C19:0 |
| ------------------------------------------------------------------ | ----: | ----: | ----: | ----: | ----: |
| famescan % total over all its FAME peaks                           |  0.67 |  0.42 |  1.47 |  0.21 |  0.74 |
| famescan % total over its FAME peaks up to the last reference peak |  0.55 |  0.38 |  1.26 |  0.19 |  0.63 |

The remaining difference comes from the peaks one side reports and the other
does not.

### Differences to the manual reports

Most rejected reference peaks are not FAMEs by the rules above. The labels below
are not from the manual reports. They are the FAME of the table whose M⁺ the
peak's spectrum shows, i.e. the FAME it would have to be. The spectra show:

- **14.1–14.5 min, M⁺ of "C14:1" (183 peaks).** An alkane (57/71/85):
  heptadecane has the same M⁺ 240 as C14:1. The ions of a mono-unsaturated FAME
  are missing. Where they are present (55/69/83/97, 74, M−32 = 208), the peak is
  accepted (7 peaks).
- **18.6–18.7 min, "C13:2" (49 peaks).** Alkene-like spectra without 67/81/95.
- **No FAME molecular ion as highest ion (47 peaks).** Mostly siloxanes
  (237/163/73/147/355/281, e.g. at 26.7, 28.7 and 30.2 min) and an alkane at
  34.3 min (base peak 57). In some manual reports these made up 25–71 % of the
  total.
- **Main ions too low.** "C19:1" at 25.2–25.3 min (32), "C16:0" at 15.8–18.2 min
  (19), "C20:0" at 27.0–27.4 min (14), and a few others. The "C19:1" at
  25.2–25.3 min has base peak 88 with 101 and M−32/M−31 (278/279), but no 74:
  neither a normal FAME nor an ethyl ester.
- **Ethyl ester (1 peak).** Ethyl oleate at 25.9 min, see below.

famescan also reports 1777 FAME peaks that are not in the manual reports:

- With them, a famescan report has a median of 17 FAME peaks instead of 11.
- Most are small: the median is below 2.3 % of the largest FAME for every FAME.
- 712 of them elute after 27.2 min, where the manual reports end. The
  spot-checked spectra were textbook: C20:0 (M⁺ 326, M−43 = 283), C20:1 (324,
  M−32 = 292), C22:0 (354, 311), C24:0 (382, 339).
- In two of the four `W_` runs, the manual reports leave out the C19:0 standard
  (26.9 min, 11 % of the total).

Ethyl esters are rejected (rule 2). Without that check, famescan would report 69
more peaks:

- ethyl oleate at 25.9 min (ECL 18.4, 43 runs) as "C19:1";
- ethyl nonadecanoate, the ethyl ester of the C19:0 standard's acid, at 27.85
  min (ECL 19.7, 19 runs) as "C20:0" (M⁺ 326);
- six "C19:0" peaks at 25.9 min and one "C18:0" at 21.9 min, where an ethyl
  ester overlaps a methyl ester.

Ethyl esters form when ethanol is present during the methylation. Peaks at 13.4,
20.8 and 25.8 min with base peak 88, 101 and 157 (M⁺ 256, 284, 312) look like
ethyl esters of C14:0, C16:0 and C18:0 as well, but lack M−45. Without 74 they
fail the main-ion rule anyway.

### Limitations

- **Isomers share the label of their M⁺.** Branched (iso/anteiso) and other
  isomers get the straight-chain label, so one label can appear more than once
  in a run, e.g.:
  - "C15:0" at ECL 14.6–14.7;
  - "C19:0" at 25.9, 26.1 and 26.45 min (ECL 18.4, 18.5 and 18.7), besides the
    standard at 26.9 min;
  - "C20:0" at 27.7 min (ECL 19.5).

  The `ecl` column in the CSV tells them apart.

- **The table ends at C24.** Longer-chain FAMEs therefore go unreported, e.g.
  C26:0/C27:0 at 33.14/33.37 min and C28:0 at 34.97 min in some of the runs.
  - To extend the table, add rows to `fame_molecular_ions.csv`.
  - Check the result with `-v`: in one run a C22:0 with a weak co-eluting ion at
    m/z 396 would become "C25:0".

- **A weak ion above M⁺ hides a FAME.** Rule 1 is strict, so an ion of only 2–7
  % above M⁺ (usually background at high oven temperature) makes a FAME
  unrecognisable. Across all runs this hides about 15 small peaks, e.g. C21:0 at
  29.4 min, C22:0 at 30.36 min and C24:0 at 32.1 min.

## Development

```sh
uv run pytest                    # FAMESCAN_DATA_ROOT=DATA_DIR also runs tests/test_data.py
ruff format . && ruff check .    # the ruff from the Nix shell, not `uv run ruff`
```
