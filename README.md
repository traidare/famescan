# famescan

Finds, identifies and integrates fatty acid methyl ester (FAME) peaks in EI
GC-MS runs, and writes FAME profiles, a peak table and ChemStation-style "Area
Percent Reports".

For each run, famescan

1. detects and integrates the peaks of the total ion chromatogram (TIC),
2. identifies FAMEs from each peak's mass spectrum: molecular ion, the main
   fragment ions of the FAME class, and the elution position,
3. computes each FAME's share of the largest FAME and of all FAMEs, optionally
   with an internal standard.

The input is Agilent ChemStation / MassHunter `.D` folders (`data.ms`, unit-mass
full scan). The peak detection thresholds, the FAME table, the rules per FAME
class and the internal standard are [settings](#settings).

## Setup

```sh
nix develop   # or direnv: cp envrc.example .envrc && direnv allow
uv sync       # creates .venv from pyproject.toml / uv.lock
```

## Usage

```sh
uv run famescan report PATH... [-o output] [settings] [--no-text-reports] [-v]
uv run famescan validate --reference FILE --data-root DIR [-o output] [settings] [--rt-tolerance 0.03]
uv run famescan config [settings]   # print the settings as TOML
```

For example:

```sh
uv run famescan report DATA_DIR              # all .D runs below DATA_DIR
uv run famescan report DATA_DIR/SUBDIR X.D   # some folders and/or runs
uv run famescan report DATA_DIR --config setup.toml --report-fames C16:0,C18:0,C18:1,C18:2
```

All paths are given explicitly. famescan only reads the runs and writes its
outputs to `output/`.

- A `PATH` is a `.D` folder or a folder searched recursively for `.D` runs.
  Blanks and alkane standards are processed like any other run.
- `-v` also logs why FAME candidates were rejected, e.g.
  `X.D 14.146 min: M+ 240 C14:1? main ions low: 69:18% 83:14% 97:14%`.

## Settings

famescan reads its settings from `famescan.toml` in the current directory, or
from a file given with `--config`.

| setting              | option                                      | meaning                                                              |
| -------------------- | ------------------------------------------- | -------------------------------------------------------------------- |
| `smoothing_window_s` |                                             | smoothing for peak detection                                         |
| `min_width_s`        |                                             | narrower maxima are spikes, not peaks                                |
| `min_separation_s`   |                                             | closer apexes are merged                                             |
| `noise_factor`       |                                             | minimum peak prominence, in multiples of the TIC noise               |
| `fame_table`         | `--fame-table CSV`                          | the FAMEs to look for                                                |
| `min_ion_pct`        |                                             | ions from this % of the base peak on are significant; M⁺ must be one |
| `classes`            |                                             | main ions and ECL window per number of double bonds                  |
| `report_fames`       | `--report-fames LABELS`                     | report only these FAMEs                                              |
| `min_area_pct`       |                                             | FAME peaks below this % of the largest FAME peak are dropped         |
| `standard`           | `--standard LABEL`, `--no-standard`         | the internal standard                                                |
| `standard_in_pct`    | `--standard-in-pct`, `--no-standard-in-pct` | whether the standard counts in % max and % total                     |

For example:

```toml
fame_table = "my_fames.csv"   # relative to this file
standard = "C17:0"
report_fames = ["C14:0", "C16:0", "C16:1", "C17:0", "C18:0", "C18:1", "C18:2"]

[classes.1]                   # one double bond; main_ions keeps its default
ecl_window = [-1.0, 1.0]
```

- Durations are in seconds. famescan converts them to whole scans per run, using
  the run's median scan interval. The defaults are set for 0.679 s per scan,
  e.g. `smoothing_window_s = 3.4` is 5 scans.
- A FAME table is a CSV file with the columns `chain_length`, `double_bonds`,
  `label`, `name` (optional) and `m_plus`, like
  `src/famescan/data/fame_molecular_ions.csv`. Each number of double bonds in it
  needs an entry in `classes`, e.g. `[classes.3]` with `main_ions` and
  `ecl_window` for C18:3.
- `report_fames` limits the report to some FAMEs; % max and % total are then
  relative to these. Identification still uses the whole FAME table.
- Without `standard_in_pct`, % max and % total are relative to the other FAME
  peaks, and so are the standard's own percentages.

## Output (`output/`)

| File                  | Content                                                                                             |
| --------------------- | --------------------------------------------------------------------------------------------------- |
| `fame_peaks.csv`      | One row per FAME peak (see below)                                                                   |
| `fame_profiles.csv`   | One row per run: % total per FAME (see below)                                                       |
| `<run>.txt`           | Percent report of one run. When two runs share a name, the name is prefixed with the parent folder. |
| `percent_reports.txt` | All text reports concatenated, readable by `famescan validate --reference`                          |
| `validation.csv`      | `famescan validate`: every reference peak with its match, plus famescan's extra FAME peaks          |

`--no-text-reports` leaves out the text reports.

`fame_peaks.csv` has these columns:

- run: `data_file`, `sample`, `acq_on`
- peak: `peak`, `rt_min`, `first_scan`, `max_scan`, `last_scan`, `type`,
  `height`, `area`
- percentages: `pct_max`, `pct_total`
- identification: `fame` (e.g. `C18:1`), `name`, `m_plus`, `m_plus_rel_pct` (M+
  in % of the base peak), `ecl` (equivalent chain length), `is_standard`

`fame_profiles.csv` has `data_file`, `sample`, `acq_on` and one column per FAME
found in any run, ordered by chain length and double bonds. A cell is the %
total of all peaks of that FAME in the run, empty if there are none. The
internal standard has a column of its own, e.g. `C19:0 (standard)`.

The text reports have the layout of ChemStation's:

- CRLF line endings, the same peak table columns and widths, and a blank line
  after every 5th peak.
- Scan numbers are 1-based; R.T. is the apex scan.
- Areas are abundance · s, heights are above the baseline.
- The parameter block lists famescan's own parameters instead of `rteint.p`.
- Peak types are `BB`/`BV`/`VB`/`VV`, without ChemStation's `r` prefix, because
  they don't come from the RTE integrator.

## How it works

The values below are the defaults; [Settings](#settings) lists those that can be
changed.

### Peaks and integration (`peaks.py`)

- Detection works on the TIC (sum over all m/z), smoothed with a Savitzky–Golay
  filter over 3.4 s (5 scans at 0.679 s, ChemStation's "Filtering: 5").
  - A maximum is a peak if its prominence is ≥ 5 × the TIC noise and its width
    is ≥ 0.88 s. Narrower maxima are spikes, such as the single scan a saturated
    peak can end with.
  - Apexes are merged if they are less than 3.4 s apart or separated by a dip of
    less than 5 %. They are also merged if the signal between them stays above
    80 % of the higher apex: the flat, dented top of a peak that saturated the
    detector.
- Boundaries are local minima of the raw TIC, reached by walking downhill from
  the apex. One-scan noise bumps are walked over, and so are small bumps on
  steep flanks. A minimum shared with the neighbouring peak is a valley (`V`);
  otherwise the boundary is a baseline point (`B`).
- Baseline: peaks joined by valleys share one straight baseline between the
  outer `B` points and are separated by vertical drop lines. If a valley lies
  below that line, the cluster is split there into two clusters, each with its
  own baseline.
- The area is the trapezoidal integral of TIC − baseline over time in seconds,
  as ChemStation's RTE integrator computes it.

### FAME identification (`fames.py`)

The spectrum of a peak is the mean over the apex ± 0.7 s (± 1 scan) minus the
background, interpolated between the peak boundaries, in % of the base peak. A
peak is a FAME if all three rules hold:

1. **Molecular ion.** The highest significant ion that co-elutes with the peak
   must be the M⁺ of a FAME in the FAME table, or one of its isotopes M+1/M+2
   (at most half as abundant as M⁺).
   - The built-in table covers the methyl esters of C10–C24 with 0, 1 or 2
     double bonds, M⁺ = 14n + 46 − 2u.
   - Significant means m/z ≥ 150 and ≥ 2 % of the base peak. M⁺ itself must
     reach ≥ 2 %.
   - Co-eluting means the ion's chromatogram peaks within ± 1.4 s of the apex,
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
   87 and stays far below 74.

3. **Retention.** The peak's equivalent chain length (ECL), interpolated between
   the saturated FAMEs of the same run, must be close to its chain length: ECL −
   n within −1 … +0.5, or −1 … +1 for two double bonds. The anchors are the
   saturated FAMEs of the largest total area whose RT increases with chain
   length. A run with fewer than two saturated FAMEs skips this rule.

FAME peaks below 0.1 % of the largest FAME peak are dropped. Then:

- % max is relative to the largest FAME peak, % total to the sum of all FAME
  peaks; by default both include the internal standard;
- the largest peak of the standard FAME (C19:0) is marked `is_standard`; other
  peaks with its label are ordinary FAME peaks.

## Validation

`famescan validate` runs the pipeline on the runs of a reference and compares
the FAME peaks:

- **Reference**: a CSV file with one row per peak, or a text file of ChemStation
  "Area Percent Report" exports.
  - CSV columns: `data_file` (the run's folder name) and `rt_min` are required;
    `area`, `pct_total`, `fame` and `group` are optional. famescan's own
    `fame_peaks.csv` qualifies, e.g. to compare two versions or settings.
  - ChemStation peaks are grouped by their peak type into `auto` (integrated by
    the RTE integrator) and `manual` (`rM`, integrated by hand).
- **Matching**: each reference peak is matched to the detected peak with the
  nearest apex RT, within `--rt-tolerance` (0.03 min). Its status is
  - `found`: also a FAME peak of famescan's report (`other FAME` if the
    reference names a different FAME),
  - `rejected`: detected, but not reported; `reason` says why,
  - `not detected`.

  famescan's FAME peaks without a reference peak are `extra`.

- **Summary**: the counts per status and group; for the found peaks, the area
  ratio, |ΔRT| and |Δ % total| per FAME; the rejected peaks by reason; the extra
  peaks per FAME. `validation.csv` has one row per reference peak and per extra
  peak.
