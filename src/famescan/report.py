"""Per-run FAME reports: a peak table, a FAME profile per run, and text reports in the layout of
ChemStation's "Area Percent Report".

% max is relative to the largest FAME peak and % total to the sum of all FAME peaks. Whether the
internal standard counts is ``Config.standard_in_pct`` (by default it does).
"""

import csv
import os
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

from . import __version__
from .config import Config
from .fames import Fame, Identification
from .io import Run
from .peaks import Peak

CRLF = "\r\n"
RUN_FIELDS = ["data_file", "sample", "acq_on"]
CSV_FIELDS = [
    *RUN_FIELDS,
    "peak",
    "rt_min",
    "first_scan",
    "max_scan",
    "last_scan",
    "type",
    "height",
    "area",
    "pct_max",
    "pct_total",
    "fame",
    "name",
    "m_plus",
    "m_plus_rel_pct",
    "ecl",
    "is_standard",
]


@dataclass
class FamePeak:
    peak: Peak
    ident: Identification
    rt: float
    height: int
    area: int
    pct_max: float
    pct_total: float
    is_standard: bool = False
    """Only the largest peak of the standard FAME in a run is the internal standard."""

    @property
    def fame(self) -> Fame:
        assert self.ident.fame is not None
        return self.ident.fame


@dataclass
class FameReport:
    run: Run
    peaks: list[FamePeak]
    config: Config
    dropped: dict[Peak, str] = field(default_factory=dict)
    """Identified FAME peaks that are not in the report, and why."""

    @property
    def sum_area(self) -> int:
        return sum(p.area for p in self.peaks)


def build_report(run: Run, peaks: list[Peak], ids: list[Identification], cfg: Config) -> FameReport:
    """Select the identified FAME peaks (only ``cfg.report_fames`` if given) that reach
    ``cfg.min_area_pct`` of the largest one, mark the internal standard and compute the
    percentages."""
    dropped: dict[Peak, str] = {}
    found = []
    for p, i in zip(peaks, ids, strict=True):
        if i.fame and cfg.report_fames and i.fame.label not in cfg.report_fames:
            dropped[p] = f"{i.fame.label} not in report_fames"
        elif i.fame:
            found.append((p, i))
    largest = max((p.area for p, _ in found), default=0.0)
    fps = []
    for p, i in found:
        if p.area > 0 and p.area >= cfg.min_area_pct / 100 * largest:
            fps.append(
                FamePeak(p, i, float(run.rt[p.apex]), round(p.height), round(p.area), 0.0, 0.0)
            )
        else:
            dropped[p] = f"{i.fame.label} below min area"
    if not fps:
        return FameReport(run, [], cfg, dropped)
    standards = [fp for fp in fps if cfg.standard and fp.fame.label == cfg.standard]
    if standards:
        max(standards, key=lambda fp: fp.area).is_standard = True
    basis = [fp.area for fp in fps if cfg.standard_in_pct or not fp.is_standard] or [
        fp.area for fp in fps
    ]
    amax, total = max(basis), sum(basis)
    for fp in fps:
        fp.pct_max, fp.pct_total = 100 * fp.area / amax, 100 * fp.area / total
    return FameReport(run, fps, cfg, dropped)


def format_row(
    n: int,
    rt: float,
    first: int,
    apex: int,
    last: int,
    peak_type: str,
    height: int,
    area: int,
    pct_max: float,
    pct_total: float,
) -> str:
    """One peak row; scans are the 1-based numbers shown in reports."""
    return (
        f"{n:3d} {rt:7.3f} {first:5d} {apex:4d} {last:4d} {peak_type:<4s} {height:7d} {area:9d} "
        f"{pct_max:6.2f}% {pct_total:7.3f}%"
    )


def _pair(label: str, value: str, right_label: str, right_value: str) -> str:
    """Two parameters on one line, the right label's colon in column 55 as in ChemStation."""
    left = f"  {label:<10}: {value}"
    return f"{left}{' ' * max(1, 54 - len(left) - len(right_label))}{right_label}: {right_value}"


def render(report: FameReport, generated: datetime | None = None) -> str:
    """The report as text with CRLF line endings, ending with an empty line like the
    ChemStation exports (so that reports can simply be concatenated)."""
    run, cfg = report.run, report.config
    tool = f"famescan {__version__}"
    lines = [
        " " * 52 + "Area Percent Report",
        "",
        f"  Data Path : {run.path.absolute().parent}{os.sep}",
        f"  Data File : {run.name:<52s}",
        f"  Acq On    : {run.acq_on}",
        "  Operator  :  ",
        f"  Sample    : {run.sample}",
        "  Misc      :  ",
        f"  ALS Vial  : {run.als_vial}   Sample Multiplier: 1",
        "",
        f"  Integration Parameters: {tool}",
        f"  Integrator: {tool}",
        _pair("Smoothing", "ON", "Filtering", f"{cfg.smoothing_window_s:g} s"),
        _pair(
            "Threshold",
            f"{cfg.noise_factor:g} x noise",
            "Min Area",
            f"{cfg.min_area_pct:g} % of largest FAME",
        ),
        _pair("Min Width", f"{cfg.min_width_s:g} s", "Peak Location", "TOP"),
        "",
        "  Baseline: one straight line per peak cluster, drop lines at valleys",
        f"  Peak separation: {cfg.min_separation_s:g} s",
        "",
        f"  Method    : {tool}",
        "  Title     : FAME peaks (molecular ion and main ions)",
        "",
        f"  Signal     : TIC: {run.name}\\data.ms",
        "",
        " peak  R.T. first  max last  PK   peak      corr.   corr.    % of",
        "   #   min   scan scan scan  TY  height     area    % max.   total",
        " ---  ----- ----- ---- ---- ---  -------   -------  ------  -------",
    ]
    for n, fp in enumerate(report.peaks, 1):
        p = fp.peak
        lines.append(
            format_row(
                n, fp.rt, *p.scan_numbers, p.type, fp.height, fp.area, fp.pct_max, fp.pct_total
            )
        )
        if n % 5 == 0:
            lines.append(" ")
    stamp = (generated or datetime.now()).strftime("%a %b %d %H:%M:%S %Y")
    lines += [
        " ",
        " ",
        " " * 24 + f"Sum of corrected areas:{report.sum_area:12d}",
        "",
        f"{tool} {stamp}   ",
        "",
        "",
    ]
    return CRLF.join(lines)


def fmt(value: float | None, digits: int) -> str:
    """A number for CSV files; ``None`` is an empty cell."""
    return "" if value is None else f"{value:.{digits}f}"


def run_fields(run: Run) -> dict[str, str]:
    """The ``RUN_FIELDS`` of a run."""
    return {"data_file": run.name, "sample": run.sample, "acq_on": run.acq_on}


def write_table(path: Path, columns: list[str], rows: Iterable[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, columns)
        writer.writeheader()
        writer.writerows(rows)


def csv_rows(report: FameReport) -> list[dict[str, str | int]]:
    rows: list[dict[str, str | int]] = []
    for n, fp in enumerate(report.peaks, 1):
        p, f = fp.peak, fp.fame
        first, apex, last = p.scan_numbers
        rows.append(
            run_fields(report.run)
            | {
                "peak": n,
                "rt_min": fmt(fp.rt, 3),
                "first_scan": first,
                "max_scan": apex,
                "last_scan": last,
                "type": p.type,
                "height": fp.height,
                "area": fp.area,
                "pct_max": fmt(fp.pct_max, 2),
                "pct_total": fmt(fp.pct_total, 3),
                "fame": f.label,
                "name": f.name,
                "m_plus": f.m_plus,
                "m_plus_rel_pct": fmt(fp.ident.m_plus_rel, 1),
                "ecl": fmt(fp.ident.ecl, 2),
                "is_standard": int(fp.is_standard),
            }
        )
    return rows


def profile_rows(reports: list[FameReport]) -> tuple[list[str], list[dict[str, str]]]:
    """The FAME profile of each run: % total per FAME label (summed over the label's peaks), the
    internal standard in a column of its own. Returns the columns and the rows."""
    labels: dict[str, Fame] = {}
    rows = []
    for r in reports:
        pct: dict[str, float] = defaultdict(float)
        for fp in r.peaks:
            label = f"{fp.fame.label} (standard)" if fp.is_standard else fp.fame.label
            labels.setdefault(label, fp.fame)
            pct[label] += fp.pct_total
        rows.append(run_fields(r.run) | {k: fmt(v, 3) for k, v in pct.items()})
    order = sorted(labels, key=lambda k: (k.endswith("(standard)"), labels[k]))
    return [*RUN_FIELDS, *order], rows


def report_names(reports: Iterable[FameReport]) -> list[str]:
    """File names for per-run reports: the run name without ``.D``, prefixed with the parent
    folder's name where two runs would otherwise share a name."""
    reports = list(reports)
    stems = [r.run.path.stem for r in reports]
    return [
        f"{r.run.path.absolute().parent.name}_{s}.txt" if stems.count(s) > 1 else f"{s}.txt"
        for r, s in zip(reports, stems, strict=True)
    ]


def write_outputs(
    reports: list[FameReport],
    outdir: Path,
    generated: datetime | None = None,
    text_reports: bool = True,
) -> list[Path]:
    """Write the peak table ``fame_peaks.csv``, the profiles ``fame_profiles.csv`` and, with
    ``text_reports``, ``<run>.txt`` per run plus all of them concatenated into
    ``percent_reports.txt``. Returns the paths written, the CSV files first."""
    outdir.mkdir(parents=True, exist_ok=True)
    peaks = outdir / "fame_peaks.csv"
    write_table(peaks, CSV_FIELDS, (row for r in reports for row in csv_rows(r)))
    profiles = outdir / "fame_profiles.csv"
    write_table(profiles, *profile_rows(reports))
    written = [peaks, profiles]
    if not text_reports:
        return written
    generated = generated or datetime.now()
    texts = [render(r, generated) for r in reports]
    for name, text in zip(report_names(reports), texts, strict=True):
        written.append(outdir / name)
        written[-1].write_bytes(text.encode("latin-1", errors="replace"))
    combined = outdir / "percent_reports.txt"
    combined.write_bytes("".join(texts).encode("latin-1", errors="replace"))
    return [*written, combined]
