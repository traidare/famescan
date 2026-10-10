"""FAME percent reports in the layout of ChemStation's "Area Percent Report", and a peak table.

As in the manual reports, % max is relative to the largest FAME peak and % total to the sum of all
FAME peaks, the internal standard (C19:0) included.
"""

import csv
import os
from collections.abc import Iterable
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from . import __version__
from .config import Config
from .fames import Fame, Identification
from .io import Run
from .peaks import Peak

CRLF = "\r\n"
CSV_FIELDS = [
    "data_file",
    "sample",
    "acq_on",
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

    @property
    def sum_area(self) -> int:
        return sum(p.area for p in self.peaks)


def build_report(
    run: Run,
    peaks: list[Peak],
    ids: list[Identification],
    cfg: Config,
    labels: frozenset[str] | None = None,
) -> FameReport:
    """Select the identified FAME peaks (optionally only the given labels, e.g. ``C16:0``) that
    reach ``cfg.min_area_pct`` of the largest one, and compute their percentages."""
    found = [
        (p, i)
        for p, i in zip(peaks, ids, strict=True)
        if i.fame and (labels is None or i.fame.label in labels)
    ]
    if not found:
        return FameReport(run, [], cfg)
    largest = max(p.area for p, _ in found)
    found = [(p, i) for p, i in found if p.area > 0 and p.area >= cfg.min_area_pct / 100 * largest]
    areas = [round(p.area) for p, _ in found]
    amax, total = max(areas), sum(areas)
    report = FameReport(
        run,
        [
            FamePeak(
                p, i, float(run.rt[p.apex]), round(p.height), a, 100 * a / amax, 100 * a / total
            )
            for (p, i), a in zip(found, areas, strict=True)
        ],
        cfg,
    )
    standards = [fp for fp in report.peaks if fp.fame.is_standard]
    if standards:
        max(standards, key=lambda fp: fp.area).is_standard = True
    return report


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
        _pair("Smoothing", "ON", "Filtering", str(cfg.smoothing_window)),
        _pair(
            "Threshold",
            f"{cfg.noise_factor:g} x noise",
            "Min Area",
            f"{cfg.min_area_pct:g} % of largest FAME",
        ),
        _pair("Min Width", f"{cfg.min_width:g} scans", "Peak Location", "TOP"),
        "",
        "  Baseline: one straight line per peak cluster, drop lines at valleys",
        f"  Peak separation: {cfg.min_separation}",
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
                n,
                fp.rt,
                p.first + 1,
                p.apex + 1,
                p.last + 1,
                p.type,
                fp.height,
                fp.area,
                fp.pct_max,
                fp.pct_total,
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


def csv_rows(report: FameReport) -> list[dict[str, str | int]]:
    run = report.run
    rows: list[dict[str, str | int]] = []
    for n, fp in enumerate(report.peaks, 1):
        p, f = fp.peak, fp.fame
        rows.append(
            {
                "data_file": run.name,
                "sample": run.sample,
                "acq_on": run.acq_on,
                "peak": n,
                "rt_min": f"{fp.rt:.3f}",
                "first_scan": p.first + 1,
                "max_scan": p.apex + 1,
                "last_scan": p.last + 1,
                "type": p.type,
                "height": fp.height,
                "area": fp.area,
                "pct_max": f"{fp.pct_max:.2f}",
                "pct_total": f"{fp.pct_total:.3f}",
                "fame": f.label,
                "name": f.name,
                "m_plus": f.m_plus,
                "m_plus_rel_pct": f"{fp.ident.m_plus_rel:.1f}",
                "ecl": "" if fp.ident.ecl is None else f"{fp.ident.ecl:.2f}",
                "is_standard": int(fp.is_standard),
            }
        )
    return rows


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
    reports: list[FameReport], outdir: Path, generated: datetime | None = None
) -> list[Path]:
    """Write ``<run>.txt`` per run, all reports concatenated into ``percent_reports.txt`` and the
    peak table ``fame_peaks.csv``. Returns the paths written."""
    outdir.mkdir(parents=True, exist_ok=True)
    generated = generated or datetime.now()
    texts = [render(r, generated) for r in reports]
    written = []
    for name, text in zip(report_names(reports), texts, strict=True):
        written.append(outdir / name)
        written[-1].write_bytes(text.encode("latin-1", errors="replace"))
    combined = outdir / "percent_reports.txt"
    combined.write_bytes("".join(texts).encode("latin-1", errors="replace"))
    table = outdir / "fame_peaks.csv"
    with table.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, CSV_FIELDS)
        writer.writeheader()
        for r in reports:
            writer.writerows(csv_rows(r))
    return [*written, combined, table]
