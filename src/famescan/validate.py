"""Compare famescan's FAME peaks with reference peaks of the same runs.

The reference is either a CSV file with one row per peak or a text file of ChemStation "Area
Percent Report" exports. CSV columns: ``data_file`` (run folder name, e.g. ``X.D``) and
``rt_min`` are required; ``area``, ``pct_total``, ``fame`` (label, e.g. ``C18:1``) and ``group``
are optional. famescan's own ``fame_peaks.csv`` is a valid reference, e.g. to compare two versions
or settings. ChemStation peaks are grouped into ``auto`` (RTE integrator) and ``manual`` (``rM``,
integrated by hand).

Every reference peak is matched to the detected peak with the nearest apex RT (within the
tolerance) and gets one status:

- ``found``: also a FAME peak of famescan's report (``other FAME`` if the reference names a
  different FAME),
- ``rejected``: detected, but not reported; ``reason`` says why,
- ``not detected``: no peak within the tolerance.

FAME peaks of famescan's reports that match no reference peak are ``extra``.
"""

import csv
import logging
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .chemstation import read_reports
from .config import Config
from .fames import Identification
from .io import find_runs
from .peaks import Peak
from .pipeline import Analysis, analyse
from .report import FamePeak, fmt, write_table

log = logging.getLogger(__name__)

STATUSES = ("found", "other FAME", "rejected", "not detected")
CSV_FIELDS = [
    "data_file",
    "group",
    "status",
    "reason",
    "ref_rt",
    "ref_fame",
    "ref_area",
    "ref_pct_total",
    "rt",
    "fame",
    "ecl",
    "type",
    "area",
    "pct_total",
    "area_ratio",
]


@dataclass
class RefPeak:
    data_file: str
    rt: float
    area: float | None = None
    pct_total: float | None = None
    fame: str = ""
    group: str = ""


@dataclass
class Match:
    status: str
    ref: RefPeak | None
    peak: Peak | None = None
    ident: Identification | None = None
    fame_peak: FamePeak | None = None
    rt: float | None = None
    reason: str = ""

    @property
    def group(self) -> str:
        return self.ref.group if self.ref else ""


@dataclass
class RunValidation:
    analysis: Analysis
    matches: list[Match]


def read_reference(path: Path) -> list[RefPeak]:
    if path.suffix.lower() == ".csv":
        with path.open(newline="", encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
        missing = {"data_file", "rt_min"} - set(rows[0] if rows else {})
        if missing:
            raise ValueError(f"{path}: missing column(s) {', '.join(sorted(missing))}")

        def num(row: dict[str, str], key: str) -> float | None:
            return float(row[key]) if row.get(key) else None

        return [
            RefPeak(
                r["data_file"],
                float(r["rt_min"]),
                num(r, "area"),
                num(r, "pct_total"),
                r.get("fame") or "",
                r.get("group") or "",
            )
            for r in rows
        ]
    return [
        RefPeak(
            report.data_file,
            p.rt,
            p.area,
            p.pct_total,
            group="manual" if "M" in p.type else "auto",
        )
        for report in read_reports(path)
        for p in report.peaks
    ]


def match_run(refs: list[RefPeak], analysis: Analysis, tolerance: float) -> list[Match]:
    run, report = analysis.run, analysis.report
    apex_rt = np.array([run.rt[p.apex] for p in analysis.peaks])
    reported = {fp.peak: fp for fp in report.peaks}
    matches, used = [], set()
    for ref in refs:
        j = int(np.argmin(np.abs(apex_rt - ref.rt))) if len(apex_rt) else -1
        if j < 0 or abs(apex_rt[j] - ref.rt) > tolerance:
            matches.append(Match("not detected", ref))
            continue
        p, ident = analysis.peaks[j], analysis.ids[j]
        fp = reported.get(p)
        reason = ""
        if fp is not None:
            used.add(p)
            status = "other FAME" if ref.fame and ref.fame != fp.fame.label else "found"
        else:
            status = "rejected"
            reason = report.dropped.get(p, ident.reason)
        matches.append(Match(status, ref, p, ident, fp, float(apex_rt[j]), reason))
    matches += [
        Match("extra", None, fp.peak, fp.ident, fp, fp.rt)
        for fp in report.peaks
        if fp.peak not in used
    ]
    return matches


def validate(
    reference: Path, data_root: Path, cfg: Config | None = None, tolerance: float = 0.03
) -> list[RunValidation]:
    """Analyse every run of the reference (found by folder name below ``data_root``) and match
    its peaks."""
    runs: dict[str, list[Path]] = defaultdict(list)
    for path in find_runs([data_root]):
        runs[path.name.lower()].append(path)
    by_run: dict[str, list[RefPeak]] = defaultdict(list)
    for ref in read_reference(reference):
        by_run[ref.data_file].append(ref)
    results = []
    for data_file, refs in by_run.items():
        paths = runs.get(data_file.lower(), [])
        if len(paths) != 1:
            found = "not found" if not paths else f"ambiguous ({len(paths)} runs)"
            log.warning("%s: %s under %s", data_file, found, data_root)
            continue
        analysis = analyse(paths[0], cfg)
        results.append(RunValidation(analysis, match_run(refs, analysis, tolerance)))
    return results


def reason_group(reason: str) -> str:
    """The reason without the peak-specific details."""
    if reason.startswith("highest ion"):
        return "highest ion is no FAME M+"
    return re.sub(r"(: .*| [\d.]+)$", "", reason)


def summary(results: list[RunValidation]) -> list[str]:
    matches = [m for r in results for m in r.matches]
    refs = [m for m in matches if m.ref is not None]
    groups = sorted({m.group for m in refs}) or [""]
    heads = [g or "all" for g in groups]
    counts = Counter((m.group, m.status) for m in refs)
    total = Counter(m.group for m in refs)

    def row(label: str, cells: list[str]) -> str:
        return f"  {label:20s}" + "".join(f"{c:>12s}" for c in cells)

    statuses = [s for s in STATUSES if s != "other FAME" or any(m.ref.fame for m in refs)]
    lines = [
        f"{len(results)} runs, {len(refs)} reference peaks",
        row("", heads),
        row("reference peaks", [str(total[g]) for g in groups]),
        *(row(s, [str(counts[g, s]) for g in groups]) for s in statuses),
        row("recall (found)", [f"{100 * counts[g, 'found'] / total[g]:.1f} %" for g in groups]),
    ]

    found = [m for m in refs if m.status in ("found", "other FAME")]
    with_area = [m for m in found if m.ref.area]
    if with_area:
        lines += ["", "area ratio famescan / reference of the found peaks:"]
        lines.append(row("", ["p5", "median", "p95", "within 5 %"]))
        for g, head in zip(groups, heads, strict=True):
            ratio = np.array([m.fame_peak.area / m.ref.area for m in with_area if m.group == g])
            if len(ratio):
                q = [f"{v:.3f}" for v in np.percentile(ratio, [5, 50, 95])]
                lines.append(row(head, [*q, f"{100 * np.mean(np.abs(ratio - 1) <= 0.05):.0f} %"]))
    if found:
        drt = np.abs([m.rt - m.ref.rt for m in found])
        lines += [
            "",
            f"|dRT| of the found peaks: median {np.median(drt):.3f}, max {drt.max():.3f} min",
        ]
    with_pct = [m for m in found if m.ref.pct_total is not None]
    if with_pct:
        per_fame = defaultdict(list)
        for m in with_pct:
            per_fame[m.fame_peak.fame].append(abs(m.fame_peak.pct_total - m.ref.pct_total))
        lines += ["", "mean |d % total| of the found peaks per FAME (percentage points):"]
        lines += [
            f"  {f.label:8s}{np.mean(v):6.2f}  (n={len(v)})" for f, v in sorted(per_fame.items())
        ]

    reasons = Counter(reason_group(m.reason) for m in refs if m.status == "rejected")
    if reasons:
        lines += ["", "rejected reference peaks by reason:"]
        lines += [f"  {n:5d}  {reason}" for reason, n in reasons.most_common()]

    extras = [m.fame_peak for m in matches if m.status == "extra"]
    lines += ["", f"extra FAME peaks (not in the reference): {len(extras)}"]
    by_fame = defaultdict(list)
    for fp in extras:
        by_fame[fp.fame].append(fp.pct_max)
    lines += [
        f"  {f.label:8s}n={len(v):5d}   % max median {np.median(v):5.2f}, max {max(v):6.2f}"
        for f, v in sorted(by_fame.items())
    ]
    return lines


def write_csv(results: list[RunValidation], path: Path) -> None:
    def row(data_file: str, m: Match) -> dict[str, str | int]:
        ref, p, ident, fp = m.ref, m.peak, m.ident, m.fame_peak
        fame = ident.fame if ident else None
        return {
            "data_file": data_file,
            "group": m.group,
            "status": m.status,
            "reason": m.reason,
            "ref_rt": fmt(ref and ref.rt, 3),
            "ref_fame": ref.fame if ref else "",
            "ref_area": fmt(ref and ref.area, 0),
            "ref_pct_total": fmt(ref and ref.pct_total, 3),
            "rt": fmt(m.rt, 3),
            "fame": fame.label if fame else "",
            "ecl": fmt(ident and ident.ecl, 2),
            "type": p.type if p else "",
            "area": round(p.area) if p else "",
            "pct_total": fmt(fp and fp.pct_total, 3),
            "area_ratio": fmt(p.area / ref.area, 3) if p and ref and ref.area else "",
        }

    write_table(path, CSV_FIELDS, (row(r.analysis.run.name, m) for r in results for m in r.matches))
