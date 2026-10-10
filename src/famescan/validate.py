"""Compare the pipeline with manual ChemStation percent reports of the same runs."""

import csv
import logging
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .chemstation import Report, ReportPeak, read_reports
from .config import Config
from .fames import Identification
from .io import find_runs
from .peaks import Peak
from .pipeline import Analysis, analyse
from .report import FamePeak

log = logging.getLogger(__name__)

CSV_FIELDS = [
    "data_file",
    "status",
    "kind",
    "ref_rt",
    "ref_type",
    "ref_area",
    "ref_pct_total",
    "rt",
    "type",
    "area",
    "pct_total",
    "area_ratio",
    "fame",
    "ecl",
    "base_mz",
    "reason",
]


@dataclass
class Match:
    data_file: str
    status: str
    """``reported`` (also a FAME peak of ours), ``rejected`` (detected, but no FAME by our
    rules), ``below min area``, ``not detected``, or ``extra`` (our FAME peak, not in the
    reference)."""
    ref: ReportPeak | None
    peak: Peak | None = None
    ident: Identification | None = None
    fame_peak: FamePeak | None = None
    rt: float | None = None

    @property
    def kind(self) -> str:
        """``auto`` for peaks integrated automatically by the RTE integrator, ``manual`` for hand
        integration (rM)."""
        if self.ref is None:
            return ""
        return "manual" if self.ref.type == "rM" else "auto"


@dataclass
class RunValidation:
    report: Report
    analysis: Analysis
    matches: list[Match]


def locate(report: Report, runs: dict[str, list[Path]]) -> Path | None:
    """The local ``.D`` folder of a report: same name, preferring the path whose folders match
    the end of the report's (Windows) data path best."""
    candidates = runs.get(report.data_file.lower(), [])
    if not candidates:
        return None
    parts = [p.lower() for p in re.split(r"[\\/]+", report.data_path) if p]

    def common_suffix(path: Path) -> int:
        n = 0
        for a, b in zip(
            reversed(parts), reversed([p.lower() for p in path.parent.parts]), strict=False
        ):
            if a != b:
                break
            n += 1
        return n

    return max(candidates, key=common_suffix)


def match_run(report: Report, analysis: Analysis, tolerance: float) -> list[Match]:
    run = analysis.run
    apex_rt = np.array([run.rt[p.apex] for p in analysis.peaks])
    reported = {id(fp.peak): fp for fp in analysis.report.peaks}
    matches, used = [], set()
    for ref in report.peaks:
        j = int(np.argmin(np.abs(apex_rt - ref.rt))) if len(apex_rt) else -1
        if j < 0 or abs(apex_rt[j] - ref.rt) > tolerance:
            matches.append(Match(report.data_file, "not detected", ref))
            continue
        p, ident = analysis.peaks[j], analysis.ids[j]
        fp = reported.get(id(p))
        if fp is not None:
            status = "reported"
            used.add(id(p))
        elif ident.fame is not None:
            status = "below min area"
        else:
            status = "rejected"
        matches.append(Match(report.data_file, status, ref, p, ident, fp, float(apex_rt[j])))
    matches += [
        Match(report.data_file, "extra", None, fp.peak, fp.ident, fp, fp.rt)
        for fp in analysis.report.peaks
        if id(fp.peak) not in used
    ]
    return matches


def validate(
    reference: Path, data_root: Path, cfg: Config | None = None, tolerance: float = 0.03
) -> list[RunValidation]:
    runs: dict[str, list[Path]] = defaultdict(list)
    for path in find_runs([data_root]):
        runs[path.name.lower()].append(path)
    results = []
    for report in read_reports(reference):
        path = locate(report, runs)
        if path is None:
            log.warning("%s: run not found under %s", report.data_file, data_root)
            continue
        analysis = analyse(path, cfg)
        results.append(RunValidation(report, analysis, match_run(report, analysis, tolerance)))
    return results


def reason_group(reason: str) -> str:
    if reason.startswith("highest ion"):
        return "highest ion is no FAME M+"
    return re.sub(r"(: .*| [\d.]+)$", "", reason)


def summary(results: list[RunValidation], tolerance: float = 0.03) -> list[str]:
    matches = [m for r in results for m in r.matches]
    refs = [m for m in matches if m.ref is not None]
    kinds = ("auto", "manual")
    statuses = ("reported", "rejected", "below min area", "not detected")
    counts = Counter((m.kind, m.status) for m in refs)
    total = Counter(m.kind for m in refs)
    lines = [
        f"{len(results)} reference reports, {len(refs)} peaks "
        f"({total['auto']} integrated automatically, {total['manual']} by hand)",
        f"{'':16s}{'auto':>8s}{'manual':>8s}",
        *(f"{s:16s}" + "".join(f"{counts[k, s]:8d}" for k in kinds) for s in statuses),
        f"{'FAME recall':16s}"
        + "".join(f"{100 * counts[k, 'reported'] / max(total[k], 1):7.1f}%" for k in kinds),
        "",
        "area ratio (ours / reference) of reported peaks, p5 p25 p50 p75 p95:",
    ]
    for k in kinds:
        ok = [m for m in refs if m.kind == k and m.status == "reported" and m.ref.area > 0]
        if not ok:
            continue
        ratio = np.array([m.fame_peak.area / m.ref.area for m in ok])
        q = " ".join(f"{v:.3f}" for v in np.percentile(ratio, [5, 25, 50, 75, 95]))
        within = 100 * np.mean(np.abs(ratio - 1) <= 0.05)
        lines.append(f"  {k:7s} {q}   ({within:.0f} % within 5 %)")
    ok = [m for m in refs if m.status == "reported"]
    if ok:
        drt = np.abs([m.rt - m.ref.rt for m in ok])
        lines.append(f"|dRT| median {np.median(drt):.3f} min, max {drt.max():.3f} min")
        # The manual reports mostly end around 27 min, ours cover the whole run. For a like-for-like
        # comparison, % total is also recomputed over our FAME peaks up to the last reference peak.
        per_fame, in_range = defaultdict(list), defaultdict(list)
        for r in results:
            last = max((p.rt for p in r.report.peaks), default=0.0) + tolerance
            total = sum(fp.area for fp in r.analysis.report.peaks if fp.rt <= last)
            for m in r.matches:
                if m.ref is None or m.status != "reported":
                    continue
                label = m.fame_peak.fame.label
                per_fame[label].append(abs(m.fame_peak.pct_total - m.ref.pct_total))
                in_range[label].append(abs(100 * m.fame_peak.area / total - m.ref.pct_total))
        lines += [
            "mean |d % total| per FAME (percentage points), ours over all FAME peaks of the run:",
            _per_label(per_fame),
            "  and over our FAME peaks up to the last reference peak of the run:",
            _per_label(in_range),
        ]
    lines += ["", "rejected reference peaks by reason:"]
    reasons = Counter(reason_group(m.ident.reason) for m in refs if m.status == "rejected")
    lines += [f"  {n:5d}  {reason}" for reason, n in reasons.most_common()]
    extras = [m for m in matches if m.status == "extra"]
    lines += ["", f"extra FAME peaks (not in the reference): {len(extras)}"]
    by_label = defaultdict(list)
    for m in extras:
        by_label[m.fame_peak.fame.label].append(m.fame_peak.pct_max)
    lines += [
        f"  {label:6s} n={len(v):4d}  % max median {np.median(v):5.2f}  max {max(v):6.2f}"
        for label, v in sorted(by_label.items(), key=lambda kv: _label_key(kv[0]))
    ]
    fields = Counter()
    for r in results:
        run, rep = r.analysis.run, r.report
        fields["Sample"] += run.sample == rep.sample
        fields["Acq On"] += run.acq_on == rep.acq_on
        fields["ALS Vial"] += str(run.als_vial) == rep.als_vial
    lines += [
        "",
        "header fields equal: " + ", ".join(f"{k} {v}/{len(results)}" for k, v in fields.items()),
    ]
    return lines


def _per_label(values: dict[str, list[float]]) -> str:
    return "  " + "  ".join(
        f"{label} {np.mean(v):.2f} (n={len(v)})"
        for label, v in sorted(values.items(), key=lambda kv: _label_key(kv[0]))
    )


def _label_key(label: str) -> tuple[int, int]:
    n, u = label[1:].split(":")
    return int(n), int(u)


def write_csv(results: list[RunValidation], path: Path) -> None:
    def fmt(value: float | None, digits: int) -> str:
        return "" if value is None else f"{value:.{digits}f}"

    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, CSV_FIELDS)
        writer.writeheader()
        for r in results:
            for m in r.matches:
                ref, p, ident, fp = m.ref, m.peak, m.ident, m.fame_peak
                fame = ident.fame if ident else None
                writer.writerow(
                    {
                        "data_file": m.data_file,
                        "status": m.status,
                        "kind": m.kind,
                        "ref_rt": fmt(ref and ref.rt, 3),
                        "ref_type": ref.type if ref else "",
                        "ref_area": ref.area if ref else "",
                        "ref_pct_total": fmt(ref and ref.pct_total, 3),
                        "rt": fmt(m.rt, 3),
                        "type": p.type if p else "",
                        "area": round(p.area) if p else "",
                        "pct_total": fmt(fp and fp.pct_total, 3),
                        "area_ratio": fmt(p.area / ref.area, 3) if p and ref and ref.area else "",
                        "fame": fame.label if fame else "",
                        "ecl": fmt(ident and ident.ecl, 2),
                        "base_mz": ident.base_mz if ident else "",
                        "reason": ident.reason if ident else "",
                    }
                )
