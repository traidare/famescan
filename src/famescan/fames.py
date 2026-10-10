"""FAME identification from the mass spectrum of a TIC peak.

A peak is a fatty acid methyl ester if

1. the highest significant ion that co-elutes with the peak (ignoring siloxane bleed) is a
   molecular ion M+ from ``data/fame_molecular_ions.csv`` or one of its isotopes M+1, M+2,
2. the main fragment ions of the FAME class implied by M+ are present (saturated: 74/87,
   mono-unsaturated: 55/69/83/97 plus the methyl ester ion 74, di-unsaturated: 67/81/95/55)
   and the ester is a methyl, not an ethyl ester (88/101 instead of 74/87), and
3. it elutes where a FAME of that chain length can elute: its equivalent chain length (ECL),
   interpolated between the saturated FAMEs found in the same run, is close to its chain length.
"""

import csv
import io
from dataclasses import dataclass
from functools import cache
from importlib.resources import files

import numpy as np

from .config import Config
from .io import Run
from .peaks import Peak


@dataclass(frozen=True)
class Fame:
    chain_length: int
    double_bonds: int
    label: str
    """e.g. ``C18:1``"""
    name: str
    m_plus: int
    is_standard: bool


@dataclass
class Identification:
    fame: Fame | None
    """The identified FAME, ``None`` if the peak is rejected."""
    candidate: Fame | None
    """FAME whose M+ matched, also for rejected peaks."""
    m_plus_rel: float
    """Molecular ion abundance in % of the base peak."""
    base_mz: int
    top_mz: int | None
    """Highest significant co-eluting ion (m/z >= ``min_ion_mz``)."""
    reason: str = ""
    """Why the peak was rejected; empty if accepted."""
    ecl: float | None = None
    """Equivalent chain length, if the run has enough saturated FAMEs to compute it."""


@cache
def fame_table() -> tuple[Fame, ...]:
    text = files("famescan").joinpath("data/fame_molecular_ions.csv").read_text(encoding="utf-8")
    return tuple(
        Fame(
            chain_length=int(r["chain_length"]),
            double_bonds=int(r["double_bonds"]),
            label=r["label"],
            name=r["name"],
            m_plus=int(r["m_plus"]),
            is_standard=r["is_standard"] == "1",
        )
        for r in csv.DictReader(io.StringIO(text))
    )


def peak_spectrum(run: Run, peak: Peak, background_scans: int = 3) -> np.ndarray:
    """Apex spectrum (mean of apex +-1 scan) minus the background interpolated between the peak
    boundaries, in % of the base peak. At a baseline boundary the background is the mean of the
    boundary scan and up to ``background_scans - 1`` scans outside; at a valley it is the valley."""
    ms = run.intensity
    lo, hi = max(peak.first, peak.apex - 1), min(peak.last, peak.apex + 1)
    top = ms[lo : hi + 1].mean(axis=0)

    def side(scan: int, kind: str, outward: int) -> np.ndarray:
        if kind == "V":
            return ms[scan]
        outer = min(max(scan + outward * (background_scans - 1), 0), len(ms) - 1)
        a, b = sorted((scan, outer))
        return ms[a : b + 1].mean(axis=0)

    left = side(peak.first, peak.start_type, -1)
    right = side(peak.last, peak.end_type, 1)
    w = (peak.apex - peak.first) / (peak.last - peak.first) if peak.last > peak.first else 0.5
    spec = np.clip(top - ((1 - w) * left + w * right), 0, None)
    base = spec.max()
    return 100 * spec / base if base > 0 else spec


def _coelutes(run: Run, peak: Peak, column: int, tolerance: int) -> bool:
    """The ion's chromatogram peaks near the TIC apex or within the upper half of the TIC peak
    (which matters for the broad flat top of a saturated peak)."""
    top = peak.first + int(np.argmax(run.intensity[peak.first : peak.last + 1, column]))
    if abs(top - peak.apex) <= tolerance:
        return True
    lo, hi = sorted((top, peak.apex))
    scans = np.arange(lo, hi + 1)
    return bool(np.all(run.tic[scans] - peak.baseline(scans) >= 0.5 * peak.height))


def _is_bleed(m: int, rel: np.ndarray, col: dict[int, int], series: tuple[int, ...]) -> bool:
    """Polysiloxane ion or one of its isotopes (see ``Config.siloxane_series``)."""

    def siloxane(x: int) -> bool:
        return (
            x in col
            and any(x >= s and (x - s) % 74 == 0 for s in series)
            and rel[col[x]] <= rel[col[73]]
        )

    return siloxane(m) or any(
        siloxane(m - k) and rel[col[m]] <= 0.5 * rel[col[m - k]] for k in (1, 2)
    )


def identify(run: Run, peak: Peak, cfg: Config) -> Identification:
    rel = peak_spectrum(run, peak)
    mz = run.mz
    col = {int(m): i for i, m in enumerate(mz)}
    base_mz = int(mz[int(np.argmax(rel))])
    by_mass = {f.m_plus: f for f in fame_table()}

    significant = [
        i
        for i in np.flatnonzero((rel >= cfg.max_higher_ion_rel) & (mz >= cfg.min_ion_mz))[::-1]
        if not _is_bleed(int(mz[i]), rel, col, cfg.siloxane_series)
        and _coelutes(run, peak, i, cfg.coelution_scans)
    ]
    if not significant:
        return Identification(
            None, None, 0.0, base_mz, None, f"no co-eluting ion >= m/z {cfg.min_ion_mz}"
        )
    top = int(mz[significant[0]])

    # The highest ion is M+ or one of its isotopes M+1, M+2, which are much weaker than M+.
    candidates = [
        by_mass[m]
        for m in range(top, top - 3, -1)
        if m in by_mass
        and m in col
        and rel[col[m]] >= cfg.min_m_plus_rel
        and all(rel[col[i]] <= 0.5 * rel[col[m]] for i in range(m + 1, top + 1))
    ]
    if not candidates:
        return Identification(None, None, 0.0, base_mz, top, f"highest ion m/z {top} is no FAME M+")
    fame = candidates[0]
    m_rel = float(rel[col[fame.m_plus]])

    missing = [
        f"{m}:{rel[col[m]]:.0f}%"
        for m, need in cfg.main_ions[fame.double_bonds].items()
        if rel[col[m]] < need
    ]
    if missing:
        return Identification(
            None,
            fame,
            m_rel,
            base_mz,
            top,
            f"M+ {fame.m_plus} {fame.label}? main ions low: {' '.join(missing)}",
        )

    ethyl = cfg.ethyl_ester_ions.items()
    if all(
        rel[col[e]] >= max(cfg.min_ethyl_ion_rel, cfg.ethyl_ester_ratio * rel[col[m]])
        for e, m in ethyl
    ):
        ions = " ".join(f"{e}:{rel[col[e]]:.0f}% ({m}:{rel[col[m]]:.0f}%)" for e, m in ethyl)
        return Identification(
            None, fame, m_rel, base_mz, top, f"M+ {fame.m_plus} {fame.label}? ethyl ester: {ions}"
        )
    return Identification(fame, fame, m_rel, base_mz, top)


def identify_run(run: Run, peaks: list[Peak], cfg: Config) -> list[Identification]:
    """Identify all peaks of a run, including the retention plausibility check (rule 3)."""
    ids = [identify(run, p, cfg) for p in peaks]
    anchors = _anchors(
        [
            (run.rt[p.apex], i.fame.chain_length, p.area)
            for p, i in zip(peaks, ids, strict=True)
            if i.fame and i.fame.double_bonds == 0
        ]
    )
    if len(anchors) < 2:
        return ids
    for p, i in zip(peaks, ids, strict=True):
        if i.candidate is None:
            continue
        i.ecl = _ecl(run.rt[p.apex], anchors)
        if i.fame:
            lo, hi = cfg.ecl_window[i.fame.double_bonds]
            if not lo <= i.ecl - i.fame.chain_length <= hi:
                i.reason = f"M+ {i.fame.m_plus} {i.fame.label}? elutes at ECL {i.ecl:.1f}"
                i.fame = None
    return ids


def _anchors(saturated: list[tuple[float, int, float]]) -> list[tuple[float, int]]:
    """Retention anchors (RT, chain length): the saturated FAMEs with the largest total area whose
    retention times increase with chain length, so that a stray misidentification cannot bend
    the scale."""
    pts = sorted(saturated, key=lambda t: (t[1], t[0]))
    best: list[tuple[float, int]] = []  # (total area of the best chain ending here, predecessor)
    for i, (rt, n, area) in enumerate(pts):
        total, prev = area, -1
        for j in range(i):
            if pts[j][1] < n and pts[j][0] < rt and best[j][0] + area > total:
                total, prev = best[j][0] + area, j
        best.append((total, prev))
    i = max(range(len(pts)), key=lambda k: best[k][0], default=-1)
    chain = []
    while i >= 0:
        chain.append((pts[i][0], pts[i][1]))
        i = best[i][1]
    return chain[::-1]


def _ecl(rt: float, anchors: list[tuple[float, int]]) -> float:
    """Piecewise linear interpolation of the chain length over RT, extrapolated at both ends."""
    rts = [a[0] for a in anchors]
    i = min(max(int(np.searchsorted(rts, rt)) - 1, 0), len(anchors) - 2)
    (r0, n0), (r1, n1) = anchors[i], anchors[i + 1]
    return n0 + (rt - r0) * (n1 - n0) / (r1 - r0)
