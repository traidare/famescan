"""FAME identification from the mass spectrum of a TIC peak.

A peak is a fatty acid methyl ester if

1. the highest significant ion that co-elutes with the peak (ignoring siloxane bleed) is a
   molecular ion M+ of the FAME table (``Config.fame_table``, by default
   ``data/fame_molecular_ions.csv``) or one of its isotopes M+1, M+2,
2. the main fragment ions of the FAME class implied by M+ are present (saturated: 74/87,
   mono-unsaturated: 55/69/83/97 plus the methyl ester ion 74, di-unsaturated: 67/81/95/55)
   and the ester is a methyl, not an ethyl ester (88/101 instead of 74/87), and
3. it elutes where a FAME of that chain length can elute: its equivalent chain length (ECL),
   interpolated between the saturated FAMEs found in the same run, is close to its chain length.

The rules of 2 and 3 are per number of double bonds (``Config.classes``), so every number in the
table needs them.
"""

import csv
import io
from dataclasses import dataclass
from functools import cache
from importlib.resources import files
from pathlib import Path

import numpy as np

from .config import Config, scans
from .io import Run
from .peaks import Peak

SPECTRUM_HALFWIDTH_S = 0.7
"""The peak spectrum is the mean over the apex +- this duration (at least the apex scan)."""
BACKGROUND_S = 2.0
"""At a baseline boundary, the background spectrum is the mean over this duration outside the
peak (at least the boundary scan)."""
COELUTION_S = 1.4
"""An ion co-elutes with a peak if its chromatogram peaks within this distance of the TIC apex or
in the upper half of the TIC peak."""
MIN_ION_MZ = 150
"""Lowest m/z considered when looking for the highest ion."""
SILOXANE_SERIES = (73, 75, 89, 207, 327, 341)
"""Polysiloxane ions (column and septum bleed), each the start of a series in steps of 74
(OSi(CH3)2), e.g. 207, 281, 355, 429. Siloxanes co-eluting with a FAME must not count as its
highest ion, so series ions are skipped if they are not more abundant than m/z 73 (Si(CH3)3+,
present in every siloxane spectrum), and so are their +1/+2 isotopes (at most half as abundant as
the series ion)."""
ETHYL_ESTER_IONS = {88: 74, 101: 87}
"""Ethyl ester ions and their methyl ester counterparts: the McLafferty ion (88 instead of 74) and
[CH2CH2COOR]+ (101 instead of 87). An ethyl ester has the M+ of the methyl ester with one more
carbon (ethyl oleate = "C19:1"). In a methyl ester, 88 and 101 are weak (88 is the isotope of
87). A peak is an ethyl ester if each ethyl ester ion reaches ``ETHYL_ESTER_RATIO`` of its
counterpart and ``MIN_ETHYL_ION_PCT`` of the base peak."""
ETHYL_ESTER_RATIO = 0.5
MIN_ETHYL_ION_PCT = 10.0


@dataclass(frozen=True, order=True)  # ordered by chain length, then double bonds
class Fame:
    chain_length: int
    double_bonds: int
    label: str
    """e.g. ``C18:1``"""
    name: str
    m_plus: int


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

    def reject(self, why: str) -> "Identification":
        """Reject the candidate FAME, stating why."""
        assert self.candidate is not None
        self.fame = None
        self.reason = f"M+ {self.candidate.m_plus} {self.candidate.label}? {why}"
        return self


@cache
def fame_table(path: str = "") -> tuple[Fame, ...]:
    """The FAMEs of a CSV table; ``""``: the built-in one."""
    if path:
        text = Path(path).read_text(encoding="utf-8")
    else:
        text = files("famescan").joinpath("data/fame_molecular_ions.csv").read_text("utf-8")
    try:
        return tuple(
            Fame(
                chain_length=int(r["chain_length"]),
                double_bonds=int(r["double_bonds"]),
                label=r["label"],
                name=r.get("name") or "",
                m_plus=int(r["m_plus"]),
            )
            for r in csv.DictReader(io.StringIO(text))
        )
    except (KeyError, TypeError, ValueError) as e:
        raise ValueError(f"FAME table {path or '(built-in)'}: bad or missing column {e}") from e


def check_config(cfg: Config) -> None:
    """Raise ``ValueError`` if the settings don't fit the FAME table."""
    table = fame_table(cfg.fame_table)
    labels = {f.label for f in table}
    unknown = sorted(set(cfg.report_fames) - labels)
    if unknown:
        raise ValueError(f"unknown FAME label(s): {', '.join(unknown)}")
    if cfg.standard and cfg.standard not in labels:
        raise ValueError(f"the standard {cfg.standard} is not in the FAME table")
    missing = sorted({f.double_bonds for f in table} - set(cfg.classes))
    if missing:
        raise ValueError(f"classes has no entry for {', '.join(map(str, missing))} double bond(s)")


def peak_spectrum(run: Run, peak: Peak) -> np.ndarray:
    """Apex spectrum (mean over apex +- ``SPECTRUM_HALFWIDTH_S``) minus the background
    interpolated between the peak boundaries, in % of the base peak. At a baseline boundary the
    background is the mean over ``BACKGROUND_S`` from the boundary outwards; at a valley it is the
    valley scan."""
    ms = run.intensity
    half = scans(SPECTRUM_HALFWIDTH_S, run.scan_interval)
    background_scans = max(1, scans(BACKGROUND_S, run.scan_interval))
    lo, hi = max(peak.first, peak.apex - half), min(peak.last, peak.apex + half)
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


def _is_bleed(m: int, rel: np.ndarray, col: dict[int, int]) -> bool:
    """Polysiloxane ion or one of its isotopes (see ``SILOXANE_SERIES``)."""

    def siloxane(x: int) -> bool:
        return (
            x in col
            and any(x >= s and (x - s) % 74 == 0 for s in SILOXANE_SERIES)
            and rel[col[x]] <= rel[col[73]]
        )

    return siloxane(m) or any(
        siloxane(m - k) and rel[col[m]] <= 0.5 * rel[col[m - k]] for k in (1, 2)
    )


def identify(run: Run, peak: Peak, cfg: Config) -> Identification:
    rel = peak_spectrum(run, peak)
    mz, col = run.mz, run.column
    base_mz = int(mz[int(np.argmax(rel))])
    by_mass = {f.m_plus: f for f in fame_table(cfg.fame_table)}
    coelution = scans(COELUTION_S, run.scan_interval)

    def ab(m: int) -> float:
        """Abundance of m/z ``m`` in % of the base peak."""
        return float(rel[col[m]])

    significant = [
        i
        for i in np.flatnonzero((rel >= cfg.min_ion_pct) & (mz >= MIN_ION_MZ))[::-1]
        if not _is_bleed(int(mz[i]), rel, col) and _coelutes(run, peak, i, coelution)
    ]
    if not significant:
        return Identification(
            None, None, 0.0, base_mz, None, f"no co-eluting ion >= m/z {MIN_ION_MZ}"
        )
    top = int(mz[significant[0]])

    # The highest ion is M+ or one of its isotopes M+1, M+2, which are much weaker than M+.
    candidates = [
        by_mass[m]
        for m in range(top, top - 3, -1)
        if m in by_mass
        and m in col
        and ab(m) >= cfg.min_ion_pct
        and all(ab(i) <= 0.5 * ab(m) for i in range(m + 1, top + 1))
    ]
    if not candidates:
        return Identification(None, None, 0.0, base_mz, top, f"highest ion m/z {top} is no FAME M+")
    fame = candidates[0]
    ident = Identification(fame, fame, ab(fame.m_plus), base_mz, top)

    main = cfg.classes[fame.double_bonds].main_ions.items()
    missing = [f"{m}:{ab(m):.0f}%" for m, need in main if ab(m) < need]
    if missing:
        return ident.reject(f"main ions low: {' '.join(missing)}")

    ethyl = ETHYL_ESTER_IONS.items()
    if all(ab(e) >= max(MIN_ETHYL_ION_PCT, ETHYL_ESTER_RATIO * ab(m)) for e, m in ethyl):
        ions = " ".join(f"{e}:{ab(e):.0f}% ({m}:{ab(m):.0f}%)" for e, m in ethyl)
        return ident.reject(f"ethyl ester: {ions}")
    return ident


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
            lo, hi = cfg.classes[i.fame.double_bonds].ecl_window
            if not lo <= i.ecl - i.fame.chain_length <= hi:
                i.reject(f"elutes at ECL {i.ecl:.1f}")
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
