"""Peak detection and integration on the total ion chromatogram.

Modelled on ChemStation's integrators: peaks are found on a smoothed TIC, a peak ends where the
signal returns to baseline (``B``) or at the valley to a neighbouring peak (``V``). Peaks joined by
valleys share one straight baseline between the outer ``B`` points and are separated by vertical
drop lines. Areas are in abundance * seconds, like ChemStation's RTE integrator.
"""

from dataclasses import dataclass

import numpy as np
from scipy.signal import find_peaks, savgol_filter

from .config import Config, scans
from .io import scan_interval

APEX_REFINE_S = 1.4
"""The apex is moved to the raw TIC maximum within this distance of the smoothed maximum."""
MIN_VALLEY_DEPTH = 0.05
"""Apexes are merged if the dip between them is below this fraction of the lower apex ..."""
FLAT_TOP_FRAC = 0.8
"""... or if the TIC between them stays above this fraction of the higher apex: the flat top of a
saturated peak."""
BOUNDARY_MAX_FRAC = 0.5
"""Local minima above this fraction of the apex never end a peak."""
BUMP_S = 0.7
"""Rises of at most this duration after a local minimum that are followed by a further fall
(within ``BUMP_LOOKAHEAD_S``) are noise and are walked over ..."""
BUMP_LOOKAHEAD_S = 1.4
DROP_S = 1.4
"""... and so are rises smaller than ``BUMP_REL`` times the fall over this duration before the
minimum (bumps on steep flanks, e.g. after a saturated peak)."""
BUMP_REL = 0.1


@dataclass(eq=False)  # compared by identity, so peaks can be dict keys
class Peak:
    first: int
    """0-based scan index of the peak start (report scan numbers are index + 1)."""
    apex: int
    last: int
    start_type: str = "B"
    end_type: str = "B"
    base_first: float = 0.0
    """Baseline value at ``first``."""
    base_last: float = 0.0
    height: float = 0.0
    area: float = 0.0

    @property
    def type(self) -> str:
        return self.start_type + self.end_type

    @property
    def scan_numbers(self) -> tuple[int, int, int]:
        """First, apex and last scan as 1-based numbers, as in reports."""
        return self.first + 1, self.apex + 1, self.last + 1

    def baseline(self, scans: np.ndarray) -> np.ndarray:
        if self.last == self.first:
            return np.full(len(scans), self.base_first)
        return np.interp(scans, [self.first, self.last], [self.base_first, self.base_last])


def robust_sigma(x: np.ndarray) -> float:
    return 1.4826 * float(np.median(np.abs(x - np.median(x))))


def smooth(y: np.ndarray, window: int) -> np.ndarray:
    if window < 3:
        return y.astype(float)
    window |= 1  # savgol needs an odd window
    return savgol_filter(y, window, 2)


def detect_peaks(tic: np.ndarray, rt_min: np.ndarray, cfg: Config) -> list[Peak]:
    dt = scan_interval(rt_min)
    y = smooth(tic, scans(cfg.smoothing_window_s, dt))
    noise = robust_sigma(np.diff(tic)) / np.sqrt(2)
    found, _ = find_peaks(
        y, prominence=cfg.noise_factor * max(noise, 1.0), width=cfg.min_width_s / dt
    )
    if len(found) == 0:
        return []

    # Apex = raw maximum near the smoothed one; boundaries are searched on the raw TIC.
    k = scans(APEX_REFINE_S, dt)
    apexes = sorted({max(a - k, 0) + int(np.argmax(tic[max(a - k, 0) : a + k + 1])) for a in found})
    apexes = _merge_apexes(apexes, tic, cfg, dt)
    valleys = [
        a + int(np.argmin(tic[a : b + 1])) for a, b in zip(apexes[:-1], apexes[1:], strict=True)
    ]
    bounds = [0, *valleys, len(tic) - 1]
    peaks = [
        Peak(
            first=_walk_down(tic, apex, -1, bounds[i], dt),
            apex=apex,
            last=_walk_down(tic, apex, 1, bounds[i + 1], dt),
        )
        for i, apex in enumerate(apexes)
    ]

    # A valley is shared ("V") only if both neighbours run into it.
    for left, right, v in zip(peaks[:-1], peaks[1:], valleys, strict=True):
        if left.last == v and right.first == v:
            left.end_type = right.start_type = "V"

    _assign_baselines(peaks, tic)
    t_sec = rt_min * 60.0
    for p in peaks:
        integrate(p, tic, t_sec)
    return peaks


def _merge_apexes(apexes: list[int], y: np.ndarray, cfg: Config, dt: float) -> list[int]:
    """Merge neighbouring apexes that are too close or only separated by a negligible dip (as on
    the flat top of a saturated peak), keeping the higher one."""
    min_separation = scans(cfg.min_separation_s, dt)
    merged = [apexes[0]]
    for a in apexes[1:]:
        b = merged[-1]
        low, high = sorted((y[a], y[b]))
        valley = y[b : a + 1].min()
        if (
            a - b < min_separation
            or low - valley < MIN_VALLEY_DEPTH * low
            or valley > FLAT_TOP_FRAC * high
        ):
            if y[a] > y[b]:
                merged[-1] = a
        else:
            merged.append(a)
    return merged


def _walk_down(y: np.ndarray, i: int, step: int, bound: int, dt: float) -> int:
    """Follow the signal downhill from the apex ``i`` towards ``bound``; return the boundary scan.

    The walk stops at a local minimum unless the rise after it is noise (``BUMP_S``, ``DROP_S``).
    Minima in the upper part of the peak (e.g. dips on the flat top of a saturated peak) never
    stop it.
    """
    bump, lookahead, drop = (scans(s, dt) for s in (BUMP_S, BUMP_LOOKAHEAD_S, DROP_S))
    floor = BOUNDARY_MAX_FRAC * y[i]
    while True:
        while i != bound and y[i + step] < y[i]:
            i += step
        if i == bound:
            return i
        if y[i] > floor:
            top = _climb(y, i, step, bound)
            if top == bound:
                return i
            i = top
            continue
        # A rise that is small against the fall just before the minimum is part of the descent.
        back = i - drop * step
        if 0 <= back < len(y):
            hill = _climb(y, i, step, bound)
            if hill != bound and y[hill] - y[i] < BUMP_REL * (y[back] - y[i]):
                i = hill
                continue
        top = _climb(y, i, step, bound, bump)
        if top == bound or y[top + step] >= y[top]:
            return i  # the rise goes on: a real neighbour
        ahead = top
        for _ in range(lookahead):
            if ahead == bound:
                return i
            ahead += step
            if y[ahead] < y[i]:
                break
        else:
            return i
        i = top


def _climb(y: np.ndarray, i: int, step: int, bound: int, limit: int | None = None) -> int:
    """Follow the signal uphill from ``i`` towards ``bound``, at most ``limit`` scans; return
    where the rise ends."""
    start = i
    while i != bound and y[i + step] >= y[i] and (limit is None or abs(i - start) < limit):
        i += step
    return i


def _clusters(peaks: list[Peak]) -> list[list[Peak]]:
    clusters: list[list[Peak]] = []
    for p in peaks:
        if clusters and p.start_type == "V":
            clusters[-1].append(p)
        else:
            clusters.append([p])
    return clusters


def _assign_baselines(peaks: list[Peak], y: np.ndarray) -> None:
    """One straight baseline per cluster; where a valley dips below it, the cluster is split."""
    split = True
    while split:
        split = False
        for cluster in _clusters(peaks):
            start, end = cluster[0].first, cluster[-1].last
            depth, left, right = max(
                (
                    (float(np.interp(a.last, [start, end], [y[start], y[end]]) - y[a.last]), a, b)
                    for a, b in zip(cluster[:-1], cluster[1:], strict=True)
                ),
                default=(0.0, None, None),
                key=lambda t: t[0],
            )
            if depth > 0:
                left.end_type = right.start_type = "B"
                split = True
    for cluster in _clusters(peaks):
        start, end = cluster[0].first, cluster[-1].last
        for p in cluster:
            p.base_first, p.base_last = np.interp(
                [p.first, p.last], [start, end], [y[start], y[end]]
            )


def integrate(p: Peak, tic: np.ndarray, t_sec: np.ndarray) -> None:
    scans = np.arange(p.first, p.last + 1)
    corrected = tic[scans] - p.baseline(scans)
    p.area = float(np.trapezoid(corrected, t_sec[scans])) if len(scans) > 1 else 0.0
    p.height = float(tic[p.apex] - p.baseline(np.array([p.apex]))[0])
