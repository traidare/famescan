import numpy as np
import pytest
from fakes import DT, retention_times

from famescan.config import Config
from famescan.peaks import detect_peaks

N = 600
SCANS = np.arange(N)
RT = retention_times(N)


def gaussian(height: float, apex: float, sigma: float) -> np.ndarray:
    return height * np.exp(-0.5 * ((SCANS - apex) / sigma) ** 2)


def gaussian_area(height: float, sigma: float) -> float:
    return height * sigma * DT * np.sqrt(2 * np.pi)  # abundance * s


def chromatogram(*components: tuple[float, float, float], seed: int = 1) -> np.ndarray:
    """Gaussian peaks (height, apex, sigma in scans) on a sloped baseline with noise."""
    rng = np.random.default_rng(seed)
    return 2e5 + 200 * SCANS + sum(gaussian(*c) for c in components) + rng.normal(0, 2e3, N)


def test_noise_only():
    assert detect_peaks(chromatogram(), RT, Config()) == []


def test_isolated_peak():
    (p,) = detect_peaks(chromatogram((5e6, 150, 3.0)), RT, Config())
    assert (p.apex, p.type) == (150, "BB")
    assert p.first < 150 - 3 * 3 and p.last > 150 + 3 * 3
    assert p.height == pytest.approx(5e6, rel=0.01)
    assert p.area == pytest.approx(gaussian_area(5e6, 3.0), rel=0.02)


def test_overlapping_peaks_share_a_drop_line():
    left, right = detect_peaks(chromatogram((1e7, 300, 3.0), (4e6, 315, 3.0)), RT, Config())
    assert (left.apex, left.type, right.apex, right.type) == (300, "BV", 315, "VB")
    assert left.last == right.first
    assert left.base_last == pytest.approx(right.base_first)
    assert left.area == pytest.approx(gaussian_area(1e7, 3.0), rel=0.02)
    assert right.area == pytest.approx(gaussian_area(4e6, 3.0), rel=0.02)


def test_separated_peaks_have_own_baselines():
    a, b = detect_peaks(chromatogram((3e6, 450, 3.0), (2e6, 490, 3.0)), RT, Config())
    assert (a.type, b.type) == ("BB", "BB")
    assert a.last < b.first
    assert a.area == pytest.approx(gaussian_area(3e6, 3.0), rel=0.02)
    assert b.area == pytest.approx(gaussian_area(2e6, 3.0), rel=0.02)


# TIC of scans 1690-1800 of a measured run: saturated C18:1 and C18:0 peaks with flat,
# dented tops, each ending with a one-scan spike (1771, 1795) and an abrupt drop.
SATURATED = [
    748473, 728962, 729572, 723550, 720499, 705171, 706447, 669538, 647155, 624431,
    624061, 654601, 710155, 859673, 1163218, 1687013, 2372863, 3228571, 4043896, 4781886,
    5224775, 5613317, 5896682, 6166690, 6558517, 7627253, 8853296, 10114278, 11875576, 13371835,
    15518369, 17601577, 19692609, 21937968, 24669107, 27824964, 30296060, 32968124, 35286537,
    36917214, 38600482, 39361538, 41327504, 43031886, 44470017, 44766909, 44665798, 44273509,
    43785144, 43755984, 43363471, 43263577, 45367959, 44420506, 46794826, 45872934, 43904157,
    43178379, 42246911, 42833368, 41089536, 38534135, 37754492, 37437520, 35225987, 34405854,
    32686347, 31723397, 31006629, 31758659, 29259145, 28208885, 26775037, 26370593, 25926815,
    25487181, 25659017, 24574858, 23884525, 23632497, 23674559, 35578653, 2910712, 1736169,
    1622567, 1524443, 1496693, 1728928, 1838631, 1530537, 1561572, 2278326, 4432018, 8298656,
    13765750, 19333490, 24813649, 30777574, 34372578, 37520155, 40106402, 41887900, 43773475,
    44831571, 44938403, 34137505, 3421401, 4771464, 2630921, 1123979, 1033172,
]  # fmt: skip


def test_saturated_peaks():
    rng = np.random.default_rng(1)
    pad = 200
    tic = np.concatenate(
        [7.3e5 + rng.normal(0, 1.5e4, pad), SATURATED, 1.0e6 + rng.normal(0, 1.5e4, pad)]
    )
    peaks = detect_peaks(tic, retention_times(len(tic)), Config())

    def scan(i: int) -> int:  # report scan number in the original run
        return i - pad + 1690

    assert not {1771, 1795} & {scan(p.apex) for p in peaks}  # spikes are no peaks
    c181, c180 = (p for p in peaks if p.height > 1e7)  # flat tops are not split
    # expected C18:1 peak: scans 1700 1744 1776, area 1266155539
    assert (scan(c181.first), scan(c181.apex), scan(c181.last)) == (1700, 1744, 1776)
    assert c181.area == pytest.approx(1266155539, rel=0.02)
    # expected C18:0 peak: scans 1779 1794 (1800), area 284085326
    assert (scan(c180.first), scan(c180.apex)) == (1779, 1794)
    assert c180.area == pytest.approx(284085326, rel=0.03)
