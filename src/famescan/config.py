"""Tunable parameters of the pipeline.

Integration defaults mirror ChemStation's ``rteint.p`` where that makes sense ("Filtering: 5",
"Peak Location: TOP"). Detection and identification thresholds were tuned with ``famescan validate``
against the manual reports.
"""

from dataclasses import dataclass, field


@dataclass(frozen=True)
class Config:
    # --- peak detection (TIC) ---
    smoothing_window: int = 5
    """Savitzky-Golay window in scans, used for detection and boundaries only."""
    noise_factor: float = 5.0
    """Minimum peak prominence as a multiple of the TIC noise."""
    min_width: float = 1.3
    """Minimum width in scans at half prominence (smoothed TIC). Narrower maxima are spikes, such
    as the one a saturated peak can end with."""
    bump_scans: int = 1
    """Peak boundaries are local minima of the raw TIC; rises of at most this many scans that are
    followed by a further fall (within ``bump_lookahead`` scans) are noise and are walked over."""
    bump_lookahead: int = 2
    boundary_max_frac: float = 0.5
    """Local minima above this fraction of the apex never end a peak."""
    bump_rel: float = 0.1
    drop_scans: int = 2
    """A rise smaller than ``bump_rel`` times the fall over the ``drop_scans`` scans before a local
    minimum is walked over as well (bumps on steep flanks, e.g. after a saturated peak)."""
    min_separation: int = 5
    """Apexes closer than this many scans are merged (ChemStation "Peak separation: 5")."""
    min_valley_depth: float = 0.05
    """Apexes are merged if the dip between them is below this fraction of the lower apex."""
    flat_top_frac: float = 0.8
    """Apexes are also merged if the TIC between them stays above this fraction of the higher
    apex: the flat top of a saturated peak."""
    apex_refine_scans: int = 2
    """The apex is moved to the raw TIC maximum within this many scans of the smoothed maximum."""

    # --- FAME identification ---
    min_m_plus_rel: float = 2.0
    """Minimum molecular ion abundance in % of the base peak."""
    max_higher_ion_rel: float = 2.0
    """Ions above M+2 must stay below this abundance (% of base peak)."""
    min_ion_mz: int = 150
    """Lowest m/z considered when looking for ions above the molecular ion."""
    siloxane_series: tuple[int, ...] = (73, 75, 89, 207, 327, 341)
    """Polysiloxane ions (column and septum bleed), each the start of a series in steps of 74
    (OSi(CH3)2), e.g. 207, 281, 355, 429. Siloxanes co-eluting with a FAME must not count as its
    highest ion, so series ions are skipped if they are not more abundant than m/z 73 (Si(CH3)3+,
    present in every siloxane spectrum), and so are their +1/+2 isotopes (at most half as
    abundant as the series ion)."""
    coelution_scans: int = 2
    """An ion co-elutes with a peak if its chromatogram peaks within this many scans of the TIC
    apex or in the upper half of the TIC peak."""
    main_ions: dict[int, dict[int, float]] = field(
        default_factory=lambda: {
            # double bonds -> {m/z: minimum % of base peak}
            0: {74: 50.0, 87: 20.0},
            1: {55: 25.0, 69: 25.0, 83: 25.0, 97: 25.0, 74: 15.0},
            2: {67: 40.0, 81: 40.0, 95: 40.0, 55: 15.0},
        }
    )
    """Main fragment ions per FAME class. The methyl ester ion m/z 74 is also required for
    mono-unsaturated FAMEs, which separates them from alkenes."""
    ethyl_ester_ions: dict[int, int] = field(default_factory=lambda: {88: 74, 101: 87})
    """Ethyl ester ions and their methyl ester counterparts: the McLafferty ion (88 instead of 74)
    and [CH2CH2COOR]+ (101 instead of 87). An ethyl ester has the M+ of the methyl ester with one
    more carbon (ethyl oleate = "C19:1"). In a methyl ester, 88 and 101 are weak (88 is the isotope
    of 87)."""
    ethyl_ester_ratio: float = 0.5
    """A peak is an ethyl ester, not a FAME, if each ethyl ester ion reaches this fraction of its
    methyl ester counterpart ..."""
    min_ethyl_ion_rel: float = 10.0
    """... and this % of the base peak."""
    ecl_window: dict[int, tuple[float, float]] = field(
        default_factory=lambda: {0: (-1.0, 0.5), 1: (-1.0, 0.5), 2: (-1.0, 1.0)}
    )
    """Allowed equivalent chain length minus chain length, per number of double bonds."""

    # --- report ---
    min_area_pct: float = 0.1
    """FAME peaks below this % of the largest FAME peak are dropped."""
