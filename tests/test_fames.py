import numpy as np
import pytest
from fakes import MZ, make_run

from famescan.config import Config, FameClass
from famescan.fames import _anchors, _coelutes, _ecl, check_config, fame_table, identify
from famescan.io import Run
from famescan.peaks import Peak, integrate

N = 41
APEX = 20

# EI spectra (% of base peak) like those of the runs
METHYL_PALMITATE = {  # C16:0, M+ 270: McLafferty ion 74, 87 (88 = its isotope), 143, M-43, M-31
    41: 30, 43: 40, 55: 30, 57: 20, 69: 10, 74: 100, 75: 15, 87: 70, 88: 3, 101: 4, 129: 10,
    143: 12, 171: 5, 185: 4, 199: 5, 213: 4, 227: 10, 239: 6, 270: 8, 271: 1.5,
}  # fmt: skip
METHYL_OLEATE = {  # C18:1, M+ 296: hydrocarbon series 55/69/83/97, 74, M-32
    41: 70, 55: 100, 69: 60, 74: 35, 83: 45, 87: 25, 96: 20, 97: 35, 110: 10, 123: 8, 137: 5,
    180: 3, 222: 5, 264: 12, 265: 4, 296: 4, 297: 0.9,
}  # fmt: skip
METHYL_LINOLEATE = {  # C18:2, M+ 294: 67/81/95
    41: 50, 55: 55, 67: 100, 68: 35, 79: 30, 81: 85, 82: 40, 95: 60, 96: 25, 109: 20, 110: 15,
    123: 10, 150: 6, 164: 4, 178: 3, 220: 2, 263: 5, 294: 12, 295: 2.4,
}  # fmt: skip
HEPTADECANE = {  # C17H36, M+ 240 = M+ of C14:1, but no ester ions (the 14.15 min peak)
    41: 35, 43: 85, 55: 20, 57: 100, 71: 60, 85: 40, 99: 12, 113: 9, 127: 6, 141: 5, 155: 4,
    169: 3, 183: 3, 197: 2, 211: 2.5, 240: 5, 241: 0.9,
}  # fmt: skip
ETHYL_OLEATE = {  # C18:1 ethyl ester, M+ 310 = M+ of C19:1: 88/101 instead of 74/87, M-46, M-45
    # (measured at 25.89 min; 74 and 87 come from co-eluting methyl esters)
    41: 50, 43: 48, 55: 100, 57: 39, 67: 38, 69: 81, 74: 34, 81: 43, 83: 75, 84: 60, 87: 26,
    88: 68, 96: 59, 97: 68, 98: 52, 101: 46, 110: 32, 111: 40, 123: 24, 180: 23, 222: 25,
    236: 23, 264: 38, 265: 31, 310: 7, 311: 1.6,
}  # fmt: skip
SILOXANE = {  # column bleed: Me3Si+ 73 and the 207/281/355/429 series with isotopes
    45: 8, 73: 100, 75: 15, 147: 12, 207: 55, 208: 11, 209: 9, 221: 8, 281: 35, 282: 8, 283: 6,
    355: 15, 356: 4, 357: 3, 429: 4,
}  # fmt: skip


def spectrum(ions: dict[int, float]) -> np.ndarray:
    s = np.zeros(len(MZ))
    for m, v in ions.items():
        s[m - MZ[0]] = v
    return s


def run_with(*components: tuple[dict[int, float], float, int]) -> tuple[Run, Peak]:
    """A run of compounds (spectrum, base peak height, apex scan) with Gaussian elution profiles
    on a small flat background, and the TIC peak around scan ``APEX``."""
    scans = np.arange(N)
    intensity = np.full((N, len(MZ)), 20.0)
    for ions, height, apex in components:
        profile = np.exp(-0.5 * ((scans - apex) / 3.0) ** 2)
        intensity += height * np.outer(profile, spectrum(ions) / 100)
    run = make_run(intensity)
    peak = Peak(
        APEX - 12, APEX, APEX + 12, base_first=run.tic[APEX - 12], base_last=run.tic[APEX + 12]
    )
    integrate(peak, run.tic, run.rt * 60)
    return run, peak


@pytest.mark.parametrize(
    ("ions", "label"),
    [(METHYL_PALMITATE, "C16:0"), (METHYL_OLEATE, "C18:1"), (METHYL_LINOLEATE, "C18:2")],
)
def test_fames_are_identified(ions, label):
    ident = identify(*run_with((ions, 1e6, APEX)), Config())
    assert ident.fame is not None and ident.fame.label == label, ident.reason
    assert ident.m_plus_rel == pytest.approx(ions[ident.fame.m_plus], rel=0.05)


@pytest.mark.parametrize(
    ("ions", "candidate", "reason"),
    [
        (HEPTADECANE, "C14:1", "M+ 240 C14:1? main ions low: 55:20% 69:0% 83:0% 97:0% 74:0%"),
        (ETHYL_OLEATE, "C19:1", "M+ 310 C19:1? ethyl ester: 88:68% (74:34%) 101:46% (87:26%)"),
        (SILOXANE, None, "no co-eluting ion >= m/z 150"),
        # a compound whose molecular ion lies above the apparent FAME M+
        ({**METHYL_PALMITATE, 288: 6.0}, None, "highest ion m/z 288 is no FAME M+"),
        # M+1 stronger than half of M+: 271 is not the isotope of 270
        ({**METHYL_PALMITATE, 271: 6.0}, None, "highest ion m/z 271 is no FAME M+"),
    ],
)
def test_rejected(ions, candidate, reason):
    ident = identify(*run_with((ions, 1e6, APEX)), Config())
    assert ident.fame is None
    assert (ident.candidate.label if ident.candidate else None) == candidate
    assert ident.reason == reason


def test_coeluting_bleed_is_ignored():
    # siloxane ions above M+ 270 (281, 355, 429, isotopes) must not hide the FAME
    run, peak = run_with((METHYL_PALMITATE, 1e6, APEX), (SILOXANE, 1.5e5, APEX))
    ident = identify(run, peak, Config())
    assert ident.fame is not None and ident.fame.label == "C16:0", ident.reason


def test_coelution():
    # an ion of a compound eluting 6 scans after the apex does not belong to the peak
    other = {57: 100, 300: 20}
    run, peak = run_with((METHYL_PALMITATE, 1e6, APEX), (other, 3e5, APEX + 6))
    assert _coelutes(run, peak, run.column[270], 2)
    assert not _coelutes(run, peak, run.column[300], 2)


def test_fame_table():
    table = fame_table()
    assert len(table) == 15 * 3
    for f in table:
        assert f.m_plus == 14 * f.chain_length + 46 - 2 * f.double_bonds
        assert f.label == f"C{f.chain_length}:{f.double_bonds}"
        assert f.name.endswith("methyl ester")
    assert {f.chain_length for f in table} == set(range(10, 25))


def test_custom_fame_table(tmp_path):
    path = tmp_path / "fames.csv"
    path.write_text("chain_length,double_bonds,label,m_plus\n16,0,C16:0,270\n18,3,C18:3,292\n")
    assert [(f.label, f.name) for f in fame_table(str(path))] == [("C16:0", ""), ("C18:3", "")]
    # the table has no C19:0 standard and a class without rules
    with pytest.raises(ValueError, match="standard C19:0"):
        check_config(Config(fame_table=str(path)))
    with pytest.raises(ValueError, match="no entry for 3 double bond"):
        check_config(Config(fame_table=str(path), standard=""))
    classes = Config().classes | {3: FameClass({79: 50.0}, (-1.0, 1.0))}
    check_config(Config(fame_table=str(path), standard="", classes=classes))
    with pytest.raises(ValueError, match="C99:0"):
        check_config(Config(report_fames=("C16:0", "C99:0")))


def test_default_config_fits_the_table():
    check_config(Config())


def test_ecl_anchors_ignore_a_stray_identification():
    # (RT, chain length, area); "C20:0" at 22.0 min elutes before C18:0 and is dropped
    saturated = [(19.92, 16, 5e7), (22.0, 20, 1e5), (25.17, 18, 2e7), (26.99, 19, 2.3e8)]
    anchors = _anchors(saturated)
    assert anchors == [(19.92, 16), (25.17, 18), (26.99, 19)]
    assert _ecl(25.17, anchors) == pytest.approx(18.0)
    assert _ecl(24.65, anchors) == pytest.approx(17.80, abs=0.01)
    assert _ecl(28.25, anchors) == pytest.approx(19.69, abs=0.01)  # extrapolated
