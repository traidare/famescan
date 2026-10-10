from pathlib import Path

from fakes import fake_report, fame

from famescan.config import Config
from famescan.fames import Identification
from famescan.peaks import Peak
from famescan.pipeline import Analysis
from famescan.report import build_report
from famescan.validate import RefPeak, RunValidation, match_run, read_reference, summary


def test_read_reference_csv(tmp_path):
    path = tmp_path / "ref.csv"
    path.write_text("data_file,rt_min,area,fame\nX.D,19.920,54868228,C16:0\nX.D,24.654,,\n")
    assert read_reference(path) == [
        RefPeak("X.D", 19.92, 54868228.0, None, "C16:0"),
        RefPeak("X.D", 24.654),
    ]


def fake_analysis(cfg: Config | None = None) -> Analysis:
    """The FAME peaks C10:0, C10:1 and the C19:0 standard of a fake report (apexes at scans 5,
    15, 25), plus a detected peak at scan 35 that is no FAME."""
    report = fake_report(Path("/data/X.D"), 3)
    peaks = [fp.peak for fp in report.peaks] + [Peak(32, 35, 39)]
    ids = [fp.ident for fp in report.peaks] + [
        Identification(None, None, 0.0, 57, 240, "M+ 240 C14:1? main ions low: 69:18%")
    ]
    return Analysis(report.run, peaks, ids, build_report(report.run, peaks, ids, cfg or Config()))


def test_match_run():
    a = fake_analysis()
    rt = a.run.rt
    refs = [
        RefPeak("X.D", float(rt[5]), area=3.3e6),  # C10:0
        RefPeak("X.D", float(rt[15]), fame="C10:0"),  # famescan: C10:1
        RefPeak("X.D", float(rt[35])),  # detected, but no FAME
        RefPeak("X.D", float(rt[-1]) + 1.0),  # nothing there
    ]
    matches = match_run(refs, a, tolerance=0.005)
    assert [m.status for m in matches] == [
        "found",
        "other FAME",
        "rejected",
        "not detected",
        "extra",
    ]
    assert matches[1].fame_peak.fame == fame("C10:1")
    assert matches[2].reason.startswith("M+ 240 C14:1? main ions low")
    assert matches[4].fame_peak.is_standard

    lines = summary([RunValidation(a, matches)])
    assert lines[0] == "1 runs, 4 reference peaks"
    assert "  recall (found)" in "\n".join(lines)
    assert "  1  M+ 240 C14:1? main ions low" in "\n".join(lines)


def test_rejected_by_the_report_settings():
    a = fake_analysis(Config(report_fames=("C10:0", "C19:0")))
    m, *extra = match_run([RefPeak("X.D", float(a.run.rt[15]))], a, tolerance=0.005)
    assert (m.status, m.reason) == ("rejected", "C10:1 not in report_fames")
    assert [x.status for x in extra] == ["extra", "extra"]
