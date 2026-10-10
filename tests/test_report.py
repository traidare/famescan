from datetime import datetime
from pathlib import Path

import numpy as np
import pytest

from famescan import __version__
from famescan.chemstation import parse_reports
from famescan.config import Config
from famescan.fames import Identification, fame_table
from famescan.io import Run
from famescan.peaks import Peak
from famescan.report import (
    _pair,
    build_report,
    format_row,
    render,
    report_names,
    write_outputs,
)

# Rows of ChemStation percent reports of the runs (RTE and hand integrations)
ROWS = [
    (
        (1, 14.146, 797, 800, 805, "rM", 2051839, 5778523, 1.02, 0.434),
        "  1  14.146   797  800  805 rM   2051839   5778523   1.02%   0.434%",
    ),
    (
        (3, 20.011, 1300, 1318, 1328, "rBV", 34325129, 205704764, 36.48, 15.467),
        "  3  20.011  1300 1318 1328 rBV  34325129 205704764  36.48%  15.467%",
    ),
    (
        (5, 24.824, 1723, 1743, 1754, "rVB", 50422279, 563910564, 100.0, 42.4),
        "  5  24.824  1723 1743 1754 rVB  50422279 563910564 100.00%  42.400%",
    ),
    (
        (6, 6.876, 153, 158, 166, "rVB2", 23047489, 33594035, 14.25, 4.569),
        "  6   6.876   153  158  166 rVB2 23047489  33594035  14.25%   4.569%",
    ),
    (
        (10, 11.417, 555, 559, 564, "rVV", 7087990, 15847000, 6.72, 2.155),
        " 10  11.417   555  559  564 rVV  7087990  15847000   6.72%   2.155%",
    ),
]


@pytest.mark.parametrize(("args", "line"), ROWS)
def test_row_layout(args, line):
    assert format_row(*args) == line


@pytest.mark.parametrize(
    ("args", "line"),
    [
        (
            ("Smoothing", "ON", "Filtering", "5"),
            "  Smoothing : ON                             Filtering: 5",
        ),
        (
            ("Sampling", "1", "Min Area", "3 % of largest Peak"),
            "  Sampling  : 1                               Min Area: 3 % of largest Peak",
        ),
        (
            ("Start Thrs", "0.2", "Max Peaks", "100"),
            "  Start Thrs: 0.2                            Max Peaks: 100",
        ),
        (
            ("Stop Thrs", "0", "Peak Location", "TOP     "),
            "  Stop Thrs : 0                          Peak Location: TOP     ",
        ),
    ],
)
def test_parameter_layout(args, line):
    assert _pair(*args) == line


def fake_run(path: Path, n_scans: int = 100) -> Run:
    rt = 5.0 + np.arange(n_scans) * 0.679 / 60
    return Run(
        path=path,
        rt=rt,
        mz=np.arange(40, 551),
        intensity=np.zeros((n_scans, 511)),
        tic=np.zeros(n_scans),
        sample="Sample 983",
        acq_on="13 Aug 2026  14:09",
        method="TEST.M",
        als_vial=4,
    )


def fake_report(path: Path, n_peaks: int, labels=None):
    """A report with ``n_peaks`` FAME peaks (C10:0, C10:1, ...) every 10 scans, the last of them
    the C19:0 standard."""
    run = fake_run(path, 10 * n_peaks + 10)
    table = [f for f in fame_table() if not f.is_standard]
    standard = next(f for f in fame_table() if f.is_standard)
    peaks, ids = [], []
    for k in range(n_peaks):
        fame = standard if k == n_peaks - 1 else table[k]
        apex = 10 * k + 5
        start, end = ("V", "B") if k % 2 else ("B", "V")
        peaks.append(
            Peak(apex - 3, apex, apex + 4, start, end, height=1e6 * (k + 1), area=3.3e6 * (k + 1))
        )
        ids.append(Identification(fame, fame, 5.0, 74, fame.m_plus))
    return build_report(run, peaks, ids, Config(), labels)


def test_build_report_percentages():
    report = fake_report(Path("/data/X.D"), 4)
    areas = [p.area for p in report.peaks]
    assert areas == [3300000, 6600000, 9900000, 13200000]
    assert report.sum_area == sum(areas)
    assert [round(p.pct_max, 2) for p in report.peaks] == [25.0, 50.0, 75.0, 100.0]
    assert sum(p.pct_total for p in report.peaks) == pytest.approx(100.0)
    assert [p.is_standard for p in report.peaks] == [False, False, False, True]


def test_build_report_label_filter_and_min_area():
    report = fake_report(Path("/data/X.D"), 4, labels=frozenset({"C10:0", "C19:0"}))
    assert [p.fame.label for p in report.peaks] == ["C10:0", "C19:0"]
    assert report.peaks[-1].pct_total == pytest.approx(100 * 13.2 / 16.5)

    report = fake_report(Path("/data/X.D"), 4)
    report.peaks[0].peak.area = 1.0  # below 0.1 % of the largest FAME
    rebuilt = build_report(
        report.run,
        [fp.peak for fp in report.peaks],
        [fp.ident for fp in report.peaks],
        Config(),
    )
    assert len(rebuilt.peaks) == 3


def test_render_layout():
    report = fake_report(Path("/data/GCMS 13.08/GSc_X.D"), 7)
    text = render(report, datetime(2026, 10, 5, 8, 55, 44))
    assert "\n" not in text.replace("\r\n", "")
    lines = text.split("\r\n")
    assert lines[0] == " " * 52 + "Area Percent Report"
    assert lines[3] == "  Data File : GSc_X.D" + " " * 45
    assert lines[4] == "  Acq On    : 13 Aug 2026  14:09"
    assert lines[8] == "  ALS Vial  : 4   Sample Multiplier: 1"
    assert "  Signal     : TIC: GSc_X.D\\data.ms" in lines
    for label in ("Filtering:", "Min Area:", "Peak Location:"):
        line = next(x for x in lines if label in x)
        assert line.index(label) + len(label) - 1 == 54
    rows = lines.index(" ---  ----- ----- ---- ---- ---  -------   -------  ------  -------") + 1
    assert lines[rows + 5] == " "  # spacer after every 5th peak
    assert lines[rows + 8 : rows + 10] == [" ", " "]
    assert lines[rows + 10] == " " * 24 + f"Sum of corrected areas:{report.sum_area:12d}"
    assert lines[-4:] == ["", f"famescan {__version__} Mon Oct 05 08:55:44 2026   ", "", ""]


def test_render_parse_roundtrip():
    report = fake_report(Path("/data/GSc_X.D"), 7)
    (parsed,) = parse_reports(render(report))
    assert parsed.data_file == "GSc_X.D"
    assert (parsed.sample, parsed.acq_on, parsed.als_vial) == (
        "Sample 983",
        "13 Aug 2026  14:09",
        "4",
    )
    assert parsed.sum_area == report.sum_area
    assert len(parsed.peaks) == len(report.peaks)
    for n, (got, fp) in enumerate(zip(parsed.peaks, report.peaks, strict=True), 1):
        assert got.number == n
        assert got.rt == pytest.approx(fp.rt, abs=5e-4)
        assert (got.first, got.max, got.last) == (
            fp.peak.first + 1,
            fp.peak.apex + 1,
            fp.peak.last + 1,
        )
        assert (got.type, got.height, got.area) == (fp.peak.type, fp.height, fp.area)
        assert got.pct_max == pytest.approx(fp.pct_max, abs=0.005)
        assert got.pct_total == pytest.approx(fp.pct_total, abs=5e-4)


def test_report_names_disambiguate():
    a = fake_report(Path("/data/GCMS 13.08/GSc_X.D"), 2)
    b = fake_report(Path("/data/GCMS 14.08/GSc_X.D"), 2)
    c = fake_report(Path("/data/GCMS 14.08/GSc_Y.D"), 2)
    assert report_names([a, b, c]) == ["GCMS 13.08_GSc_X.txt", "GCMS 14.08_GSc_X.txt", "GSc_Y.txt"]


def test_write_outputs(tmp_path):
    reports = [fake_report(Path("/data/A.D"), 3), fake_report(Path("/data/B.D"), 6)]
    written = write_outputs(reports, tmp_path / "out")
    assert [p.name for p in written] == ["A.txt", "B.txt", "percent_reports.txt", "fame_peaks.csv"]
    combined = parse_reports(written[2].read_bytes().decode("latin-1"))
    assert [r.data_file for r in combined] == ["A.D", "B.D"]
    csv_lines = written[3].read_text(encoding="utf-8").splitlines()
    assert csv_lines[0].startswith("data_file,sample,acq_on,peak,rt_min")
    assert len(csv_lines) == 1 + 3 + 6
