from datetime import datetime
from pathlib import Path

import pytest
from fakes import fake_report

from famescan import __version__
from famescan.chemstation import parse_reports
from famescan.config import Config
from famescan.report import (
    _pair,
    build_report,
    format_row,
    profile_rows,
    render,
    report_names,
    write_outputs,
)

# Rows of ChemStation percent reports, with different number widths and peak types
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


def test_build_report_percentages():
    report = fake_report(Path("/data/X.D"), 4)
    areas = [p.area for p in report.peaks]
    assert areas == [3300000, 6600000, 9900000, 13200000]
    assert report.sum_area == sum(areas)
    assert [round(p.pct_max, 2) for p in report.peaks] == [25.0, 50.0, 75.0, 100.0]
    assert sum(p.pct_total for p in report.peaks) == pytest.approx(100.0)
    assert [p.is_standard for p in report.peaks] == [False, False, False, True]


def test_build_report_label_filter_and_min_area():
    report = fake_report(Path("/data/X.D"), 4, Config(report_fames=("C10:0", "C19:0")))
    assert [p.fame.label for p in report.peaks] == ["C10:0", "C19:0"]
    assert report.peaks[-1].pct_total == pytest.approx(100 * 13.2 / 16.5)
    assert sorted(report.dropped.values()) == [
        "C10:1 not in report_fames",
        "C10:2 not in report_fames",
    ]

    report = fake_report(Path("/data/X.D"), 4)
    peaks = [fp.peak for fp in report.peaks]
    peaks[0].area = 1.0  # below 0.1 % of the largest FAME
    rebuilt = build_report(report.run, peaks, [fp.ident for fp in report.peaks], Config())
    assert len(rebuilt.peaks) == 3
    assert rebuilt.dropped == {peaks[0]: "C10:0 below min area"}


def test_standard_options():
    # areas 3.3, 6.6, 9.9 and the standard 13.2 (x 1e6)
    report = fake_report(Path("/data/X.D"), 4, Config(standard_in_pct=False))
    assert [round(p.pct_max, 2) for p in report.peaks] == [33.33, 66.67, 100.0, 133.33]
    assert [round(p.pct_total, 2) for p in report.peaks] == [16.67, 33.33, 50.0, 66.67]
    assert report.peaks[-1].is_standard

    report = fake_report(Path("/data/X.D"), 4, Config(standard=""))
    assert not any(p.is_standard for p in report.peaks)
    assert sum(p.pct_total for p in report.peaks) == pytest.approx(100.0)

    # C10:1 as standard: it counts in the percentages, C19:0 is an ordinary FAME
    report = fake_report(Path("/data/X.D"), 4, Config(standard="C10:1"))
    assert [p.is_standard for p in report.peaks] == [False, True, False, False]


def test_profile_rows():
    a, b = fake_report(Path("/data/A.D"), 3), fake_report(Path("/data/B.D"), 2)
    b.peaks[0].ident.fame = a.peaks[1].fame  # B: C10:1 + standard, A: C10:0 + C10:1 + standard
    columns, rows = profile_rows([a, b])
    assert columns == ["data_file", "sample", "acq_on", "C10:0", "C10:1", "C19:0 (standard)"]
    assert rows[0] | {"sample": "", "acq_on": ""} == {
        "data_file": "A.D",
        "sample": "",
        "acq_on": "",
        "C10:0": "16.667",
        "C10:1": "33.333",
        "C19:0 (standard)": "50.000",
    }
    assert "C10:0" not in rows[1]


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
        assert (got.first, got.max, got.last) == fp.peak.scan_numbers
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
    assert [p.name for p in written] == [
        "fame_peaks.csv",
        "fame_profiles.csv",
        "A.txt",
        "B.txt",
        "percent_reports.txt",
    ]
    combined = parse_reports(written[-1].read_bytes().decode("latin-1"))
    assert [r.data_file for r in combined] == ["A.D", "B.D"]
    csv_lines = written[0].read_text(encoding="utf-8").splitlines()
    assert csv_lines[0].startswith("data_file,sample,acq_on,peak,rt_min")
    assert len(csv_lines) == 1 + 3 + 6
    assert len(written[1].read_text(encoding="utf-8").splitlines()) == 1 + 2

    written = write_outputs(reports, tmp_path / "csv", text_reports=False)
    assert sorted(p.name for p in (tmp_path / "csv").iterdir()) == [
        "fame_peaks.csv",
        "fame_profiles.csv",
    ]
