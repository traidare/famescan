"""Checks against a real run; only if ``FAMESCAN_DATA_ROOT`` names the folder with the runs."""

import os
from pathlib import Path

import pytest

from famescan.chemstation import parse_reports
from famescan.pipeline import analyse
from famescan.report import render

DATA_ROOT = os.environ.get("FAMESCAN_DATA_ROOT")
RUN = Path(DATA_ROOT or ".") / "GCMS 13.08" / "GSc_EBaS1004_20260813.D"

pytestmark = pytest.mark.skipif(
    not DATA_ROOT or not RUN.is_dir(), reason="FAMESCAN_DATA_ROOT not set or run not found"
)

# The manual "ONLY FAME PEAKS" report of the run: RT, FAME, area, and the
# allowed area deviation. Peak 1 was integrated by hand (rM); the small C18:1 isomer on the tail
# of the large one gets a slightly different boundary.
EXPECTED = [
    (19.343, "C16:1", 39449827, 0.15),
    (19.920, "C16:0", 54868228, 0.05),
    (24.461, "C18:2", 30226151, 0.05),
    (24.654, "C18:1", 108953040, 0.05),
    (24.733, "C18:1", 7134360, 0.10),
    (25.174, "C18:0", 19193509, 0.05),
    (26.986, "C19:0", 235807479, 0.05),
]


@pytest.fixture(scope="module")
def analysis():
    return analyse(RUN)


def test_header(analysis):
    run = analysis.run
    assert (run.sample, run.acq_on, run.als_vial) == ("Sample 983", "13 Aug 2026  14:09", 4)


@pytest.mark.parametrize(("rt", "label", "area", "tolerance"), EXPECTED)
def test_reference_fame_peaks(analysis, rt, label, area, tolerance):
    fp = min(analysis.report.peaks, key=lambda fp: abs(fp.rt - rt))
    assert abs(fp.rt - rt) <= 0.02
    assert fp.fame.label == label
    assert fp.area == pytest.approx(area, rel=tolerance)


def test_standard(analysis):
    (standard,) = [fp for fp in analysis.report.peaks if fp.is_standard]
    assert standard.fame.label == "C19:0"
    assert abs(standard.rt - 26.986) <= 0.02


def test_report_roundtrip(analysis):
    (parsed,) = parse_reports(render(analysis.report))
    assert parsed.data_file == RUN.name
    assert [p.area for p in parsed.peaks] == [fp.area for fp in analysis.report.peaks]
    assert parsed.sum_area == analysis.report.sum_area
