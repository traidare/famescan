"""Synthetic runs and reports for the tests."""

from pathlib import Path

import numpy as np

from famescan.config import Config
from famescan.fames import Fame, Identification, fame_table
from famescan.io import Run
from famescan.peaks import Peak
from famescan.report import FameReport, build_report

DT = 0.679  # s per scan
MZ = np.arange(40, 551)


def retention_times(n_scans: int) -> np.ndarray:
    """RT in minutes of ``n_scans`` scans from 5 min on."""
    return 5.0 + np.arange(n_scans) * DT / 60


def make_run(intensity: np.ndarray, path: Path = Path("/data/X.D")) -> Run:
    """A run with the spectra ``intensity`` (scans x ``MZ``)."""
    rt, tic = retention_times(len(intensity)), intensity.sum(axis=1)
    return Run(path, rt, MZ, intensity, tic, "Sample 983", "13 Aug 2026  14:09", "TEST.M", 4)


def fame(label: str) -> Fame:
    return next(f for f in fame_table() if f.label == label)


def fake_report(path: Path, n_peaks: int, cfg: Config | None = None) -> FameReport:
    """A report with ``n_peaks`` FAME peaks (C10:0, C10:1, ...) every 10 scans, the last of them
    the C19:0 standard."""
    run = make_run(np.zeros((10 * n_peaks + 10, len(MZ))), path)
    others = [f for f in fame_table() if f.label != "C19:0"]
    peaks, ids = [], []
    for k in range(n_peaks):
        f = fame("C19:0") if k == n_peaks - 1 else others[k]
        apex = 10 * k + 5
        start, end = ("V", "B") if k % 2 else ("B", "V")
        peaks.append(
            Peak(apex - 3, apex, apex + 4, start, end, height=1e6 * (k + 1), area=3.3e6 * (k + 1))
        )
        ids.append(Identification(f, f, 5.0, 74, f.m_plus))
    return build_report(run, peaks, ids, cfg or Config())
