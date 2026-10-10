"""The whole pipeline for one run: load, detect and integrate peaks, identify FAMEs, report."""

from dataclasses import dataclass
from pathlib import Path

from .config import Config
from .fames import Identification, identify_run
from .io import Run, load_run
from .peaks import Peak, detect_peaks
from .report import FameReport, build_report


@dataclass
class Analysis:
    run: Run
    peaks: list[Peak]
    """All detected TIC peaks."""
    ids: list[Identification]
    """The identification of each peak."""
    report: FameReport


def analyse(path: str | Path, cfg: Config | None = None) -> Analysis:
    cfg = cfg or Config()
    run = load_run(path)
    peaks = detect_peaks(run.tic, run.rt, cfg)
    ids = identify_run(run, peaks, cfg)
    return Analysis(run, peaks, ids, build_report(run, peaks, ids, cfg))
