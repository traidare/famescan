"""Reading Agilent ChemStation / MassHunter GC-MS runs (``*.D/data.ms``)."""

import struct
from collections.abc import Iterable
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np
import rainbow


@dataclass
class Run:
    path: Path
    rt: np.ndarray
    """Retention time per scan in minutes. Scan numbers in reports are 1-based indices into this."""
    mz: np.ndarray
    """Unit m/z of the intensity columns."""
    intensity: np.ndarray
    """Abundance, scans x m/z."""
    tic: np.ndarray
    sample: str
    acq_on: str
    """Acquisition time formatted like ChemStation reports, e.g. ``10 Aug 2026  13:02``."""
    method: str
    als_vial: int

    @property
    def name(self) -> str:
        return self.path.name


def _pascal_string(raw: bytes, offset: int) -> str:
    return raw[offset + 1 : offset + 1 + raw[offset]].decode("latin-1").strip()


def _format_acq_on(header_date: str) -> str:
    # data.ms stores e.g. "10 Aug 26  01:02 pm"
    try:
        dt = datetime.strptime(" ".join(header_date.split()), "%d %b %y %I:%M %p")
    except ValueError:
        return header_date
    return dt.strftime("%d %b %Y  %H:%M")


def _data_ms(path: Path) -> Path:
    for f in path.iterdir():
        if f.name.lower() == "data.ms":
            return f
    raise FileNotFoundError(f"no data.ms in {path}")


def load_run(path: str | Path) -> Run:
    path = Path(path)
    ms_path = _data_ms(path)
    ms = rainbow.read(str(path)).get_file(ms_path.name)
    intensity = np.asarray(ms.data, dtype=np.float64)

    # Header fields rainbow does not expose (ChemStation "GC / MS Data File" layout).
    with ms_path.open("rb") as f:
        raw = f.read(0x110)
    return Run(
        path=path,
        rt=np.asarray(ms.xlabels, dtype=np.float64),
        mz=np.rint(np.asarray(ms.ylabels, dtype=np.float64)).astype(int),
        intensity=intensity,
        tic=intensity.sum(axis=1),
        sample=_pascal_string(raw, 0x18),
        acq_on=_format_acq_on(_pascal_string(raw, 0xB2)),
        method=_pascal_string(raw, 0xE4),
        als_vial=struct.unpack_from(">H", raw, 0xFE)[0],
    )


def find_runs(paths: Iterable[str | Path]) -> list[Path]:
    """Expand ``.D`` folders and folders containing them (recursively), sorted and deduplicated."""
    runs: set[Path] = set()
    for p in map(Path, paths):
        if p.suffix.upper() == ".D" and p.is_dir():
            runs.add(p)
        elif p.is_dir():
            runs.update(d for d in p.rglob("*") if d.is_dir() and d.suffix.upper() == ".D")
        else:
            raise FileNotFoundError(f"not a .D folder or directory: {p}")
    return sorted(runs)
