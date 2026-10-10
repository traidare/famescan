"""Parser for ChemStation "Area Percent Report" text exports (one or many reports per file)."""

import re
from dataclasses import dataclass, field
from pathlib import Path

_ROW = re.compile(
    r"^\s*(\d+)\s+(\d+\.\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)\s+(\d+)\s+"
    r"([\d.]+)%\s+([\d.]+)%\s*$"
)
_FIELD = re.compile(r"^\s{2}(Data Path|Data File|Acq On|Sample|ALS Vial)\s*:\s?(.*?)\s*$")
_SUM = re.compile(r"Sum of corrected areas:\s*(\d+)")


@dataclass
class ReportPeak:
    number: int
    rt: float
    first: int
    max: int
    last: int
    type: str
    height: int
    area: int
    pct_max: float
    pct_total: float


@dataclass
class Report:
    data_path: str = ""
    data_file: str = ""
    acq_on: str = ""
    sample: str = ""
    als_vial: str = ""
    peaks: list[ReportPeak] = field(default_factory=list)
    sum_area: int | None = None


def parse_reports(text: str) -> list[Report]:
    reports: list[Report] = []
    for line in text.replace("\r", "").split("\n"):
        if m := _FIELD.match(line):
            key, value = m.groups()
            if key == "Data Path":  # each report starts here; title lines are not reliable
                reports.append(Report())
            if not reports:
                continue
            r = reports[-1]
            if key == "ALS Vial":
                value = value.split()[0] if value else ""
            setattr(r, key.lower().replace(" ", "_"), value)
        elif reports and (m := _ROW.match(line)):
            g = m.groups()
            reports[-1].peaks.append(
                ReportPeak(
                    int(g[0]),
                    float(g[1]),
                    int(g[2]),
                    int(g[3]),
                    int(g[4]),
                    g[5],
                    int(g[6]),
                    int(g[7]),
                    float(g[8]),
                    float(g[9]),
                )
            )
        elif reports and (m := _SUM.search(line)):
            reports[-1].sum_area = int(m.group(1))
    return reports


def read_reports(path: str | Path) -> list[Report]:
    return parse_reports(Path(path).read_text(encoding="latin-1"))
