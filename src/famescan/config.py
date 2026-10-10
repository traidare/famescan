"""The settings, and loading them from a TOML file.

The settings are what depends on the instrument, the column and the samples. The defaults are for
EI full scan runs (Agilent 5977C) at 0.679 s per scan, with C19:0 as internal standard, on a column
on which unsaturated FAMEs elute shortly before the saturated FAME of the same chain length. The
fixed parts of the rules are constants in ``peaks.py`` and ``fames.py``.

Durations are in seconds and converted to whole scans per run (``scans``), so the settings carry
over to other scan rates.
"""

import tomllib
from dataclasses import asdict, dataclass, field, fields, replace
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class FameClass:
    """The rules for the FAMEs with one number of double bonds."""

    main_ions: dict[int, float]
    """Fragment ions that must be present: m/z -> minimum % of the base peak."""
    ecl_window: tuple[float, float]
    """Allowed equivalent chain length minus chain length."""


@dataclass(frozen=True)
class Config:
    # --- peak detection (TIC); defaults like ChemStation's "Filtering: 5", "Peak separation: 5"
    smoothing_window_s: float = 3.4
    """Savitzky-Golay window, used for detection and boundaries only. Below 3 scans: no
    smoothing."""
    min_width_s: float = 0.88
    """Minimum width at half prominence (smoothed TIC). Narrower maxima are spikes, such as the
    single scan a saturated peak can end with."""
    min_separation_s: float = 3.4
    """Closer apexes are merged."""
    noise_factor: float = 5.0
    """Minimum peak prominence as a multiple of the TIC noise."""

    # --- FAME identification ---
    fame_table: str = ""
    """CSV file with the FAMEs to look for (columns ``chain_length, double_bonds, label, name,
    m_plus``); empty: the built-in table (C10-C24 methyl esters with 0-2 double bonds)."""
    min_ion_pct: float = 2.0
    """Ions from this % of the base peak on are significant; M+ must be one of them."""
    classes: dict[int, FameClass] = field(
        default_factory=lambda: {
            0: FameClass({74: 50.0, 87: 20.0}, (-1.0, 0.5)),
            # 74, the methyl ester ion, separates mono-unsaturated FAMEs from alkenes
            1: FameClass({55: 25.0, 69: 25.0, 83: 25.0, 97: 25.0, 74: 15.0}, (-1.0, 0.5)),
            2: FameClass({67: 40.0, 81: 40.0, 95: 40.0, 55: 15.0}, (-1.0, 1.0)),
        }
    )
    """The rules per number of double bonds; every number in the FAME table needs one. The ECL
    windows suit columns on which unsaturated FAMEs elute shortly before the saturated ones; on
    polar (cyanopropyl) columns they elute after them."""

    # --- report ---
    report_fames: tuple[str, ...] = ()
    """Report only these FAMEs, e.g. ``("C16:0", "C18:1")``; empty: all FAMEs of the table."""
    min_area_pct: float = 0.1
    """FAME peaks below this % of the largest FAME peak are dropped."""
    standard: str = "C19:0"
    """Label of the internal standard; its largest peak in a run is marked as the standard.
    Empty: no internal standard."""
    standard_in_pct: bool = True
    """Whether the standard counts in % max and % total. If not, the percentages are relative to
    the other FAME peaks (the standard's own percentages too)."""


def scans(seconds: float, interval: float) -> int:
    """A duration in whole scans of ``interval`` seconds."""
    return max(0, round(seconds / interval))


def _classes(value: dict, default: dict[int, FameClass]) -> dict[int, FameClass]:
    """``[classes.N]`` tables; a class of the defaults keeps the rules that are not given."""
    classes = dict(default)
    for key, rules in value.items():
        u = int(key)
        unknown = set(rules) - {"main_ions", "ecl_window"}
        if unknown:
            raise ValueError(f"classes.{u}: unknown setting(s) {', '.join(sorted(unknown))}")
        given = {}
        if "main_ions" in rules:
            given["main_ions"] = {int(m): float(x) for m, x in rules["main_ions"].items()}
        if "ecl_window" in rules:
            lo, hi = rules["ecl_window"]
            given["ecl_window"] = (float(lo), float(hi))
        if u in classes:
            classes[u] = replace(classes[u], **given)
        elif len(given) == 2:
            classes[u] = FameClass(**given)
        else:
            raise ValueError(f"classes.{u} needs main_ions and ecl_window")
    return classes


def _convert(name: str, value: Any, default: Any) -> Any:
    """A TOML or command line value in the type of the field's default."""
    if name == "classes":
        return _classes(value, default)
    if isinstance(default, tuple):
        if not isinstance(value, list | tuple):
            raise ValueError(f"{name} must be a list")
        return tuple(str(v) for v in value)
    if isinstance(default, bool):
        if not isinstance(value, bool):
            raise ValueError(f"{name} must be true or false")
        return value
    return type(default)(value)


def load_config(path: Path | None = None, **overrides: Any) -> Config:
    """The defaults, updated from a TOML file (keys = field names) and then from ``overrides``
    (``None`` values are ignored). A relative ``fame_table`` in the file is relative to it."""
    values: dict[str, Any] = {}
    if path is not None:
        with path.open("rb") as f:
            values = tomllib.load(f)
        if values.get("fame_table") and not Path(values["fame_table"]).is_absolute():
            values["fame_table"] = str(path.parent / values["fame_table"])
    values |= {k: v for k, v in overrides.items() if v is not None}
    defaults = Config()
    unknown = sorted(set(values) - {f.name for f in fields(Config)})
    if unknown:
        raise ValueError(f"unknown setting(s): {', '.join(unknown)}")
    return replace(defaults, **{k: _convert(k, v, getattr(defaults, k)) for k, v in values.items()})


def to_toml(cfg: Config) -> str:
    """The settings as TOML, readable by ``load_config``."""

    def value(v: Any) -> str:
        if isinstance(v, bool):
            return "true" if v else "false"
        if isinstance(v, str):
            return '"' + v.replace("\\", "\\\\").replace('"', '\\"') + '"'
        if isinstance(v, tuple | list):
            return "[" + ", ".join(value(x) for x in v) + "]"
        if isinstance(v, dict):
            return "{ " + ", ".join(f"{k} = {value(x)}" for k, x in v.items()) + " }"
        return repr(v)

    settings = asdict(cfg)
    classes = settings.pop("classes")
    lines = [f"{k} = {value(v)}" for k, v in settings.items()]
    for u, rules in classes.items():
        lines += ["", f"[classes.{u}]", *(f"{k} = {value(v)}" for k, v in rules.items())]
    return "\n".join(lines) + "\n"
