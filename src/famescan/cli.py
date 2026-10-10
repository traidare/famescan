"""Command line interface: ``famescan report``, ``famescan validate`` and ``famescan config``."""

import argparse
import logging
import sys
from pathlib import Path

from . import __version__
from .config import Config, load_config, to_toml
from .fames import check_config
from .io import find_runs
from .pipeline import analyse
from .report import write_outputs
from .validate import summary, validate, write_csv

log = logging.getLogger("famescan")

LOCAL_CONFIG = Path("famescan.toml")
"""The settings file used without ``--config``, if the current directory has one."""


def _labels(value: str) -> tuple[str, ...]:
    return tuple(v.strip() for v in value.split(",") if v.strip())


def _config(args: argparse.Namespace) -> Config:
    """The settings: defaults < settings file < command line options."""
    path = args.config
    if path is None and LOCAL_CONFIG.is_file():
        path = LOCAL_CONFIG
    if path is not None:
        log.info("settings: %s", path)
    cfg = load_config(
        path,
        fame_table=args.fame_table and str(args.fame_table),
        report_fames=args.report_fames,
        standard=args.standard,
        standard_in_pct=args.standard_in_pct,
    )
    check_config(cfg)
    return cfg


def _report(args: argparse.Namespace) -> int:
    runs = find_runs(args.paths)
    if not runs:
        log.error("no .D runs found")
        return 1
    cfg = _config(args)
    reports, failed = [], 0
    for path in runs:
        try:
            a = analyse(path, cfg)
        except Exception as e:  # a broken run must not stop the batch
            log.error("%s: %s", path, e)
            failed += 1
            continue
        reports.append(a.report)
        for p, i in zip(a.peaks, a.ids, strict=True):
            if i.candidate and not i.fame:
                log.debug("%s %.3f min: %s", a.run.name, a.run.rt[p.apex], i.reason)
        labels = " ".join(fp.fame.label for fp in a.report.peaks)
        log.info("%s: %d FAME peaks  %s", a.run.name, len(a.report.peaks), labels)
        if not a.report.peaks:
            log.warning("%s: no FAME peaks found", a.run.name)
    written = write_outputs(reports, args.output, text_reports=args.text_reports)
    log.info("wrote %d files to %s", len(written), args.output)
    return 1 if failed else 0


def _validate(args: argparse.Namespace) -> int:
    results = validate(args.reference, args.data_root, _config(args), args.rt_tolerance)
    if not results:
        log.error("no reference run was found below %s", args.data_root)
        return 1
    print("\n".join(summary(results)))
    args.output.mkdir(parents=True, exist_ok=True)
    path = args.output / "validation.csv"
    write_csv(results, path)
    log.info("wrote %s", path)
    return 0


def _print_config(args: argparse.Namespace) -> int:
    print(to_toml(_config(args)), end="")
    return 0


class _Formatter(logging.Formatter):
    def format(self, record: logging.LogRecord) -> str:
        msg = super().format(record)
        return msg if record.levelno == logging.INFO else f"{record.levelname.lower()}: {msg}"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="famescan",
        description="Find, identify and integrate fatty acid methyl ester (FAME) peaks in GC-MS "
        "runs (Agilent .D folders) and write FAME profiles and percent reports.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument(
        "-v", "--verbose", action="store_true", help="also log why FAME candidates were rejected"
    )
    settings = common.add_argument_group("settings (override the settings file)")
    settings.add_argument(
        "--config",
        type=Path,
        metavar="TOML",
        help="settings file (default: famescan.toml in the current directory, if there is one)",
    )
    settings.add_argument(
        "--fame-table", type=Path, metavar="CSV", help="FAMEs to look for (default: built-in)"
    )
    settings.add_argument(
        "--report-fames",
        type=_labels,
        metavar="LABELS",
        help="report only these FAMEs, e.g. C14:0,C16:0,C16:1,C18:0,C18:1,C18:2,C19:0 "
        "(default: all FAMEs of the table)",
    )
    standard = settings.add_mutually_exclusive_group()
    standard.add_argument(
        "--standard", metavar="LABEL", help="internal standard FAME (default: C19:0)"
    )
    standard.add_argument(
        "--no-standard", dest="standard", action="store_const", const="", help="no standard"
    )
    settings.add_argument(
        "--standard-in-pct",
        action=argparse.BooleanOptionalAction,
        help="count the standard in %% max and %% total (default: yes)",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser(
        "report",
        parents=[common],
        help="detect, identify and integrate FAME peaks, write reports",
    )
    p.add_argument("paths", nargs="+", type=Path, help=".D folders or folders containing them")
    p.add_argument("-o", "--output", type=Path, default=Path("output"), help="output folder")
    p.add_argument(
        "--text-reports",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="also write ChemStation-style text reports (default: yes)",
    )
    p.set_defaults(func=_report)

    p = sub.add_parser(
        "validate", parents=[common], help="compare the FAME peaks with reference peaks"
    )
    p.add_argument(
        "--reference",
        type=Path,
        required=True,
        help="reference peaks: a CSV file or ChemStation percent reports",
    )
    p.add_argument("--data-root", type=Path, required=True, help="folder with the .D runs")
    p.add_argument("-o", "--output", type=Path, default=Path("output"), help="output folder")
    p.add_argument(
        "--rt-tolerance", type=float, default=0.03, help="peak matching tolerance in minutes"
    )
    p.set_defaults(func=_validate)

    p = sub.add_parser("config", parents=[common], help="print the settings as TOML")
    p.set_defaults(func=_print_config)

    args = parser.parse_args(argv)
    if not log.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(_Formatter("%(message)s"))
        log.addHandler(handler)
    log.setLevel(logging.DEBUG if args.verbose else logging.INFO)
    try:
        return args.func(args)
    except (OSError, ValueError) as e:
        log.error("%s", e)
        return 2


if __name__ == "__main__":
    sys.exit(main())
