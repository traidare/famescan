"""Command line interface: ``famescan report`` and ``famescan validate``."""

import argparse
import logging
import sys
from pathlib import Path

from . import __version__
from .config import Config
from .fames import fame_table
from .io import find_runs
from .pipeline import analyse
from .report import write_outputs
from .validate import summary, validate, write_csv

log = logging.getLogger("famescan")


def _labels(value: str) -> frozenset[str]:
    labels = frozenset(v.strip() for v in value.split(",") if v.strip())
    unknown = labels - {f.label for f in fame_table()}
    if unknown:
        raise argparse.ArgumentTypeError(f"unknown FAME label(s): {', '.join(sorted(unknown))}")
    return labels


def _report(args: argparse.Namespace) -> int:
    runs = find_runs(args.paths)
    if not runs:
        log.error("no .D runs found")
        return 1
    cfg = Config()
    reports, failed = [], 0
    for path in runs:
        try:
            a = analyse(path, cfg, args.fames)
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
    written = write_outputs(reports, args.output)
    log.info("wrote %d per-run reports, %s and %s", len(written) - 2, *written[-2:])
    return 1 if failed else 0


def _validate(args: argparse.Namespace) -> int:
    results = validate(args.reference, args.data_root, Config(), args.rt_tolerance)
    if not results:
        log.error("no reference report could be matched to a run")
        return 1
    print("\n".join(summary(results, args.rt_tolerance)))
    args.output.mkdir(parents=True, exist_ok=True)
    path = args.output / "validation.csv"
    write_csv(results, path)
    log.info("wrote %s", path)
    return 0


class _Formatter(logging.Formatter):
    def format(self, record: logging.LogRecord) -> str:
        msg = super().format(record)
        return msg if record.levelno == logging.INFO else f"{record.levelname.lower()}: {msg}"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="famescan",
        description="FAME area percent reports from Agilent GC/MSD .D runs.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument(
        "-v", "--verbose", action="store_true", help="also log why FAME candidates were rejected"
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
        "--fames",
        type=_labels,
        metavar="LABELS",
        help="report only these FAMEs, e.g. C14:0,C16:0,C16:1,C18:0,C18:1,C18:2,C19:0 "
        "(default: all FAMEs of the table)",
    )
    p.set_defaults(func=_report)

    p = sub.add_parser(
        "validate", parents=[common], help="compare the pipeline with manual ChemStation reports"
    )
    p.add_argument("--reference", type=Path, required=True, help="ChemStation percent reports")
    p.add_argument("--data-root", type=Path, required=True, help="folder with the .D runs")
    p.add_argument("-o", "--output", type=Path, default=Path("output"), help="output folder")
    p.add_argument(
        "--rt-tolerance", type=float, default=0.03, help="peak matching tolerance in minutes"
    )
    p.set_defaults(func=_validate)

    args = parser.parse_args(argv)
    if not log.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(_Formatter("%(message)s"))
        log.addHandler(handler)
    log.setLevel(logging.DEBUG if args.verbose else logging.INFO)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
