"""Command line interface: ``inbrain-dementia <command> ...``.

Commands
--------
check-config     load and validate the configuration files
init-config      copy the active configuration to a folder (to edit it)
list-variables   list the configured variables (optionally of one dataset)
validate         check folders and raw tables under a root path
preprocess       run the per-dataset pipelines
integrate        merge the per-dataset tables and build the expanded master
cohort           run a cohort tool described in a YAML file (see configs/cohorts/)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import yaml

from .config import ConfigError, load_config
from .logging_utils import enable_logging

COHORT_TOOLS = {
    "anchor": "create_anchor_cohort",
    "cross_sectional": "create_cross_sectional_cohort",
    "multimodal": "create_multimodal_cohort",
    "longitudinal": "create_longitudinal_cohort",
    "longitudinal_multimodal": "create_longitudinal_multimodal_cohort",
}


def _cfg(args):
    return load_config(args.config_dir)


def cmd_check_config(args) -> int:
    cfg = _cfg(args)
    print(f"Configuration OK: {cfg.config_dir}")
    print(f"  datasets : {', '.join(cfg.datasets)}")
    print(f"  variables: {len(cfg.variables)} ({sum(v.expand for v in cfg.variables.values())} expanded)")
    print(f"  imaging types: {len(cfg.imaging['types'])}")
    print(f"  sha256: {cfg.fingerprint()[:16]}...")
    return 0


def cmd_init_config(args) -> int:
    import shutil

    src = _cfg(args).config_dir
    dst = Path(args.folder)
    if dst.exists() and any(dst.iterdir()):
        raise SystemExit(f"{dst} is not empty")
    shutil.copytree(src, dst, dirs_exist_ok=True)
    print(f"Copied {src} -> {dst}")
    print(f"Use it with --config-dir {dst} (CLI), load_config('{dst}') (Python) "
          f"or the environment variable INBRAIN_CONFIG_DIR.")
    return 0


def cmd_list_variables(args) -> int:
    cfg = _cfg(args)
    rows = []
    for name, v in cfg.variables.items():
        for ds, srcs in v.sources.items():
            if args.dataset and cfg.dataset(args.dataset).name != ds:
                continue
            for s in srcs:
                where = s.table or cfg.datasets[ds].default_table
                what = (s.transform + "(" + ",".join(s.items) + ")" if s.transform
                        else s.instrument + "_<tp>_" + (s.suffix or "") if s.instrument
                        else " | ".join(s.column))
                what = what if len(what) <= 60 else what[:57] + "..."
                rows.append((name, ds, where, what, "yes" if v.expand else ""))
        if not v.sources and v.builtin and not args.dataset:
            rows.append((name, "(builtin)", "", "", "yes" if v.expand else ""))
    w = [max(len(str(r[i])) for r in rows + [("variable", "dataset", "table", "source", "expand")])
         for i in range(5)]
    head = ("variable", "dataset", "table", "source", "expand")
    print("  ".join(h.ljust(w[i]) for i, h in enumerate(head)))
    for r in rows:
        print("  ".join(str(x).ljust(w[i]) for i, x in enumerate(r)))
    return 0


def cmd_validate(args) -> int:
    from .validation import validate_all

    reps = validate_all(args.root, args.datasets, config=_cfg(args))
    for rep in reps.values():
        print(rep.summary())
    if args.save:
        import pandas as pd

        pd.concat([r.to_frame() for r in reps.values()]).to_csv(args.save, index=False)
        print(f"Saved -> {args.save}")
    return 0 if all(r.ok for r in reps.values()) else 1


def cmd_preprocess(args) -> int:
    from .preprocess import preprocess_all

    res = preprocess_all(args.root, args.datasets, config=_cfg(args), allow_missing=args.allow_missing,
                         stop_on_error=args.stop_on_error)
    print(res.to_string(index=False))
    return 0 if (res["status"] == "ok").all() else 1


def cmd_integrate(args) -> int:
    from .integration import integrate

    paths = integrate(args.root, args.datasets, config=_cfg(args), expand=not args.no_expand)
    for k, p in paths.items():
        print(f"{k:<9} {p}")
    return 0


def cmd_cohort(args) -> int:
    from . import cohorts
    from .integration import integrated_paths

    cfg = _cfg(args)
    spec = yaml.safe_load(Path(args.spec).read_text(encoding="utf-8")) or {}
    tool = spec.pop("tool", None)
    if tool not in COHORT_TOOLS:
        raise SystemExit(f"{args.spec}: 'tool' must be one of {sorted(COHORT_TOOLS)}")
    input_path = args.input or spec.pop("input_path", None)
    if input_path is None:
        if not args.root:
            raise SystemExit("give the root path (or --input <expanded master>)")
        input_path = integrated_paths(args.root, cfg)["expanded"]
    spec.pop("input_path", None)
    if args.out_dir:
        spec["out_dir"] = args.out_dir
    fn = getattr(cohorts, COHORT_TOOLS[tool])
    res = fn(input_path, config=cfg, **spec)
    print(f"{tool}: {len(res.table):,} rows, {res.n_subjects:,} subjects")
    for k, p in res.files.items():
        print(f"  {k:<9} {p}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="inbrain-dementia", description=__doc__.split("\n")[0])
    p.add_argument("--config-dir", default=None, help="configuration folder (default: repository configs/)")
    p.add_argument("--log-level", default="INFO", help="DEBUG, INFO, WARNING (default INFO)")
    sub = p.add_subparsers(dest="command", required=True)

    s = sub.add_parser("check-config", help="validate the configuration files")
    s.set_defaults(func=cmd_check_config)

    s = sub.add_parser("init-config", help="copy the configuration to a folder for editing")
    s.add_argument("folder")
    s.set_defaults(func=cmd_init_config)

    s = sub.add_parser("list-variables", help="list configured variables")
    s.add_argument("--dataset", default=None)
    s.set_defaults(func=cmd_list_variables)

    for name, func, h in (("validate", cmd_validate, "check folders and raw tables"),
                          ("preprocess", cmd_preprocess, "run the per-dataset pipelines"),
                          ("integrate", cmd_integrate, "merge + expanded master")):
        s = sub.add_parser(name, help=h)
        s.add_argument("root", help="root folder holding the Tables_<DATASET> folders")
        s.add_argument("--datasets", nargs="*", default=None, help="subset of datasets")
        if name == "validate":
            s.add_argument("--save", default=None, help="write the report as CSV")
        if name == "preprocess":
            s.add_argument("--allow-missing", action="store_true",
                           help="run even when required tables are missing")
            s.add_argument("--stop-on-error", action="store_true")
        if name == "integrate":
            s.add_argument("--no-expand", action="store_true", help="only merge")
        s.set_defaults(func=func)

    s = sub.add_parser("cohort", help="run a cohort tool from a YAML spec (configs/cohorts/*.yaml)")
    s.add_argument("spec", help="cohort YAML file")
    s.add_argument("root", nargs="?", default=None, help="root folder (to locate the expanded master)")
    s.add_argument("--input", default=None, help="expanded master CSV (overrides root)")
    s.add_argument("--out-dir", default=None)
    s.set_defaults(func=cmd_cohort)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    enable_logging(args.log_level)
    try:
        return int(args.func(args) or 0)
    except ConfigError as exc:
        print(f"Configuration error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
