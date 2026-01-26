# /// script
# requires-python = ">=3.12"
# dependencies = []
# ///

from __future__ import annotations

import argparse
import platform
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Sequence

from run_utils import resolve_base_name, run_root, update_summary


SCRIPT_ORDER: list[tuple[str, str, list[str]]] = [
    ("QC", "qc_analysis.py", []),
    ("Core traits", "query_snps.py", []),
    ("Variant verification", "verify_variants.py", []),
    ("Aging/lifestyle", "life_aging_analysis.py", []),
    ("Hidden risks", "check_extra_snps.py", []),
    ("Expanded panels", "additional_panels.py", []),
    ("Research template", "build_research_findings.py", ["--write-empty"]),
    ("Clinical trials", "search_trials_for_findings.py", []),
    ("Report", "generate_report.py", []),
]


def _run_command(args: Sequence[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        args,
        text=True,
        check=True,
        capture_output=False,
    )


def _safe_version(cmd: Sequence[str]) -> str | None:
    try:
        result = subprocess.run(cmd, text=True, check=True, capture_output=True)
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None
    return result.stdout.strip() or result.stderr.strip() or None


def _collect_manifest() -> dict[str, str | None]:
    return {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "uv_version": _safe_version(["uv", "--version"]),
        "git_commit": _safe_version(["git", "rev-parse", "--short", "HEAD"]),
    }


def _ensure_input_file(base_name: str) -> Path:
    input_file = Path(f"{base_name}.txt")
    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")
    return input_file


def run_pipeline(
    base_name: str,
    *,
    sex: str | None,
    age: int | None,
    skip_trials: bool,
    build_gwas: str | None,
) -> None:
    _ensure_input_file(base_name)
    run_dir = run_root(base_name)
    research_path = run_dir / "research_findings.json"
    update_summary(
        run_dir,
        {
            "run_folder": str(run_dir),
            "run_manifest": _collect_manifest(),
            "reported_sex": sex,
            "reported_age": age,
        },
    )

    if build_gwas:
        _run_command(["uv", "run", "--script", "build_gwas_risk_table.py", build_gwas])

    for label, script, extra_args in SCRIPT_ORDER:
        if skip_trials and script == "search_trials_for_findings.py":
            print("Skipping clinical trials search (--skip-trials)")
            continue
        if script == "build_research_findings.py" and research_path.exists():
            print("research_findings.json already exists; skipping template generation.")
            continue
        print(f"\n==> {label}: {script}")
        _run_command(["uv", "run", "--script", script, base_name, *extra_args])


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the DNA analysis pipeline end-to-end.")
    parser.add_argument("base_name", help="Base filename without extension")
    parser.add_argument("--sex", choices=["female", "male"], help="Reported sex for hormone notes")
    parser.add_argument("--age", type=int, help="Reported age in years (optional)")
    parser.add_argument("--skip-trials", action="store_true", help="Skip clinical trials search")
    parser.add_argument(
        "--build-gwas",
        help="Optional: path to GWAS TSV/CSV to build data/gwas_risk_alleles.json",
    )
    args = parser.parse_args()

    base_name = resolve_base_name(args.base_name)
    try:
        run_pipeline(
            base_name,
            sex=args.sex,
            age=args.age,
            skip_trials=args.skip_trials,
            build_gwas=args.build_gwas,
        )
    except FileNotFoundError as exc:
        print(str(exc))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
