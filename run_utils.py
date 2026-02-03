from __future__ import annotations

import json
from datetime import date
from pathlib import Path
from typing import Any, TypedDict

_VALID_BASES = {"A", "C", "G", "T"}


class GenotypeCall(TypedDict):
    kind: str
    genotype: str | None
    raw: str | None


def resolve_base_name(arg: str | None, default: str = "ancestrydna-test-file") -> str:
    if not arg:
        return default
    base = arg
    for suffix in (".normalized.parquet", ".parquet", ".txt"):
        if base.endswith(suffix):
            base = base[: -len(suffix)]
    return base


def run_root(base_name: str) -> Path:
    run_date = date.today().strftime("%Y%m%d")
    root = Path("runs") / run_date / base_name
    root.mkdir(parents=True, exist_ok=True)
    return root


def resolve_parquet_path(base_name: str) -> Path:
    root = run_root(base_name)
    return root / f"{base_name}.normalized.parquet"


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def normalize_genotype(allele1: str | None, allele2: str | None) -> str | None:
    """Normalize a diploid genotype to sorted A/C/G/T pairs.

    Returns None for missing or non-ACGT calls (including indels), so callers
    can label the marker as not assessed rather than misinterpreting.
    """
    a1 = (allele1 or "").strip().upper()
    a2 = (allele2 or "").strip().upper()
    if not a1 or not a2:
        return None
    if a1 in {"0", "--"} or a2 in {"0", "--"}:
        return None
    if a1 not in _VALID_BASES or a2 not in _VALID_BASES:
        return None
    return "".join(sorted([a1, a2]))


def classify_genotype(allele1: str | None, allele2: str | None) -> GenotypeCall:
    """Classify genotype as ACGT SNP, missing, or non-SNP call (indel/repeat)."""
    a1 = (allele1 or "").strip().upper()
    a2 = (allele2 or "").strip().upper()
    if not a1 or not a2 or a1 in {"0", "--"} or a2 in {"0", "--"}:
        return {"kind": "missing", "genotype": None, "raw": None}
    if a1 in _VALID_BASES and a2 in _VALID_BASES:
        return {
            "kind": "acgt",
            "genotype": "".join(sorted([a1, a2])),
            "raw": None,
        }
    return {"kind": "non_snp", "genotype": None, "raw": f"{a1}/{a2}"}


def load_summary(root: Path) -> dict[str, Any]:
    summary_path = root / "summary.json"
    if not summary_path.exists():
        return {}
    return json.loads(summary_path.read_text(encoding="utf-8"))


def update_summary(root: Path, updates: dict[str, Any]) -> None:
    summary = load_summary(root)
    summary.update(updates)
    write_json(root / "summary.json", summary)
