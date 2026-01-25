from __future__ import annotations

import json
from datetime import date
from pathlib import Path
from typing import Any


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


def load_summary(root: Path) -> dict[str, Any]:
    summary_path = root / "summary.json"
    if not summary_path.exists():
        return {}
    return json.loads(summary_path.read_text(encoding="utf-8"))


def update_summary(root: Path, updates: dict[str, Any]) -> None:
    summary = load_summary(root)
    summary.update(updates)
    write_json(root / "summary.json", summary)
