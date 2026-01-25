from __future__ import annotations

from pathlib import Path
from typing import Iterable

import polars as pl

REFERENCE_PATH = Path(__file__).resolve().parent / "data" / "snp_reference.csv"


def load_reference() -> pl.DataFrame:
    if not REFERENCE_PATH.exists():
        raise FileNotFoundError(f"Missing SNP reference file: {REFERENCE_PATH}")
    return pl.read_csv(REFERENCE_PATH)


def panel_records(reference: pl.DataFrame, panel_name: str) -> list[dict[str, str]]:
    filtered = reference.filter(pl.col("panel") == panel_name)
    return filtered.to_dicts()


def panels_to_records(
    reference: pl.DataFrame, panel_names: Iterable[str]
) -> dict[str, list[dict[str, str]]]:
    return {panel: panel_records(reference, panel) for panel in panel_names}
