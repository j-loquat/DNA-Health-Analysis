# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "polars",
# ]
# ///

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Iterable

import polars as pl

from snp_reference import load_reference


def _parse_risk_allele(value: str | None) -> str | None:
    if not value:
        return None
    token = value.strip().upper()
    if "-" in token:
        token = token.split("-")[-1]
    if token in {"A", "C", "G", "T"}:
        return token
    return None


def _split_snps(value: str | None) -> list[str]:
    if not value:
        return []
    items = re.split(r"[;,\s]+", value.strip())
    return [item for item in items if item.startswith("rs")]


def _load_input(path: Path) -> pl.DataFrame:
    if path.suffix.lower() in {".tsv", ".tab"}:
        return pl.read_csv(
            path,
            separator="\t",
            infer_schema_length=0,
            quote_char=None,
            ignore_errors=True,
        )
    return pl.read_csv(
        path,
        infer_schema_length=0,
        quote_char=None,
        ignore_errors=True,
    )


def _build_mapping(df: pl.DataFrame) -> dict[str, str | None]:
    columns = {col.lower(): col for col in df.columns}
    if "rsid" in columns and "risk_allele" in columns:
        rsid_col = columns["rsid"]
        risk_col = columns["risk_allele"]
        mapping: dict[str, str | None] = {}
        for row in df.select([rsid_col, risk_col]).iter_rows(named=True):
            rsid = str(row[rsid_col]).strip()
            risk = _parse_risk_allele(str(row[risk_col]) if row[risk_col] is not None else None)
            if rsid:
                mapping[rsid] = risk
        return mapping

    # GWAS Catalog TSV style
    snps_col = columns.get("snps")
    strongest_col = columns.get("strongest snp-risk allele")
    if snps_col and strongest_col:
        mapping: dict[str, str | None] = {}
        for row in df.select([snps_col, strongest_col]).iter_rows(named=True):
            rsids = _split_snps(str(row[snps_col]) if row[snps_col] is not None else None)
            risk = _parse_risk_allele(str(row[strongest_col]) if row[strongest_col] is not None else None)
            for rsid in rsids:
                mapping[rsid] = risk
        return mapping

    raise ValueError(
        "Input must contain either columns 'rsid' and 'risk_allele' or GWAS Catalog columns "
        "'SNPS' and 'STRONGEST SNP-RISK ALLELE'."
    )


def _filter_to_reference(mapping: dict[str, str | None], rsids: Iterable[str]) -> dict[str, str | None]:
    allowed = set(rsids)
    return {rsid: allele for rsid, allele in mapping.items() if rsid in allowed}


def main() -> int:
    parser = argparse.ArgumentParser(description="Build local GWAS risk allele table.")
    parser.add_argument("input_path", help="TSV/CSV file with risk allele data")
    parser.add_argument(
        "--output",
        default="data/gwas_risk_alleles.json",
        help="Output JSON path",
    )
    parser.add_argument(
        "--mode",
        choices=["overwrite", "merge"],
        default="merge",
        help="Overwrite or merge with existing JSON",
    )
    args = parser.parse_args()

    input_path = Path(args.input_path)
    output_path = Path(args.output)

    reference = load_reference()
    rsids = reference.select("rsid").to_series().to_list()

    df = _load_input(input_path)
    mapping = _build_mapping(df)
    mapping = _filter_to_reference(mapping, rsids)

    if args.mode == "merge" and output_path.exists():
        existing = json.loads(output_path.read_text(encoding="utf-8"))
    else:
        existing = {rsid: None for rsid in rsids}

    updated = 0
    for rsid, allele in mapping.items():
        if allele is not None:
            existing[rsid] = allele
            updated += 1

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(existing, indent=2), encoding="utf-8")

    print(f"Updated {updated} rsIDs. Output: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
