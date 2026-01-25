# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "polars",
# ]
# ///

import json
import sys
from pathlib import Path
from typing import Final

import polars as pl

from run_utils import resolve_base_name, resolve_parquet_path, run_root, update_summary, write_json
from snp_reference import load_reference, panels_to_records

CLINICAL_PATH = Path(__file__).resolve().parent / "data" / "clinical_interpretations.json"


def _load_apoe_map() -> dict[str, str]:
    if not CLINICAL_PATH.exists():
        return {}
    data = json.loads(CLINICAL_PATH.read_text(encoding="utf-8"))
    return data.get("apoe_haplotype_map", {})


def query_core_traits(parquet_path: str, base_name: str) -> None:
    print(f"Querying Core Traits from {parquet_path}...")

    df = pl.read_parquet(parquet_path)

    reference = load_reference()
    panel_names: Final[list[str]] = [
        "Core Wellness",
        "Structural & Functional",
        "Dietary Sensitivities",
        "Celiac / Gluten Tags",
        "Sleep & Chronotype",
    ]
    panels = panels_to_records(reference, panel_names)

    target_snps = [row["rsid"] for rows in panels.values() for row in rows]

    results = df.filter(pl.col("rsid").is_in(target_snps))

    found_snps: dict[str, str] = {}
    for row in results.iter_rows(named=True):
        rsid = row["rsid"]
        a1 = row["allele1"]
        a2 = row["allele2"]
        genotype = "".join(sorted([a1, a2]))
        found_snps[rsid] = genotype

    print("\n--- CORE WELLNESS AND LIFESTYLE REPORT ---")
    # Note: Genotypes are sorted alphabetical (e.g., AG, not GA)
    for panel_name, entries in panels.items():
        print(f"\n>> {panel_name}")
        for entry in entries:
            label = entry["label"]
            rsid = entry["rsid"]
            print(f"{label} ({rsid}): {found_snps.get(rsid, 'Not Found')}")

    apoe_map = _load_apoe_map()
    rs429358 = found_snps.get("rs429358")
    rs7412 = found_snps.get("rs7412")
    if rs429358 and rs7412 and apoe_map:
        key = f"{rs429358}|{rs7412}"
        haplotype = apoe_map.get(key, "Unknown")
        print(f"\nAPOE haplotype (rs429358|rs7412): {key} -> {haplotype}")

    print("----------------------------\n")

    run_dir = run_root(base_name)
    payload = {"panels": panels, "genotypes": found_snps}
    write_json(run_dir / "core_traits.json", payload)
    update_summary(run_dir, {"core_traits_path": str(run_dir / "core_traits.json")})

if __name__ == "__main__":
    base_name = resolve_base_name(sys.argv[1] if len(sys.argv) > 1 else None)
    parquet_file = resolve_parquet_path(base_name)
    if not parquet_file.exists():
        print("Parquet file not found. Run qc_analysis.py first.")
        sys.exit(1)
    query_core_traits(str(parquet_file), base_name)
