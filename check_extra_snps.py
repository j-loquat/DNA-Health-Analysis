# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "polars",
# ]
# ///

import sys

import polars as pl

from run_utils import (
    normalize_genotype,
    resolve_base_name,
    resolve_parquet_path,
    run_root,
    update_summary,
    write_json,
)
from snp_reference import load_reference, panel_records

def check_risks(parquet_path: str, base_name: str) -> None:
    print(f"Checking Hidden Risks from {parquet_path}...")
    df = pl.read_parquet(parquet_path)
    
    reference = load_reference()
    records = panel_records(reference, "Hidden Actionable Risks")
    targets = [row["rsid"] for row in records]
    
    results = df.filter(pl.col("rsid").is_in(targets))
    found: dict[str, str] = {}
    for row in results.iter_rows(named=True):
        genotype = normalize_genotype(row["allele1"], row["allele2"])
        if genotype:
            found[row["rsid"]] = genotype

    print("\n--- HIDDEN RISKS REPORT ---")
    for entry in records:
        rsid = entry["rsid"]
        label = entry["label"]
        print(f"{label} ({rsid}): {found.get(rsid, 'Not Found')}")
    print("----------------------------\n")

    run_dir = run_root(base_name)
    payload = {"panel": "Hidden Actionable Risks", "records": records, "genotypes": found}
    write_json(run_dir / "hidden_risks.json", payload)
    update_summary(run_dir, {"hidden_risks_path": str(run_dir / "hidden_risks.json")})

if __name__ == "__main__":
    base_name = resolve_base_name(sys.argv[1] if len(sys.argv) > 1 else None)
    parquet_path = resolve_parquet_path(base_name)
    if not parquet_path.exists():
        sys.exit(1)
    check_risks(str(parquet_path), base_name)
