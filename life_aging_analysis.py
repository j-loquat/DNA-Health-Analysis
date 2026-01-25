# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "polars",
# ]
# ///

import sys

import polars as pl

from run_utils import resolve_base_name, resolve_parquet_path, run_root, update_summary, write_json
from snp_reference import load_reference, panel_records

def analyze_aging(parquet_path: str, base_name: str) -> None:
    print(f"Analyzing Aging & Lifestyle from {parquet_path}...")
    df = pl.read_parquet(parquet_path)
    
    reference = load_reference()
    records = panel_records(reference, "Healthy Aging")
    targets = [row["rsid"] for row in records]
    
    results = df.filter(pl.col("rsid").is_in(targets))
    found = {}
    for row in results.iter_rows(named=True):
        found[row['rsid']] = "".join(sorted([row['allele1'], row['allele2']]))

    print("\n--- AGING & LIFESTYLE REPORT ---")
    for entry in records:
        rsid = entry["rsid"]
        label = entry["label"]
        print(f"{label} ({rsid}): {found.get(rsid, 'Not Found')}")
    print("----------------------------\n")

    run_dir = run_root(base_name)
    payload = {"panel": "Healthy Aging", "records": records, "genotypes": found}
    write_json(run_dir / "healthy_aging.json", payload)
    update_summary(run_dir, {"healthy_aging_path": str(run_dir / "healthy_aging.json")})

if __name__ == "__main__":
    base_name = resolve_base_name(sys.argv[1] if len(sys.argv) > 1 else None)
    parquet_path = resolve_parquet_path(base_name)
    if not parquet_path.exists():
        sys.exit(1)
    analyze_aging(str(parquet_path), base_name)
