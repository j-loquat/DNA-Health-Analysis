#!/usr/bin/env -S uv --quiet run --active --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "polars",
# ]
# ///

from __future__ import annotations

import sys
from typing import Final

import polars as pl

from run_utils import (
    classify_genotype,
    resolve_base_name,
    resolve_parquet_path,
    run_root,
    update_summary,
    write_json,
)
from snp_reference import load_reference, panels_to_records, panel_records


def _format_fun_trait(entry: dict[str, str], genotype: str | None) -> str:
    label = entry["label"]
    rsid = entry["rsid"]
    effect_allele = entry.get("effect_allele") or ""
    effect_trait = entry.get("effect_trait") or ""
    non_effect_trait = entry.get("non_effect_trait") or ""
    evidence_note = entry.get("evidence_note") or ""

    if not genotype:
        return f"{label} ({rsid}): Not Found"
    if not effect_allele:
        return f"{label} ({rsid}): Genotype {genotype}"

    allele_count = sum(1 for allele in genotype if allele == effect_allele)
    if allele_count == 0:
        if non_effect_trait:
            message = (
                f"{label} ({rsid}): Genotype {genotype}. "
                f"This marker is more often linked to {non_effect_trait}."
            )
        else:
            message = (
                f"{label} ({rsid}): Genotype {genotype}. "
                f"This marker does not carry the allele most often linked to {effect_trait}."
            )
    elif allele_count == 1:
        message = (
            f"{label} ({rsid}): Genotype {genotype}. "
            f"One copy of the allele has been associated with {effect_trait}."
        )
    else:
        message = (
            f"{label} ({rsid}): Genotype {genotype}. "
            f"Two copies of the allele have been associated with {effect_trait}, "
            "often with a stronger effect."
        )
    if evidence_note:
        return f"{message} {evidence_note}"
    return message


def check_panels(parquet_path: str, base_name: str) -> None:
    print(f"Checking Expanded Panels from {parquet_path}...")
    df = pl.read_parquet(parquet_path)

    reference = load_reference()
    panel_names: Final[list[str]] = [
        "Cardiometabolic",
        "Neuro/Psych",
        "Cancer",
        "Pharmacogenomics",
        "CYP2D6 Risk Factors",
        "Lifestyle",
        "Functional Health - Histamine",
        "Functional Health - Detox/Acetylation",
        "Functional Health - Inflammation",
        "Functional Health - VDR/Bone",
        "Functional Health - Autoimmune",
        "Functional Health - Hormone",
        "Functional Health - Methylation",
        "Functional Health - Longevity",
        "Functional Health - Neuroplasticity",
        "Functional Health - Oxidative Stress",
        "Functional Health - Metabolic",
        "Functional Health - Iron Metabolism",
    ]
    panels = panels_to_records(reference, panel_names)

    fun_panel_names: Final[list[str]] = [
        "Fun Traits - Appearance",
        "Fun Traits - Sensory",
        "Fun Traits - Body & Appetite",
        "Fun Traits - Sleep",
        "Fun Traits - Behavior",
    ]
    fun_panels = {name: panel_records(reference, name) for name in fun_panel_names}

    all_targets: list[str] = [row["rsid"] for rows in panels.values() for row in rows]
    all_targets.extend([row["rsid"] for rows in fun_panels.values() for row in rows])

    results = df.filter(pl.col("rsid").is_in(all_targets))
    found: dict[str, str] = {}
    non_snp: dict[str, str] = {}
    for row in results.iter_rows(named=True):
        call = classify_genotype(row["allele1"], row["allele2"])
        if call["kind"] == "acgt" and call["genotype"]:
            found[row["rsid"]] = call["genotype"]
        elif call["kind"] == "non_snp" and call["raw"]:
            non_snp[row["rsid"]] = call["raw"]

    print("\n--- EXPANDED PANELS REPORT ---")
    for panel_name in panel_names:
        print(f">> {panel_name}")
        for entry in panels[panel_name]:
            rsid = entry["rsid"]
            label = entry["label"]
            print(f"{label} ({rsid}): {found.get(rsid, 'Not Found')}")
        print("")

    print(">> Fun Traits & Appearance")
    for section_name in fun_panel_names:
        print(f"\n{section_name.replace('Fun Traits - ', '')}")
        for entry in fun_panels[section_name]:
            rsid = entry["rsid"]
            print(_format_fun_trait(entry, found.get(rsid)))
    print(
        "These are statistical associations that can vary by ancestry and environment. "
        "This is not medical advice."
    )

    print("----------------------------\n")

    run_dir = run_root(base_name)
    payload = {
        "panels": panels,
        "fun_panels": fun_panels,
        "genotypes": found,
        "non_snp_genotypes": non_snp,
    }
    write_json(run_dir / "expanded_panels.json", payload)
    update_summary(run_dir, {"expanded_panels_path": str(run_dir / "expanded_panels.json")})


if __name__ == "__main__":
    base_name = resolve_base_name(sys.argv[1] if len(sys.argv) > 1 else None)
    parquet_path = resolve_parquet_path(base_name)
    if not parquet_path.exists():
        sys.exit(1)
    check_panels(str(parquet_path), base_name)

