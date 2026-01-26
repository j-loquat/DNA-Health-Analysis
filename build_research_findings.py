#!/usr/bin/env -S uv --quiet run --active --script
# /// script
# requires-python = ">=3.12"
# dependencies = []
# ///

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import generate_report as report
from run_utils import resolve_base_name, write_json


LEVEL_ORDER = {"high": 0, "med": 1, "low": 2, "neutral": 3}

RESEARCH_TEMPLATES: dict[str, dict[str, str]] = {
    "Clotting Risk": {
        "topic": "Prothrombin G20210A / Factor V Leiden: VTE recurrence and guidance",
        "suggested_query": "Prothrombin G20210A heterozygous recurrence VTE guideline 2024",
    },
    "Clopidogrel Response": {
        "topic": "CYP2C19 and clopidogrel: CPIC and clinical outcomes",
        "suggested_query": "CPIC clopidogrel CYP2C19 2022 update",
    },
    "Statin Myopathy Risk": {
        "topic": "SLCO1B1 and statin-associated myopathy: CPIC guidance",
        "suggested_query": "CPIC statin SLCO1B1 2022 guideline",
    },
    "Tacrolimus Metabolism": {
        "topic": "CYP3A5 and tacrolimus dosing: CPIC guidance",
        "suggested_query": "CPIC tacrolimus CYP3A5 guideline",
    },
    "Efavirenz Metabolism": {
        "topic": "CYP2B6 and efavirenz: CPIC guidance",
        "suggested_query": "CPIC efavirenz CYP2B6 guideline 2019",
    },
    "Fluoropyrimidine Toxicity": {
        "topic": "DPYD and fluoropyrimidines: CPIC guidance and updates",
        "suggested_query": "CPIC DPYD fluoropyrimidine guideline 2024 clarification",
    },
    "Thiopurine Toxicity": {
        "topic": "TPMT/NUDT15 and thiopurines: CPIC guidance",
        "suggested_query": "CPIC thiopurine TPMT NUDT15 guideline",
    },
    "Abacavir Hypersensitivity": {
        "topic": "HLA-B*57:01 and abacavir: CPIC guidance",
        "suggested_query": "CPIC abacavir HLA-B*57:01 guideline",
    },
    "Warfarin Sensitivity": {
        "topic": "VKORC1/CYP2C9 and warfarin: CPIC dosing algorithms",
        "suggested_query": "CPIC warfarin CYP2C9 VKORC1 2017 guideline",
    },
    "CYP2C9 Reduced Function": {
        "topic": "CYP2C9 reduced function: warfarin/NSAID dosing guidance",
        "suggested_query": "CPIC CYP2C9 dosing warfarin NSAID guidance",
    },
    "Detox/Drug Metabolism": {
        "topic": "NAT2 and isoniazid: genotype-guided dosing evidence",
        "suggested_query": "NAT2 genotype guided isoniazid dosing randomized trial",
    },
}


def _load_genotypes(run_dir: Path) -> dict[str, str]:
    core_traits = report._load_json(run_dir / "core_traits.json")
    healthy = report._load_json(run_dir / "healthy_aging.json")
    hidden = report._load_json(run_dir / "hidden_risks.json")
    expanded = report._load_json(run_dir / "expanded_panels.json")
    return report._merge_genotypes(core_traits, healthy, hidden, expanded)


def _risk_cards(run_dir: Path, min_level: str) -> list[dict[str, Any]]:
    genotypes = _load_genotypes(run_dir)
    variant_verification = report._load_json(run_dir / "variant_verification.json")
    lookup = report._variant_lookup(variant_verification if isinstance(variant_verification, list) else [])
    cards = report._build_risk_cards(genotypes, lookup)
    min_rank = LEVEL_ORDER[min_level]
    return [
        card for card in cards
        if LEVEL_ORDER.get(card.get("level", "neutral"), 99) <= min_rank
    ]


def build_research_findings(
    base_name: str,
    *,
    run_date: str | None,
    min_level: str,
    write_empty: bool,
) -> int:
    run_dir = report._find_run_dir(base_name, run_date)
    cards = _risk_cards(run_dir, min_level)
    if not cards and not write_empty:
        print("No high-priority findings detected; skipping research_findings.json.")
        return 0

    findings: list[dict[str, str]] = []
    for card in cards:
        label = card.get("label", "Research Topic")
        template = RESEARCH_TEMPLATES.get(label, {})
        topic = template.get("topic", label)
        query = template.get("suggested_query", "")
        findings.append(
            {
                "topic": topic,
                "content": "",
                "source": "",
                "suggested_query": query,
            }
        )

    output_path = run_dir / "research_findings.json"
    write_json(output_path, findings)
    print(f"Wrote research template to {output_path}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create a research_findings.json template for high-priority findings."
    )
    parser.add_argument("base_name", help="Base filename without extension")
    parser.add_argument("--run-date", help="Run date in YYYYMMDD (optional)")
    parser.add_argument(
        "--min-level",
        choices=["high", "med", "low"],
        default="high",
        help="Minimum severity to include in the research template.",
    )
    parser.add_argument(
        "--write-empty",
        action="store_true",
        help="Write an empty file even if no findings meet the threshold.",
    )
    args = parser.parse_args()

    base_name = resolve_base_name(args.base_name)
    return build_research_findings(
        base_name,
        run_date=args.run_date,
        min_level=args.min_level,
        write_empty=args.write_empty,
    )


if __name__ == "__main__":
    raise SystemExit(main())
