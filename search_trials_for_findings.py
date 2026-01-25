#!/usr/bin/env -S uv --quiet run --active --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "requests",
# ]
# ///

from __future__ import annotations

from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Any

import argparse
import requests

import generate_report as report
from run_utils import resolve_base_name, run_root, update_summary, write_json


@dataclass(frozen=True)
class TrialQuery:
    label: str
    condition: str
    query_term: str | None = None


FINDING_QUERY_MAP: dict[str, TrialQuery] = {
    "Vision": TrialQuery("Vision", "age-related macular degeneration"),
    "Heart Health": TrialQuery("Heart Health", "lipoprotein(a)"),
    "Clotting Risk": TrialQuery("Clotting Risk", "venous thromboembolism"),
    "Addiction Risk": TrialQuery("Addiction Risk", "smoking cessation"),
    "Methylation": TrialQuery("Methylation", "hyperhomocysteinemia"),
    "Detox/Drug Metabolism": TrialQuery(
        "Detox/Drug Metabolism",
        "tuberculosis",
        query_term="isoniazid pharmacogenomics",
    ),
    "Fluoropyrimidine Toxicity": TrialQuery(
        "Fluoropyrimidine Toxicity",
        "cancer",
        query_term="fluoropyrimidine pharmacogenomics",
    ),
}


def _fetch_trials(
    query: TrialQuery,
    *,
    location: str | None,
    geo_filter: str | None,
    page_size: int = 50,
) -> dict[str, Any]:
    url = "https://clinicaltrials.gov/api/v2/studies"
    params = {
        "query.cond": query.condition,
        "filter.overallStatus": "RECRUITING",
        "pageSize": page_size,
        "sort": "LastUpdatePostDate:desc",
    }
    if query.query_term:
        params["query.term"] = query.query_term
    if location:
        params["query.locn"] = location
    if geo_filter:
        params["filter.geo"] = geo_filter

    response = requests.get(url, params=params, timeout=20)
    response.raise_for_status()
    data = response.json()

    studies = data.get("studies", [])
    total_count = data.get("totalCount", len(studies))

    recruiting: list[dict[str, str]] = []
    for study in studies:
        protocol = study.get("protocolSection", {})
        status = protocol.get("statusModule", {})
        overall = status.get("overallStatus", "").upper()
        if overall != "RECRUITING":
            continue
        ident = protocol.get("identificationModule", {})
        design = protocol.get("designModule", {})
        nct_id = ident.get("nctId", "N/A")
        recruiting.append(
            {
                "nct_id": nct_id,
                "title": ident.get("briefTitle", "N/A"),
                "phase": ", ".join(design.get("phases", []) or ["N/A"]),
                "overall_status": overall,
                "url": f"https://clinicaltrials.gov/study/{nct_id}" if nct_id != "N/A" else "",
            }
        )

    return {
        "total_count": total_count,
        "recruiting_count": len(recruiting),
        "recruiting_studies": recruiting[:25],
    }


def _load_genotypes(run_dir: Path) -> dict[str, str]:
    core_traits = report._load_json(run_dir / "core_traits.json")
    healthy = report._load_json(run_dir / "healthy_aging.json")
    hidden = report._load_json(run_dir / "hidden_risks.json")
    expanded = report._load_json(run_dir / "expanded_panels.json")
    return report._merge_genotypes(core_traits, healthy, hidden, expanded)


def search_trials_for_findings(
    base_name: str,
    *,
    location: str | None,
    geo_filter: str | None,
) -> None:
    run_dir = run_root(base_name)
    genotypes = _load_genotypes(run_dir)
    risk_cards = [
        card for card in report._build_risk_cards(genotypes)
        if card.get("category") == "clinical"
    ]

    if not risk_cards:
        print("No critical findings detected; skipping clinical trials search.")
        return

    findings_payload: list[dict[str, Any]] = []
    for card in risk_cards:
        query = FINDING_QUERY_MAP.get(card["label"])
        if not query:
            continue
        result = _fetch_trials(query, location=location, geo_filter=geo_filter)
        findings_payload.append(
            {
                "finding_label": card["label"],
                "finding_level": card["level"],
                "query_term": query.query_term or query.condition,
                "total_count": result["total_count"],
                "recruiting_count": result["recruiting_count"],
                "recruiting_studies": result["recruiting_studies"],
            }
        )

    if not findings_payload:
        print("No trial queries mapped to detected findings; skipping output.")
        return

    payload = {
        "query_date": date.today().isoformat(),
        "findings": findings_payload,
    }
    output_path = run_dir / "trials_by_finding.json"
    write_json(output_path, payload)
    update_summary(run_dir, {"trials_by_finding_path": str(output_path)})
    print(f"Wrote trials to {output_path}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Search ClinicalTrials.gov based on detected critical findings."
    )
    parser.add_argument("base_name", nargs="?", help="Base filename without extension")
    parser.add_argument(
        "--location",
        help="Optional location term (city/state/country) to filter recruiting studies.",
    )
    parser.add_argument(
        "--geo",
        dest="geo_filter",
        help="Optional geo filter, e.g. distance(39.00,-77.10,50mi).",
    )
    args = parser.parse_args()

    base_name = resolve_base_name(args.base_name)
    search_trials_for_findings(base_name, location=args.location, geo_filter=args.geo_filter)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
