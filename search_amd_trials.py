#!/usr/bin/env -S uv --quiet run --active --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "requests",
# ]
# ///

from __future__ import annotations

from datetime import date

import requests

import sys

from run_utils import resolve_base_name, run_root, update_summary, write_json


def search_amd_trials(base_name: str) -> None:
    url = "https://clinicaltrials.gov/api/v2/studies"
    params = {
        "query.term": "age-related macular degeneration",
        "pageSize": 50,
    }

    try:
        response = requests.get(url, params=params, timeout=20)
        response.raise_for_status()
        data = response.json()

        studies = data.get("studies", [])
        total_count = data.get("totalCount")
        if total_count is None:
            total_count = len(studies)

        recruiting: list[dict[str, str]] = []
        for study in studies:
            protocol = study.get("protocolSection", {})
            status = protocol.get("statusModule", {})
            overall = status.get("overallStatus", "").upper()
            if overall != "RECRUITING":
                continue
            ident = protocol.get("identificationModule", {})
            design = protocol.get("designModule", {})
            recruiting.append(
                {
                    "nct_id": ident.get("nctId", "N/A"),
                    "title": ident.get("briefTitle", "N/A"),
                    "phase": ", ".join(design.get("phases", []) or ["N/A"]),
                    "overall_status": overall,
                }
            )

        print(f"Total AMD studies (API reported): {total_count}")
        print(f"Recruiting studies in first page: {len(recruiting)}")
        print("\n--- TOP RECRUITING TRIALS ---")
        for study in recruiting[:5]:
            print(f"[{study['nct_id']}] {study['title']}")
            print(f"Phase: {study['phase']}")
            print("-" * 20)

        run_dir = run_root(base_name)
        payload = {
            "query_date": date.today().isoformat(),
            "query": params,
            "total_count": total_count,
            "recruiting_count": len(recruiting),
            "recruiting_studies": recruiting[:25],
        }
        write_json(run_dir / "amd_trials.json", payload)
        update_summary(run_dir, {"amd_trials_path": str(run_dir / "amd_trials.json")})

    except requests.RequestException as exc:
        print(f"Error: {exc}")


if __name__ == "__main__":
    base_name = resolve_base_name(sys.argv[1] if len(sys.argv) > 1 else None)
    search_amd_trials(base_name)

