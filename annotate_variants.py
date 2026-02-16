# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "requests",
# ]
# ///

from __future__ import annotations

import csv
import json
import os
import time
from datetime import datetime
from pathlib import Path
from typing import Any

import requests

from run_utils import resolve_base_name, run_root, update_summary

NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GNOMAD_API = "https://gnomad.broadinstitute.org/api"
REQUEST_TIMEOUT = (5, 20)
MAX_RETRIES = 4


def _load_json(path: Path) -> Any:
    if not path.exists():
        return {}
    raw = path.read_bytes()
    for encoding in ("utf-8-sig", "utf-8", "cp1252"):
        try:
            return json.loads(raw.decode(encoding))
        except (UnicodeDecodeError, json.JSONDecodeError):
            continue
    return {}


def _load_reference_rsids() -> set[str]:
    path = Path("data") / "snp_reference.csv"
    if not path.exists():
        return set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return {str(row.get("rsid", "")).strip() for row in reader if str(row.get("rsid", "")).strip()}


def _hidden_actionable_rsids() -> set[str]:
    path = Path("data") / "snp_reference.csv"
    if not path.exists():
        return set()
    rsids: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if str(row.get("panel", "")).strip() == "Hidden Actionable Risks":
                rsid = str(row.get("rsid", "")).strip()
                if rsid:
                    rsids.add(rsid)
    return rsids


def _request_json(
    session: requests.Session,
    url: str,
    *,
    params: dict[str, Any] | None = None,
    payload: dict[str, Any] | None = None,
) -> dict[str, Any] | None:
    for attempt in range(MAX_RETRIES):
        try:
            if payload is None:
                resp = session.get(url, params=params, timeout=REQUEST_TIMEOUT)
            else:
                resp = session.post(url, json=payload, timeout=REQUEST_TIMEOUT)
        except requests.RequestException:
            time.sleep(1 + attempt)
            continue
        if resp.status_code == 200:
            try:
                return resp.json()
            except ValueError:
                return None
        if resp.status_code in {429, 500, 502, 503, 504}:
            retry_after = resp.headers.get("Retry-After")
            if retry_after and str(retry_after).isdigit():
                time.sleep(max(int(retry_after), 1))
            else:
                time.sleep(1 + attempt * 2)
            continue
        return None
    return None


def _extract_clinvar_best(entry_map: dict[str, Any]) -> dict[str, Any] | None:
    """Pick the most useful ClinVar summary entry from an esummary result block."""
    if not entry_map:
        return None

    # The esummary response puts records at result.<uid>; this function expects a uid->record mapping.
    candidates: list[dict[str, Any]] = []
    for record in entry_map.values():
        if not isinstance(record, dict):
            continue
        germline = record.get("germline_classification")
        if not isinstance(germline, dict):
            germline = {}
        description = str(germline.get("description") or "").strip().lower()
        review = str(germline.get("review_status") or "").strip().lower()
        obj_type = str(record.get("obj_type") or "").strip().lower()

        score = 0
        if obj_type and obj_type != "haplotype":
            score += 2
        if description and description not in {"other", "not provided", "unspecified"}:
            score += 3
        if review and "no assertion" not in review:
            score += 1
        candidates.append({"score": score, "record": record})

    if not candidates:
        return None
    candidates.sort(key=lambda item: item["score"], reverse=True)
    return candidates[0]["record"]


def _clinvar_summary(session: requests.Session, rsid: str) -> dict[str, Any] | None:
    api_key = os.environ.get("NCBI_API_KEY")
    email = os.environ.get("NCBI_EMAIL")
    common_params = {"retmode": "json", "tool": "DNA-Health-Analysis"}
    if api_key:
        common_params["api_key"] = api_key
    if email:
        common_params["email"] = email

    search = _request_json(
        session,
        f"{NCBI_EUTILS}/esearch.fcgi",
        params={**common_params, "db": "clinvar", "term": f"{rsid}[variant]"},
    )
    idlist = (((search or {}).get("esearchresult") or {}).get("idlist") or []) if isinstance(search, dict) else []
    ids = [str(item) for item in idlist if str(item).strip()][:6]
    if not ids:
        return None

    summary = _request_json(
        session,
        f"{NCBI_EUTILS}/esummary.fcgi",
        params={**common_params, "db": "clinvar", "id": ",".join(ids)},
    )
    result = (summary or {}).get("result") if isinstance(summary, dict) else None
    if not isinstance(result, dict):
        return None

    uid_map: dict[str, Any] = {}
    for uid in ids:
        # Esumm uses numeric property names; JSON parsing keeps them as strings.
        record = result.get(uid)
        if isinstance(record, dict):
            uid_map[uid] = record
    best = _extract_clinvar_best(uid_map)
    if not best:
        return None

    germline = best.get("germline_classification") if isinstance(best.get("germline_classification"), dict) else {}
    traits = []
    for trait in germline.get("trait_set", []) if isinstance(germline, dict) else []:
        if isinstance(trait, dict):
            name = str(trait.get("trait_name") or "").strip()
            if name:
                traits.append(name)

    return {
        "uid": str(best.get("uid") or ""),
        "accession": str(best.get("accession") or ""),
        "accession_version": str(best.get("accession_version") or ""),
        "title": str(best.get("title") or ""),
        "obj_type": str(best.get("obj_type") or ""),
        "description": str(germline.get("description") or ""),
        "review_status": str(germline.get("review_status") or ""),
        "last_evaluated": str(germline.get("last_evaluated") or ""),
        "traits": traits[:6],
    }


def _parse_ensembl_location(location: str) -> tuple[str, int] | None:
    # "19:19379549-19379549"
    if not location or ":" not in location:
        return None
    chrom, rest = location.split(":", 1)
    start = rest.split("-", 1)[0]
    if not start.isdigit():
        return None
    return chrom, int(start)


def _gnomad_variant_af(
    session: requests.Session,
    *,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    dataset: str = "gnomad_r2_1",
) -> dict[str, Any] | None:
    variant_id = f"{chrom}-{pos}-{ref}-{alt}"
    query = (
        "query($variantId: String!, $dataset: DatasetId!) {"
        "  variant(variantId: $variantId, dataset: $dataset) {"
        "    variantId rsid exome { af } genome { af }"
        "  }"
        "}"
    )
    payload = {"query": query, "variables": {"variantId": variant_id, "dataset": dataset}}
    data = _request_json(session, GNOMAD_API, payload=payload)
    if not isinstance(data, dict):
        return None
    if data.get("errors"):
        return None
    variant = (data.get("data") or {}).get("variant")
    if not isinstance(variant, dict):
        return None
    exome = variant.get("exome") if isinstance(variant.get("exome"), dict) else {}
    genome = variant.get("genome") if isinstance(variant.get("genome"), dict) else {}
    return {
        "dataset": dataset,
        "variant_id": str(variant.get("variantId") or variant_id),
        "rsid": str(variant.get("rsid") or ""),
        "exome_af": exome.get("af"),
        "genome_af": genome.get("af"),
    }


def annotate(base_name: str) -> None:
    run_dir = run_root(base_name)

    core_traits = _load_json(run_dir / "core_traits.json")
    healthy = _load_json(run_dir / "healthy_aging.json")
    hidden = _load_json(run_dir / "hidden_risks.json")
    expanded = _load_json(run_dir / "expanded_panels.json")
    variant_verification = _load_json(run_dir / "variant_verification.json")

    verification_lookup: dict[str, dict[str, Any]] = {}
    if isinstance(variant_verification, list):
        for entry in variant_verification:
            if isinstance(entry, dict) and isinstance(entry.get("rsid"), str):
                verification_lookup[entry["rsid"]] = entry

    # Use the same merge logic as the report (without importing generate_report.py).
    def merge_genotypes(*payloads: Any) -> dict[str, str]:
        merged: dict[str, str] = {}
        for payload in payloads:
            if not isinstance(payload, dict):
                continue
            genotypes = payload.get("genotypes", {})
            if isinstance(genotypes, dict):
                for rsid, genotype in genotypes.items():
                    if isinstance(rsid, str) and isinstance(genotype, str) and genotype.strip():
                        merged[rsid] = genotype.strip()
        return merged

    genotypes = merge_genotypes(core_traits, healthy, hidden, expanded)

    target_rsids = _hidden_actionable_rsids()
    # Ensure newly-added high-value markers are annotated even if panel placement changes.
    target_rsids.update({"rs58542926", "rs73885319", "rs60910145", "rs71785313", "rs10484555"})
    # Avoid wasting calls on rsids we don't recognize locally.
    reference_rsids = _load_reference_rsids()
    target_rsids = {rsid for rsid in target_rsids if rsid in reference_rsids}

    session = requests.Session()
    annotations: dict[str, dict[str, Any]] = {}
    for rsid in sorted(target_rsids):
        verification = verification_lookup.get(rsid, {})
        allele_string = str(verification.get("ensembl_alleles") or "").strip().upper()
        ref = str(verification.get("ensembl_ref_allele") or "").strip().upper()
        location = str(verification.get("ensembl_location") or "").strip()
        assembly = str(verification.get("ensembl_assembly") or "").strip() or "GRCh37"

        payload: dict[str, Any] = {"updated_at": datetime.now().isoformat(timespec="seconds")}
        clinvar = _clinvar_summary(session, rsid)
        if clinvar:
            payload["clinvar"] = clinvar

        # Best-effort gnomAD AF (may be rate-limited / overloaded).
        loc = _parse_ensembl_location(location)
        alleles = [a.strip() for a in allele_string.replace("|", "/").split("/") if a.strip()]
        alleles = [a for a in alleles if a in {"A", "C", "G", "T"}]
        if not ref and alleles:
            ref = alleles[0]
        if loc and ref in {"A", "C", "G", "T"} and assembly.upper().startswith("GRCH37"):
            chrom, pos = loc
            alts = [a for a in alleles if a != ref]
            best: dict[str, Any] | None = None
            best_af = -1.0
            for alt in alts[:3]:
                result = _gnomad_variant_af(session, chrom=chrom, pos=pos, ref=ref, alt=alt)
                if not result:
                    continue
                af_candidates = []
                for key in ("exome_af", "genome_af"):
                    value = result.get(key)
                    if isinstance(value, (int, float)):
                        af_candidates.append(float(value))
                af_value = max(af_candidates) if af_candidates else -1.0
                if af_value > best_af:
                    best_af = af_value
                    best = result
                time.sleep(0.2)
            if best:
                payload["gnomad"] = best

        if any(key in payload for key in ("clinvar", "gnomad")):
            annotations[rsid] = payload

        # Be polite to NCBI; keep request rate low when no API key is configured.
        time.sleep(0.34 if not os.environ.get("NCBI_API_KEY") else 0.12)

    out_path = run_dir / "variant_annotations.json"
    out_path.write_text(json.dumps(annotations, indent=2), encoding="utf-8")
    update_summary(run_dir, {"variant_annotations_path": str(out_path)})
    print(f"Wrote {len(annotations)} variant annotations to {out_path}")


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(description="Annotate high-value variants with ClinVar and gnomAD metadata.")
    parser.add_argument("base_name", help="Base filename without extension")
    args = parser.parse_args()

    base_name = resolve_base_name(args.base_name)
    annotate(base_name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

