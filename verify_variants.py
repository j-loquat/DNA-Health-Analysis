# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "polars",
#     "requests",
# ]
# ///

from __future__ import annotations

from dataclasses import dataclass
import json
import sys
import time
from pathlib import Path
from typing import Any, Iterable, TypedDict

import polars as pl
import requests

from snp_reference import load_reference
from run_utils import classify_genotype

ENSEMBL_GRCH37 = "https://grch37.rest.ensembl.org"
REQUEST_TIMEOUT = (5, 10)
MAX_RETRIES = 3
GWAS_RISK_PATH = Path(__file__).resolve().parent / "data" / "gwas_risk_alleles.json"
_GWAS_RISK_CACHE: dict[str, str | None] | None = None
_REFERENCE_RISK_CACHE: dict[str, str] | None = None
_PROXY_NOTES: dict[str, str] = {
    "rs2395029": "Proxy marker for HLA-B*57:01.",
    "rs4349859": "Proxy marker for HLA-B*27.",
    "rs2844682": "Proxy marker for HLA-B*15:02.",
    "rs3909184": "Proxy marker for HLA-B*15:02.",
    "rs1061235": "Proxy marker for HLA-A*31:01.",
    "rs9263726": "Proxy marker for HLA-B*58:01.",
    "rs887829": "Proxy marker for UGT1A1*28.",
}


def _proxy_note(rsid: str) -> str | None:
    return _PROXY_NOTES.get(rsid)


@dataclass(frozen=True)
class VariantVerification:
    rsid: str
    observed_genotype: str | None
    observed_alleles: str | None
    ensembl_alleles: str | None
    ensembl_strand: int | None
    match_status: str
    gwas_risk_allele: str | None
    note: str | None
    proxy_note: str | None


class EnsemblCacheEntry(TypedDict):
    allele_string: str | None
    strand: int | None


class VariantCache(TypedDict):
    ensembl: dict[str, EnsemblCacheEntry]


def _sleep_on_rate_limit(response: requests.Response) -> bool:
    if response.status_code != 429:
        return False
    retry_after = int(response.headers.get("Retry-After", "1"))
    time.sleep(max(retry_after, 1))
    return True


def _get_json(
    session: requests.Session,
    url: str,
    params: dict[str, Any] | None = None,
) -> dict[str, Any] | None:
    headers = {"Accept": "application/json"}
    for attempt in range(MAX_RETRIES):
        try:
            response = session.get(url, params=params, headers=headers, timeout=REQUEST_TIMEOUT)
        except requests.RequestException:
            if attempt == MAX_RETRIES - 1:
                return None
            time.sleep(1 + attempt)
            continue
        if response.status_code == 200:
            return response.json()
        if response.status_code in {400, 404}:
            return None
        if _sleep_on_rate_limit(response):
            continue
        if 500 <= response.status_code < 600:
            time.sleep(1 + attempt)
            continue
        return None
    return None


def _allele_set_from_string(allele_string: str | None) -> set[str]:
    if not allele_string:
        return set()
    cleaned = allele_string.strip().upper()
    if "/" not in cleaned and "|" not in cleaned and len(cleaned) == 2:
        return {base for base in cleaned if base in {"A", "C", "G", "T"}}
    allele_string = allele_string.replace("|", "/")
    parts = [part.strip().upper() for part in allele_string.split("/")]
    return {part for part in parts if part in {"A", "C", "G", "T"}}


def _complement(alleles: Iterable[str]) -> set[str]:
    mapping = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return {mapping[a] for a in alleles if a in mapping}


def _normalize_genotype(allele1: str, allele2: str) -> str | None:
    alleles = [allele1.upper(), allele2.upper()]
    if "0" in alleles or "--" in alleles:
        return None
    cleaned = [a for a in alleles if a in {"A", "C", "G", "T"}]
    if len(cleaned) != 2:
        return None
    return "".join(sorted(cleaned))


def fetch_ensembl_alleles(
    session: requests.Session,
    cache: VariantCache,
    rsid: str,
) -> tuple[str | None, int | None]:
    cached = cache["ensembl"].get(rsid)
    if cached is not None:
        return cached["allele_string"], cached["strand"]
    url = f"{ENSEMBL_GRCH37}/variation/homo_sapiens/{rsid}"
    data = _get_json(session, url)
    if not data:
        cache["ensembl"][rsid] = {"allele_string": None, "strand": None}
        return None, None
    mappings = data.get("mappings") or []
    if not mappings:
        cache["ensembl"][rsid] = {"allele_string": None, "strand": None}
        return None, None
    mapping = mappings[0]
    allele_string = mapping.get("allele_string")
    strand = mapping.get("strand")
    cache["ensembl"][rsid] = {"allele_string": allele_string, "strand": strand}
    return allele_string, strand


def _load_gwas_risk_cache() -> dict[str, str | None]:
    global _GWAS_RISK_CACHE
    if _GWAS_RISK_CACHE is not None:
        return _GWAS_RISK_CACHE
    if not GWAS_RISK_PATH.exists():
        _GWAS_RISK_CACHE = {}
        return _GWAS_RISK_CACHE
    _GWAS_RISK_CACHE = json.loads(GWAS_RISK_PATH.read_text(encoding="utf-8"))
    return _GWAS_RISK_CACHE


def _load_reference_risk_cache() -> dict[str, str]:
    global _REFERENCE_RISK_CACHE
    if _REFERENCE_RISK_CACHE is not None:
        return _REFERENCE_RISK_CACHE
    try:
        reference = load_reference()
    except FileNotFoundError:
        _REFERENCE_RISK_CACHE = {}
        return _REFERENCE_RISK_CACHE
    rows = reference.select(["rsid", "effect_allele"]).to_dicts()
    _REFERENCE_RISK_CACHE = {
        row["rsid"]: row["effect_allele"]
        for row in rows
        if row.get("effect_allele")
    }
    return _REFERENCE_RISK_CACHE


def fetch_gwas_risk_allele(
    rsid: str,
) -> str | None:
    risk_table = _load_gwas_risk_cache()
    risk = risk_table.get(rsid)
    if risk:
        return risk
    reference_risk = _load_reference_risk_cache()
    return reference_risk.get(rsid)


def _load_cache(path: str) -> VariantCache:
    if not Path(path).exists():
        return {"ensembl": {}}
    with open(path, "r", encoding="utf-8") as handle:
        data = json.load(handle)
    return {
        "ensembl": {key: value for key, value in data.get("ensembl", {}).items()},
    }


def _save_cache(path: str, cache: VariantCache) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(cache, handle, indent=2)


def verify_variants(
    parquet_path: str,
    rsids: list[str],
    cache_path: str,
) -> list[VariantVerification]:
    df = pl.read_parquet(parquet_path)
    results = df.filter(pl.col("rsid").is_in(rsids))
    observed_map: dict[str, str | None] = {}
    non_snp_map: dict[str, str] = {}
    for row in results.iter_rows(named=True):
        call = classify_genotype(row["allele1"], row["allele2"])
        if call["kind"] == "acgt":
            observed_map[row["rsid"]] = call["genotype"]
        elif call["kind"] == "non_snp" and call["raw"]:
            non_snp_map[row["rsid"]] = call["raw"]

    verifications: list[VariantVerification] = []
    cache = _load_cache(cache_path)
    session = requests.Session()
    for rsid in rsids:
        proxy_note = _proxy_note(rsid)
        if rsid in non_snp_map:
            allele_string, strand = fetch_ensembl_alleles(session, cache, rsid)
            note = "Non-SNP allele call (indel/repeat); not validated by allele orientation."
            verifications.append(
                VariantVerification(
                    rsid=rsid,
                    observed_genotype=non_snp_map.get(rsid),
                    observed_alleles=non_snp_map.get(rsid),
                    ensembl_alleles=allele_string,
                    ensembl_strand=strand,
                    match_status="non_snp",
                    gwas_risk_allele=None,
                    note=note,
                    proxy_note=proxy_note,
                )
            )
            continue
        genotype = observed_map.get(rsid)
        if not genotype:
            verifications.append(
                VariantVerification(
                    rsid=rsid,
                    observed_genotype=None,
                    observed_alleles=None,
                    ensembl_alleles=None,
                    ensembl_strand=None,
                    match_status="missing_in_file",
                    gwas_risk_allele=None,
                    note=None,
                    proxy_note=proxy_note,
                )
            )
            continue
        observed_set = _allele_set_from_string(genotype)
        ensembl_allele_string, strand = fetch_ensembl_alleles(session, cache, rsid)
        ensembl_set = _allele_set_from_string(ensembl_allele_string)
        comp_set = _complement(ensembl_set)

        match_status = "unknown"
        if observed_set and ensembl_set:
            if observed_set.issubset(ensembl_set):
                match_status = "match"
            elif observed_set.issubset(comp_set):
                match_status = "reverse_complement"
            else:
                match_status = "mismatch"

        gwas_risk = fetch_gwas_risk_allele(rsid)
        note = None
        if match_status == "reverse_complement":
            note = "Observed alleles match the reverse complement of the reference; treat with caution."
        if rsid == "rs4349859" and match_status in {"reverse_complement", "mismatch"}:
            extra = (
                "HLA-B27 proxy (rs4349859) shows strand disagreement. "
                "If file shows T/C while reference shows A/G, flip before interpreting."
            )
            note = f"{note} {extra}".strip() if note else extra
        verifications.append(
            VariantVerification(
                rsid=rsid,
                observed_genotype=genotype,
                observed_alleles="".join(sorted(observed_set)) if observed_set else None,
                ensembl_alleles=ensembl_allele_string,
                ensembl_strand=strand,
                match_status=match_status,
                gwas_risk_allele=gwas_risk,
                note=note,
                proxy_note=proxy_note,
            )
        )
        time.sleep(0.1)
    _save_cache(cache_path, cache)

    return verifications


def _resolve_base_name(arg: str | None) -> str:
    if not arg:
        return "ancestrydna-test-file"
    base = arg
    for suffix in (".normalized.parquet", ".parquet", ".txt"):
        if base.endswith(suffix):
            base = base[: -len(suffix)]
    return base


def main() -> int:
    base_name = _resolve_base_name(sys.argv[1] if len(sys.argv) > 1 else None)
    from run_utils import resolve_parquet_path, run_root, update_summary
    parquet_path = resolve_parquet_path(base_name)
    if not parquet_path.exists():
        print(f"Missing normalized file: {parquet_path}")
        return 1

    rsids = [
        "rs4988235", "rs671", "rs762551", "rs5751876", "rs1815739", "rs4680",
        "rs429358", "rs7412", "rs713598", "rs72921001",
        "rs1801133", "rs1801131", "rs2802292", "rs2282679", "rs602662",
        "rs174546", "rs2274924", "rs1061170", "rs6265", "rs9939609",
        "rs4149056", "rs4244285", "rs1800562", "rs1799945", "rs6025",
        "rs1799963", "rs1333049", "rs3798220",
        "rs7903146", "rs738409", "rs10455872", "rs662799", "rs708272",
        "rs34637584", "rs1052553", "rs17070145", "rs25531",
        "rs6983267", "rs2736100", "rs2981582", "rs10993994",
        "rs9923231", "rs1799853", "rs1057910", "rs1800460", "rs1142345",
        "rs116855232", "rs2395029", "rs3918290", "rs67376798", "rs55886062",
        "rs56038477", "rs75017182", "rs4986893", "rs12248560", "rs776746",
        "rs3745274", "rs2279343", "rs4148323",
        "rs2108622", "rs12777823", "rs28371686", "rs9332131", "rs7900194",
        "rs28371685", "rs2231142", "rs2306283", "rs1800462", "rs887829",
        "rs8175347", "rs2844682", "rs3909184", "rs1061235", "rs9263726",
        "rs16947", "rs1135840", "rs28371725", "rs35742686", "rs5030655",
        "rs334", "rs113993960", "rs28929474", "rs17580", "rs1050828",
        "rs1050829", "rs5742904", "rs80357906", "rs80359550", "rs1801155",
        "rs17879961", "rs10490924",
        "rs1229984", "rs16969968",
        "rs10156191", "rs2052129", "rs11558538",
        "rs1801280", "rs1799930", "rs1799931",
        "rs3177928", "rs7197", "rs4349859", "rs2234693",
    ]

    run_dir = run_root(base_name)
    cache_path = run_dir / "variant_api_cache.json"
    verifications = verify_variants(str(parquet_path), rsids, str(cache_path))

    output_path = run_dir / "variant_verification.json"
    with open(output_path, "w", encoding="utf-8") as handle:
        json.dump([v.__dict__ for v in verifications], handle, indent=2)

    print("\n--- STRAND & ALLELE VERIFICATION ---")
    for v in verifications:
        print(
            f"{v.rsid}: observed={v.observed_genotype or 'NA'} "
            f"| ensembl={v.ensembl_alleles or 'NA'} "
            f"| strand={v.ensembl_strand or 'NA'} "
            f"| match={v.match_status} "
            f"| gwas_risk={v.gwas_risk_allele or 'NA'}"
        )
        if v.note:
            print(f"  NOTE: {v.note}")
    print(f"\nSaved verification results to {output_path}")
    print(f"Saved API cache to {cache_path}")
    hla_warning = any(v.rsid == "rs4349859" and v.note for v in verifications)
    reverse_complement_rsids = [v.rsid for v in verifications if v.match_status == "reverse_complement"]
    mismatch_rsids = [v.rsid for v in verifications if v.match_status == "mismatch"]
    update_summary(
        run_dir,
        {
            "variant_verification_path": str(output_path),
            "hla_b27_strand_warning": hla_warning,
            "reverse_complement_count": len(reverse_complement_rsids),
            "reverse_complement_rsids": reverse_complement_rsids,
            "mismatch_rsids": mismatch_rsids,
        },
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
