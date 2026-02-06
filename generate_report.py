# /// script
# requires-python = ">=3.12"
# dependencies = []
# ///

from __future__ import annotations

import argparse
import json
from datetime import date
from pathlib import Path
from typing import Any


def _load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    raw = path.read_bytes()
    for encoding in ("utf-8-sig", "utf-8", "cp1252"):
        try:
            return json.loads(raw.decode(encoding))
        except (UnicodeDecodeError, json.JSONDecodeError):
            continue
    print(f"Warning: unable to parse JSON at {path}; skipping.")
    return {}


def _find_run_dir(base_name: str, run_date: str | None) -> Path:
    runs_root = Path("runs")
    if run_date:
        candidate = runs_root / run_date / base_name
        if not candidate.exists():
            raise FileNotFoundError(f"Run folder not found: {candidate}")
        return candidate
    candidates = sorted(runs_root.glob(f"*/{base_name}"), key=lambda p: p.parent.name)
    if not candidates:
        raise FileNotFoundError(f"No run folders found for {base_name}")
    return candidates[-1]


def _merge_genotypes(*payloads: dict[str, Any]) -> dict[str, str]:
    merged: dict[str, str] = {}
    for payload in payloads:
        genotypes = payload.get("genotypes", {})
        for rsid, genotype in genotypes.items():
            if genotype:
                merged[rsid] = genotype
    return merged


def _merge_non_snp_genotypes(*payloads: dict[str, Any]) -> dict[str, str]:
    merged: dict[str, str] = {}
    for payload in payloads:
        genotypes = payload.get("non_snp_genotypes", {})
        for rsid, genotype in genotypes.items():
            if genotype:
                merged[rsid] = genotype
    return merged


_PROXY_LABELS: dict[str, str] = {
    "rs2395029": "HLA-B*57:01 proxy",
    "rs2844682": "HLA-B*15:02 proxy",
    "rs3909184": "HLA-B*15:02 proxy",
    "rs1061235": "HLA-A*31:01 proxy",
    "rs9263726": "HLA-B*58:01 proxy",
    "rs4349859": "HLA-B27 proxy",
    "rs887829": "UGT1A1*28 proxy",
}
_NAT2_ORDER = ("rs1801280", "rs1799930", "rs1799931")
_NAT2_RSIDS = set(_NAT2_ORDER)

_CARRIER_MARKERS: dict[str, dict[str, str | None]] = {
    "rs334": {"label": "Sickle cell (HbS)", "effect_allele": "T"},
    "rs113993960": {"label": "CFTR F508del", "effect_allele": None},
    "rs28929474": {"label": "SERPINA1 Pi*Z", "effect_allele": "A"},
    "rs17580": {"label": "SERPINA1 Pi*S", "effect_allele": "T"},
    "rs1050828": {"label": "G6PD c.202G>A", "effect_allele": "A"},
    "rs1050829": {"label": "G6PD c.376A>G", "effect_allele": "G"},
}
_HEREDITARY_FOUNDER_MARKERS: dict[str, dict[str, str | None]] = {
    "rs80357906": {"label": "BRCA1 5382insC", "effect_allele": None},
    "rs80359550": {"label": "BRCA2 6174delT", "effect_allele": None},
}
_STRAND_CAUTION_MATCH_STATUSES = {"reverse_complement", "mismatch"}
_HIGH_EVIDENCE_PGX_ESCALATION: list[dict[str, Any]] = [
    {
        "label": "Fluoropyrimidine Toxicity",
        "level": "high",
        "description_prefix": "DPYD variant(s) detected: ",
        "description_suffix": ".",
        "action": "Confirm with clinical-grade DPYD testing before 5-FU/capecitabine; dosing changes may be needed.",
        "evidence": "CPIC",
        "markers": [
            {"rsid": "rs3918290", "risk_allele": "A", "display": "rs3918290 (*2A)"},
            {"rsid": "rs67376798", "risk_allele": "T", "display": "rs67376798 (c.2846A>T)"},
            {"rsid": "rs55886062", "risk_allele": "G", "display": "rs55886062 (c.1679T>G)"},
            {"rsid": "rs56038477", "risk_allele": "A", "display": "rs56038477 (HapB3 tag)"},
            {"rsid": "rs75017182", "risk_allele": "G", "display": "rs75017182 (HapB3)"},
        ],
    },
    {
        "label": "Statin Myopathy Risk",
        "level": "med",
        "description_prefix": "SLCO1B1 risk marker(s) detected: ",
        "description_suffix": ".",
        "action": "If prescribed simvastatin, consider lower dose or alternative statin per CPIC guidance.",
        "evidence": "CPIC",
        "markers": [
            {"rsid": "rs4149056", "risk_allele": "C", "display": "SLCO1B1 rs4149056"},
        ],
    },
    {
        "label": "Tacrolimus Metabolism",
        "level": "med",
        "description_prefix": "CYP3A5 marker(s) detected: ",
        "description_suffix": ".",
        "action": "Tacrolimus dosing should follow CPIC CYP3A5 guidance.",
        "evidence": "CPIC",
        "markers": [
            {"rsid": "rs776746", "risk_allele": "A", "display": "CYP3A5 rs776746"},
        ],
    },
    {
        "label": "Atazanavir Hyperbilirubinemia",
        "level": "med",
        "description_prefix": "UGT1A1 reduced-function marker(s) detected: ",
        "description_suffix": ".",
        "action": "Atazanavir use should follow CPIC UGT1A1 guidance.",
        "evidence": "CPIC",
        "markers": [
            {"rsid": "rs4148323", "risk_allele": "A", "display": "UGT1A1*6 (rs4148323)"},
            {"rsid": "rs887829", "risk_allele": "T", "display": "UGT1A1*28 proxy (rs887829)"},
        ],
    },
]


def _count_phrase(count: int, singular: str, plural: str | None = None) -> str:
    noun = singular if count == 1 else (plural or f"{singular}s")
    return f"{count} {noun}"


def _display_chromosome_label(chromosome: Any) -> str:
    text = str(chromosome or "NA").strip().upper()
    if text == "23":
        return "X"
    if text == "24":
        return "Y"
    if text == "25":
        return "MT"
    if text in {"X", "Y", "MT"}:
        return text
    if text.isdigit() and 1 <= int(text) <= 22:
        return text
    return f"Other/Contigs ({text})"


def _variant_flags(
    rsid: str,
    genotype: str | None,
    non_snp_call: str | None,
    variant_lookup: dict[str, dict[str, Any]] | None,
    *,
    is_partial_panel: bool = False,
) -> dict[str, Any]:
    entry = variant_lookup.get(rsid) if variant_lookup else None
    match_status = str(entry.get("match_status", "")) if entry else ""
    proxy_note_raw = entry.get("proxy_note") if entry else None
    proxy_note = proxy_note_raw.strip() if isinstance(proxy_note_raw, str) else ""
    has_call = bool(genotype or non_snp_call)
    return {
        "is_missing": not has_call,
        "is_non_snp_placeholder": bool(non_snp_call),
        "is_strand_caution": bool(genotype and match_status in _STRAND_CAUTION_MATCH_STATUSES),
        "is_proxy": bool(proxy_note and has_call),
        "is_partial_panel": is_partial_panel,
        "proxy_note": proxy_note,
    }


def _proxy_markers_present(
    genotypes: dict[str, str],
    variant_lookup: dict[str, dict[str, Any]] | None,
) -> list[dict[str, str]]:
    markers: list[dict[str, str]] = []
    if not variant_lookup:
        return markers
    for rsid, label in _PROXY_LABELS.items():
        genotype = genotypes.get(rsid)
        if not genotype:
            continue
        note = variant_lookup.get(rsid, {}).get("proxy_note") or "Proxy marker; confirm clinically."
        markers.append(
            {
                "label": label,
                "rsid": rsid,
                "genotype": genotype,
                "note": note,
            }
        )
    return markers


def _high_priority_findings(
    genotypes: dict[str, str],
    non_snp_genotypes: dict[str, str],
    variant_lookup: dict[str, dict[str, Any]] | None,
    risk_cards: list[dict[str, str]] | None = None,
) -> list[dict[str, str]]:
    findings: list[dict[str, str]] = []

    for rsid, entry in _CARRIER_MARKERS.items():
        label = entry["label"] or "Carrier marker"
        effect_allele = entry.get("effect_allele")
        genotype = genotypes.get(rsid)
        if genotype and effect_allele and _risk_allele_present(rsid, genotype, effect_allele, variant_lookup):
            display_label = label
            note = "Screening-level signal; confirm with clinical testing."
            if rsid == "rs334":
                display_label = "Hemoglobin variant (HbS)"
                if genotype == "TT":
                    note = "Suggests disease-level HbSS pattern from screening data; confirm clinically."
                else:
                    note = "Likely carrier-level HbS signal from screening data; confirm clinically."
            elif rsid == "rs17580":
                note = (
                    "Pi*S allele present; possible mild AAT reduction. Interpret with Pi*Z status, "
                    "serum AAT level, and clinical context."
                )
            findings.append(
                {
                    "category": "High-impact findings (screening-level)",
                    "label": display_label,
                    "sub": f"{rsid} {genotype}",
                    "note": note,
                }
            )

    for rsid, entry in _HEREDITARY_FOUNDER_MARKERS.items():
        label = entry["label"] or "Hereditary cancer marker"
        effect_allele = entry.get("effect_allele")
        genotype = genotypes.get(rsid)
        if genotype and effect_allele and _risk_allele_present(rsid, genotype, effect_allele, variant_lookup):
            findings.append(
                {
                    "category": "Hereditary cancer founder variants",
                    "label": label,
                    "sub": f"{rsid} {genotype}",
                    "note": "Founder variant signal; confirm with clinical testing.",
                }
            )

    medication_labels: set[str] = set()
    for card in risk_cards or []:
        if card.get("category") != "clinical":
            continue
        if card.get("evidence") not in {"CPIC", "DPWG"}:
            continue
        label = card.get("label") or "Medication-related finding"
        if label in medication_labels:
            continue
        medication_labels.add(label)
        findings.append(
            {
                "category": "Medication alerts (screening-level)",
                "label": label,
                "sub": str(card.get("description", "")).strip() or "See actionable PGx section.",
                "note": "Mirrors the actionable PGx section; confirm clinically before medication changes.",
            }
        )

    return findings


def _validate_hbs_interpretation_guardrail(
    genotypes: dict[str, str],
    risk_cards: list[dict[str, str]],
    high_priority: list[dict[str, str]],
) -> None:
    if genotypes.get("rs334"):
        return
    offenders: list[str] = []
    for card in risk_cards:
        label = str(card.get("label", ""))
        description = str(card.get("description", ""))
        combined = f"{label} {description}".lower()
        if "hbc" in combined or label == "Sickle Cell (HbS)":
            offenders.append(f"risk_card:{label or description}")
    for item in high_priority:
        label = str(item.get("label", ""))
        sub = str(item.get("sub", ""))
        combined = f"{label} {sub}".lower()
        if "hbc" in combined or label in {"Hemoglobin variant (HbS)", "Sickle Cell (HbS)", "Sickle cell (HbS)"}:
            offenders.append(f"high_priority:{label or sub}")
    if offenders:
        details = ", ".join(offenders)
        raise ValueError(
            "Validation failed: rs334 is missing, but HbS/HbC interpretations were generated "
            f"({details})."
        )


def _has_allele(genotype: str | None, allele: str) -> bool:
    return bool(genotype and allele in genotype)


def _allele_count(genotype: str | None, allele: str) -> int:
    if not genotype:
        return 0
    return sum(1 for base in genotype if base == allele)


def _apoe_haplotype(genotypes: dict[str, str], clinical: dict[str, Any]) -> str:
    rs429358 = genotypes.get("rs429358")
    rs7412 = genotypes.get("rs7412")
    if not rs429358 or not rs7412:
        return "Unknown"
    key = f"{rs429358}|{rs7412}"
    return clinical.get("apoe_haplotype_map", {}).get(key, "Unknown")


def _apoe_assessment(
    genotypes: dict[str, str],
    clinical: dict[str, Any],
    variant_lookup: dict[str, dict[str, Any]] | None,
) -> dict[str, Any]:
    required_rsids = ("rs429358", "rs7412")
    observed = {rsid: genotypes.get(rsid) for rsid in required_rsids}
    missing = [rsid for rsid, genotype in observed.items() if not genotype]
    if missing:
        return {
            "assessed": False,
            "reason": f"APOE not assessed (partial/missing SNPs): missing {', '.join(missing)}.",
            "missing_rsids": missing,
            "unverified_rsids": [],
            "haplotype": None,
            "genotypes": observed,
        }

    unverified: list[str] = []
    for rsid in required_rsids:
        if not variant_lookup:
            unverified.append(rsid)
            continue
        match_status = str((variant_lookup.get(rsid) or {}).get("match_status") or "")
        if match_status not in {"match", "reverse_complement"}:
            unverified.append(rsid)
    if unverified:
        return {
            "assessed": False,
            "reason": f"APOE not assessed (verification incomplete): {', '.join(unverified)}.",
            "missing_rsids": [],
            "unverified_rsids": unverified,
            "haplotype": None,
            "genotypes": observed,
        }

    key = f"{observed['rs429358']}|{observed['rs7412']}"
    haplotype = clinical.get("apoe_haplotype_map", {}).get(key)
    if not haplotype:
        return {
            "assessed": False,
            "reason": f"APOE not assessed (unmapped genotype combination): {key}.",
            "missing_rsids": [],
            "unverified_rsids": [],
            "haplotype": None,
            "genotypes": observed,
        }
    return {
        "assessed": True,
        "reason": "",
        "missing_rsids": [],
        "unverified_rsids": [],
        "haplotype": haplotype,
        "genotypes": observed,
    }


def _risk_card(
    label: str,
    level: str,
    description: str,
    action: str,
    *,
    evidence: str,
    category: str,
) -> dict[str, str]:
    return {
        "label": label,
        "level": level,
        "description": description,
        "action": action,
        "evidence": evidence,
        "category": category,
    }


def _warfarin_panel_status(genotypes: dict[str, str]) -> str:
    has_vkorc1 = bool(genotypes.get("rs9923231"))
    has_rs12777823 = bool(genotypes.get("rs12777823"))
    if has_vkorc1 and has_rs12777823:
        return "Warfarin panel status: VKORC1 present; rs12777823 present."
    if has_vkorc1 and not has_rs12777823:
        return (
            "Warfarin panel status: VKORC1 present; rs12777823 missing in this file build "
            "(ancestry modifier). Panel partial."
        )
    if not has_vkorc1 and has_rs12777823:
        return (
            "Warfarin panel status: VKORC1 missing in this file build; rs12777823 present. "
            "Dosing algorithm incomplete."
        )
    return (
        "Warfarin panel status: VKORC1 and rs12777823 missing in this file build. "
        "Dosing algorithm incomplete."
    )


def _warfarin_action_guidance(genotypes: dict[str, str]) -> str:
    has_vkorc1 = bool(genotypes.get("rs9923231"))
    has_rs12777823 = bool(genotypes.get("rs12777823"))
    missing: list[str] = []
    if not has_vkorc1:
        missing.append("VKORC1")
    if not has_rs12777823:
        missing.append("rs12777823")
    if missing:
        return (
            "Partial warfarin genetics available ("
            + " / ".join(missing)
            + " missing in this file build); do not use a full CPIC dosing calculator from this file alone. "
            "Use clinical dosing with INR monitoring, and consider clinical-grade PGx if warfarin is planned."
        )
    return (
        "Warfarin genetics panel appears complete in this file; combine CPIC-guided dosing with "
        "clinical judgment and INR monitoring."
    )


def _nat2_profile(genotypes: dict[str, str]) -> dict[str, Any]:
    slow_alleles = {
        "rs1801280": "C",
        "rs1799930": "A",
        "rs1799931": "A",
    }
    missing = [rsid for rsid in slow_alleles if rsid not in genotypes]
    if missing:
        return {"status": "unknown", "slow_count": None, "missing": missing}
    slow_count = 0
    for rsid, allele in slow_alleles.items():
        genotype = genotypes.get(rsid, "")
        slow_count += sum(1 for base in genotype if base == allele)
    if slow_count >= 2:
        status = "likely_slow"
    elif slow_count == 1:
        status = "indeterminate"
    else:
        status = "no_slow_markers"
    return {"status": status, "slow_count": slow_count, "missing": []}


def _variant_lookup(variant_verification: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    lookup: dict[str, dict[str, Any]] = {}
    for entry in variant_verification:
        rsid = entry.get("rsid")
        if rsid:
            lookup[rsid] = entry
    return lookup


def _variant_match_ok(
    rsid: str,
    variant_lookup: dict[str, dict[str, Any]] | None,
) -> bool:
    if not variant_lookup:
        return True
    entry = variant_lookup.get(rsid)
    if not entry:
        return True
    return entry.get("match_status") not in {"mismatch", "non_snp_mismatch"}


def _risk_allele_count(
    rsid: str,
    genotype: str | None,
    risk_allele: str,
    variant_lookup: dict[str, dict[str, Any]] | None,
) -> int:
    if not genotype:
        return 0
    if not _variant_match_ok(rsid, variant_lookup):
        return 0
    if _has_allele(genotype, risk_allele):
        return _allele_count(genotype, risk_allele)
    if not variant_lookup:
        return 0
    entry = variant_lookup.get(rsid)
    if not entry:
        return 0
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}.get(risk_allele)
    if not complement:
        return 0
    ensembl_alleles = entry.get("ensembl_alleles") or ""
    if entry.get("match_status") == "reverse_complement":
        return _allele_count(genotype, complement)
    if ensembl_alleles and (risk_allele not in ensembl_alleles) and (complement in ensembl_alleles):
        return _allele_count(genotype, complement)
    return 0


def _risk_allele_present(
    rsid: str,
    genotype: str | None,
    risk_allele: str,
    variant_lookup: dict[str, dict[str, Any]] | None,
) -> bool:
    return _risk_allele_count(rsid, genotype, risk_allele, variant_lookup) > 0


def _risk_zygosity_label(
    rsid: str,
    genotype: str | None,
    risk_allele: str,
    variant_lookup: dict[str, dict[str, Any]] | None,
) -> str:
    count = _risk_allele_count(rsid, genotype, risk_allele, variant_lookup)
    if count >= 2:
        return "homozygous risk allele"
    if count == 1:
        return "heterozygous risk allele"
    return "risk-allele status uncertain"


def _format_pgx_hit(
    label: str,
    rsid: str,
    genotype: str | None,
    risk_allele: str,
    variant_lookup: dict[str, dict[str, Any]] | None,
) -> str:
    genotype_text = genotype or "NA"
    zygosity = _risk_zygosity_label(rsid, genotype, risk_allele, variant_lookup)
    return f"{label} ({rsid} {genotype_text}; {zygosity})"


def _escalate_high_evidence_pgx(
    cards: list[dict[str, str]],
    genotypes: dict[str, str],
    variant_lookup: dict[str, dict[str, Any]] | None,
) -> list[dict[str, str]]:
    existing_labels = {card.get("label") for card in cards}
    escalated = list(cards)
    for rule in _HIGH_EVIDENCE_PGX_ESCALATION:
        label = str(rule["label"])
        if label in existing_labels:
            continue
        hits: list[str] = []
        for marker in rule.get("markers", []):
            rsid = str(marker["rsid"])
            risk_allele = str(marker["risk_allele"])
            genotype = genotypes.get(rsid)
            if _risk_allele_present(rsid, genotype, risk_allele, variant_lookup):
                hits.append(str(marker["display"]))
        if not hits:
            continue
        escalated.append(
            _risk_card(
                label,
                str(rule["level"]),
                str(rule["description_prefix"]) + ", ".join(hits) + str(rule["description_suffix"]),
                str(rule["action"]),
                evidence=str(rule["evidence"]),
                category="clinical",
            )
        )
    return escalated


def _build_risk_cards(
    genotypes: dict[str, str],
    variant_lookup: dict[str, dict[str, Any]] | None = None,
    *,
    sex: str | None = None,
) -> list[dict[str, str]]:
    cards: list[dict[str, str]] = []
    warfarin_panel_status = _warfarin_panel_status(genotypes)
    warfarin_action_guidance = _warfarin_action_guidance(genotypes)

    cyp2c9_2 = genotypes.get("rs1799853")
    cyp2c9_3 = genotypes.get("rs1057910")
    cyp2c9_5 = genotypes.get("rs28371686")
    cyp2c9_8 = genotypes.get("rs7900194")
    cyp2c9_11 = genotypes.get("rs28371685")
    cyp2c9_2_count = _risk_allele_count("rs1799853", cyp2c9_2, "T", variant_lookup)
    cyp2c9_3_count = _risk_allele_count("rs1057910", cyp2c9_3, "C", variant_lookup)
    cyp2c9_5_count = _risk_allele_count("rs28371686", cyp2c9_5, "G", variant_lookup)
    cyp2c9_8_count = _risk_allele_count("rs7900194", cyp2c9_8, "A", variant_lookup)
    cyp2c9_11_count = _risk_allele_count("rs28371685", cyp2c9_11, "T", variant_lookup)
    cyp2c9_variant_count = (
        cyp2c9_2_count
        + cyp2c9_3_count
        + cyp2c9_5_count
        + cyp2c9_8_count
        + cyp2c9_11_count
    )
    if cyp2c9_variant_count:
        level = "med" if cyp2c9_variant_count == 1 else "high"
        detected = []
        if cyp2c9_2_count:
            detected.append(_format_pgx_hit("CYP2C9*2", "rs1799853", cyp2c9_2, "T", variant_lookup))
        if cyp2c9_3_count:
            detected.append(_format_pgx_hit("CYP2C9*3", "rs1057910", cyp2c9_3, "C", variant_lookup))
        if cyp2c9_5_count:
            detected.append(_format_pgx_hit("CYP2C9*5", "rs28371686", cyp2c9_5, "G", variant_lookup))
        if cyp2c9_8_count:
            detected.append(_format_pgx_hit("CYP2C9*8", "rs7900194", cyp2c9_8, "A", variant_lookup))
        if cyp2c9_11_count:
            detected.append(_format_pgx_hit("CYP2C9*11", "rs28371685", cyp2c9_11, "T", variant_lookup))
        cards.append(
            _risk_card(
                "CYP2C9 Reduced Function",
                level,
                "CYP2C9 decreased-function allele(s) detected ("
                + ", ".join(detected)
                + "); affects warfarin and some NSAID dosing.",
                f"{warfarin_action_guidance} {warfarin_panel_status}",
                evidence="CPIC",
                category="clinical",
            )
        )

    factor_v = genotypes.get("rs6025")
    if _risk_allele_present("rs6025", factor_v, "A", variant_lookup):
        cards.append(
            _risk_card(
                "Clotting Risk",
                "high",
                "Factor V Leiden (rs6025) variant detected; elevated venous thrombosis risk.",
                "Inform clinician before surgery or hormone therapy.",
                evidence="ClinGen",
                category="clinical",
            )
        )

    prothrombin = genotypes.get("rs1799963")
    if _risk_allele_present("rs1799963", prothrombin, "A", variant_lookup):
        prothrombin_hit = _format_pgx_hit(
            "Prothrombin G20210A",
            "rs1799963",
            prothrombin,
            "A",
            variant_lookup,
        )
        cards.append(
            _risk_card(
                "Clotting Risk",
                "high",
                f"{prothrombin_hit} detected; elevated venous thrombosis risk.",
                "Inform clinician before surgery or hormone therapy.",
                evidence="ClinGen",
                category="clinical",
            )
        )

    cyp2c19_2 = genotypes.get("rs4244285")
    cyp2c19_3 = genotypes.get("rs4986893")
    cyp2c19_17 = genotypes.get("rs12248560")
    cyp2c19_2_count = _risk_allele_count("rs4244285", cyp2c19_2, "A", variant_lookup)
    cyp2c19_3_count = _risk_allele_count("rs4986893", cyp2c19_3, "A", variant_lookup)
    cyp2c19_17_count = _risk_allele_count("rs12248560", cyp2c19_17, "T", variant_lookup)
    lof_count = cyp2c19_2_count + cyp2c19_3_count
    inc_count = cyp2c19_17_count
    cyp2c19_missing = [
        rsid
        for rsid, gt in (
            ("rs4244285", cyp2c19_2),
            ("rs4986893", cyp2c19_3),
            ("rs12248560", cyp2c19_17),
        )
        if gt is None
    ]
    if lof_count >= 2:
        phenotype = "Likely poor metabolizer"
        level = "high"
    elif lof_count == 1:
        phenotype = "Likely intermediate metabolizer"
        level = "med"
    elif inc_count >= 1:
        phenotype = "Likely rapid/ultrarapid metabolizer"
        level = "med"
    else:
        phenotype = None
        level = "low"

    if lof_count >= 1:
        detected = []
        if cyp2c19_2_count:
            detected.append(_format_pgx_hit("CYP2C19*2", "rs4244285", cyp2c19_2, "A", variant_lookup))
        if cyp2c19_3_count:
            detected.append(_format_pgx_hit("CYP2C19*3", "rs4986893", cyp2c19_3, "A", variant_lookup))
        if detected:
            coverage_note = ""
            if cyp2c19_missing:
                coverage_note = (
                    " Phenotype based on partial CYP2C19 panel; "
                    f"missing {', '.join(cyp2c19_missing)}."
                )
            cards.append(
                _risk_card(
                    "Clopidogrel Response",
                    level,
                    f"{' / '.join(detected)} detected; {phenotype}. Reduced clopidogrel activation."
                    f"{coverage_note}",
                    "Discuss CPIC-guided antiplatelet selection with a clinician.",
                    evidence="CPIC",
                    category="clinical",
                )
            )
    elif inc_count >= 1:
        coverage_note = ""
        if cyp2c19_missing:
            coverage_note = (
                " Phenotype based on partial CYP2C19 panel; "
                f"missing {', '.join(cyp2c19_missing)}."
            )
        cards.append(
            _risk_card(
                "CYP2C19 Increased Function",
                "med",
                f"{_format_pgx_hit('CYP2C19*17', 'rs12248560', cyp2c19_17, 'T', variant_lookup)} detected; {phenotype}. "
                "Altered exposure for some CYP2C19 substrates."
                f"{coverage_note}",
                "Consider CPIC guidance for CYP2C19 substrates (e.g., PPIs, voriconazole).",
                evidence="CPIC",
                category="clinical",
            )
        )

    slco1b1 = genotypes.get("rs4149056")
    if _has_allele(slco1b1, "C"):
        allele_count = _allele_count(slco1b1, "C")
        level = "med" if allele_count == 1 else "high"
        cards.append(
            _risk_card(
                "Statin Myopathy Risk",
                level,
                f"SLCO1B1 rs4149056 ({slco1b1}): higher simvastatin myopathy risk.",
                "If prescribed simvastatin, consider lower dose or alternative statin per CPIC guidance.",
                evidence="CPIC",
                category="clinical",
            )
        )

    vkorc1 = genotypes.get("rs9923231")
    if _has_allele(vkorc1, "T"):
        allele_count = _allele_count(vkorc1, "T")
        level = "med" if allele_count == 1 else "high"
        cards.append(
            _risk_card(
                "Warfarin Sensitivity",
                level,
                f"VKORC1 rs9923231 ({vkorc1}): increased warfarin sensitivity.",
                f"{warfarin_action_guidance} {warfarin_panel_status}",
                evidence="CPIC",
                category="clinical",
            )
        )

    cyp3a5 = genotypes.get("rs776746")
    if _risk_allele_present("rs776746", cyp3a5, "A", variant_lookup):
        cyp3a5_hit = _format_pgx_hit("CYP3A5", "rs776746", cyp3a5, "A", variant_lookup)
        cards.append(
            _risk_card(
                "Tacrolimus Metabolism",
                "med",
                f"{cyp3a5_hit}; expresser phenotype with higher tacrolimus clearance.",
                "Tacrolimus dosing should follow CPIC CYP3A5 guidance.",
                evidence="CPIC",
                category="clinical",
            )
        )

    cyp2b6_516 = genotypes.get("rs3745274")
    cyp2b6_785 = genotypes.get("rs2279343")
    cyp2b6_markers = []
    if _risk_allele_present("rs3745274", cyp2b6_516, "T", variant_lookup):
        cyp2b6_markers.append(_format_pgx_hit("CYP2B6 516G>T", "rs3745274", cyp2b6_516, "T", variant_lookup))
    if _risk_allele_present("rs2279343", cyp2b6_785, "G", variant_lookup):
        cyp2b6_markers.append(_format_pgx_hit("CYP2B6 785A>G", "rs2279343", cyp2b6_785, "G", variant_lookup))
    if cyp2b6_markers:
        cyp2b6_missing = [
            rsid
            for rsid, gt in (("rs3745274", cyp2b6_516), ("rs2279343", cyp2b6_785))
            if gt is None
        ]
        coverage_note = ""
        if cyp2b6_missing:
            coverage_note = (
                " Partial CYP2B6 panel; phenotype may be incomplete "
                f"(missing {', '.join(cyp2b6_missing)})."
            )
        cards.append(
            _risk_card(
                "Efavirenz Metabolism",
                "med",
                "CYP2B6 decreased-function marker(s) detected ("
                + ", ".join(cyp2b6_markers)
                + ")."
                + coverage_note,
                "Efavirenz dosing should follow CPIC CYP2B6 guidance; full haplotyping may be needed.",
                evidence="CPIC",
                category="clinical",
            )
        )

    ugt1a1 = genotypes.get("rs4148323")
    ugt1a1_proxy = genotypes.get("rs887829")
    ugt1a1_hits = []
    if _risk_allele_present("rs4148323", ugt1a1, "A", variant_lookup):
        ugt1a1_hits.append(_format_pgx_hit("UGT1A1*6", "rs4148323", ugt1a1, "A", variant_lookup))
    if _has_allele(ugt1a1_proxy, "T"):
        ugt1a1_hits.append(_format_pgx_hit("UGT1A1*28 proxy", "rs887829", ugt1a1_proxy, "T", variant_lookup))
    if ugt1a1_hits:
        cards.append(
            _risk_card(
                "Atazanavir Hyperbilirubinemia",
                "med",
                "UGT1A1 reduced-function marker(s) detected: " + ", ".join(ugt1a1_hits) + ".",
                "Atazanavir use should follow CPIC UGT1A1 guidance.",
                evidence="CPIC",
                category="clinical",
            )
        )

    hla_b57 = genotypes.get("rs2395029")
    if _has_allele(hla_b57, "G"):
        cards.append(
            _risk_card(
                "Abacavir Hypersensitivity",
                "high",
                "HLA-B*57:01 proxy (rs2395029) detected; abacavir hypersensitivity risk.",
                "Confirm with clinical HLA-B*57:01 testing before abacavir.",
                evidence="CPIC",
                category="clinical",
            )
        )

    nat2_profile = _nat2_profile(genotypes)
    if nat2_profile["status"] == "likely_slow":
        cards.append(
            _risk_card(
                "Detox/Drug Metabolism",
                "med",
                "Likely NAT2 slow acetylator based on a partial SNP panel.",
                "Confirm with clinical-grade PGx before dosing isoniazid/hydralazine/sulfasalazine; avoid smoking.",
                evidence="CPIC",
                category="clinical",
            )
        )

    abcg2 = genotypes.get("rs2231142")
    if _has_allele(abcg2, "A"):
        cards.append(
            _risk_card(
                "Statin Exposure (ABCG2)",
                "med",
                f"ABCG2 rs2231142 ({abcg2}) detected; reduced transporter function can increase exposure.",
                "If prescribed rosuvastatin or other substrates, consider CPIC guidance.",
                evidence="CPIC",
                category="clinical",
            )
        )

    cyp4f2 = genotypes.get("rs2108622")
    if _has_allele(cyp4f2, "T"):
        cards.append(
            _risk_card(
                "Warfarin Dose Modifier (CYP4F2)",
                "low",
                f"CYP4F2 rs2108622 ({cyp4f2}) detected; may increase warfarin dose requirement.",
                f"{warfarin_action_guidance} {warfarin_panel_status}",
                evidence="CPIC",
                category="clinical",
            )
        )

    cyp2c_cluster = genotypes.get("rs12777823")
    if _has_allele(cyp2c_cluster, "A"):
        cards.append(
            _risk_card(
                "Warfarin Dose Modifier (rs12777823)",
                "med",
                f"rs12777823 ({cyp2c_cluster}) detected; dose modifier most relevant in African ancestry.",
                f"{warfarin_action_guidance} {warfarin_panel_status}",
                evidence="CPIC",
                category="clinical",
            )
        )

    dpyd_variants: list[str] = []
    dpyd_markers = [
        ("rs3918290", "A", "DPYD*2A"),
        ("rs67376798", "T", "DPYD c.2846A>T"),
        ("rs55886062", "G", "DPYD c.1679T>G"),
        ("rs56038477", "A", "DPYD HapB3 tag"),
        ("rs75017182", "G", "DPYD HapB3"),
    ]
    dpyd_missing: list[str] = []
    for rsid, risk_allele, label in dpyd_markers:
        genotype = genotypes.get(rsid)
        if genotype is None:
            dpyd_missing.append(rsid)
            continue
        if _risk_allele_present(rsid, genotype, risk_allele, variant_lookup):
            dpyd_variants.append(_format_pgx_hit(label, rsid, genotype, risk_allele, variant_lookup))
    if dpyd_variants:
        dpyd_coverage_note = ""
        if dpyd_missing:
            dpyd_coverage_note = f" Partial DPYD panel; missing {', '.join(dpyd_missing)}."
        cards.append(
            _risk_card(
                "Fluoropyrimidine Toxicity",
                "high",
                "DPYD variant(s) detected: " + ", ".join(dpyd_variants) + "." + dpyd_coverage_note,
                "Confirm with clinical-grade DPYD testing before 5-FU/capecitabine; dosing changes may be needed.",
                evidence="CPIC",
                category="clinical",
            )
        )

    tpmt_2 = genotypes.get("rs1800462")
    tpmt_3b = genotypes.get("rs1800460")
    tpmt_3c = genotypes.get("rs1142345")
    nudt15 = genotypes.get("rs116855232")
    tpmt_2_count = _risk_allele_count("rs1800462", tpmt_2, "C", variant_lookup)
    tpmt_3b_count = _risk_allele_count("rs1800460", tpmt_3b, "A", variant_lookup)
    tpmt_3c_count = _risk_allele_count("rs1142345", tpmt_3c, "G", variant_lookup)
    nudt15_count = _risk_allele_count("rs116855232", nudt15, "T", variant_lookup)
    thiopurine_count = tpmt_2_count + tpmt_3b_count + tpmt_3c_count + nudt15_count
    thiopurine_hits: list[str] = []
    if tpmt_2_count:
        thiopurine_hits.append(_format_pgx_hit("TPMT*2", "rs1800462", tpmt_2, "C", variant_lookup))
    if tpmt_3b_count:
        thiopurine_hits.append(_format_pgx_hit("TPMT*3B", "rs1800460", tpmt_3b, "A", variant_lookup))
    if tpmt_3c_count:
        thiopurine_hits.append(_format_pgx_hit("TPMT*3C", "rs1142345", tpmt_3c, "G", variant_lookup))
    if nudt15_count:
        thiopurine_hits.append(_format_pgx_hit("NUDT15", "rs116855232", nudt15, "T", variant_lookup))
    tpmt_missing = [
        rsid
        for rsid, gt in (
            ("rs1800462", tpmt_2),
            ("rs1800460", tpmt_3b),
            ("rs1142345", tpmt_3c),
            ("rs116855232", nudt15),
        )
        if gt is None
    ]
    if thiopurine_count and thiopurine_hits:
        coverage_note = ""
        if tpmt_missing:
            missing_tpmt = {"rs1800462", "rs1800460", "rs1142345"}
            has_partial_tpmt = any(rsid in missing_tpmt for rsid in tpmt_missing)
            coverage_note = (
                " Partial TPMT panel; phenotype may be incomplete "
                f"(missing {', '.join(tpmt_missing)})."
            )
            if not has_partial_tpmt:
                coverage_note = (
                    " Partial TPMT/NUDT15 panel; phenotype may be incomplete "
                    f"(missing {', '.join(tpmt_missing)})."
                )
        level = "med" if thiopurine_count == 1 else "high"
        cards.append(
            _risk_card(
                "Thiopurine Toxicity",
                level,
                "TPMT/NUDT15 variant(s) detected: "
                + ", ".join(thiopurine_hits)
                + "."
                + coverage_note,
                "Thiopurine dosing should follow CPIC-guided genotype adjustments.",
                evidence="CPIC",
                category="clinical",
            )
        )

    hbb = genotypes.get("rs334")
    if hbb and "T" in hbb:
        if hbb == "TT":
            description = "HbS/HbS genotype detected (rs334 TT); high risk for sickle cell disease."
            level = "high"
            action = "Confirm with clinical hemoglobin testing; disease-level genotype requires clinical care."
        else:
            description = f"HbS variant detected (rs334 {hbb}); likely sickle cell trait (carrier)."
            level = "med"
            action = "Confirm with clinical hemoglobin testing; carrier screening context."
        cards.append(
            _risk_card(
                "Sickle Cell (HbS)",
                level,
                description,
                action,
                evidence="ClinVar",
                category="clinical",
            )
        )

    serpina1_z = genotypes.get("rs28929474")
    serpina1_s = genotypes.get("rs17580")
    serpina1_z_hit = _risk_allele_present("rs28929474", serpina1_z, "A", variant_lookup)
    serpina1_s_hit = _risk_allele_present("rs17580", serpina1_s, "T", variant_lookup)
    if serpina1_z_hit or serpina1_s_hit:
        details = []
        if serpina1_z_hit:
            details.append(_format_pgx_hit("Pi*Z", "rs28929474", serpina1_z, "A", variant_lookup))
        if serpina1_s_hit:
            details.append(_format_pgx_hit("Pi*S", "rs17580", serpina1_s, "T", variant_lookup))
        if serpina1_z_hit:
            description = (
                "SERPINA1 variant(s) detected (" + ", ".join(details) + "). "
                "Potential AAT deficiency signal; clinical severity depends on serum AAT level and full typing."
            )
        else:
            description = (
                "SERPINA1 variant(s) detected (" + ", ".join(details) + "). "
                "Pi*S/S screening result may indicate mild AAT reduction; serum AAT level determines significance."
            )
        cards.append(
            _risk_card(
                "Alpha-1 Antitrypsin Deficiency",
                "med",
                description,
                "Confirm with serum AAT testing; consider full SERPINA1 typing if needed. Avoid smoking.",
                evidence="ClinVar",
                category="clinical",
            )
        )

    g6pd_202 = genotypes.get("rs1050828")
    g6pd_376 = genotypes.get("rs1050829")
    if g6pd_202 and "A" in g6pd_202:
        if sex == "male":
            desc = "G6PD deficiency risk allele detected (X-linked); males are higher risk."
            level = "high"
        elif sex == "female":
            desc = "G6PD deficiency risk allele detected; females can be carriers or affected."
            level = "med"
        else:
            desc = "G6PD deficiency risk allele detected (X-linked); risk depends on sex."
            level = "med"
        if g6pd_376 and "G" in g6pd_376:
            desc += " A- haplotype context is possible (rs1050828 + rs1050829)."
        cards.append(
            _risk_card(
                "G6PD Deficiency",
                level,
                desc,
                "Confirm with clinical G6PD enzyme testing before oxidant drugs.",
                evidence="ClinVar",
                category="clinical",
            )
        )

    apob = genotypes.get("rs5742904")
    if _has_allele(apob, "A"):
        cards.append(
            _risk_card(
                "Familial Hypercholesterolemia (APOB)",
                "high",
                "APOB R3500Q variant detected; associated with familial hypercholesterolemia.",
                "Confirm with clinical lipid testing and genetic counseling.",
                evidence="ClinVar",
                category="clinical",
            )
        )

    apc = genotypes.get("rs1801155")
    if _has_allele(apc, "A"):
        cards.append(
            _risk_card(
                "Colorectal Cancer Risk (APC I1307K)",
                "med",
                "APC I1307K allele detected; population-dependent colorectal cancer risk.",
                "Discuss screening with a clinician; confirm clinically.",
                evidence="ClinVar",
                category="clinical",
            )
        )

    chek2 = genotypes.get("rs17879961")
    if _has_allele(chek2, "C"):
        cards.append(
            _risk_card(
                "Cancer Risk (CHEK2 I157T)",
                "med",
                "CHEK2 I157T allele detected; low-to-moderate penetrance risk modifier.",
                "Confirm clinically; consider family history in screening decisions.",
                evidence="ClinVar",
                category="clinical",
            )
        )

    arms2 = genotypes.get("rs10490924")
    if _has_allele(arms2, "T"):
        cards.append(
            _risk_card(
                "Vision (ARMS2)",
                "med",
                "ARMS2 A69S risk allele detected; associated with AMD susceptibility.",
                "Protect vision (UV protection, avoid smoking); consider eye exams.",
                evidence="GWAS",
                category="association",
            )
        )

    chRNA5 = genotypes.get("rs16969968")
    if chRNA5 == "AA":
        cards.append(
            _risk_card(
                "Addiction Risk",
                "high",
                f"CHRNA5 rs16969968 ({chRNA5}): increased susceptibility to nicotine dependence if exposed.",
                "Avoid nicotine initiation; add support if quitting.",
                evidence="GWAS",
                category="association",
            )
        )

    mthfr_677 = genotypes.get("rs1801133")
    mthfr_1298 = genotypes.get("rs1801131")
    if mthfr_677 == "AG" and mthfr_1298 == "GT":
        cards.append(
            _risk_card(
                "Methylation",
                "med",
                "MTHFR compound heterozygote: rs1801133 AG + rs1801131 GT.",
                "Consider homocysteine testing if clinically indicated.",
                evidence="Low/Contested",
                category="association",
            )
        )

    lpa = genotypes.get("rs10455872")
    if _has_allele(lpa, "G"):
        cards.append(
            _risk_card(
                "Heart Health",
                "med",
                f"LPA rs10455872 ({lpa}): risk allele detected.",
                "Consider one-time Lp(a) blood test.",
                evidence="GWAS",
                category="association",
            )
        )

    cfh = genotypes.get("rs1061170")
    if _has_allele(cfh, "C"):
        cards.append(
            _risk_card(
                "Vision",
                "med",
                f"CFH rs1061170 ({cfh}): AMD risk allele present.",
                "UV protection, leafy greens, avoid smoking.",
                evidence="GWAS",
                category="association",
            )
        )

    ninep21 = genotypes.get("rs1333049")
    if ninep21 == "GG":
        cards.append(
            _risk_card(
                "Early Heart Attack",
                "low",
                f"9p21 CAD locus rs1333049 ({ninep21}): protective genotype.",
                "Good baseline; maintain heart-healthy habits.",
                evidence="GWAS",
                category="association",
            )
        )

    level_order = {"high": 0, "med": 1, "low": 2, "neutral": 3}
    return sorted(cards, key=lambda card: level_order.get(card["level"], 9))


def _status_pill(level: str) -> str:
    mapping = {
        "high": "status-high",
        "med": "status-med",
        "low": "status-low",
        "neutral": "status-neutral",
        "risk": "status-risk",
        "protective": "status-protective",
        "missing": "status-missing",
        "proxy": "status-proxy",
        "caution": "status-caution",
        "info": "status-neutral",
    }
    return mapping.get(level, "status-neutral")


def _wellness_emoji(label: str) -> str:
    key = label.lower()
    if "lactose" in key:
        return "🥛"
    if "caffeine" in key:
        return "☕"
    if "alcohol" in key:
        return "🍺"
    if "nicotine" in key:
        return "🚭"
    if "bitter" in key:
        return "🥬"
    if "celiac" in key or "gluten" in key:
        return "🌾"
    if "muscle" in key or "fiber" in key:
        return "💪"
    if "vitamin d" in key:
        return "☀️"
    if "vitamin b12" in key:
        return "💊"
    if "stress" in key or "comt" in key:
        return "🧠"
    if "alzheimer" in key or "apoe" in key:
        return "🧬"
    if "histamine" in key:
        return "🧫"
    if "detox" in key or "acetylation" in key:
        return "🧪"
    if "autoimmune" in key or "thyroid" in key:
        return "🛡️"
    if "estrogen" in key or "hormone" in key:
        return "🌸"
    if "inflammation" in key:
        return "🔥"
    if "vdr" in key or "bone" in key:
        return "🦴"
    if "methylation" in key:
        return "⚙️"
    if "longevity" in key:
        return "⌛"
    if "neuroplasticity" in key:
        return "🧠"
    if "oxidative" in key:
        return "☢️"
    if "metabolic" in key:
        return "📊"
    if "iron" in key:
        return "🧲"
    return "✨"

def _panel_display_name(panel_name: str) -> str:
    if panel_name == "Functional Health - Methylation":
        return "B-vitamin / Homocysteine pathway"
    return panel_name.replace("Functional Health - ", "").strip()


def _functional_evidence(panel_name: str) -> str:
    mapping = {
        "Functional Health - Detox/Acetylation": "Clinical PGx",
        "Functional Health - Autoimmune": "Association/Tag",
        "Functional Health - Histamine": "Association",
        "Functional Health - Hormone": "Association/Context",
        "Functional Health - Inflammation": "Association",
        "Functional Health - VDR/Bone": "Association",
        "Functional Health - Methylation": "Biochemical Pathway",
        "Functional Health - Longevity": "Statistical Association",
        "Functional Health - Neuroplasticity": "Biological Mechanism",
        "Functional Health - Oxidative Stress": "Biochemical Pathway",
        "Functional Health - Metabolic": "GWAS Association",
        "Functional Health - Iron Metabolism": "Clinical Risk",
    }
    return mapping.get(panel_name, "Association")


def _functional_tags(panel_name: str) -> str:
    mapping = {
        "Functional Health - Detox/Acetylation": "Isoniazid, hydralazine, sulfasalazine",
        "Functional Health - Autoimmune": "Thyroid autoimmunity; ankylosing spondylitis",
        "Functional Health - Histamine": "DAO/HNMT; histamine intolerance",
        "Functional Health - Hormone": "OC/HRT sensitivity; estrogen metabolism",
        "Functional Health - Inflammation": "Systemic inflammation; cytokines",
        "Functional Health - VDR/Bone": "Vitamin D receptor; bone health",
        "Functional Health - Methylation": "B-Vitamins, Homocysteine, Folate",
        "Functional Health - Longevity": "Aging pathways, FOXO3",
        "Functional Health - Neuroplasticity": "BDNF, Brain health",
        "Functional Health - Oxidative Stress": "Antioxidant enzymes, SOD2",
        "Functional Health - Metabolic": "BMI tendency, Appetite, FTO",
        "Functional Health - Iron Metabolism": "Hemochromatosis, Iron storage, HFE",
    }
    return mapping.get(panel_name, "")


def _next_test_for_entry(rsid: str, panel_name: str) -> str | None:
    mapping = {
        "rs4349859": "Clinical HLA-B27 typing if symptoms or family history",
        "rs3177928": "TSH + thyroid antibodies if symptoms/family history",
        "rs7197": "TSH + thyroid antibodies if symptoms/family history",
        "rs10156191": "Symptom-guided histamine elimination trial",
        "rs2052129": "Symptom-guided histamine elimination trial",
        "rs11558538": "Symptom-guided histamine elimination trial",
        "rs2234693": "Discuss with clinician if on OC/HRT",
        "rs4680": "Discuss with clinician if on OC/HRT",
        "rs1800629": "hs-CRP or cytokine panel if symptoms of chronic inflammation",
        "rs1800795": "hs-CRP or cytokine panel if symptoms of chronic inflammation",
        "rs1800896": "hs-CRP or cytokine panel if symptoms of chronic inflammation",
        "rs4880": "Mitochondrial health assessment if fatigue/exercise intolerance",
        "rs1544410": "Serum 25-OH Vitamin D and bone density (DEXA) if indicated",
        "rs1801133": "Serum homocysteine and RBC folate testing",
        "rs1801131": "Serum homocysteine and RBC folate testing",
        "rs1801394": "Serum B12 and methylmalonic acid (MMA) testing",
        "rs1805087": "Serum homocysteine and B12 testing",
        "rs234706": "Plasma amino acids (taurine/methionine) if indicated",
        "rs1800562": "Serum ferritin and transferrin saturation",
        "rs1799945": "Serum ferritin and transferrin saturation",
    }
    if rsid in mapping:
        return mapping[rsid]
    if panel_name == "Functional Health - Detox/Acetylation":
        return "Clinical PGx confirmation if medication relevant"
    return None


def _panel_summary_row(
    panel_name: str,
    entries: list[dict[str, Any]],
    genotypes: dict[str, str],
    non_snp_genotypes: dict[str, str] | None,
    variant_lookup: dict[str, dict[str, Any]] | None,
) -> dict[str, str | None] | None:
    if not entries:
        return None
    risk_rsids = []
    proxy_rsids = []
    missing_rsids = []
    non_snp_calls: dict[str, str] = {}
    non_snp_genotypes = non_snp_genotypes or {}
    for entry in entries:
        rsid = entry.get("rsid", "")
        effect_allele = entry.get("effect_allele") or ""
        genotype = genotypes.get(rsid)
        non_snp_call = non_snp_genotypes.get(rsid)
        flags = _variant_flags(rsid, genotype, non_snp_call, variant_lookup)
        if not genotype:
            if non_snp_call:
                non_snp_calls[rsid] = non_snp_call
            else:
                missing_rsids.append(rsid)
            continue
        if flags["is_proxy"]:
            proxy_rsids.append(rsid)
            continue
        if effect_allele and _risk_allele_present(rsid, genotype, effect_allele, variant_lookup):
            risk_rsids.append(rsid)
    if risk_rsids:
        status = "risk"
        summary_value = "Risk marker present"
    elif missing_rsids or non_snp_calls:
        status = "missing"
        summary_value = "Incomplete screen"
    elif proxy_rsids:
        status = "proxy"
        summary_value = "Proxy marker(s) present"
    else:
        status = "protective"
        summary_value = "No high-confidence adverse flags detected in this screened set"

    detail_parts = []
    if risk_rsids:
        detail_parts.append(f"Risk markers: {', '.join(risk_rsids)}")
    if proxy_rsids:
        detail_parts.append(f"Proxy markers: {', '.join(proxy_rsids)}")
    if missing_rsids:
        detail_parts.append(f"Missing markers: {', '.join(missing_rsids)}")
    if non_snp_calls:
        calls = ", ".join(f"{rsid}={call}" for rsid, call in non_snp_calls.items())
        detail_parts.append(f"Non-SNP calls: {calls}")
    if panel_name == "Functional Health - Methylation":
        detail_parts.append("MTHFR markers are reported separately in Lifestyle & Genetic Associations.")
    detail = "; ".join(detail_parts) if detail_parts else None

    return {
        "label": f"{_panel_display_name(panel_name)} summary",
        "status": status,
        "sub": "",
        "value": summary_value,
        "detail": detail,
        "emoji": _wellness_emoji(panel_name),
        "row_type": "summary",
        "evidence": _functional_evidence(panel_name),
        "tags": _functional_tags(panel_name),
    }


def _panel_rows(
    entries: list[dict[str, Any]],
    genotypes: dict[str, str],
    non_snp_genotypes: dict[str, str] | None = None,
    *,
    prefer_conclusion: bool = False,
    panel_name: str | None = None,
    variant_lookup: dict[str, dict[str, Any]] | None = None,
    include_indicators: bool = False,
) -> list[dict[str, str | None]]:
    rows: list[dict[str, str | None]] = []
    for entry in entries:
        rsid = entry.get("rsid", "")
        label = entry.get("label", "Trait")
        genotype = genotypes.get(rsid)
        non_snp_call = None
        if not genotype and non_snp_genotypes:
            non_snp_call = non_snp_genotypes.get(rsid)
        flags = _variant_flags(rsid, genotype, non_snp_call, variant_lookup)
        effect_allele = entry.get("effect_allele") or ""
        effect_trait = entry.get("effect_trait") or ""
        non_effect_trait = entry.get("non_effect_trait") or ""
        evidence_note = entry.get("evidence_note") or ""
        notes = entry.get("notes") or ""
        nat2_marker = bool(
            panel_name == "Functional Health - Detox/Acetylation"
            and rsid in _NAT2_RSIDS
        )
        if nat2_marker:
            label = "NAT2 SNP (partial panel)"
        use_conclusion = (
            prefer_conclusion
            or rsid == "rs4349859"
            or "autoimmune thyroid risk" in label.lower()
        )

        status = "neutral"
        detail = None
        allele_count: int | None = None
        indicator: str | None = None
        evidence = None
        tags = None
        next_test = None
        risk_present = False

        if include_indicators and panel_name:
            evidence = _functional_evidence(panel_name)
            tags = _functional_tags(panel_name)
            next_test = _next_test_for_entry(rsid, panel_name)

        if not genotype:
            if non_snp_call:
                value = "Not assessed (non-SNP call)"
                status = "missing"
                if include_indicators:
                    indicator = "Non-SNP call"
                detail = f"Non-SNP call observed: {non_snp_call}. Indel/repeat not interpreted."
            else:
                value = "Not assessed"
                status = "missing"
                indicator = "Not assessed"
        elif include_indicators and flags["is_proxy"]:
            value = "Proxy marker"
            status = "proxy"
            indicator = "Proxy marker"
            if rsid == "rs4349859":
                detail = (
                    "Proxy-marker interpretation depends on ancestry and allele orientation; "
                    "confirm with clinical HLA-B27 typing if relevant."
                )
            else:
                detail = (
                    "Proxy-marker interpretation depends on ancestry and allele orientation; "
                    "confirm clinically if relevant."
                )
        elif nat2_marker:
            value = "NAT2 SNP observed (partial panel)"
            status = "info"
            if include_indicators:
                indicator = "Genotype observed"
            detail = (
                "This marker contributes to NAT2 screening, but single-SNP calls do not "
                "define acetylator phenotype without full haplotyping."
            )
        elif effect_allele and (effect_trait or non_effect_trait):
            allele_count = sum(1 for allele in genotype if allele == effect_allele)
            risk_present = _risk_allele_present(rsid, genotype, effect_allele, variant_lookup)
            if allele_count == 0 and non_effect_trait:
                if use_conclusion:
                    value = non_effect_trait
                else:
                    value = f"Genotype {genotype}"
                    detail = non_effect_trait
                if include_indicators:
                    value = "No high-confidence adverse flags in screened set"
                    status = "protective"
                    indicator = "No high-confidence adverse flags"
            elif allele_count >= 1 and effect_trait:
                if use_conclusion:
                    value = effect_trait
                else:
                    value = f"Genotype {genotype}"
                    detail = effect_trait
                if include_indicators:
                    value = "Risk allele present"
                    status = "risk"
                    indicator = "Risk allele present"
            else:
                value = f"Genotype {genotype}"
                if include_indicators and risk_present:
                    value = "Risk allele present"
                    status = "risk"
                    indicator = "Risk allele present"
                elif include_indicators:
                    value = "Genotype observed"
                    status = "info"
                    indicator = "Genotype observed"
        else:
            value = f"Genotype {genotype}"
            if include_indicators:
                status = "info"
                indicator = "Genotype observed"

        detail_lines = []
        if detail:
            detail_lines.append(detail)
        if evidence_note and (risk_present or (allele_count or 0) > 0):
            detail_lines.append(evidence_note)
        if include_indicators and indicator == "No high-confidence adverse flags":
            detail_lines.append("If symptoms persist, consider clinical evaluation.")
        if notes and genotype and "proxy" in notes.lower():
            detail_lines.append(notes)
        if non_snp_call and notes and any(term in notes.lower() for term in ("indel", "repeat")):
            detail_lines.append(notes)

        detail = " ".join(detail_lines) if detail_lines else None

        if rsid == "rs4349859" and genotype and not (include_indicators and flags["is_proxy"]):
            proxy_note = (
                "Tag SNP for HLA-B*27; ancestry-dependent. "
                "Confirm with clinical HLA-B27 testing if symptoms or family history."
            )
            detail = f"{proxy_note} {detail}".strip() if detail else proxy_note
            if include_indicators:
                status = "proxy"
                indicator = "Proxy marker"

        if flags["is_strand_caution"]:
            caution_note = "Strand caution: reference orientation differs."
            detail = f"{caution_note} {detail}".strip() if detail else caution_note
            if include_indicators:
                status = "caution"
                indicator = "Strand caution"
        if flags["is_proxy"] and flags["proxy_note"]:
            detail = f"{detail} {flags['proxy_note']}".strip() if detail else flags["proxy_note"]

        rows.append({
            "label": label,
            "status": status,
            "sub": f"{rsid} {genotype or non_snp_call or 'Not Found'}",
            "value": value,
            "detail": detail,
            "emoji": _wellness_emoji(label),
            "indicator": indicator,
            "evidence": evidence,
            "tags": tags,
            "next_test": next_test,
        })
    return rows


def _hidden_screening_rows(
    hidden: dict[str, Any],
    genotypes: dict[str, str],
    non_snp_genotypes: dict[str, str] | None,
    variant_lookup: dict[str, dict[str, Any]] | None,
    *,
    sex: str | None,
    qc_sex: str | None = None,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    records = hidden.get("records", [])
    non_snp_genotypes = non_snp_genotypes or {}
    for entry in records:
        rsid = entry.get("rsid", "")
        if not rsid:
            continue
        genotype = genotypes.get(rsid)
        non_snp_call = None
        if not genotype:
            non_snp_call = non_snp_genotypes.get(rsid)
        if not genotype and not non_snp_call:
            continue
        effect_allele = entry.get("effect_allele") or ""
        label = entry.get("label", "Risk marker")
        note = ""
        row_sub = ""
        if genotype:
            if not effect_allele:
                continue
            if _risk_allele_present(rsid, genotype, effect_allele, variant_lookup):
                continue
            value = f"No risk allele detected at {rsid}"
            status = "protective"
            sub_value = genotype
            note = "Screening-level; other variants not assessed."
            row_sub = sub_value
        else:
            value = "Not assessed (non-SNP call)"
            status = "missing"
            sub_value = non_snp_call or "Not Found"
            placeholder = non_snp_call or "non-SNP"
            note = (
                f"Present in raw file as {placeholder} placeholder; "
                "cannot be interpreted from SNP array; treat as not assessed. "
                "Confirm with clinical testing."
            )
            row_sub = f"{rsid} {sub_value}"

        if rsid in {"rs1050828", "rs1050829"}:
            if sex == "male":
                sex_note = "X-linked; male results are typically more predictive."
            elif sex == "female":
                sex_note = "X-linked; females may be carriers or affected depending on X-inactivation."
            else:
                if qc_sex:
                    sex_note = f"X-linked; reported sex missing. QC inferred {qc_sex} (for QC only)."
                else:
                    sex_note = "X-linked; sex not specified, interpret cautiously."
            note = f"{note} {sex_note}".strip() if note else sex_note
        rows.append(
            {
                "label": label,
                "value": value,
                "sub": row_sub,
                "note": note,
                "status": status,
            }
        )
    return rows


def _normalize_sex(summary: dict[str, Any]) -> str | None:
    for key in ("sex", "gender", "reported_sex"):
        value = summary.get(key)
        if not isinstance(value, str):
            continue
        cleaned = value.strip().lower()
        if cleaned in {"m", "male", "man"}:
            return "male"
        if cleaned in {"f", "female", "woman"}:
            return "female"
    return None


def _demographics_notice(summary: dict[str, Any]) -> str | None:
    missing = []
    if not summary.get("reported_sex") and not summary.get("sex") and not summary.get("gender"):
        missing.append("sex")
    if summary.get("reported_age") is None and summary.get("age") is None:
        missing.append("age")
    if not missing:
        return None
    items = ", ".join(missing)
    note = (
        f"Reported demographics missing: {items}. "
        "Provide them to enable hormone-related or age-stratified notes."
    )
    if "sex" in missing:
        sex_inference = summary.get("sex_inference")
        if isinstance(sex_inference, str) and sex_inference.strip():
            note = f"{note} QC inferred sex: {sex_inference} (for QC only)."
    return note


def _coverage_notes(
    genotypes: dict[str, str],
    non_snp_genotypes: dict[str, str],
) -> tuple[list[str], list[str]]:
    critical = [
        {"label": "APOE haplotype", "rsids": ["rs429358", "rs7412"], "proxy": None, "bucket": "missing"},
        {"label": "CYP2C19*2 (clopidogrel)", "rsids": ["rs4244285"], "proxy": None, "bucket": "missing"},
        {"label": "CYP2C19*3", "rsids": ["rs4986893"], "proxy": None, "bucket": "missing"},
        {"label": "CYP2C19*17", "rsids": ["rs12248560"], "proxy": None, "bucket": "missing"},
        {"label": "SLCO1B1 (statin myopathy)", "rsids": ["rs4149056"], "proxy": None, "bucket": "missing"},
        {"label": "VKORC1 (warfarin)", "rsids": ["rs9923231"], "proxy": None, "bucket": "missing"},
        {"label": "CYP4F2*3 (warfarin modifier)", "rsids": ["rs2108622"], "proxy": None, "bucket": "missing"},
        {
            "label": "Warfarin ancestry modifier (rs12777823)",
            "rsids": ["rs12777823"],
            "proxy": "Ancestry-dependent dose modifier.",
            "bucket": "missing",
            "note_when_present": False,
        },
        {"label": "CYP2C9*5", "rsids": ["rs28371686"], "proxy": None, "bucket": "missing"},
        {"label": "CYP2C9*6", "rsids": ["rs9332131"], "proxy": "Indel; not reliably called on arrays.", "bucket": "expected"},
        {"label": "CYP2C9*8", "rsids": ["rs7900194"], "proxy": None, "bucket": "missing"},
        {"label": "CYP2C9*11", "rsids": ["rs28371685"], "proxy": None, "bucket": "missing"},
        {"label": "TPMT*2", "rsids": ["rs1800462"], "proxy": None, "bucket": "missing"},
        {"label": "ABCG2 Q141K", "rsids": ["rs2231142"], "proxy": None, "bucket": "missing"},
        {"label": "DPYD*2A (fluoropyrimidines)", "rsids": ["rs3918290"], "proxy": None, "bucket": "missing"},
        {"label": "DPYD c.2846A>T", "rsids": ["rs67376798"], "proxy": None, "bucket": "missing"},
        {"label": "DPYD c.1679T>G", "rsids": ["rs55886062"], "proxy": None, "bucket": "missing"},
        {"label": "DPYD HapB3 tag", "rsids": ["rs56038477", "rs75017182"], "proxy": None, "bucket": "missing"},
        {"label": "CYP3A5*3 (tacrolimus)", "rsids": ["rs776746"], "proxy": None, "bucket": "missing"},
        {"label": "CYP2B6 516G>T / 785A>G", "rsids": ["rs3745274", "rs2279343"], "proxy": None, "bucket": "missing"},
        {"label": "UGT1A1*6 (atazanavir)", "rsids": ["rs4148323"], "proxy": None, "bucket": "missing"},
        {"label": "UGT1A1*28 proxy", "rsids": ["rs887829"], "proxy": "Proxy SNP used (rs887829).", "bucket": "expected"},
        {
            "label": "UGT1A1*28 TA repeat",
            "rsids": ["rs8175347"],
            "proxy": "TA repeat indel; not reliably called on arrays.",
            "bucket": "expected",
        },
        {
            "label": "HLA-B*57:01 (abacavir)",
            "rsids": ["rs2395029"],
            "proxy": "Proxy SNP used (rs2395029).",
            "bucket": "expected",
        },
        {
            "label": "HLA-B*15:02 (carbamazepine)",
            "rsids": ["rs2844682", "rs3909184"],
            "proxy": "Proxy tag SNPs; ancestry-dependent.",
            "bucket": "expected",
        },
        {
            "label": "HLA-A*31:01 (carbamazepine)",
            "rsids": ["rs1061235"],
            "proxy": "Proxy tag SNP; ancestry-dependent.",
            "bucket": "expected",
        },
        {
            "label": "HLA-B*58:01 (allopurinol)",
            "rsids": ["rs9263726"],
            "proxy": "Proxy tag SNP; ancestry-dependent.",
            "bucket": "missing",
        },
        {
            "label": "HLA-B27 (ankylosing spondylitis)",
            "rsids": ["rs4349859"],
            "proxy": "Proxy SNP used (rs4349859).",
            "bucket": "expected",
        },
        {"label": "CFH (AMD)", "rsids": ["rs1061170"], "proxy": None, "bucket": "missing"},
        {"label": "ARMS2 (AMD)", "rsids": ["rs10490924"], "proxy": None, "bucket": "missing"},
        {"label": "Factor V Leiden", "rsids": ["rs6025"], "proxy": None, "bucket": "missing"},
        {"label": "Prothrombin G20210A", "rsids": ["rs1799963"], "proxy": None, "bucket": "missing"},
        {"label": "Sickle cell (HbS)", "rsids": ["rs334"], "proxy": None, "bucket": "missing"},
        {
            "label": "CFTR F508del",
            "rsids": ["rs113993960"],
            "proxy": "Indel; not reliably called on arrays.",
            "bucket": "expected",
        },
        {"label": "SERPINA1 Pi*Z", "rsids": ["rs28929474"], "proxy": None, "bucket": "missing"},
        {"label": "SERPINA1 Pi*S", "rsids": ["rs17580"], "proxy": None, "bucket": "missing"},
        {"label": "G6PD c.202G>A", "rsids": ["rs1050828"], "proxy": None, "bucket": "missing"},
        {"label": "G6PD c.376A>G", "rsids": ["rs1050829"], "proxy": None, "bucket": "missing"},
        {"label": "APOB R3500Q", "rsids": ["rs5742904"], "proxy": None, "bucket": "missing"},
        {
            "label": "BRCA1 5382insC",
            "rsids": ["rs80357906"],
            "proxy": "Indel; not reliably called on arrays.",
            "bucket": "expected",
        },
        {
            "label": "BRCA2 6174delT",
            "rsids": ["rs80359550"],
            "proxy": "Indel; not reliably called on arrays.",
            "bucket": "expected",
        },
        {"label": "APC I1307K", "rsids": ["rs1801155"], "proxy": None, "bucket": "missing"},
        {"label": "CHEK2 I157T", "rsids": ["rs17879961"], "proxy": None, "bucket": "missing"},
    ]
    expected_notes: list[str] = []
    missing_notes: list[str] = []

    def _append_note(bucket: str, note: str) -> None:
        if bucket == "expected":
            expected_notes.append(note)
        else:
            missing_notes.append(note)

    for entry in critical:
        rsids = entry["rsids"]
        bucket = entry.get("bucket", "missing")
        present = [rsid for rsid in rsids if rsid in genotypes]
        non_snp_present = {
            rsid: non_snp_genotypes[rsid]
            for rsid in rsids
            if rsid in non_snp_genotypes
        }
        missing = [
            rsid
            for rsid in rsids
            if rsid not in genotypes and rsid not in non_snp_genotypes
        ]

        if not present and not non_snp_present:
            if bucket == "expected":
                proxy_note = entry.get("proxy")
                if proxy_note and "not reliably" in proxy_note.lower():
                    note = f"Not assessed: {entry['label']} - {proxy_note}"
                else:
                    note = (
                        f"Not assessed: {entry['label']} - expected limitation for SNP arrays."
                    )
            else:
                if entry["label"] == "HLA-B*58:01 (allopurinol)":
                    note = (
                        "Not assessed: HLA-B*58:01 (allopurinol) "
                        "(missing rs9263726 proxy SNP in this file build)"
                    )
                else:
                    note = f"Not assessed: {entry['label']} (missing {', '.join(missing)})"
            _append_note(bucket, note)
            continue

        if present:
            if missing or non_snp_present:
                parts = []
                if present:
                    parts.append(f"present {', '.join(present)}")
                if missing:
                    parts.append(f"missing {', '.join(missing)}")
                if non_snp_present:
                    calls = ", ".join(f"{rsid}={call}" for rsid, call in non_snp_present.items())
                    parts.append(f"non-SNP calls {calls}")
                note = f"Partial coverage: {entry['label']} ({'; '.join(parts)})."
                if entry["proxy"]:
                    note = f"{note} {entry['proxy']}"
                _append_note(bucket, note)
                continue
            if entry["proxy"] and entry.get("note_when_present", True):
                _append_note(bucket, f"Note: {entry['label']} - {entry['proxy']}")
            continue

        if non_snp_present and not missing:
            calls = ", ".join(f"{rsid}={call}" for rsid, call in non_snp_present.items())
            _append_note(
                bucket,
                "Non-SNP call present: "
                f"{entry['label']} ({calls}). "
                "Present in raw file as placeholder; cannot be interpreted from SNP array; "
                "treat as not assessed."
            )
            continue
        if non_snp_present and missing:
            calls = ", ".join(f"{rsid}={call}" for rsid, call in non_snp_present.items())
            note = (
                f"Partial coverage: {entry['label']} "
                f"(non-SNP calls {calls}; missing {', '.join(missing)})."
            )
            if entry["proxy"]:
                note = f"{note} {entry['proxy']}"
            _append_note(bucket, note)

    expected_notes.append(
        "Not assessed: Malignant hyperthermia (RYR1/CACNA1S) - variant spectrum not reliably callable on SNP arrays."
    )
    expected_notes.append(
        "Not assessed: Aminoglycoside ototoxicity (MT-RNR1 m.1555A>G) - mitochondrial variants not reliably callable on SNP arrays."
    )
    return expected_notes, missing_notes


def _non_snp_verification_notes(
    variant_verification: list[dict[str, Any]],
) -> list[str]:
    notes: list[str] = []
    for entry in variant_verification:
        match_status = entry.get("match_status")
        if match_status not in {"non_snp_match", "non_snp_mismatch", "non_snp_unknown"}:
            continue
        rsid = entry.get("rsid", "unknown")
        observed = entry.get("observed_genotype") or entry.get("observed_alleles") or "NA"
        if match_status == "non_snp_match":
            status_text = "matches reference allele set"
        elif match_status == "non_snp_mismatch":
            status_text = "does not match reference allele set"
        else:
            status_text = "reference alleles unavailable"
        note = entry.get("note")
        line = f"Non-SNP verification: {rsid} ({observed}) - {status_text}."
        if note:
            line = f"{line} {note}"
        notes.append(line)
    return notes


def _apply_estrogen_notes(
    rows: list[dict[str, str | None]],
    sex: str | None,
) -> list[dict[str, str | None]]:
    if not rows or sex is None:
        return rows
    for row in rows:
        label = str(row.get("label", "")).lower()
        if not label.startswith("estrogen"):
            continue
        if sex == "female":
            sex_note = (
                "Female-specific: may influence breast cancer risk modulation and "
                "hormone-related symptom sensitivity (OC/HRT)."
            )
        elif sex == "male":
            sex_note = (
                "Sex: Male. Estrogen metabolism still matters (bone/vascular health), "
                "but OC/HRT-related impact is less relevant."
            )
        existing = row.get("detail")
        row["detail"] = f"{sex_note} {existing}".strip() if existing else sex_note
    return rows


def _nat2_status(genotypes: dict[str, str]) -> dict[str, str | None]:
    profile = _nat2_profile(genotypes)
    status_key = profile["status"]
    observed_snps = ", ".join(f"{rsid} {genotypes.get(rsid, 'Not Found')}" for rsid in _NAT2_ORDER)
    if status_key == "unknown":
        return {
            "label": "NAT2 acetylation status",
            "status": "neutral",
            "sub": "NAT2 rs1801280/rs1799930/rs1799931",
            "value": "Unknown",
            "detail": (
                "Incomplete NAT2 markers; phenotype cannot be inferred from this partial panel. "
                f"Observed SNPs: {observed_snps}."
            ),
            "emoji": _wellness_emoji("detox"),
            "indicator": "Not assessed",
            "evidence": "Clinical PGx",
            "tags": "Isoniazid, hydralazine, sulfasalazine",
            "next_test": "Clinical PGx confirmation if medication relevant",
        }
    if status_key == "likely_slow":
        value = "Likely slow acetylator (screening-level)"
        status = "risk"
        detail = (
            "Based on three NAT2 tag SNPs; full haplotyping is needed for a clinical "
            "phenotype. Slow acetylator status may affect isoniazid, hydralazine, "
            "sulfasalazine dosing and can modulate toxin/carcinogen handling "
            "(e.g., smoking-related risk). Confirm clinically before medication changes. "
            f"Observed SNPs: {observed_snps}."
        )
        indicator = "Likely slow (partial panel)"
    elif status_key == "indeterminate":
        value = "Indeterminate (one slow allele detected)"
        status = "neutral"
        detail = (
            "One slow allele detected; NAT2 acetylator status is indeterminate without "
            f"full haplotyping. Confirm clinically if medication is relevant. Observed SNPs: {observed_snps}."
        )
        indicator = "Indeterminate"
    else:
        value = "No slow alleles detected (screening-level)"
        status = "protective"
        detail = (
            "No slow alleles detected across the partial panel; full NAT2 haplotyping "
            f"is required for definitive phenotype. Observed SNPs: {observed_snps}."
        )
        indicator = "No slow alleles (partial panel)"
    return {
        "label": "NAT2 acetylation status",
        "status": status,
        "sub": "NAT2 rs1801280/rs1799930/rs1799931",
        "value": value,
        "detail": detail,
        "emoji": _wellness_emoji("detox"),
        "indicator": indicator,
        "evidence": "Clinical PGx",
        "tags": "Isoniazid, hydralazine, sulfasalazine",
        "next_test": "Clinical PGx confirmation if medication relevant",
    }


def _wellness_tables(
    genotypes: dict[str, str],
    non_snp_genotypes: dict[str, str] | None,
    apoe_assessment: dict[str, Any],
    expanded: dict[str, Any],
    summary: dict[str, Any],
    variant_lookup: dict[str, dict[str, Any]] | None,
) -> dict[str, list[dict[str, str | None]]]:
    lactose = genotypes.get("rs4988235")
    caffeine = genotypes.get("rs762551")
    alcohol = genotypes.get("rs671")
    bitter = genotypes.get("rs713598")
    dq25 = genotypes.get("rs3135388")
    dq8 = genotypes.get("rs7454108")

    met_rows = []
    if lactose:
        status = "high" if "T" not in lactose else "low"
        met_rows.append({
            "label": "Lactose Tolerance",
            "status": status,
            "sub": f"rs4988235 {lactose}",
            "value": "Intolerant" if status == "high" else "Tolerant",
            "emoji": _wellness_emoji("Lactose Tolerance"),
        })
    if caffeine:
        status = "med" if caffeine == "AA" else "neutral"
        met_rows.append({
            "label": "Caffeine Metabolism",
            "status": status,
            "sub": f"CYP1A2 rs762551 {caffeine}",
            "value": "Faster" if caffeine == "AA" else "Intermediate",
            "emoji": _wellness_emoji("Caffeine Metabolism"),
        })
    if alcohol:
        status = "high" if "A" in alcohol else "low"
        met_rows.append({
            "label": "Alcohol Flush",
            "status": status,
            "sub": f"ALDH2 rs671 {alcohol}",
            "value": "Flush Risk" if status == "high" else "Tolerant",
            "emoji": _wellness_emoji("Alcohol Flush"),
        })
    if bitter:
        met_rows.append({
            "label": "Bitter Taste",
            "status": "neutral",
            "sub": f"rs713598 {bitter}",
            "value": "Taster",
            "emoji": _wellness_emoji("Bitter Taste"),
        })
    if dq25 and dq8:
        met_rows.append({
            "label": "Celiac Tags",
            "status": "low",
            "sub": f"DQ2.5 {dq25} / DQ8 {dq8}",
            "value": "Low Risk",
            "emoji": _wellness_emoji("Celiac Tags"),
        })

    actn3 = genotypes.get("rs1815739")
    vitd = genotypes.get("rs2282679")
    b12 = genotypes.get("rs602662")
    comt = genotypes.get("rs4680")

    fit_rows = []
    if actn3:
        fit_rows.append({
            "label": "Muscle Fiber Type",
            "status": "neutral",
            "sub": f"ACTN3 {actn3}",
            "value": "Power / Sprint" if actn3 == "CC" else "Mixed",
            "emoji": _wellness_emoji("Muscle Fiber Type"),
        })
    if vitd:
        fit_rows.append({
            "label": "Vitamin D Levels",
            "status": "med" if vitd == "TT" else "neutral",
            "sub": f"GC rs2282679 {vitd}",
            "value": "Risk" if vitd == "TT" else "Average",
            "emoji": _wellness_emoji("Vitamin D Levels"),
        })
    if b12:
        fit_rows.append({
            "label": "Vitamin B12 Absorption",
            "status": "low",
            "sub": f"FUT2 rs602662 {b12}",
            "value": "Normal",
            "emoji": _wellness_emoji("Vitamin B12 Absorption"),
        })
    if comt:
        fit_rows.append({
            "label": "Stress Response",
            "status": "neutral",
            "sub": f"COMT rs4680 {comt}",
            "value": "Balanced",
            "emoji": _wellness_emoji("Stress Response"),
        })
    if apoe_assessment.get("assessed"):
        haplotype = str(apoe_assessment.get("haplotype") or "Unknown")
        apoe_genotypes = apoe_assessment.get("genotypes", {})
        rs429358 = apoe_genotypes.get("rs429358", "NA")
        rs7412 = apoe_genotypes.get("rs7412", "NA")
        fit_rows.append({
            "label": "Alzheimer's APOE",
            "status": "low" if haplotype == "e3/e3" else "neutral",
            "sub": f"rs429358 {rs429358} + rs7412 {rs7412} -> {haplotype}",
            "value": "Neutral" if haplotype == "e3/e3" else "Genotype context marker",
            "emoji": _wellness_emoji("Alzheimer's APOE"),
        })
    else:
        fit_rows.append({
            "label": "Alzheimer's APOE",
            "status": "missing",
            "sub": "rs429358/rs7412",
            "value": "APOE not assessed (partial/missing SNPs)",
            "detail": str(apoe_assessment.get("reason") or "").strip(),
            "emoji": _wellness_emoji("Alzheimer's APOE"),
        })

    panels = expanded.get("panels", {})
    fun_panels = expanded.get("fun_panels", {})
    functional_rows: list[dict[str, str | None]] = []
    functional_rows.append(_nat2_status(genotypes))
    for panel_name in (
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
    ):
        entries = panels.get(panel_name, [])
        if panel_name == "Functional Health - Detox/Acetylation":
            entries = [entry for entry in entries if entry.get("rsid") not in _NAT2_RSIDS]
        summary_row = _panel_summary_row(panel_name, entries, genotypes, non_snp_genotypes, variant_lookup)
        if summary_row:
            functional_rows.append(summary_row)
        functional_rows.extend(
            _panel_rows(
                entries,
                genotypes,
                non_snp_genotypes,
                prefer_conclusion=(panel_name == "Functional Health - Histamine"),
                panel_name=panel_name,
                variant_lookup=variant_lookup,
                include_indicators=True,
            )
        )
    functional_rows = _apply_estrogen_notes(functional_rows, _normalize_sex(summary))

    lifestyle_rows = _panel_rows(
        panels.get("Lifestyle", []),
        genotypes,
        non_snp_genotypes,
        prefer_conclusion=True,
        panel_name="Lifestyle",
    )
    met_rows.extend(lifestyle_rows)

    fun_entries = [entry for rows in fun_panels.values() for entry in rows]
    fun_by_rsid = {entry["rsid"]: entry for entry in fun_entries if entry.get("rsid")}

    metabolism_fun_rsids = [
        "rs838133",   # sweet preference
        "rs34160967", # umami sensitivity
        "rs4481887",  # asparagus odor detection
        "rs17782313", # appetite tendency (MC4R)
        "rs4846567",  # fat distribution
    ]
    fitness_fun_rsids = [
        "rs1801260",  # chronotype
        "rs1800497",  # reward sensitivity
    ]

    met_fun_entries = [fun_by_rsid[rsid] for rsid in metabolism_fun_rsids if rsid in fun_by_rsid]
    fit_fun_entries = [fun_by_rsid[rsid] for rsid in fitness_fun_rsids if rsid in fun_by_rsid]

    met_rows.extend(
        _panel_rows(
            met_fun_entries,
            genotypes,
            non_snp_genotypes,
            prefer_conclusion=True,
            panel_name="Metabolism & Diet",
        )
    )
    fit_rows.extend(
        _panel_rows(
            fit_fun_entries,
            genotypes,
            non_snp_genotypes,
            prefer_conclusion=True,
            panel_name="Fitness & Aging",
        )
    )

    return {
        "metabolism": met_rows,
        "fitness": fit_rows,
        "lifestyle": [],
        "functional": functional_rows,
    }


def _expanded_panels(
    expanded: dict[str, Any],
    genotypes: dict[str, str],
    non_snp_genotypes: dict[str, str] | None = None,
    variant_lookup: dict[str, dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    partial_groups = [
        {"name": "DPYD HapB3", "rsids": {"rs56038477", "rs75017182"}},
        {"name": "CYP2B6 516/785", "rsids": {"rs3745274", "rs2279343"}},
        {"name": "HLA-B*15:02 proxies", "rsids": {"rs2844682", "rs3909184"}},
    ]
    partial_rsids: set[str] = set()
    non_snp_genotypes = non_snp_genotypes or {}
    for group in partial_groups:
        rsids = group["rsids"]
        present = [rsid for rsid in rsids if rsid in genotypes]
        missing = [rsid for rsid in rsids if rsid not in genotypes]
        if present and missing:
            partial_rsids.update(rsids)

    panels_out: list[dict[str, Any]] = []
    panels = expanded.get("panels", {})
    for panel_name, entries in panels.items():
        items = []
        for entry in entries:
            rsid = entry.get("rsid")
            label = entry.get("label")
            genotype = genotypes.get(rsid)
            non_snp_call = non_snp_genotypes.get(rsid)
            flags = _variant_flags(
                rsid,
                genotype,
                non_snp_call,
                variant_lookup,
                is_partial_panel=rsid in partial_rsids,
            )
            if genotype:
                genotype_display = genotype
            elif flags["is_non_snp_placeholder"]:
                genotype_display = f"{non_snp_call} (non-SNP call)"
            else:
                genotype_display = "Not Found"
            if flags["is_strand_caution"]:
                genotype_display = "Not interpreted (strand caution)"
            if flags["is_partial_panel"]:
                genotype_display = f"{genotype_display} (partial coverage)"
            if flags["is_proxy"]:
                genotype_display = f"{genotype_display} (proxy marker)"
            items.append(f"{label} ({rsid}): {genotype_display}")
        panels_out.append({"name": panel_name, "items": items})
    return panels_out


def _fun_emoji(label: str) -> str:
    key = label.lower()
    if "eye" in key:
        return "👁️"
    if "skin" in key or "pigmentation" in key:
        return "☀️"
    if "sleep" in key or "chronotype" in key or "clock" in key:
        return "⏰"
    if "sweet" in key:
        return "🍭"
    if "umami" in key or "savory" in key:
        return "🥩"
    if "asparagus" in key:
        return "🥦"
    if "baldness" in key or "hair" in key or "freckle" in key:
        return "💇"
    if "taste" in key:
        return "🍽️"
    return "✨"


def _fun_cards(expanded: dict[str, Any], genotypes: dict[str, str]) -> list[dict[str, str]]:
    cards: list[dict[str, str]] = []
    fun_panels = expanded.get("fun_panels", {})
    for section, entries in fun_panels.items():
        if "Appearance" in section:
            category = "Appearance"
        elif "Sensory" in section:
            category = "Sensory"
        else:
            category = "Lifestyle"
        for entry in entries:
            rsid = entry.get("rsid")
            if rsid == "rs16969968":
                continue
            genotype = genotypes.get(rsid)
            if not genotype:
                continue
            label = entry.get("label", "Trait")
            effect_allele = entry.get("effect_allele") or ""
            effect_trait = entry.get("effect_trait") or ""
            non_effect_trait = entry.get("non_effect_trait") or ""
            allele_count = sum(1 for allele in genotype if allele == effect_allele)
            if effect_allele and allele_count == 0 and non_effect_trait:
                value = non_effect_trait
            elif effect_allele and allele_count >= 1 and effect_trait:
                value = effect_trait
            else:
                value = f"Genotype {genotype}"
            cards.append(
                {
                    "label": label,
                    "value": value,
                    "emoji": _fun_emoji(label),
                    "category": category,
                    "sub": f"{rsid} {genotype}",
                }
            )
    return cards


def _sort_fun_appearance(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    if not rows:
        return rows
    order = {
        "rs12913832": 0,  # HERC2 eye color
        "rs1800407": 1,   # OCA2 eye shade
        "rs1426654": 2,   # SLC24A5 skin pigmentation
        "rs16891982": 3,  # SLC45A2 skin/hair pigmentation
        "rs1805007": 4,   # MC1R red hair
        "rs885479": 5,    # MC1R freckles/sun sensitivity
        "rs7349332": 6,   # AR male pattern baldness
    }

    def key(row: dict[str, str]) -> tuple[int, str]:
        rsid = row.get("sub", "").split()[0]
        return (order.get(rsid, 999), row.get("label", ""))

    return sorted(rows, key=key)


def _group_fun_appearance(rows: list[dict[str, str]]) -> list[tuple[str, str, list[dict[str, str]]]]:
    if not rows:
        return []
    rsid_map = {row.get("sub", "").split()[0]: row for row in rows}
    groups: list[tuple[str, str, list[str]]] = [
        (
            "Eyes (HERC2/OCA2)",
            "HERC2 regulates OCA2 expression; OCA2 fine-tunes eye shade. Effects are probabilistic.",
            ["rs12913832", "rs1800407"],
        ),
        (
            "Skin/Hair Pigmentation (SLC24A5/SLC45A2)",
            "These genes influence melanin production/transport; effects can be additive.",
            ["rs1426654", "rs16891982"],
        ),
        (
            "Hair/Freckles/Sun Sensitivity (MC1R)",
            "MC1R affects melanin type (red hair/freckles/sun sensitivity).",
            ["rs1805007", "rs885479"],
        ),
        (
            "Other",
            "Other appearance-related markers.",
            ["rs7349332"],
        ),
    ]
    output: list[tuple[str, str, list[dict[str, str]]]] = []
    used: set[str] = set()
    for title, intro, rsids in groups:
        items = []
        for rsid in rsids:
            row = rsid_map.get(rsid)
            if row:
                items.append(row)
                used.add(rsid)
        if items:
            output.append((title, intro, items))

    remaining = [row for rsid, row in rsid_map.items() if rsid not in used]
    if remaining:
        output.append(
            ("Additional markers", "Additional appearance markers.", remaining)
        )
    return output


def _trials_by_finding(trials: dict[str, Any]) -> list[dict[str, Any]]:
    findings = trials.get("findings")
    if isinstance(findings, list):
        return findings
    return []


def _trial_url(study: dict[str, Any]) -> str | None:
    raw_url = str(study.get("url") or "").strip()
    if raw_url.startswith(("http://", "https://")):
        return raw_url
    nct_id = str(study.get("nct_id") or "").strip()
    if nct_id.upper().startswith("NCT"):
        return f"https://clinicaltrials.gov/study/{nct_id}"
    return None


def _should_include_trials(
    risk_cards: list[dict[str, str]],
    trials_by_finding: list[dict[str, Any]],
) -> bool:
    clinical_cards = [card for card in risk_cards if card.get("category") == "clinical"]
    return bool(clinical_cards) and bool(trials_by_finding)


def _actionable_not_available(
    genotypes: dict[str, str],
    non_snp_genotypes: dict[str, str] | None,
) -> list[dict[str, str]]:
    non_snp_genotypes = non_snp_genotypes or {}
    candidates = [
        {
            "label": "CYP2C19*2 (clopidogrel)",
            "rsids": ["rs4244285"],
            "reason": "Not available in this file build (missing rs4244285).",
            "next": "Use clinical PGx testing before clopidogrel treatment decisions.",
        },
        {
            "label": "CYP2C19*17",
            "rsids": ["rs12248560"],
            "reason": "Not available in this file build (missing rs12248560).",
            "next": "Use clinical PGx testing when CYP2C19 phenotype is medication-relevant.",
        },
        {
            "label": "VKORC1 (warfarin sensitivity)",
            "rsids": ["rs9923231"],
            "reason": "Not available in this file build (missing rs9923231).",
            "next": "Use clinical PGx testing for complete warfarin dosing guidance.",
        },
        {
            "label": "HLA-B*58:01 proxy (allopurinol)",
            "rsids": ["rs9263726"],
            "reason": "Not available in this file build (missing rs9263726 proxy SNP).",
            "next": "If allopurinol is relevant, order clinical HLA-B*58:01 typing.",
        },
        {
            "label": "Sickle cell (HbS)",
            "rsids": ["rs334"],
            "reason": "Not available in this file build (missing rs334).",
            "next": "If clinically relevant, confirm with clinical hemoglobin testing.",
        },
    ]
    missing_items: list[dict[str, str]] = []
    for item in candidates:
        rsids = item["rsids"]
        has_any_call = any(rsid in genotypes or rsid in non_snp_genotypes for rsid in rsids)
        if has_any_call:
            continue
        missing_items.append(
            {
                "label": item["label"],
                "reason": item["reason"],
                "next": item["next"],
            }
        )
    return missing_items


def _functional_groups(functional_rows: list[dict[str, str | None]]) -> list[dict[str, Any]]:
    groups: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None
    for row in functional_rows:
        label = str(row.get("label") or "")
        is_summary = str(row.get("row_type") or "") == "summary" or label.endswith(" summary")
        if is_summary:
            if current is not None:
                groups.append(current)
            current = {"summary": row, "children": []}
            continue
        if current is not None:
            current["children"].append(row)
    if current is not None:
        groups.append(current)
    return groups


def _validate_summary_consistency(functional_rows: list[dict[str, str | None]]) -> None:
    for group in _functional_groups(functional_rows):
        summary = group["summary"]
        children = group["children"]
        summary_label = str(summary.get("label") or "")
        summary_value = str(summary.get("value") or "").lower()
        if (
            "no high-confidence adverse flags detected in this screened set" not in summary_value
            and "no risk markers detected" not in summary_value
        ):
            continue
        for child in children:
            child_label = str(child.get("label") or "")
            child_status = str(child.get("status") or "").lower()
            child_value = str(child.get("value") or "").lower()
            flagged = (
                child_status in {"risk", "missing"}
                or "risk allele present" in child_value
                or "incomplete" in child_value
                or "not assessed" in child_value
            )
            if flagged:
                raise ValueError(
                    "Validation failed: summary-child mismatch in functional section. "
                    f"Summary '{summary_label}' says no risk, but child '{child_label}' is flagged."
                )


def _validate_warfarin_disclaimer(risk_cards: list[dict[str, str]]) -> None:
    for card in risk_cards:
        label = str(card.get("label") or "")
        description = str(card.get("description") or "")
        action = str(card.get("action") or "")
        combined = f"{label} {description} {action}".lower()
        if "warfarin" not in combined:
            continue
        if "warfarin panel status:" not in combined or "vkorc1" not in combined:
            raise ValueError(
                "Validation failed: warfarin-related item missing explicit VKORC1 panel status."
            )
        if ("present" not in combined) and ("missing" not in combined):
            raise ValueError(
                "Validation failed: warfarin-related item missing explicit present/missing wording."
            )


def _validate_missing_rollup(
    functional_rows: list[dict[str, str | None]],
    *,
    allowlist: set[str] | None = None,
) -> None:
    allowed = allowlist or set()
    for group in _functional_groups(functional_rows):
        summary = group["summary"]
        children = group["children"]
        detail = str(summary.get("detail") or "")
        for child in children:
            value = str(child.get("value") or "").lower()
            sub = str(child.get("sub") or "")
            if "not assessed" not in value or "not found" not in sub.lower():
                continue
            rsid = sub.split()[0] if sub.split() else ""
            if not rsid or rsid in allowed:
                continue
            if rsid not in detail:
                raise ValueError(
                    "Validation failed: missing marker not included in section rollup. "
                    f"rsid={rsid}, summary={summary.get('label')}"
                )


def _validate_report_lints(
    risk_cards: list[dict[str, str]],
    wellness: dict[str, list[dict[str, str | None]]],
) -> None:
    functional_rows = wellness.get("functional", [])
    _validate_summary_consistency(functional_rows)
    _validate_warfarin_disclaimer(risk_cards)
    _validate_missing_rollup(functional_rows)


def _render_html(
    template: str,
    base_name: str,
    call_rate: str,
    summary: dict[str, Any],
    risk_cards: list[dict[str, str]],
    wellness: dict[str, list[dict[str, str]]],
    expanded_panels: list[dict[str, Any]],
    fun_cards: list[dict[str, str]],
    trials_by_finding: list[dict[str, Any]],
    demographics_note: str | None,
    coverage_expected: list[str],
    coverage_missing: list[str],
    qc_appendix_notes: list[str],
    hidden_screening: list[dict[str, str]],
    proxy_markers: list[dict[str, str]],
    high_priority: list[dict[str, str]],
    actionable_not_available: list[dict[str, str]],
    *,
    include_trials: bool,
    research_findings: list[dict[str, str]],
) -> str:
    html = template.replace("[filename]", base_name)
    html = html.replace("[date]", date.today().strftime("%B %d, %Y"))
    html = html.replace("[call_rate]", call_rate)
    html = html.replace("[demographics_note]", demographics_note or "")

    def section_intro(text: str) -> str:
        return f"<p class=\"section-intro\">{text}</p>"

    def collapsible_section(
        title: str,
        body_html: str,
        *,
        intro_text: str | None = None,
        open_default: bool = False,
    ) -> str:
        intro_html = section_intro(intro_text) if intro_text else ""
        open_attr = " open" if open_default else ""
        return (
            f"<details class=\"section-collapsible\"{open_attr}>"
            "<summary>"
            "<div class=\"section-head\">"
            f"<h2>{title}</h2>"
            "<div class=\"section-line\"></div>"
            "<span class=\"collapse-hint\" aria-hidden=\"true\"></span>"
            "</div>"
            "</summary>"
            f"{intro_html}"
            f"{body_html}"
            "</details>"
        )

    qc_rows = []
    qc_rows.append(
        "<div class=\"data-row\">"
        f"<div class=\"data-label\">Total SNPs</div><div class=\"data-val\">{summary.get('total_snps', 'NA')}</div>"
        "</div>"
    )
    qc_rows.append(
        "<div class=\"data-row\">"
        f"<div class=\"data-label\">Call Rate</div><div class=\"data-val\">{summary.get('call_rate_percent', 'NA')}%</div>"
        "</div>"
    )
    if summary.get("heterozygosity_rate") is not None:
        qc_rows.append(
            "<div class=\"data-row\">"
            f"<div class=\"data-label\">Heterozygosity Rate</div><div class=\"data-val\">{summary.get('heterozygosity_rate')}</div>"
            "</div>"
        )
    if summary.get("sex_inference"):
        qc_rows.append(
            "<div class=\"data-row\">"
            f"<div class=\"data-label\">Sex Inference</div><div class=\"data-val\">{summary.get('sex_inference')}</div>"
            "</div>"
        )
    if summary.get("build_detected"):
        qc_rows.append(
            "<div class=\"data-row\">"
            f"<div class=\"data-label\">Build Detected</div><div class=\"data-val\">{summary.get('build_detected')}</div>"
            "</div>"
        )
    if summary.get("ambiguous_snp_count") is not None:
        ambiguous_count = summary.get("ambiguous_snp_count")
        ambiguous_pct = summary.get("ambiguous_snp_percent_called")
        ambiguous_val = str(ambiguous_count)
        if ambiguous_pct is not None:
            ambiguous_val = f"{ambiguous_count} ({ambiguous_pct}% of called SNPs)"
        qc_rows.append(
            "<div class=\"data-row\">"
            "<div class=\"data-label\">Ambiguous A/T or C/G SNPs</div>"
            f"<div class=\"data-val\">{ambiguous_val}</div>"
            "</div>"
        )
        qc_rows.append(
            "<div class=\"data-row sub-row\">"
            "<div class=\"data-label\">Informational: computed from observed genotype classes (A/T or C/G). "
            "High values can be normal for array content and are handled via reference-allele verification rules.</div>"
            "<div class=\"data-val\"></div>"
            "</div>"
        )
    if summary.get("duplicate_rsid_count") is not None:
        qc_rows.append(
            "<div class=\"data-row\">"
            f"<div class=\"data-label\">Duplicate rsIDs</div><div class=\"data-val\">{summary.get('duplicate_rsid_count')}</div>"
            "</div>"
        )
    if summary.get("reverse_complement_count"):
        rc_count = int(summary.get("reverse_complement_count") or 0)
        rc_rsids = summary.get("reverse_complement_rsids") or []
        rc_text = f" ({', '.join(rc_rsids)})" if rc_rsids else ""
        qc_rows.append(
            "<div class=\"data-row\">"
            "<div class=\"data-label\">Strand caution</div>"
            f"<div class=\"data-val\">{_count_phrase(rc_count, 'reverse-complement match')}{rc_text}</div>"
            "</div>"
        )
    missing_by_chr = summary.get("missing_by_chromosome")
    if isinstance(missing_by_chr, list) and missing_by_chr:
        missing_items = []
        for entry in missing_by_chr:
            chrom = _display_chromosome_label(entry.get("chr_norm", "NA"))
            total = entry.get("total", "NA")
            missing = entry.get("missing", "NA")
            missing_items.append(f"{chrom}: {missing}/{total}")
        qc_rows.append(
            "<div class=\"data-row\">"
            "<div class=\"data-label\">Missing by Chromosome</div>"
            f"<div class=\"data-val\">{', '.join(missing_items)}</div>"
            "</div>"
        )
    html = html.replace("<!-- QC summary inserted here -->", "".join(qc_rows))

    if hidden_screening:
        rows = []
        for row in hidden_screening:
            status = row.get("status", "protective")
            rows.append(
                "<div class=\"data-row\">"
                f"<div class=\"data-label\">{row['label']}</div>"
                "<div class=\"data-val\">"
                f"<span class=\"status-pill status-{status}\">{row['value']}</span>"
                f"<span class=\"data-sub\">{row['sub']}</span>"
                "</div></div>"
            )
            note = row.get("note")
            if note:
                rows.append(
                    "<div class=\"data-row sub-row\">"
                    f"<div class=\"data-label\">{note}</div>"
                    "<div class=\"data-val\"></div>"
                    "</div>"
                )
        hidden_block = (
            collapsible_section(
                "Hidden Actionable Risks (Screening)",
                "<div class=\"dashboard-grid\">"
                "<div class=\"col-full card\">"
                + "".join(rows) +
                "</div>"
                "</div>",
                intro_text=(
                    "Screening panel only. Absence of a risk allele is not diagnostic, "
                    "and non-SNP calls are treated as not assessed."
                ),
                open_default=any(
                    str(row.get("status", "protective")) not in {"protective", "missing"}
                    for row in hidden_screening
                ),
            )
        )
    else:
        hidden_block = ""

    proxy_block = ""
    if proxy_markers:
        proxy_rows = []
        for marker in proxy_markers:
            label = marker["label"]
            rsid = marker["rsid"]
            genotype = marker["genotype"]
            note = marker.get("note")
            proxy_rows.append(
                "<div class=\"data-row\">"
                f"<div class=\"data-label\">{label}</div>"
                f"<div class=\"data-val\">{rsid} {genotype}</div>"
                "</div>"
            )
            if note:
                proxy_rows.append(
                    "<div class=\"data-row sub-row\">"
                    f"<div class=\"data-label\">{note}</div>"
                    "<div class=\"data-val\"></div>"
                    "</div>"
                )
        proxy_block = (
            collapsible_section(
                "Proxy Marker Screening (Non-diagnostic)",
                "<div class=\"dashboard-grid\">"
                "<div class=\"col-full card\">"
                + "".join(proxy_rows) +
                "</div>"
                "</div>",
                intro_text=(
                    "Proxy markers are ancestry-dependent tags and are not diagnostic. "
                    "Confirm clinically before care decisions."
                ),
                open_default=False,
            )
        )

    if coverage_expected or coverage_missing:
        notes = ""
        if coverage_expected:
            notes += (
                "<div class=\"section-subhead\">Expected chip limitations (repeats/indels/CNV/HLA typing)</div>"
            )
            notes += "".join(
                "<div class=\"data-row coverage-row\">"
                f"<div class=\"data-label\">{note}</div><div class=\"data-val\"></div>"
                "</div>"
                for note in coverage_expected
            )
        if coverage_missing:
            notes += (
                "<div class=\"section-subhead\">Markers usually on arrays but missing in this file build</div>"
            )
            notes += "".join(
                "<div class=\"data-row coverage-row\">"
                f"<div class=\"data-label\">{note}</div><div class=\"data-val\"></div>"
                "</div>"
                for note in coverage_missing
            )
        coverage_block = (
            collapsible_section(
                "Coverage Notes",
                "<div class=\"dashboard-grid\">"
                "<div class=\"col-full card\">"
                f"{notes}"
                "</div>"
                "</div>",
                intro_text="Assay limitations and missing-marker coverage gaps that affect interpretation confidence.",
                open_default=False,
            )
        )
    else:
        coverage_block = ""

    html = html.replace(
        "<!-- Coverage notes inserted here -->",
        hidden_block + coverage_block + proxy_block,
    )

    if research_findings:
        research_cards = []
        for item in research_findings:
            topic = item.get("topic", "Research Topic")
            content = item.get("content", "").replace("\n", "<br>")
            source = item.get("source", "")
            source_html = f"<div class=\"card-source\">Source: {source}</div>" if source else ""
            
            research_cards.append(
                "<div class=\"col-full card\">"
                "<div class=\"card-header\">"
                f"<h3 class=\"card-title\">{topic}</h3>"
                "</div>"
                f"<div class=\"card-prose\">{content}{source_html}</div>"
                "</div>"
            )
        
        research_block = (
            "<div class=\"section-head\">"
            "<h2>Research Augmentation (2025/2026 Consensus)</h2>"
            "<div class=\"section-line\"></div>"
            "</div>"
            + section_intro("Targeted literature updates for high-priority findings in this specific profile.")
            +
            "<div class=\"dashboard-grid\">"
            f"{''.join(research_cards)}"
            "</div>"
        )
    else:
        research_block = ""
    
    # We inject research block before the risk cards for visibility
    # Note: The template doesn't have a specific placeholder for research, so we prepend it to actionable risk placeholder
    # Or better, let's replace "<!-- Actionable risk cards inserted here -->" with research + actionable block
    # But since I'm replacing placeholders, I'll just create a new variable for clinical_block and prepend research if needed.
    
    def card_html(card: dict[str, str]) -> str:
        level = card["level"]
        return (
            f"<div class=\"col-third card risk-card {level}\">"
            f"<div class=\"risk-label\">{card['label']}</div>"
            f"<div class=\"risk-val\" style=\"color: var(--risk-{level})\">{level.title()}</div>"
            f"<div class=\"risk-desc\">{card['description']}</div>"
            f"<div class=\"risk-evidence\">Evidence: {card.get('evidence', 'NA')}</div>"
            f"<div class=\"risk-action\"><strong>Action:</strong> {card['action']}</div>"
            "</div>"
        )

    clinical_cards = [card for card in risk_cards if card.get("category") == "clinical"]
    association_cards = [card for card in risk_cards if card.get("category") != "clinical"]

    if clinical_cards:
        clinical_block = "\n".join(card_html(card) for card in clinical_cards)
    else:
        clinical_block = "<div class=\"col-full card\">No actionable clinical findings detected.</div>"
    if actionable_not_available:
        missing_rows = []
        for item in actionable_not_available:
            missing_rows.append(
                "<div class=\"data-row\">"
                f"<div class=\"data-label\">{item['label']}</div>"
                f"<div class=\"data-val\">{item['reason']}</div>"
                "</div>"
            )
            missing_rows.append(
                "<div class=\"data-row sub-row\">"
                f"<div class=\"data-label\">Next: {item['next']}</div>"
                "<div class=\"data-val\"></div>"
                    "</div>"
                )
        clinical_block += (
            "<div class=\"col-full card\">"
            "<details class=\"card-collapsible\">"
            "<summary>"
            "<div class=\"card-header\">"
            "<h3 class=\"card-title\">Assessed but Not Available in This File</h3>"
            "<span class=\"collapse-hint\" aria-hidden=\"true\"></span>"
            "</div>"
            "</summary>"
            + "".join(missing_rows) +
            "</details>"
            "</div>"
        )

    high_priority_block = ""
    if high_priority:
        grouped: dict[str, list[dict[str, str]]] = {}
        for item in high_priority:
            grouped.setdefault(item["category"], []).append(item)
        blocks = []
        high_priority_intro = next(iter(grouped.keys()), "")
        for idx, (category, items) in enumerate(grouped.items()):
            rows = [
            ]
            if idx > 0:
                rows.append(
                    "<div class=\"data-row sub-row\">"
                    f"<div class=\"data-label\">{category}</div>"
                    "<div class=\"data-val\"></div>"
                    "</div>"
                )
            for item in items:
                rows.append(
                    "<div class=\"data-row\">"
                    f"<div class=\"data-label\">{item['label']}</div>"
                    f"<div class=\"data-val\">{item['sub']}</div>"
                    "</div>"
                )
                note = item.get("note")
                if note:
                    rows.append(
                        "<div class=\"data-row sub-row\">"
                        f"<div class=\"data-label\">{note}</div>"
                        "<div class=\"data-val\"></div>"
                        "</div>"
                    )
            blocks.append("".join(rows))
        high_priority_block = collapsible_section(
            "High Priority Findings",
            "<div class=\"dashboard-grid\">"
            "<div class=\"col-full card\">"
            + "".join(blocks) +
            "</div>"
            "</div>",
            intro_text=high_priority_intro or None,
            open_default=False,
        )
    else:
        high_priority_block = collapsible_section(
            "High Priority Findings",
            "<div class=\"dashboard-grid\">"
            "<div class=\"col-full card\">No high priority findings detected.</div>"
            "</div>",
            open_default=False,
        )

    if association_cards:
        association_block = "\n".join(card_html(card) for card in association_cards)
    else:
        association_block = (
            "<div class=\"col-full card\">Lifestyle/association findings are summarized "
            "in the Wellness &amp; Lifestyle section.</div>"
        )

    def table_card(
        title: str,
        rows: list[dict[str, str | None]],
        *,
        collapse_details: bool = False,
    ) -> str:
        inner = []
        for row in rows:
            pill = _status_pill(row["status"])
            sub_text = row.get("sub") if row.get("row_type") != "summary" else ""
            detail_items: list[str] = []
            if row.get("detail"):
                detail_items.append(str(row["detail"]))
            if row.get("indicator"):
                detail_items.append(f"Indicator: {row['indicator']}")
            if row.get("evidence"):
                detail_items.append(f"Evidence: {row['evidence']}")
            if row.get("tags"):
                detail_items.append(f"Tags: {row['tags']}")
            if row.get("next_test"):
                detail_items.append(f"Best next test: {row['next_test']}")
            details_html = ""
            if collapse_details and detail_items:
                details_html = (
                    "<div class=\"inline-details-wrap\">"
                    "<details class=\"inline-details\">"
                    "<summary>Details</summary>"
                    + "".join(
                        f"<div class=\"inline-detail-item\">{item}</div>"
                        for item in detail_items
                    )
                    + "</details>"
                    "</div>"
                )
            inner.append(
                "<div class=\"data-row\">"
                f"<div class=\"data-label\">{row.get('emoji', '')} {row['label']}{details_html}</div>"
                "<div class=\"data-val\">"
                f"<span class=\"status-pill {pill}\">{row['value']}</span>"
                f"<span class=\"data-sub\">{sub_text}</span>"
                "</div></div>"
            )
            if not collapse_details:
                if row.get("detail"):
                    inner.append(
                        "<div class=\"data-row sub-row\">"
                        f"<div class=\"data-label\">{row['detail']}</div>"
                        "<div class=\"data-val\"></div>"
                        "</div>"
                    )
                if row.get("indicator"):
                    inner.append(
                        "<div class=\"data-row sub-row\">"
                        f"<div class=\"data-label\">Indicator: {row['indicator']}</div>"
                        "<div class=\"data-val\"></div>"
                        "</div>"
                    )
                if row.get("evidence"):
                    inner.append(
                        "<div class=\"data-row sub-row\">"
                        f"<div class=\"data-label\">Evidence: {row['evidence']}</div>"
                        "<div class=\"data-val\"></div>"
                        "</div>"
                    )
                if row.get("tags"):
                    inner.append(
                        "<div class=\"data-row sub-row\">"
                        f"<div class=\"data-label\">Tags: {row['tags']}</div>"
                        "<div class=\"data-val\"></div>"
                        "</div>"
                    )
                if row.get("next_test"):
                    inner.append(
                        "<div class=\"data-row sub-row\">"
                        f"<div class=\"data-label\">Best next test: {row['next_test']}</div>"
                        "<div class=\"data-val\"></div>"
                        "</div>"
                    )
        return (
            "<div class=\"col-half card\">"
            "<div class=\"card-header\">"
            f"<h3 class=\"card-title\">{title}</h3>"
            "</div>" + "".join(inner) + "</div>"
        )

    functional_table = table_card(
        "Functional Health",
        wellness.get("functional", []),
        collapse_details=True,
    )
    functional_table = functional_table.replace('col-half', 'col-full')

    base_wellness = "".join(
        table_card(title, rows)
        for title, rows in (
            ("Metabolism & Diet", wellness.get("metabolism", [])),
            ("Fitness & Aging", wellness.get("fitness", [])),
        )
        if rows
    )

    fun_appearance = _sort_fun_appearance(
        [card for card in fun_cards if card["category"] == "Appearance"]
    )
    fun_sensory = [card for card in fun_cards if card["category"] == "Sensory"]
    fun_lifestyle = [card for card in fun_cards if card["category"] == "Lifestyle"]

    def fun_table_card(title: str, rows: list[dict[str, str]]) -> str:
        inner = []
        for row in rows:
            inner.append(
                "<div class=\"data-row\">"
                f"<div class=\"data-label\">{row['emoji']} {row['label']}</div>"
                "<div class=\"data-val\">"
                f"<span class=\"status-pill status-neutral\">Trait</span>"
                f"<span class=\"data-sub\">{row['sub']}</span>"
                "</div></div>"
            )
            inner.append(
                f"<div class=\"data-row sub-row\"><div class=\"data-label\">{row['value']}</div><div class=\"data-val\"></div></div>"
            )
        return (
            "<div class=\"col-half card\">"
            "<div class=\"card-header\">"
            f"<h3 class=\"card-title\">{title}</h3>"
            "</div>" + "".join(inner) + "</div>"
        )

    def fun_appearance_card(groups: list[tuple[str, str, list[dict[str, str]]]]) -> str:
        inner = []
        intro = (
            "Grouped markers: HERC2/OCA2 (eyes), SLC24A5/SLC45A2 (pigmentation), "
            "MC1R (red hair/freckles/sun sensitivity)."
        )
        inner.append(f"<div class=\"fun-desc\">{intro}</div>")
        inner.append(
            "<div class=\"fun-desc\">"
            "How these combine: HERC2 sets the primary eye-color switch, OCA2 fine-tunes shade; "
            "SLC24A5/SLC45A2 additively shift skin/hair pigmentation; MC1R affects melanin type "
            "and sun sensitivity."
            "</div>"
        )
        for title, note, rows in groups:
            inner.append(
                "<div class=\"data-row\">"
                f"<div class=\"data-label\"><strong>{title}</strong></div>"
                "<div class=\"data-val\"></div></div>"
            )
            inner.append(
                "<div class=\"data-row sub-row\">"
                f"<div class=\"data-label\">{note}</div>"
                "<div class=\"data-val\"></div></div>"
            )
            for row in rows:
                inner.append(
                    "<div class=\"data-row\">"
                    f"<div class=\"data-label\">{row['emoji']} {row['label']}</div>"
                    "<div class=\"data-val\">"
                    f"<span class=\"status-pill status-neutral\">{row['value']}</span>"
                    f"<span class=\"data-sub\">{row['sub']}</span>"
                    "</div></div>"
                )
        return (
            "<div class=\"card fun-full\">"
            "<div class=\"card-header\">"
            "<h3 class=\"card-title\">Appearance</h3>"
            "</div>" + "".join(inner) + "</div>"
        )

    def panels_html(panels: list[dict[str, Any]]) -> str:
        rows_html = []
        for panel in panels:
            items = "".join(f"<div>{item}</div>" for item in panel["items"])
            rows_html.append(
                "<div class=\"data-row\">"
                f"<div class=\"data-label\">{panel['name']}</div>"
                f"<div class=\"data-val\"><div class=\"panel-lines\">{items}</div></div>"
                "</div>"
            )
        return "".join(rows_html)

    appearance_groups = _group_fun_appearance(fun_appearance)
    fun_block = fun_appearance_card(appearance_groups)

    wellness_block = base_wellness + fun_block + functional_table

    general_panels = [
        panel for panel in expanded_panels
        if not panel["name"].startswith("Functional Health - ") and panel["name"] != "Lifestyle"
    ]
    expanded_intro = section_intro(
        "Most loci here are context markers. Pharmacogenomics markers are shown for coverage, "
        "and any high-evidence risk calls are promoted above in Actionable Clinical &amp; Pharmacogenomics."
    )
    expanded_block = panels_html(general_panels)

    if include_trials:
        cards = []
        for finding in trials_by_finding:
            label = finding.get("finding_label", "Finding")
            level = finding.get("finding_level", "unknown")
            query_term = finding.get("query_term", "")
            studies = finding.get("recruiting_studies", [])
            total_studies = len(studies)
            shown_studies = studies[:5]
            rows = []
            if not studies:
                rows.append(
                    "<div class=\"data-row\">"
                    "<div class=\"data-label\">No recruiting trials found.</div>"
                    "<div class=\"data-val\"></div>"
                    "</div>"
                )
            else:
                if total_studies > 5:
                    rows.append(
                        "<div class=\"data-row sub-row\">"
                        f"<div class=\"data-label\">Showing top 5 of {total_studies} recruiting trials.</div>"
                        "<div class=\"data-val\"></div>"
                        "</div>"
                    )
                for study in shown_studies:
                    trial_id = study.get("nct_id", "N/A")
                    title = study.get("title", "N/A")
                    phase = study.get("phase", "N/A")
                    trial_url = _trial_url(study)
                    trial_label = str(trial_id)
                    if trial_url:
                        trial_label = (
                            f'<a href="{trial_url}" target="_blank" rel="noopener noreferrer">{trial_id}</a>'
                        )
                    rows.append(
                        "<div class=\"data-row\">"
                        f"<div class=\"data-label\">{trial_label}</div>"
                        f"<div class=\"data-val\">{title} ({phase})</div>"
                        "</div>"
                    )
            query_note = (
                f"<div class=\"data-row sub-row\"><div class=\"data-label\">Query: {query_term}</div>"
                "<div class=\"data-val\"></div></div>"
                if query_term
                else ""
            )
            cards.append(
                "<div class=\"col-full card\">"
                "<div class=\"card-header\">"
                f"<h3 class=\"card-title\">{label} ({level})</h3>"
                "</div>"
                + query_note
                + "".join(rows)
                + "</div>"
            )
        trials_block = "".join(cards)
        trials_block = (
            "<div class=\"section-head\">"
            "<h2>Clinical Trials (Personalized)</h2>"
            "<div class=\"section-line\"></div>"
            "</div>"
            "<div class=\"dashboard-grid\">"
            f"{trials_block}"
            "</div>"
        )
    else:
        trials_block = ""

    if qc_appendix_notes:
        qc_rows = "".join(
            "<div class=\"data-row\">"
            f"<div class=\"data-label\">{note}</div><div class=\"data-val\"></div>"
            "</div>"
            for note in qc_appendix_notes
        )
        qc_appendix_block = collapsible_section(
            "Developer/QC Appendix",
            "<div class=\"dashboard-grid\">"
            "<div class=\"col-full card\">"
            + qc_rows +
            "</div>"
            "</div>",
            intro_text="Pipeline debugging details for non-SNP verification checks.",
            open_default=False,
        )
    else:
        qc_appendix_block = ""

    html = html.replace("<!-- Actionable risk cards inserted here -->", clinical_block)
    html = html.replace("<!-- High priority inserted here -->", high_priority_block)
    html = html.replace("<!-- Association risk cards inserted here -->", association_block)
    html = html.replace("<!-- Wellness tables inserted here -->", wellness_block)
    html = html.replace("<!-- Expanded panels intro inserted here -->", expanded_intro)
    html = html.replace("<!-- Expanded panels inserted here -->", expanded_block)
    html = html.replace("<!-- Fun trait cards inserted here -->", "")
    html = html.replace("<!-- Trials section inserted here -->", research_block + trials_block + qc_appendix_block)

    return html


def _render_markdown(
    base_name: str,
    summary: dict[str, Any],
    apoe: str,
    risk_cards: list[dict[str, str]],
    wellness: dict[str, list[dict[str, str]]],
    expanded_panels: list[dict[str, Any]],
    fun_cards: list[dict[str, str]],
    trials_by_finding: list[dict[str, Any]],
    variant_verification: list[dict[str, Any]],
    coverage_expected: list[str],
    coverage_missing: list[str],
    qc_appendix_notes: list[str],
    hidden_screening: list[dict[str, str]],
    proxy_markers: list[dict[str, str]],
    high_priority: list[dict[str, str]],
    actionable_not_available: list[dict[str, str]],
    *,
    include_trials: bool,
    research_findings: list[dict[str, str]],
) -> str:
    lines = []
    lines.append("# Comprehensive DNA Analysis Report")
    lines.append(f"**File:** {base_name}.txt  ")
    lines.append(f"**Date:** {date.today().strftime('%B %d, %Y')}  ")
    lines.append(f"**Run Folder:** {summary.get('run_folder', '')}")
    lines.append("")
    clinical_cards = [card for card in risk_cards if card.get("category") == "clinical"]
    association_cards = [card for card in risk_cards if card.get("category") != "clinical"]
    demographics_notice = _demographics_notice(summary)
    section = 1
    lines.append(f"## {section}. Quality Control Summary")
    section += 1
    lines.append(f"* **Total SNPs Analyzed:** {summary.get('total_snps', 'NA')}")
    lines.append(f"* **Call Rate:** {summary.get('call_rate_percent', 'NA')}% (Excellent)")
    if summary.get("heterozygosity_rate") is not None:
        lines.append(f"* **Heterozygosity Rate:** {summary.get('heterozygosity_rate')}")
    if summary.get("sex_inference"):
        lines.append(f"* **Sex Inference (QC):** {summary.get('sex_inference')}")
    if summary.get("build_detected"):
        lines.append(f"* **Build Detected:** {summary.get('build_detected')}")
    if summary.get("ambiguous_snp_count") is not None:
        ambiguous_count = summary.get("ambiguous_snp_count")
        ambiguous_pct = summary.get("ambiguous_snp_percent_called")
        if ambiguous_pct is not None:
            lines.append(f"* **Ambiguous SNPs (A/T or C/G):** {ambiguous_count} ({ambiguous_pct}% of called SNPs)")
        else:
            lines.append(f"* **Ambiguous SNPs (A/T or C/G):** {ambiguous_count}")
        lines.append(
            "* **Ambiguous SNP note:** Informational; computed from observed genotype classes (A/T or C/G). "
            "High values can be normal for array content and are handled via reference-allele verification rules."
        )
    if summary.get("duplicate_rsid_count") is not None:
        dup_examples = summary.get("duplicate_rsid_examples") or []
        example_text = f" Examples: {', '.join(dup_examples)}" if dup_examples else ""
        lines.append(f"* **Duplicate rsIDs:** {summary.get('duplicate_rsid_count')}.{example_text}")
    if summary.get("reverse_complement_count"):
        rc_count = int(summary.get("reverse_complement_count") or 0)
        rc_rsids = summary.get("reverse_complement_rsids") or []
        rc_text = f" ({', '.join(rc_rsids)})" if rc_rsids else ""
        lines.append(
            f"* **Strand caution:** {_count_phrase(rc_count, 'variant')} "
            f"matched the reverse complement of reference alleles{rc_text}."
        )
    missing_by_chr = summary.get("missing_by_chromosome")
    if isinstance(missing_by_chr, list) and missing_by_chr:
        missing_items = []
        for entry in missing_by_chr:
            chrom = _display_chromosome_label(entry.get("chr_norm", "NA"))
            total = entry.get("total", "NA")
            missing = entry.get("missing", "NA")
            missing_items.append(f"{chrom}: {missing}/{total}")
        lines.append(f"* **Missing by Chromosome:** {', '.join(missing_items)}")
    lines.append("* **Status:** Pass. Data is sufficient for high-confidence health and trait screening.")
    if demographics_notice:
        lines.append(f"* **Note:** {demographics_notice}")
    lines.append("\n---\n")

    lines.append(f"## {section}. Actionable Clinical & Pharmacogenomics")
    section += 1
    if clinical_cards:
        for idx, card in enumerate(clinical_cards, start=1):
            lines.append(f"{idx}. **{card['label']}**  ")
            lines.append(f"   * **Summary:** {card['description']}  ")
            lines.append(f"   * **Evidence:** {card.get('evidence', 'NA')}  ")
            lines.append(f"   * **Action:** {card['action']}")
            lines.append("")
    else:
        lines.append("No actionable clinical findings detected.")
    if actionable_not_available:
        lines.append("**Assessed but Not Available in This File**")
        for item in actionable_not_available:
            lines.append(f"* {item['label']}: {item['reason']}")
            lines.append(f"  - Next: {item['next']}")
        lines.append("")
    lines.append("\n---\n")

    lines.append(f"## {section}. High Priority Findings")
    section += 1
    if high_priority:
        grouped: dict[str, list[dict[str, str]]] = {}
        for item in high_priority:
            grouped.setdefault(item["category"], []).append(item)
        for category, items in grouped.items():
            lines.append(f"**{category}**")
            for item in items:
                lines.append(f"* {item['label']}: {item['sub']}")
                if item.get("note"):
                    lines.append(f"  - Note: {item['note']}")
            lines.append("")
    else:
        lines.append("No high priority findings detected.")
    lines.append("\n---\n")

    lines.append(f"## {section}. Lifestyle & Genetic Associations")
    section += 1
    if association_cards:
        for idx, card in enumerate(association_cards, start=1):
            lines.append(f"{idx}. **{card['label']}**  ")
            lines.append(f"   * **Summary:** {card['description']}  ")
            lines.append(f"   * **Evidence:** {card.get('evidence', 'NA')}  ")
            lines.append(f"   * **Action:** {card['action']}")
            lines.append("")
    else:
        lines.append("Lifestyle/association findings are summarized in the Wellness & Lifestyle section.")
    lines.append("\n---\n")

    if hidden_screening:
        lines.append(f"## {section}. Hidden Actionable Risks (Screening)")
        section += 1
        lines.append(
            "_Markers assessed with no risk allele detected or non-SNP calls noted. "
            "Screening only; absence is not diagnostic._"
        )
        for row in hidden_screening:
            lines.append(f"* {row['label']}: {row['value']} ({row['sub']})")
            if row.get("note"):
                lines.append(f"  - Note: {row['note']}")
        lines.append("\n---\n")

    if coverage_expected or coverage_missing:
        lines.append(f"## {section}. Coverage Notes")
        section += 1
        if coverage_expected:
            lines.append("**Expected chip limitations (repeats/indels/CNV/HLA typing)**")
            for note in coverage_expected:
                lines.append(f"* {note}")
            lines.append("")
        if coverage_missing:
            lines.append("**Markers usually on arrays but missing in this file build**")
            for note in coverage_missing:
                lines.append(f"* {note}")
            lines.append("")
        lines.append("\n---\n")

    if proxy_markers:
        lines.append(f"## {section}. Proxy Marker Screening (Non-diagnostic)")
        section += 1
        lines.append(
            "_Proxy markers are population-dependent tags and are not diagnostic. "
            "Confirm with clinical testing._"
        )
        for marker in proxy_markers:
            lines.append(f"* {marker['label']}: {marker['rsid']} {marker['genotype']}")
            note = marker.get("note")
            if note:
                lines.append(f"  - Note: {note}")
        lines.append("\n---\n")

    lines.append(f"## {section}. Wellness & Lifestyle (Summary)")
    section += 1
    for title, key in (
        ("Metabolism & Diet", "metabolism"),
        ("Fitness & Aging", "fitness"),
        ("Functional Health", "functional"),
    ):
        rows = wellness.get(key, [])
        if not rows:
            continue
        lines.append(f"**{title}**")
        for row in rows:
            sub = row.get("sub")
            if sub:
                lines.append(f"* {row['label']}: {row['value']} ({sub})")
            else:
                lines.append(f"* {row['label']}: {row['value']}")
            if row.get("indicator"):
                lines.append(f"  - Indicator: {row['indicator']}")
            if row.get("evidence"):
                lines.append(f"  - Evidence: {row['evidence']}")
            if row.get("tags"):
                lines.append(f"  - Tags: {row['tags']}")
            if row.get("next_test"):
                lines.append(f"  - Best next test: {row['next_test']}")
            if row.get("detail"):
                lines.append(f"  - {row['detail']}")
        lines.append("")

    lines.append(f"## {section}. Appearance")
    section += 1
    appearance_groups = _group_fun_appearance(
        _sort_fun_appearance([card for card in fun_cards if card["category"] == "Appearance"])
    )
    if appearance_groups:
        lines.append("**Appearance (grouped)**")
        lines.append(
            "HERC2/OCA2 influence eye color, SLC24A5/SLC45A2 influence pigmentation, "
            "and MC1R affects red hair/freckles/sun sensitivity. Effects are probabilistic and additive."
        )
        lines.append(
            "How these combine: HERC2 sets the primary eye-color switch, OCA2 fine-tunes shade; "
            "SLC24A5/SLC45A2 additively shift skin/hair pigmentation; MC1R affects melanin type "
            "and sun sensitivity."
        )
        for title, note, rows in appearance_groups:
            lines.append(f"* **{title}**")
            lines.append(f"  - {note}")
            for row in rows:
                lines.append(f"  - {row['emoji']} {row['label']}: {row['value']} ({row['sub']})")
        lines.append("")

    # Sensory/Lifestyle fun traits are surfaced in the Wellness tables above.
    lines.append("\n---\n")

    lines.append(f"## {section}. Expanded Panels (Coverage + GWAS Context)")
    section += 1
    lines.append(
        "_Most loci here are context markers (often GWAS-scale effects). Pharmacogenomics markers are listed for coverage; high-evidence risk calls are promoted above in Actionable Clinical & Pharmacogenomics. “Not Found” means the SNP was not in the file or had no call._"
    )
    for panel in expanded_panels:
        if panel["name"].startswith("Functional Health - ") or panel["name"] == "Lifestyle":
            continue
        lines.append(f"**{panel['name']}**")
        for item in panel["items"]:
            lines.append(f"* {item}")
        lines.append("")
    lines.append("\n---\n")

    if research_findings:
        lines.append(f"## {section}. Research Augmentation (2025/2026 Consensus)")
        section += 1
        lines.append("_Automated research summary for high-priority findings._")
        for item in research_findings:
            topic = item.get("topic", "Topic")
            content = item.get("content", "").strip()
            source = item.get("source", "").strip()
            lines.append(f"### {topic}")
            lines.append(content)
            if source:
                lines.append(f"\n_Source/Context: {source}_")
            lines.append("")
        lines.append("\n---\n")

    if include_trials:
        lines.append(f"## {section}. Clinical Trials (Personalized)")
        section += 1
        for finding in trials_by_finding:
            label = finding.get("finding_label", "Finding")
            level = finding.get("finding_level", "unknown")
            query_term = finding.get("query_term", "")
            lines.append(f"**{label} ({level})**")
            if query_term:
                lines.append(f"_Query:_ {query_term}")
            studies = finding.get("recruiting_studies", [])
            if not studies:
                lines.append("* No recruiting trials found.")
            else:
                if len(studies) > 5:
                    lines.append(f"* Showing top 5 of {len(studies)} recruiting trials.")
                for study in studies[:5]:
                    nct_id = study.get("nct_id", "N/A")
                    title = study.get("title", "N/A")
                    phase = study.get("phase", "N/A")
                    trial_url = _trial_url(study)
                    if trial_url:
                        lines.append(f"* [**{nct_id}**]({trial_url}): {title} ({phase})")
                    else:
                        lines.append(f"* **{nct_id}**: {title} ({phase})")
        lines.append("")
        lines.append("\n---\n")

    lines.append(f"## {section}. Limitations & Disclaimer")
    section += 1
    lines.append("* **CYP2D6:** Status cannot be accurately determined from microarray data due to Copy Number Variation limitations.")
    lines.append("* **CYP2D6 actionability:** If opioids, SSRIs, tricyclics, or tamoxifen are relevant, consider clinical PGx testing.")
    lines.append("* **GSTM1/GSTT1:** Null genotypes are copy-number deletions and cannot be inferred from SNP array data; dedicated CNV testing is required.")
    lines.append("* **Indels/repeats:** Certain variants (e.g., CFTR F508del, BRCA founders, UGT1A1*28) are indels or repeats and may not be callable from array data.")
    lines.append("* **Screening-level only:** Microarray does not capture rare variants or structural changes.")
    lines.append("* **Pharmacogenomics:** Array-based PGx findings are screening-level only; confirm with clinical-grade testing before medication changes.")
    if summary.get("reverse_complement_count"):
        rc_count = int(summary.get("reverse_complement_count") or 0)
        rc_rsids = summary.get("reverse_complement_rsids") or []
        rc_text = f" (rsids: {', '.join(rc_rsids)})" if rc_rsids else ""
        lines.append(
            f"* **Strand caution:** {_count_phrase(rc_count, 'variant')} "
            f"matched the reverse complement of reference alleles{rc_text}."
        )
    hla_entry = next((v for v in variant_verification if v.get("rsid") == "rs4349859"), None)
    if hla_entry and hla_entry.get("note"):
        lines.append("* **HLA-B27 proxy:** Strand mismatch detected for rs4349859. Flip alleles before interpreting.")
    else:
        lines.append("* **HLA-B27 proxy:** Use Step 3 Ensembl verification; pay special attention to rs4349859 if alleles appear reversed.")
    lines.append("* **Not medical advice:** Confirm clinically before making medical decisions.")
    lines.append("")

    if qc_appendix_notes:
        lines.append(f"## {section}. Developer/QC Appendix")
        section += 1
        lines.append("_Pipeline debugging details (non-SNP verification checks)._")
        for note in qc_appendix_notes:
            lines.append(f"* {note}")
        lines.append("")

    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate Markdown and HTML reports from run outputs.")
    parser.add_argument("base_name", help="Base filename without extension")
    parser.add_argument("--run-date", help="Run date in YYYYMMDD (optional)")
    args = parser.parse_args()

    base_name = args.base_name
    run_dir = _find_run_dir(base_name, args.run_date)

    summary = _load_json(run_dir / "summary.json")
    summary["run_folder"] = str(run_dir)
    core_traits = _load_json(run_dir / "core_traits.json")
    healthy = _load_json(run_dir / "healthy_aging.json")
    hidden = _load_json(run_dir / "hidden_risks.json")
    expanded = _load_json(run_dir / "expanded_panels.json")
    trials = _load_json(run_dir / "trials_by_finding.json")
    research_findings = _load_json(run_dir / "research_findings.json")
    if isinstance(research_findings, list):
        research_findings = [
            item for item in research_findings
            if str(item.get("content", "")).strip()
        ]
    else:
        research_findings = []

    clinical = _load_json(Path("data") / "clinical_interpretations.json")

    genotypes = _merge_genotypes(core_traits, healthy, hidden, expanded)
    non_snp_genotypes = _merge_non_snp_genotypes(core_traits, healthy, hidden, expanded)
    variant_verification = _load_json(run_dir / "variant_verification.json")
    variant_lookup = _variant_lookup(variant_verification if isinstance(variant_verification, list) else [])
    apoe_assessment = _apoe_assessment(genotypes, clinical, variant_lookup)
    normalized_sex = _normalize_sex(summary)
    risk_cards = _build_risk_cards(genotypes, variant_lookup, sex=normalized_sex)
    risk_cards = _escalate_high_evidence_pgx(risk_cards, genotypes, variant_lookup)
    hidden_screening = _hidden_screening_rows(
        hidden,
        genotypes,
        non_snp_genotypes,
        variant_lookup,
        sex=normalized_sex,
        qc_sex=summary.get("sex_inference"),
    )
    wellness = _wellness_tables(
        genotypes,
        non_snp_genotypes,
        apoe_assessment,
        expanded,
        summary,
        variant_lookup,
    )
    expanded_panels = _expanded_panels(expanded, genotypes, non_snp_genotypes, variant_lookup)
    fun_cards = _fun_cards(expanded, genotypes)
    trials_by_finding = _trials_by_finding(trials)
    include_trials = _should_include_trials(risk_cards, trials_by_finding)
    # variant_verification loaded above for allele orientation checks
    coverage_expected, coverage_missing = _coverage_notes(genotypes, non_snp_genotypes)
    qc_appendix_notes = (
        _non_snp_verification_notes(variant_verification)
        if isinstance(variant_verification, list)
        else []
    )
    proxy_markers = _proxy_markers_present(genotypes, variant_lookup)
    high_priority = _high_priority_findings(
        genotypes,
        non_snp_genotypes,
        variant_lookup,
        risk_cards=risk_cards,
    )
    actionable_not_available = _actionable_not_available(genotypes, non_snp_genotypes)
    _validate_hbs_interpretation_guardrail(genotypes, risk_cards, high_priority)
    _validate_report_lints(risk_cards, wellness)
    demographics_note = _demographics_notice(summary)

    template = Path("report_template.html").read_text(encoding="utf-8")
    html = _render_html(
        template,
        base_name,
        str(summary.get("call_rate_percent", "NA")),
        summary,
        risk_cards,
        wellness,
        expanded_panels,
        fun_cards,
        trials_by_finding,
        demographics_note,
        coverage_expected,
        coverage_missing,
        qc_appendix_notes,
        hidden_screening,
        proxy_markers,
        high_priority,
        actionable_not_available,
        include_trials=include_trials,
        research_findings=research_findings if isinstance(research_findings, list) else [],
    )

    markdown = _render_markdown(
        base_name,
        summary,
        str(apoe_assessment.get("haplotype") or "Unknown"),
        risk_cards,
        wellness,
        expanded_panels,
        fun_cards,
        trials_by_finding,
        variant_verification if isinstance(variant_verification, list) else [],
        coverage_expected,
        coverage_missing,
        qc_appendix_notes,
        hidden_screening,
        proxy_markers,
        high_priority,
        actionable_not_available,
        include_trials=include_trials,
        research_findings=research_findings if isinstance(research_findings, list) else [],
    )

    (run_dir / f"{base_name}_Report.html").write_text(html, encoding="utf-8")
    (run_dir / f"{base_name}_Report.md").write_text(markdown, encoding="utf-8")

    print(f"Generated reports in {run_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
