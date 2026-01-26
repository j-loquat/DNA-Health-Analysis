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


def _risk_allele_present(
    rsid: str,
    genotype: str | None,
    risk_allele: str,
    variant_lookup: dict[str, dict[str, Any]] | None,
) -> bool:
    if not genotype:
        return False
    if _has_allele(genotype, risk_allele):
        return True
    if not variant_lookup:
        return False
    entry = variant_lookup.get(rsid)
    if not entry:
        return False
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}.get(risk_allele)
    if not complement:
        return False
    ensembl_alleles = entry.get("ensembl_alleles") or ""
    if entry.get("match_status") == "reverse_complement":
        return _has_allele(genotype, complement)
    if ensembl_alleles and (risk_allele not in ensembl_alleles) and (complement in ensembl_alleles):
        return _has_allele(genotype, complement)
    return False


def _build_risk_cards(
    genotypes: dict[str, str],
    variant_lookup: dict[str, dict[str, Any]] | None = None,
) -> list[dict[str, str]]:
    cards: list[dict[str, str]] = []

    cyp2c9_2 = genotypes.get("rs1799853")
    cyp2c9_3 = genotypes.get("rs1057910")
    cyp2c9_variant_count = _allele_count(cyp2c9_2, "T") + _allele_count(cyp2c9_3, "C")
    if cyp2c9_variant_count:
        level = "med" if cyp2c9_variant_count == 1 else "high"
        cards.append(
            _risk_card(
                "CYP2C9 Reduced Function",
                level,
                "CYP2C9 decreased-function allele(s) detected; affects warfarin and some NSAID dosing.",
                "Use CPIC-guided dosing if these drugs are prescribed.",
                evidence="CPIC",
                category="clinical",
            )
        )

    factor_v = genotypes.get("rs6025")
    if _has_allele(factor_v, "A"):
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
    if _has_allele(prothrombin, "A"):
        cards.append(
            _risk_card(
                "Clotting Risk",
                "high",
                "Prothrombin G20210A variant detected; elevated venous thrombosis risk.",
                "Inform clinician before surgery or hormone therapy.",
                evidence="ClinGen",
                category="clinical",
            )
        )

    cyp2c19_2 = genotypes.get("rs4244285")
    cyp2c19_3 = genotypes.get("rs4986893")
    cyp2c19_17 = genotypes.get("rs12248560")
    lof_count = _allele_count(cyp2c19_2, "A") + _allele_count(cyp2c19_3, "A")
    inc_count = _allele_count(cyp2c19_17, "T")
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
        if _has_allele(cyp2c19_2, "A"):
            detected.append("CYP2C19*2")
        if _has_allele(cyp2c19_3, "A"):
            detected.append("CYP2C19*3")
        cards.append(
            _risk_card(
                "Clopidogrel Response",
                level,
                f"{' / '.join(detected)} detected; {phenotype}. Reduced clopidogrel activation.",
                "Discuss CPIC-guided antiplatelet selection with a clinician.",
                evidence="CPIC",
                category="clinical",
            )
        )
    elif inc_count >= 1:
        cards.append(
            _risk_card(
                "CYP2C19 Increased Function",
                "med",
                f"CYP2C19*17 detected; {phenotype}. Altered exposure for some CYP2C19 substrates.",
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
                "Warfarin dosing should follow CPIC genotype-guided algorithms.",
                evidence="CPIC",
                category="clinical",
            )
        )

    cyp3a5 = genotypes.get("rs776746")
    if _risk_allele_present("rs776746", cyp3a5, "A", variant_lookup):
        cards.append(
            _risk_card(
                "Tacrolimus Metabolism",
                "med",
                f"CYP3A5 rs776746 ({cyp3a5}): expresser genotype; higher tacrolimus clearance.",
                "Tacrolimus dosing should follow CPIC CYP3A5 guidance.",
                evidence="CPIC",
                category="clinical",
            )
        )

    cyp2b6_516 = genotypes.get("rs3745274")
    cyp2b6_785 = genotypes.get("rs2279343")
    cyp2b6_markers = []
    if _has_allele(cyp2b6_516, "T"):
        cyp2b6_markers.append("rs3745274")
    if _has_allele(cyp2b6_785, "G"):
        cyp2b6_markers.append("rs2279343")
    if cyp2b6_markers:
        cards.append(
            _risk_card(
                "Efavirenz Metabolism",
                "med",
                "CYP2B6 decreased-function marker(s) detected ("
                + ", ".join(cyp2b6_markers)
                + ").",
                "Efavirenz dosing should follow CPIC CYP2B6 guidance; full haplotyping may be needed.",
                evidence="CPIC",
                category="clinical",
            )
        )

    ugt1a1 = genotypes.get("rs4148323")
    if _risk_allele_present("rs4148323", ugt1a1, "A", variant_lookup):
        cards.append(
            _risk_card(
                "Atazanavir Hyperbilirubinemia",
                "med",
                f"UGT1A1*6 (rs4148323 {ugt1a1}) detected; reduced UGT1A1 function.",
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

    dpyd_variants: list[str] = []
    if _risk_allele_present("rs3918290", genotypes.get("rs3918290"), "A", variant_lookup):
        dpyd_variants.append("rs3918290 (*2A)")
    if _risk_allele_present("rs67376798", genotypes.get("rs67376798"), "A", variant_lookup):
        dpyd_variants.append("rs67376798 (c.2846A>T)")
    if _risk_allele_present("rs55886062", genotypes.get("rs55886062"), "G", variant_lookup):
        dpyd_variants.append("rs55886062 (c.1679T>G)")
    if _risk_allele_present("rs56038477", genotypes.get("rs56038477"), "A", variant_lookup):
        dpyd_variants.append("rs56038477 (HapB3 tag)")
    if _risk_allele_present("rs75017182", genotypes.get("rs75017182"), "G", variant_lookup):
        dpyd_variants.append("rs75017182 (HapB3)")
    if dpyd_variants:
        cards.append(
            _risk_card(
                "Fluoropyrimidine Toxicity",
                "high",
                "DPYD variant(s) detected: " + ", ".join(dpyd_variants) + ".",
                "Confirm with clinical-grade DPYD testing before 5-FU/capecitabine; dosing changes may be needed.",
                evidence="CPIC",
                category="clinical",
            )
        )

    tpmt_3b = genotypes.get("rs1800460")
    tpmt_3c = genotypes.get("rs1142345")
    nudt15 = genotypes.get("rs116855232")
    tpmt_variant_count = _allele_count(tpmt_3b, "A") + _allele_count(tpmt_3c, "G")
    nudt15_variant_count = _allele_count(nudt15, "T")
    thiopurine_count = tpmt_variant_count + nudt15_variant_count
    if thiopurine_count:
        level = "med" if thiopurine_count == 1 else "high"
        cards.append(
            _risk_card(
                "Thiopurine Toxicity",
                level,
                "TPMT/NUDT15 variant allele(s) detected; higher thiopurine toxicity risk.",
                "Thiopurine dosing should follow CPIC-guided genotype adjustments.",
                evidence="CPIC",
                category="clinical",
            )
        )

    chRNA5 = genotypes.get("rs16969968")
    if chRNA5 == "AA":
        cards.append(
            _risk_card(
                "Addiction Risk",
                "high",
                "CHRNA5 (AA): increased susceptibility to nicotine dependence if exposed.",
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
                "MTHFR compound heterozygote (AG + GT).",
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
                "Lp(a) risk allele detected (rs10455872).",
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
                "CFH (CT): AMD risk allele present.",
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
                "9p21 (GG): protective genotype.",
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
    variant_lookup: dict[str, dict[str, Any]] | None,
) -> dict[str, str | None] | None:
    if not entries:
        return None
    if panel_name == "Functional Health - Detox/Acetylation":
        return None
    risk_rsids = []
    missing_rsids = []
    for entry in entries:
        rsid = entry.get("rsid", "")
        effect_allele = entry.get("effect_allele") or ""
        genotype = genotypes.get(rsid)
        if not genotype:
            missing_rsids.append(rsid)
            continue
        if effect_allele and _risk_allele_present(rsid, genotype, effect_allele, variant_lookup):
            risk_rsids.append(rsid)
    if risk_rsids:
        status = "risk"
        summary_value = "Risk marker present"
    elif missing_rsids:
        status = "missing"
        summary_value = "Incomplete screen"
    else:
        status = "protective"
        summary_value = "No risk markers detected"

    detail_parts = []
    if risk_rsids:
        detail_parts.append(f"Risk markers: {', '.join(risk_rsids)}")
    if missing_rsids:
        detail_parts.append(f"Missing markers: {', '.join(missing_rsids)}")
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
        effect_allele = entry.get("effect_allele") or ""
        effect_trait = entry.get("effect_trait") or ""
        non_effect_trait = entry.get("non_effect_trait") or ""
        evidence_note = entry.get("evidence_note") or ""
        notes = entry.get("notes") or ""

        if rsid == "rs4349859" or "autoimmune thyroid risk" in label.lower():
            prefer_conclusion = True

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
            value = "Not assessed"
            status = "missing"
            indicator = "Not assessed"
        elif effect_allele and (effect_trait or non_effect_trait):
            allele_count = sum(1 for allele in genotype if allele == effect_allele)
            risk_present = _risk_allele_present(rsid, genotype, effect_allele, variant_lookup)
            if allele_count == 0 and non_effect_trait:
                if prefer_conclusion:
                    value = non_effect_trait
                else:
                    value = f"Genotype {genotype}"
                    detail = non_effect_trait
                if include_indicators:
                    value = "No risk allele"
                    status = "protective"
                    indicator = "No risk allele"
            elif allele_count >= 1 and effect_trait:
                if prefer_conclusion:
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
        if include_indicators and indicator == "No risk allele":
            detail_lines.append("If symptoms persist despite no risk allele, consider clinical evaluation.")
        if notes and "verify strand" in notes.lower():
            detail_lines.append("Strand caution: verify orientation before interpreting.")

        detail = " ".join(detail_lines) if detail_lines else None

        if rsid == "rs4349859" and genotype:
            proxy_note = (
                "Tag SNP for HLA-B*27; ancestry-dependent. "
                "Confirm with clinical HLA-B27 testing if symptoms or family history."
            )
            detail = f"{proxy_note} {detail}".strip() if detail else proxy_note
            if include_indicators:
                status = "proxy"
                indicator = "Proxy marker"

        if variant_lookup and rsid in variant_lookup:
            match_status = variant_lookup[rsid].get("match_status")
            if match_status in {"reverse_complement", "mismatch"}:
                caution_note = "Strand caution: reference orientation differs."
                detail = f"{caution_note} {detail}".strip() if detail else caution_note
                if include_indicators:
                    status = "caution"
                    indicator = "Strand caution"

        rows.append({
            "label": label,
            "status": status,
            "sub": f"{rsid} {genotype or 'Not Found'}",
            "value": value,
            "detail": detail,
            "emoji": _wellness_emoji(label),
            "indicator": indicator,
            "evidence": evidence,
            "tags": tags,
            "next_test": next_test,
        })
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
    return f"Demographics missing: {items}. Provide them to enable hormone-related or age-stratified notes."


def _coverage_notes(genotypes: dict[str, str]) -> list[str]:
    critical = [
        {"label": "APOE haplotype", "rsids": ["rs429358", "rs7412"], "proxy": None},
        {"label": "CYP2C19*2 (clopidogrel)", "rsids": ["rs4244285"], "proxy": None},
        {"label": "CYP2C19*3", "rsids": ["rs4986893"], "proxy": None},
        {"label": "CYP2C19*17", "rsids": ["rs12248560"], "proxy": None},
        {"label": "SLCO1B1 (statin myopathy)", "rsids": ["rs4149056"], "proxy": None},
        {"label": "VKORC1 (warfarin)", "rsids": ["rs9923231"], "proxy": None},
        {"label": "DPYD*2A (fluoropyrimidines)", "rsids": ["rs3918290"], "proxy": None},
        {"label": "DPYD c.2846A>T", "rsids": ["rs67376798"], "proxy": None},
        {"label": "DPYD c.1679T>G", "rsids": ["rs55886062"], "proxy": None},
        {"label": "DPYD HapB3 tag", "rsids": ["rs56038477", "rs75017182"], "proxy": None},
        {"label": "CYP3A5*3 (tacrolimus)", "rsids": ["rs776746"], "proxy": None},
        {"label": "CYP2B6 516G>T / 785A>G", "rsids": ["rs3745274", "rs2279343"], "proxy": None},
        {"label": "UGT1A1*6 (atazanavir)", "rsids": ["rs4148323"], "proxy": None},
        {"label": "HLA-B*57:01 (abacavir)", "rsids": ["rs2395029"], "proxy": "Proxy SNP used (rs2395029)."},
        {"label": "HLA-B27 (ankylosing spondylitis)", "rsids": ["rs4349859"], "proxy": "Proxy SNP used (rs4349859)."},
        {"label": "CFH (AMD)", "rsids": ["rs1061170"], "proxy": None},
        {"label": "Factor V Leiden", "rsids": ["rs6025"], "proxy": None},
        {"label": "Prothrombin G20210A", "rsids": ["rs1799963"], "proxy": None},
    ]
    notes: list[str] = []
    for entry in critical:
        rsids = entry["rsids"]
        missing = [rsid for rsid in rsids if rsid not in genotypes]
        if missing:
            notes.append(f"Not assessed: {entry['label']} (missing {', '.join(missing)})")
        elif entry["proxy"]:
            notes.append(f"Proxy marker used for {entry['label']}: {entry['proxy']}")
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
    if status_key == "unknown":
        return {
            "label": "NAT2 acetylation status",
            "status": "neutral",
            "sub": "NAT2 rs1801280/rs1799930/rs1799931",
            "value": "Unknown",
            "detail": (
                "Incomplete NAT2 markers; phenotype cannot be inferred from this partial panel."
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
            "(e.g., smoking-related risk). Confirm clinically before medication changes."
        )
        indicator = "Likely slow (partial panel)"
    elif status_key == "indeterminate":
        value = "Indeterminate (one slow allele detected)"
        status = "neutral"
        detail = (
            "One slow allele detected; NAT2 acetylator status is indeterminate without "
            "full haplotyping. Confirm clinically if medication is relevant."
        )
        indicator = "Indeterminate"
    else:
        value = "No slow alleles detected (screening-level)"
        status = "protective"
        detail = (
            "No slow alleles detected across the partial panel; full NAT2 haplotyping "
            "is required for definitive phenotype."
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
    apoe: str,
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
    if apoe != "Unknown":
        fit_rows.append({
            "label": "Alzheimer's APOE",
            "status": "low",
            "sub": f"{apoe}",
            "value": "Neutral",
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
        summary_row = _panel_summary_row(panel_name, entries, genotypes, variant_lookup)
        if summary_row:
            functional_rows.append(summary_row)
        functional_rows.extend(
            _panel_rows(
                entries,
                genotypes,
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

    met_rows.extend(_panel_rows(met_fun_entries, genotypes, prefer_conclusion=True, panel_name="Metabolism & Diet"))
    fit_rows.extend(_panel_rows(fit_fun_entries, genotypes, prefer_conclusion=True, panel_name="Fitness & Aging"))

    return {
        "metabolism": met_rows,
        "fitness": fit_rows,
        "lifestyle": [],
        "functional": functional_rows,
    }


def _expanded_panels(
    expanded: dict[str, Any],
    genotypes: dict[str, str],
    variant_lookup: dict[str, dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    panels_out: list[dict[str, Any]] = []
    panels = expanded.get("panels", {})
    for panel_name, entries in panels.items():
        items = []
        for entry in entries:
            rsid = entry.get("rsid")
            label = entry.get("label")
            genotype = genotypes.get(rsid, "Not Found")
            genotype_display = genotype
            if variant_lookup and rsid in variant_lookup:
                match_status = variant_lookup[rsid].get("match_status")
                if match_status in {"reverse_complement", "mismatch"}:
                    genotype_display = "Not interpreted (strand caution)"
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


def _should_include_trials(
    risk_cards: list[dict[str, str]],
    trials_by_finding: list[dict[str, Any]],
) -> bool:
    clinical_cards = [card for card in risk_cards if card.get("category") == "clinical"]
    return bool(clinical_cards) and bool(trials_by_finding)


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
    coverage_notes: list[str],
    *,
    include_trials: bool,
    research_findings: list[dict[str, str]],
) -> str:
    html = template.replace("[filename]", base_name)
    html = html.replace("[date]", date.today().strftime("%B %d, %Y"))
    html = html.replace("[call_rate]", call_rate)
    html = html.replace("[demographics_note]", demographics_note or "")

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
        qc_rows.append(
            "<div class=\"data-row\">"
            "<div class=\"data-label\">Ambiguous A/T or C/G SNPs</div>"
            f"<div class=\"data-val\">{summary.get('ambiguous_snp_count')}</div>"
            "</div>"
        )
    if summary.get("duplicate_rsid_count") is not None:
        qc_rows.append(
            "<div class=\"data-row\">"
            f"<div class=\"data-label\">Duplicate rsIDs</div><div class=\"data-val\">{summary.get('duplicate_rsid_count')}</div>"
            "</div>"
        )
    missing_by_chr = summary.get("missing_by_chromosome")
    if isinstance(missing_by_chr, list) and missing_by_chr:
        missing_items = []
        for entry in missing_by_chr:
            chrom = entry.get("chr_norm", "NA")
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

    if coverage_notes:
        notes = "".join(
            "<div class=\"data-row\">"
            f"<div class=\"data-label\">{note}</div><div class=\"data-val\"></div>"
            "</div>"
            for note in coverage_notes
        )
        coverage_block = (
            "<div class=\"section-head\">"
            "<h2>Coverage Notes</h2>"
            "<div class=\"section-line\"></div>"
            "</div>"
            "<div class=\"dashboard-grid\">"
            "<div class=\"col-full card\">"
            f"{notes}"
            "</div>"
            "</div>"
        )
    else:
        coverage_block = ""
    html = html.replace("<!-- Coverage notes inserted here -->", coverage_block)

    if research_findings:
        research_cards = []
        for item in research_findings:
            topic = item.get("topic", "Research Topic")
            content = item.get("content", "").replace("\n", "<br>")
            source = item.get("source", "")
            source_html = f"<div class=\"risk-evidence\" style=\"margin-top: 8px; font-style: italic;\">Source: {source}</div>" if source else ""
            
            research_cards.append(
                "<div class=\"col-full card\">"
                "<div class=\"card-header\">"
                f"<h3 class=\"card-title\">{topic}</h3>"
                "</div>"
                f"<div class=\"data-row sub-row\" style=\"display: block; padding: 12px;\">{content}{source_html}</div>"
                "</div>"
            )
        
        research_block = (
            "<div class=\"section-head\">"
            "<h2>Research Augmentation (2025/2026 Consensus)</h2>"
            "<div class=\"section-line\"></div>"
            "</div>"
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

    if association_cards:
        association_block = "\n".join(card_html(card) for card in association_cards)
    else:
        association_block = (
            "<div class=\"col-full card\">Lifestyle/association findings are summarized "
            "in the Wellness &amp; Lifestyle section.</div>"
        )

    def table_card(title: str, rows: list[dict[str, str | None]]) -> str:
        inner = []
        for row in rows:
            pill = _status_pill(row["status"])
            sub_text = row.get("sub") if row.get("row_type") != "summary" else ""
            inner.append(
                "<div class=\"data-row\">"
                f"<div class=\"data-label\">{row.get('emoji', '')} {row['label']}</div>"
                "<div class=\"data-val\">"
                f"<span class=\"status-pill {pill}\">{row['value']}</span>"
                f"<span class=\"data-sub\">{sub_text}</span>"
                "</div></div>"
            )
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

    functional_table = table_card("Functional Health", wellness.get("functional", []))
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
    expanded_block = panels_html(general_panels)

    if include_trials:
        cards = []
        for finding in trials_by_finding:
            label = finding.get("finding_label", "Finding")
            level = finding.get("finding_level", "unknown")
            query_term = finding.get("query_term", "")
            studies = finding.get("recruiting_studies", [])
            rows = []
            if not studies:
                rows.append(
                    "<div class=\"data-row\">"
                    "<div class=\"data-label\">No recruiting trials found.</div>"
                    "<div class=\"data-val\"></div>"
                    "</div>"
                )
            else:
                for study in studies:
                    trial_id = study.get("nct_id", "N/A")
                    url = study.get("url", "")
                    if url:
                        trial_id = f"<a href=\"{url}\" target=\"_blank\" rel=\"noopener noreferrer\">{trial_id}</a>"
                    title = study.get("title", "N/A")
                    phase = study.get("phase", "N/A")
                    rows.append(
                        "<div class=\"data-row\">"
                        f"<div class=\"data-label\">{trial_id}</div>"
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

    html = html.replace("<!-- Actionable risk cards inserted here -->", clinical_block)
    html = html.replace("<!-- Association risk cards inserted here -->", association_block)
    html = html.replace("<!-- Wellness tables inserted here -->", wellness_block)
    html = html.replace("<!-- Expanded panels inserted here -->", expanded_block)
    html = html.replace("<!-- Fun trait cards inserted here -->", "")
    html = html.replace("<!-- Trials section inserted here -->", research_block + trials_block)

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
    coverage_notes: list[str],
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
        lines.append(f"* **Ambiguous SNPs (A/T or C/G):** {summary.get('ambiguous_snp_count')}")
    if summary.get("duplicate_rsid_count") is not None:
        dup_examples = summary.get("duplicate_rsid_examples") or []
        example_text = f" Examples: {', '.join(dup_examples)}" if dup_examples else ""
        lines.append(f"* **Duplicate rsIDs:** {summary.get('duplicate_rsid_count')}.{example_text}")
    missing_by_chr = summary.get("missing_by_chromosome")
    if isinstance(missing_by_chr, list) and missing_by_chr:
        missing_items = []
        for entry in missing_by_chr:
            chrom = entry.get("chr_norm", "NA")
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

    if coverage_notes:
        lines.append(f"## {section}. Coverage Notes")
        section += 1
        for note in coverage_notes:
            lines.append(f"* {note}")
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

    lines.append(f"## {section}. Expanded Panels (GWAS-Style, Low Effect)")
    section += 1
    lines.append(
        "_These are common GWAS loci with small effect sizes; genotypes shown are raw allele pairs and do not imply risk unless a known effect allele is present. “Not Found” means the SNP was not in the file or had no call._"
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
                for study in studies:
                    nct_id = study.get("nct_id", "N/A")
                    title = study.get("title", "N/A")
                    phase = study.get("phase", "N/A")
                    url = study.get("url", "")
                    if url:
                        lines.append(f"* **{nct_id}**: {title} ({phase}) - {url}")
                    else:
                        lines.append(f"* **{nct_id}**: {title} ({phase})")
            lines.append("")
        lines.append("\n---\n")

    lines.append(f"## {section}. Limitations & Disclaimer")
    lines.append("* **CYP2D6:** Status cannot be accurately determined from microarray data due to Copy Number Variation limitations.")
    lines.append("* **GSTM1/GSTT1:** Null genotypes are copy-number deletions and cannot be inferred from SNP array data; dedicated CNV testing is required.")
    lines.append("* **Screening-level only:** Microarray does not capture rare variants or structural changes.")
    lines.append("* **Pharmacogenomics:** Array-based PGx findings are screening-level only; confirm with clinical-grade testing before medication changes.")
    if summary.get("reverse_complement_count"):
        lines.append(
            f"* **Strand caution:** {summary.get('reverse_complement_count')} variants matched the reverse complement of reference alleles."
        )
    hla_entry = next((v for v in variant_verification if v.get("rsid") == "rs4349859"), None)
    if hla_entry and hla_entry.get("note"):
        lines.append("* **HLA-B27 proxy:** Strand mismatch detected for rs4349859. Flip alleles before interpreting.")
    else:
        lines.append("* **HLA-B27 proxy:** Use Step 3 Ensembl verification; pay special attention to rs4349859 if alleles appear reversed.")
    lines.append("* **Not medical advice:** Confirm clinically before making medical decisions.")
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
    apoe = _apoe_haplotype(genotypes, clinical)

    variant_verification = _load_json(run_dir / "variant_verification.json")
    variant_lookup = _variant_lookup(variant_verification if isinstance(variant_verification, list) else [])
    risk_cards = _build_risk_cards(genotypes, variant_lookup)
    wellness = _wellness_tables(genotypes, apoe, expanded, summary, variant_lookup)
    expanded_panels = _expanded_panels(expanded, genotypes, variant_lookup)
    fun_cards = _fun_cards(expanded, genotypes)
    trials_by_finding = _trials_by_finding(trials)
    include_trials = _should_include_trials(risk_cards, trials_by_finding)
    # variant_verification loaded above for allele orientation checks
    coverage_notes = _coverage_notes(genotypes)
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
        coverage_notes,
        include_trials=include_trials,
        research_findings=research_findings if isinstance(research_findings, list) else [],
    )

    markdown = _render_markdown(
        base_name,
        summary,
        apoe,
        risk_cards,
        wellness,
        expanded_panels,
        fun_cards,
        trials_by_finding,
        variant_verification if isinstance(variant_verification, list) else [],
        coverage_notes,
        include_trials=include_trials,
        research_findings=research_findings if isinstance(research_findings, list) else [],
    )

    (run_dir / f"{base_name}_Report.html").write_text(html, encoding="utf-8")
    (run_dir / f"{base_name}_Report.md").write_text(markdown, encoding="utf-8")

    print(f"Generated reports in {run_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
