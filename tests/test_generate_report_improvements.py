import unittest

from generate_report import (
    _actionable_not_available,
    _apoe_assessment,
    _build_risk_cards,
    _count_phrase,
    _coverage_notes,
    _display_chromosome_label,
    _escalate_high_evidence_pgx,
    _expanded_panels,
    _hidden_screening_rows,
    _high_priority_findings,
    _panel_rows,
    _render_markdown,
    _validate_hbs_interpretation_guardrail,
    _validate_report_lints,
    _wellness_tables,
)


class GenerateReportImprovementsTests(unittest.TestCase):
    def test_strand_caution_only_applies_to_flagged_variants(self) -> None:
        rows = _panel_rows(
            [
                {
                    "rsid": "rs2052129",
                    "label": "Histamine intolerance risk",
                    "effect_allele": "T",
                    "effect_trait": "Reduced DAO activity.",
                    "non_effect_trait": "No risk detected.",
                    "notes": "Verify strand in Step 3.",
                }
            ],
            {"rs2052129": "GG"},
            {},
            panel_name="Functional Health - Histamine",
            variant_lookup={"rs2052129": {"match_status": "match", "proxy_note": None}},
            include_indicators=True,
        )
        self.assertEqual(len(rows), 1)
        detail = rows[0].get("detail") or ""
        self.assertNotIn("Strand caution", detail)
        self.assertEqual(rows[0].get("indicator"), "No high-confidence adverse flags")

    def test_strand_caution_is_added_when_variant_is_flagged(self) -> None:
        rows = _panel_rows(
            [
                {
                    "rsid": "rs10993994",
                    "label": "MSMB prostate risk",
                    "effect_allele": "C",
                    "effect_trait": "Risk marker present.",
                    "non_effect_trait": "No risk marker.",
                    "notes": "",
                }
            ],
            {"rs10993994": "TT"},
            {},
            panel_name="Cancer",
            variant_lookup={"rs10993994": {"match_status": "reverse_complement", "proxy_note": None}},
            include_indicators=True,
        )
        self.assertEqual(len(rows), 1)
        detail = rows[0].get("detail") or ""
        self.assertIn("Strand caution", detail)
        self.assertEqual(rows[0].get("indicator"), "Strand caution")

    def test_hla_b58_missing_stays_in_missing_bucket(self) -> None:
        expected_notes, missing_notes = _coverage_notes({}, {})
        self.assertFalse(any("HLA-B*58:01" in note for note in expected_notes))
        self.assertTrue(any("HLA-B*58:01" in note for note in missing_notes))
        self.assertTrue(
            any("missing rs9263726 proxy SNP in this file build" in note for note in missing_notes)
        )

    def test_dpyd_rs67376798_t_allele_generates_actionable_card(self) -> None:
        cards = _build_risk_cards({"rs67376798": "TT"})
        labels = {card.get("label") for card in cards}
        self.assertIn("Fluoropyrimidine Toxicity", labels)

    def test_high_evidence_escalation_adds_missing_pgx_card(self) -> None:
        cards = _escalate_high_evidence_pgx([], {"rs67376798": "TT"}, None)
        labels = {card.get("label") for card in cards}
        self.assertIn("Fluoropyrimidine Toxicity", labels)

    def test_high_priority_includes_medication_alerts(self) -> None:
        findings = _high_priority_findings(
            {},
            {},
            None,
            risk_cards=[
                {
                    "label": "Clopidogrel Response",
                    "category": "clinical",
                    "evidence": "CPIC",
                    "description": "CYP2C19*2 detected; reduced clopidogrel activation.",
                }
            ],
        )
        categories = {item.get("category") for item in findings}
        self.assertIn("Medication alerts (screening-level)", categories)
        labels = {item.get("label") for item in findings}
        self.assertIn("Clopidogrel Response", labels)

    def test_hla_b27_proxy_not_rendered_as_no_risk(self) -> None:
        rows = _panel_rows(
            [
                {
                    "rsid": "rs4349859",
                    "label": "HLA-B27 proxy",
                    "effect_allele": "A",
                    "effect_trait": "Risk detected",
                    "non_effect_trait": "No risk",
                    "notes": "Proxy tag SNP; ancestry-dependent.",
                }
            ],
            {"rs4349859": "GG"},
            {},
            panel_name="Functional Health - Autoimmune",
            variant_lookup={"rs4349859": {"match_status": "match", "proxy_note": "Proxy marker for HLA-B*27."}},
            include_indicators=True,
        )
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].get("value"), "Proxy marker")
        self.assertEqual(rows[0].get("indicator"), "Proxy marker")
        detail = rows[0].get("detail") or ""
        self.assertIn("ancestry and allele orientation", detail)

    def test_high_priority_hbs_tt_is_not_labeled_carrier(self) -> None:
        findings = _high_priority_findings(
            {"rs334": "TT"},
            {},
            None,
            risk_cards=[],
        )
        self.assertTrue(findings)
        self.assertEqual(findings[0].get("category"), "High-impact findings (screening-level)")
        self.assertIn("disease-level HbSS pattern", findings[0].get("note", ""))

    def test_dpyd_actionable_has_genotype_zygosity_and_completeness(self) -> None:
        cards = _build_risk_cards({"rs67376798": "TT"})
        dpyd = next(card for card in cards if card["label"] == "Fluoropyrimidine Toxicity")
        self.assertIn("rs67376798 TT", dpyd["description"])
        self.assertIn("homozygous risk allele", dpyd["description"])
        self.assertIn("Partial DPYD panel", dpyd["description"])

    def test_hidden_screening_uses_standardized_no_risk_phrase(self) -> None:
        rows = _hidden_screening_rows(
            {
                "records": [
                    {
                        "rsid": "rs4149056",
                        "label": "Statin muscle risk",
                        "effect_allele": "C",
                    }
                ]
            },
            {"rs4149056": "TT"},
            {},
            None,
            sex=None,
            qc_sex=None,
        )
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["value"], "No risk allele detected at rs4149056")
        self.assertIn("Screening-level; other variants not assessed.", rows[0]["note"])

    def test_association_summaries_include_rsid_for_cfh_and_9p21(self) -> None:
        cards = _build_risk_cards({"rs1061170": "CT", "rs1333049": "GG"})
        by_label = {card["label"]: card for card in cards}
        self.assertIn("Vision", by_label)
        self.assertIn("rs1061170", by_label["Vision"]["description"])
        self.assertIn("Early Heart Attack", by_label)
        self.assertIn("rs1333049", by_label["Early Heart Attack"]["description"])

    def test_count_phrase_grammar(self) -> None:
        self.assertEqual(_count_phrase(1, "variant"), "1 variant")
        self.assertEqual(_count_phrase(2, "variant"), "2 variants")

    def test_expanded_panel_proxy_marker_requires_real_proxy_note(self) -> None:
        panels = _expanded_panels(
            {
                "panels": {
                    "Pharmacogenomics": [
                        {"rsid": "rs67376798", "label": "DPYD c.2846A>T"},
                    ]
                }
            },
            {"rs67376798": "TT"},
            {},
            {"rs67376798": {"match_status": "match", "proxy_note": None}},
        )
        self.assertEqual(len(panels), 1)
        self.assertEqual(len(panels[0]["items"]), 1)
        self.assertNotIn("(proxy marker)", panels[0]["items"][0])

    def test_apoe_assessment_requires_complete_and_verified_markers(self) -> None:
        partial = _apoe_assessment(
            {"rs429358": "TT"},
            {"apoe_haplotype_map": {"TT|CC": "e3/e3"}},
            {"rs429358": {"match_status": "match"}},
        )
        self.assertFalse(partial["assessed"])
        self.assertIn("missing rs7412", partial["reason"])

        unverified = _apoe_assessment(
            {"rs429358": "TT", "rs7412": "CC"},
            {"apoe_haplotype_map": {"TT|CC": "e3/e3"}},
            {"rs429358": {"match_status": "match"}, "rs7412": {"match_status": "missing_in_file"}},
        )
        self.assertFalse(unverified["assessed"])
        self.assertIn("verification incomplete", unverified["reason"])

        verified = _apoe_assessment(
            {"rs429358": "TT", "rs7412": "CC"},
            {"apoe_haplotype_map": {"TT|CC": "e3/e3"}},
            {"rs429358": {"match_status": "match"}, "rs7412": {"match_status": "match"}},
        )
        self.assertTrue(verified["assessed"])
        self.assertEqual(verified["haplotype"], "e3/e3")

    def test_wellness_apoe_row_includes_rsids_or_not_assessed_status(self) -> None:
        assessed_rows = _wellness_tables(
            {"rs429358": "TT", "rs7412": "CC"},
            {},
            {
                "assessed": True,
                "haplotype": "e3/e3",
                "genotypes": {"rs429358": "TT", "rs7412": "CC"},
            },
            {"panels": {}, "fun_panels": {}},
            {},
            {"rs429358": {"match_status": "match"}, "rs7412": {"match_status": "match"}},
        )["fitness"]
        apoe_assessed = next(row for row in assessed_rows if row.get("label") == "Alzheimer's APOE")
        self.assertIn("rs429358 TT + rs7412 CC -> e3/e3", apoe_assessed.get("sub", ""))

        missing_rows = _wellness_tables(
            {},
            {},
            {
                "assessed": False,
                "reason": "APOE not assessed (partial/missing SNPs): missing rs7412.",
            },
            {"panels": {}, "fun_panels": {}},
            {},
            None,
        )["fitness"]
        apoe_missing = next(row for row in missing_rows if row.get("label") == "Alzheimer's APOE")
        self.assertEqual(apoe_missing.get("value"), "APOE not assessed (partial/missing SNPs)")

    def test_prothrombin_actionable_includes_genotype_and_zygosity(self) -> None:
        cards = _build_risk_cards({"rs1799963": "AG"})
        clot = next(card for card in cards if card["label"] == "Clotting Risk")
        self.assertIn("rs1799963 AG", clot["description"])
        self.assertIn("heterozygous risk allele", clot["description"])

    def test_hbs_guardrail_blocks_interpretation_when_rs334_missing(self) -> None:
        with self.assertRaises(ValueError):
            _validate_hbs_interpretation_guardrail(
                {},
                [{"label": "Sickle Cell (HbS)", "description": "HbS variant detected."}],
                [],
            )

    def test_chromosome_display_labels_nonstandard_values(self) -> None:
        self.assertEqual(_display_chromosome_label("1"), "1")
        self.assertEqual(_display_chromosome_label("23"), "X")
        self.assertEqual(_display_chromosome_label("26"), "Other/Contigs (26)")

    def test_association_summaries_include_rsid_and_genotype_for_lifestyle_items(self) -> None:
        cards = _build_risk_cards(
            {
                "rs16969968": "AA",
                "rs1801133": "AG",
                "rs1801131": "GT",
                "rs10455872": "AG",
            }
        )
        by_label = {card["label"]: card for card in cards}
        self.assertIn("Addiction Risk", by_label)
        self.assertIn("rs16969968", by_label["Addiction Risk"]["description"])
        self.assertIn("(AA)", by_label["Addiction Risk"]["description"])
        self.assertIn("Methylation", by_label)
        self.assertIn("rs1801133 AG + rs1801131 GT", by_label["Methylation"]["description"])
        self.assertIn("Heart Health", by_label)
        self.assertIn("rs10455872 (AG)", by_label["Heart Health"]["description"])

    def test_actionable_serpina1_includes_genotype_and_zygosity(self) -> None:
        cards = _build_risk_cards({"rs17580": "TT"})
        alpha1 = next(card for card in cards if card["label"] == "Alpha-1 Antitrypsin Deficiency")
        self.assertIn("rs17580 TT", alpha1["description"])
        self.assertIn("homozygous risk allele", alpha1["description"])

    def test_nat2_marker_rows_use_non_contradictory_partial_panel_wording(self) -> None:
        rows = _panel_rows(
            [
                {
                    "rsid": "rs1801280",
                    "label": "Slow acetylation",
                    "effect_allele": "C",
                    "effect_trait": "Slow acetylation tendency",
                    "non_effect_trait": "No slow marker",
                }
            ],
            {"rs1801280": "TT"},
            {},
            panel_name="Functional Health - Detox/Acetylation",
            include_indicators=True,
        )
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].get("label"), "NAT2 SNP (partial panel)")
        self.assertEqual(rows[0].get("value"), "NAT2 SNP observed (partial panel)")
        self.assertEqual(rows[0].get("indicator"), "Genotype observed")

    def test_warfarin_cards_include_vkorc1_panel_status(self) -> None:
        cards = _build_risk_cards({"rs1057910": "AC", "rs2108622": "CT"})
        relevant = [
            card for card in cards
            if "warfarin" in (card.get("description", "") + " " + card.get("action", "")).lower()
        ]
        self.assertTrue(relevant)
        for card in relevant:
            text = f"{card['description']} {card['action']}".lower()
            self.assertIn("warfarin panel status:", text)
            self.assertIn("vkorc1", text)
            self.assertIn("missing", text)
            self.assertIn("do not use a full cpic dosing calculator", text)

    def test_actionable_not_available_lists_major_missing_items(self) -> None:
        items = _actionable_not_available({}, {})
        labels = {item["label"] for item in items}
        self.assertIn("VKORC1 (warfarin sensitivity)", labels)
        self.assertIn("CYP2C19*2 (clopidogrel)", labels)
        self.assertIn("CYP2C19*17", labels)

    def test_wellness_tables_compacts_nat2_rows(self) -> None:
        wellness = _wellness_tables(
            {"rs1801280": "TT", "rs1799930": "GG", "rs1799931": "GG"},
            {},
            {"assessed": False, "reason": "APOE not assessed (partial/missing SNPs): missing rs429358, rs7412."},
            {
                "panels": {
                    "Functional Health - Detox/Acetylation": [
                        {"rsid": "rs1801280", "label": "Slow acetylation", "effect_allele": "C"},
                        {"rsid": "rs1799930", "label": "Slow acetylation", "effect_allele": "A"},
                        {"rsid": "rs1799931", "label": "Slow acetylation", "effect_allele": "A"},
                        {"rsid": "rs1695", "label": "GSTP1 Ile105Val", "effect_allele": "G"},
                    ]
                },
                "fun_panels": {},
            },
            {},
            None,
        )
        labels = [str(row.get("label")) for row in wellness["functional"]]
        self.assertIn("NAT2 acetylation status", labels)
        self.assertNotIn("NAT2 SNP (partial panel)", labels)

    def test_report_lints_fail_for_summary_child_mismatch(self) -> None:
        with self.assertRaises(ValueError):
            _validate_report_lints(
                [{"label": "CYP2C9 Reduced Function", "description": "warfarin", "action": "missing vkorc1"}],
                {
                    "functional": [
                        {
                            "label": "Inflammation summary",
                            "value": "No high-confidence adverse flags detected in this screened set",
                            "row_type": "summary",
                        },
                        {"label": "IL-6", "status": "risk", "value": "Risk allele present", "sub": "rs1800795 GC"},
                    ]
                },
            )

    def test_serpina1_pis_message_is_softened(self) -> None:
        cards = _build_risk_cards({"rs17580": "TT"})
        alpha1 = next(card for card in cards if card["label"] == "Alpha-1 Antitrypsin Deficiency")
        self.assertIn("mild AAT reduction", alpha1["description"])
        self.assertIn("serum AAT testing", alpha1["action"])

    def test_markdown_trials_are_capped_and_codes_are_clickable(self) -> None:
        trials = [
            {
                "nct_id": f"NCT0000000{i}",
                "title": f"Trial {i}",
                "phase": "Phase 2",
                "url": f"https://clinicaltrials.gov/study/NCT0000000{i}",
            }
            for i in range(1, 8)
        ]
        output = _render_markdown(
            "sample",
            {"total_snps": 1, "call_rate_percent": 99.0},
            "Unknown",
            [],
            {"metabolism": [], "fitness": [], "functional": []},
            [],
            [],
            [{"finding_label": "Clotting Risk", "finding_level": "high", "query_term": "clot", "recruiting_studies": trials}],
            [],
            [],
            [],
            [],
            [],
            [],
            [],
            [],
            include_trials=True,
            research_findings=[],
        )
        self.assertIn("Showing top 5 of 7 recruiting trials.", output)
        self.assertIn("[**NCT00000001**](https://clinicaltrials.gov/study/NCT00000001)", output)
        self.assertIn("[**NCT00000005**](https://clinicaltrials.gov/study/NCT00000005)", output)
        self.assertNotIn("NCT00000006", output)
        self.assertIn("https://clinicaltrials.gov/study/NCT00000001", output)


if __name__ == "__main__":
    unittest.main()
