# DNA Analysis Pipeline: Reproducible Workflow

This document captures the exact steps, scripts (you may need to create), and skills used to generate a comprehensive DNA analysis reports. Use this guide to replicate the analysis for any new AncestryDNA (or similar 23andMe/MyHeritage) raw data file.

## Prerequisites
*   **Use uv and Python 3.12+** and use PEP 723 inline metadata to handle dependencies in all Python scripts that run with `uv run`.
*   **Libraries:** `polars`, `requests` (add to python scripts in PEP 723 inline metadata as needed).
*   **Input File:** Your raw DNA text file (e.g., `ancestrydna-test-file-1.txt`).
*   **Local Reference Tables (required):**
    * `data/snp_reference.csv` (all panels and labels)
    * `data/clinical_interpretations.json` (APOE, CYP2C9, SLCO1B1, LPA mappings)
    * `data/gwas_risk_alleles.json` (optional local risk-allele table; can be empty)

**Naming Convention:** For all steps below, replace `[filename]` with the base name of your input file (e.g., `ancestrydna-test-file-1`) to avoid overwriting previous results.
**Run Directory Rule:** All derived outputs (parquet, JSON, reports, caches) are written to `runs/YYYYMMDD/[filename]/` and scripts will not read legacy outputs from the project root. If legacy outputs exist, re-run Step 1 or move them into the run folder.

---

## Recommended: One-Command Pipeline
**Goal:** Run Steps 1–9 end-to-end with consistent metadata capture.

1. **Script:** `uv run --script pipeline.py [filename]`
2. **Optional inputs:**
   * `--sex female|male` and `--age N` to enable hormone/age notes.
   * `--skip-trials` to disable clinical-trials lookup.
   * `--build-gwas path/to/gwas_file.tsv` to refresh `data/gwas_risk_alleles.json` before Step 3.
3. **Metadata:** Writes a `run_manifest` (uv/python/git) into `summary.json`.

## Optional Step 0: Build Local GWAS Risk Allele Table
**Goal:** Populate `data/gwas_risk_alleles.json` from a curated TSV/CSV so Step 3 can use local risk alleles without web calls.

1. **Script:** `uv run build_gwas_risk_table.py path/to/gwas_file.tsv`
2. **Accepted Inputs:**
   * CSV/TSV with columns `rsid` and `risk_allele`, **or**
   * GWAS Catalog TSV with columns `SNPS` and `STRONGEST SNP-RISK ALLELE`.
3. **Output:** `data/gwas_risk_alleles.json` (merged by default).
   * **If the file is missing:** Download it from: `https://www.ebi.ac.uk/gwas/api/search/downloads/associations/v1.0.2?split=false` (this triggers a download of `gwas_catalog_v1.0.2-associations_e115_r2026-01-19_full.zip` or similar). Unzip it and use the resulting TSV as the input for the script.
   * **Manual Download Alternative:** If the direct link fails, browse to `https://www.ebi.ac.uk/gwas/docs/file-downloads` and select Description: "All associations v1.0.2 - with added ontology annotations, GWAS Catalog study accession numbers and genotyping technology" and click "Download (full list)".
   * If you already have the GWAS Catalog TSV locally (e.g., `data/gwas-catalog-download-associations-alt-full.tsv`), reuse it; no re-download needed.
4. **Skip if already built:** If `data/gwas_risk_alleles.json` exists and has non-null entries for your rsIDs, you can skip rebuilding.
   * Quick check: `uv run python -c "import json; d=json.load(open('data/gwas_risk_alleles.json')); print(sum(v is not None for v in d.values()))"`
5. **Troubleshooting:** If you see `ModuleNotFoundError: polars`, run `uv run --with polars build_gwas_risk_table.py ...`

## Step 1: Ingestion & Quality Control
**Goal:** Validate file integrity, check for missing data, and normalize to a fast-reading format.

1.  **Script:** `uv run --script qc_analysis.py [filename]`
2.  **Action:** Reads the raw text file (skipping comments), checks for valid alleles (A, C, G, T), counts missing rows, and saves a normalized Parquet file.
    * QC adds: heterozygosity rate, sex inference (X/Y metrics), missingness by chromosome, duplicate rsID detection, ambiguous A/T or C/G SNP count, and build detection from header (supports GRCh37/HG19 and “build 37.x” wording).
3.  **Output:** `runs/YYYYMMDD/[filename]/[filename].normalized.parquet` (Used by all subsequent steps).
4.  **Additional Output:** `runs/YYYYMMDD/[filename]/summary.json` (QC stats + downstream file paths).
5.  **Optional:** If the user provides sex/gender, add `sex` to `summary.json` (values: `female`/`male`) so estrogen-related notes render correctly.
6.  **Skill Used:** `polars` (Fast dataframe processing).

## Step 2: Core Wellness & Trait Query
**Goal:** Check "Famous" SNPs for immediate health/wellness insights.

1.  **Script:** `uv run --script query_snps.py [filename]`
2.  **Variants Checked (from local table):**
    *   `rs4988235` (Lactose Intolerance)
    *   `rs671` (Alcohol Flush)
    *   `rs762551` (Caffeine Metabolism - CYP1A2)
    *   `rs5751876` (Caffeine Anxiety - ADORA2A) **[NEW]**
    *   `rs1815739` (Muscle Performance - ACTN3)
    *   `rs4680` (Stress/Cognition - COMT)
    *   `rs429358` + `rs7412` (Alzheimer's Risk - APOE Haplotype)
    *   `rs713598` (Bitter Taste)
    *   `rs72921001` (Cilantro Aversion - OR6A2) **[NEW]**
3.  **Action:** Reads `data/snp_reference.csv` for panels, filters the parquet file for those rsIDs, and prints genotypes. Writes `runs/YYYYMMDD/[filename]/core_traits.json`.

## Step 3: Strand & Allele Verification
**Goal:** Ensure we aren't misinterpreting "Minus Strand" reporting (common in microarray data).

1.  **Script:** `uv run --script verify_variants.py [filename]`
2.  **Action:** Queries the **Ensembl REST API** (GRCh37) for allele strings; caches results to `runs/YYYYMMDD/[filename]/variant_api_cache.json`.
3.  **Key Check:** Compares the "Observed" alleles in our file with Ensembl alleles to confirm strand orientation.
4.  **GWAS Risk Alleles:** Use local `data/gwas_risk_alleles.json` if populated. If an rsID is missing there, fall back to `data/snp_reference.csv` effect alleles (still local). Do not call GWAS REST API in this step.
5.  **Skills Used:** `ensembl-database` (local cache), `gwas-database` only if manual lookup needed.
6.  **Gotcha Check:** Confirm `rs4349859` (HLA-B27 proxy) strand alignment. If the file shows T/C while Ensembl reports A/G, it is a strand flip; treat the proxy cautiously and do not over-interpret.

## Step 4: Daily Life & Aging Analysis
**Goal:** Dig deeper into longevity, vitamins, and methylation.

1.  **Script:** `uv run --script life_aging_analysis.py [filename]`
2.  **Variants Checked (from local table):**
    *   `rs1801133` + `rs1801131` (MTHFR - Methylation status)
    *   `rs2802292` (FOXO3 - Longevity/Centenarian gene)
    *   `rs2282679` (Vitamin D levels)
    *   `rs602662` (Vitamin B12 absorption - FUT2) **[NEW]**
    *   `rs174546` (Omega-3 Conversion - FADS1) **[NEW]**
    *   `rs2274924` (Magnesium Absorption - TRPM6) **[NEW]**
    *   `rs1061170` (CFH - Macular Degeneration risk)
    *   `rs10490924` (ARMS2 - Macular Degeneration risk) **[NEW]**
    *   `rs6265` (BDNF - Brain plasticity/Exercise response)
    *   `rs9939609` (FTO - Appetite)

## Step 5: "Hidden" Actionable Risks
**Goal:** Check specific high-stakes pharmacogenomic and disease risks that are actionable.

1.  **Script:** `uv run --script check_extra_snps.py [filename]`
2.  **Variants Checked (from local table):**
    *   `rs4149056` (SLCO1B1 - Statin muscle pain risk)
    *   `rs4244285` (CYP2C19 - Clopidogrel response) **[NEW]**
    *   `rs1800562` + `rs1799945` (HFE - Hemochromatosis/Iron Overload)
    *   `rs6025` (Factor V Leiden - Clotting Risk) **[NEW]**
    *   `rs1799963` (Prothrombin G20210A - Clotting Risk) **[NEW]**
    *   `rs1333049` (9p21 - Early Heart Attack locus)
    *   `rs3798220` (LPA - Lipoprotein(a) heart risk) **[NEW]**
    *   `rs334` (HBB - Sickle cell trait/disease) **[NEW]**
    *   `rs113993960` (CFTR F508del - Cystic fibrosis carrier) **[NEW]**
    *   `rs28929474` + `rs17580` (SERPINA1 - Alpha-1 antitrypsin deficiency) **[NEW]**
    *   `rs1050828` + `rs1050829` (G6PD deficiency markers) **[NEW]**
    *   `rs5742904` (APOB R3500Q - Familial hypercholesterolemia) **[NEW]**
    *   `rs1801155` (APC I1307K - Colorectal cancer risk allele) **[NEW]**
    *   `rs17879961` (CHEK2 I157T - Moderate cancer risk allele) **[NEW]**
    *   `rs80357906` (BRCA1 5382insC founder) **[NEW]**
    *   `rs80359550` (BRCA2 6174delT founder) **[NEW]**

## Step 5b: Structural & Functional Lifestyle
**Goal:** Add structural health, injury risk, targeted diet sensitivities, celiac risk tags, and sleep/chronotype.

1.  **Script:** `uv run --script query_snps.py [filename]`
2.  **Panels & Variants (from local table):**
    * **Injury & Structural Health (Athlete Panel)**
        * `rs12722` (COL5A1 - Tendon/Ligament injury risk)
        * `rs1800012` (COL1A1 - Bone density/soft tissue strength)
        * `rs1799971` (OPRM1 - Pain sensitivity/opioid response)
    * **Specific Dietary Sensitivities**
        * `rs5082` (APOA2 - Saturated fat response)
        * `rs699` (AGT - Salt sensitivity)
        * `rs4343` (ACE proxy - Salt sensitivity)
        * `rs16890979` (SLC2A9 - Uric acid clearance/gout risk)
    * **Celiac & Gluten Risk (Tag SNPs)**
        * `rs3135388` (HLA-DQ2.5 tag)
        * `rs7454108` (HLA-DQ8 tag)
    * **Sleep & Chronotype**
        * `rs1801260` (CLOCK - Eveningness tendency)
        * `rs73598374` (ADA - Deep sleep intensity)
3.  **Notes:**
    * Tag SNPs (DQ2.5/DQ8) are screening tools, not diagnostic.
    * Absence of risk alleles at both DQ2.5 and DQ8 makes celiac disease extremely unlikely (very high negative predictive value).

## Step 6: Expanded SNP Panels (Cardiometabolic, Neuro, Cancer, PGx, Lifestyle, Functional Health)
**Goal:** Broaden coverage with additional high-value SNP panels.

1.  **Script:** `uv run --script additional_panels.py [filename]`
2.  **Panels & Variants (from local table):**
    * **Cardiometabolic**
        * `rs7903146` (TCF7L2 - Type 2 Diabetes)
        * `rs738409` (PNPLA3 I148M - NAFLD)
        * `rs10455872` (LPA - Lipoprotein(a) levels)
        * `rs662799` (APOA5 - Triglycerides)
        * `rs708272` (CETP - HDL-C)
    * **Neuro/Psych**
        * `rs34637584` (LRRK2 G2019S - Parkinson's)
        * `rs1052553` (MAPT H1/H2 - Parkinson's/PSP)
        * `rs17070145` (KIBRA - Memory performance)
        * `rs25531` (SLC6A4 proxy - Serotonin transporter)
    * **Cancer (common risk loci)**
        * `rs6983267` (8q24 - Colorectal/Prostate)
        * `rs2736100` (TERT - Multiple cancers)
        * `rs2981582` (FGFR2 - Breast cancer)
        * `rs10993994` (MSMB - Prostate cancer)
    * **Pharmacogenomics**
        * `rs9923231` (VKORC1 - Warfarin sensitivity)
        * `rs1799853` (CYP2C9*2 - Warfarin/NSAIDs)
        * `rs1057910` (CYP2C9*3 - Warfarin/NSAIDs)
        * `rs2108622` (CYP4F2*3 - Warfarin dose modifier) **[NEW]**
        * `rs12777823` (Warfarin ancestry modifier) **[NEW]**
        * `rs28371686` (CYP2C9*5) **[NEW]**
        * `rs9332131` (CYP2C9*6 indel) **[NEW]**
        * `rs7900194` (CYP2C9*8) **[NEW]**
        * `rs28371685` (CYP2C9*11) **[NEW]**
        * `rs1800460` (TPMT*3B - Thiopurines)
        * `rs1800462` (TPMT*2 - Thiopurines) **[NEW]**
        * `rs1142345` (TPMT*3C - Thiopurines)
        * `rs116855232` (NUDT15 - Thiopurines)
        * `rs2395029` (HLA-B*57:01 proxy - Abacavir)
        * `rs3918290` (DPYD*2A - Fluoropyrimidines)
        * `rs2231142` (ABCG2 Q141K - Statin exposure) **[NEW]**
        * `rs2306283` (SLCO1B1 388A>G - Haplotype support) **[NEW]**
        * `rs887829` (UGT1A1*28 proxy - Bilirubin risk) **[NEW]**
        * `rs8175347` (UGT1A1*28 TA repeat) **[NEW]**
        * `rs2844682` + `rs3909184` (HLA-B*15:02 proxy - Carbamazepine) **[NEW]**
        * `rs1061235` (HLA-A*31:01 proxy - Carbamazepine) **[NEW]**
        * `rs9263726` (HLA-B*58:01 proxy - Allopurinol) **[NEW]**
    * **CYP2D6 Risk Factors**
        * `rs1065852` (CYP2D6*10 - Risk factor only)
        * `rs3892097` (CYP2D6*4 - Risk factor only)
        * `rs16947` (CYP2D6*2 tag) **[NEW]**
        * `rs1135840` (CYP2D6*2 tag) **[NEW]**
        * `rs28371725` (CYP2D6*41 tag) **[NEW]**
        * `rs35742686` (CYP2D6*3 indel) **[NEW]**
        * `rs5030655` (CYP2D6*6 indel) **[NEW]**
    * **Lifestyle**
        * `rs1229984` (ADH1B - Alcohol metabolism)
        * `rs16969968` (CHRNA5 - Nicotine dependence)
    * **Functional Health - Histamine Intolerance**
        * `rs10156191` (AOC1/DAO - Histamine breakdown)
        * `rs2052129` (AOC1/DAO - Histamine expression)
        * `rs11558538` (HNMT - CNS histamine clearance)
    * **Functional Health - Detox/Acetylation**
        * `rs1801280` (NAT2*5 - Slow acetylation)
        * `rs1799930` (NAT2*6 - Slow acetylation)
        * `rs1799931` (NAT2*7 - Slow acetylation)
    * **Functional Health - Autoimmune Screening**
        * `rs3177928` (HLA-DRA - Thyroid autoimmune risk)
        * `rs7197` (HLA-DRA - Thyroid autoimmune risk, linked to Graves')
        * `rs4349859` (HLA-B27 proxy - Ankylosing spondylitis tag)
    * **Functional Health - Hormone Metabolism**
        * `rs2234693` (ESR1 - Estrogen sensitivity)
        * `rs4680` (COMT - Estrogen clearance context)
        * `rs1799853` (CYP2C9*2 - Estrogen clearance context; OC/HRT)
        * `rs1057910` (CYP2C9*3 - Estrogen clearance context; OC/HRT)
    * **Fun Traits & Appearance (user-friendly text + evidence note)**
        * **Appearance**
            * `rs12913832` (HERC2 - Eye/Hair Lightness)
            * `rs1805007` (MC1R - Red Hair/Freckles)
            * `rs1800407` (OCA2 - Eye Color Shade)
            * `rs1426654` (SLC24A5 - Skin Pigmentation)
            * `rs16891982` (SLC45A2 - Skin/Hair Pigmentation)
            * `rs885479` (MC1R - Freckles/Sun Sensitivity)
            * `rs7349332` (AR region - Male Pattern Baldness)
        * **Sensory**
            * `rs838133` (FGF21 - Sweet Preference)
            * `rs34160967` (TAS1R1 - Umami/Savory Sensitivity)
            * `rs4481887` (OR2M7 region - Asparagus Odor Detection)
        * **Body & Appetite**
            * `rs9930506` (FTO - BMI/Overeating Tendency)
            * `rs17782313` (MC4R - Appetite/Overeating)
            * `rs4846567` (LYPLAL1 - Fat Distribution)
        * **Sleep**
            * `rs1801260` (CLOCK - Morning/Evening Type)
            * `rs73598374` (ADA - Sleep Depth/Need)
        * **Behavior**
            * `rs1800497` (DRD2/ANKK1 - Reward/Risk Taking)
            * `rs16969968` (CHRNA5 - Nicotine Dependence Tendency)

**Interpretation Notes (use in reports):**
- **GWAS/common risk loci (e.g., cancer, many cardiometabolic markers):** Low effect sizes; not diagnostic. Use to contextualize risk, not to make clinical decisions.
- **Pharmacogenomics:** Actionable only when **clinical-grade PGx** testing confirms genotype and phenotype (e.g., warfarin dosing, thiopurines, DPYD, HLA-B*57:01).
- **HLA proxies:** Tag SNPs (e.g., HLA-B*15:02, HLA-A*31:01, HLA-B*58:01) are ancestry-dependent proxies; confirm with clinical HLA typing.
- **CYP2D6:** Do not assign metabolizer phenotype from array data; CNVs drive poor/ultra-rapid status and are not captured.
- **NAT2:** Use rs1801280/rs1799930/rs1799931 as a **partial panel**. If any marker is missing, report “unknown.” If all three are present, **≥2 slow alleles = likely slow acetylator (screening-level)**; 1 slow allele = indeterminate without full haplotyping. Highlight drug dosing/toxin-exposure relevance but require clinical confirmation.
- **Lifestyle/behavior markers:** Often ancestry-dependent and environment-modulated; treat as informational, not prescriptive.
- **Fun Traits & Appearance:** Use the report's friendly association text + short evidence note; always include "This is not medical advice."
- **Indels/repeats:** Some markers (e.g., BRCA founders, CFTR F508del, UGT1A1*28) are indels or repeats and may not be callable from array data; treat as not assessed.

## Step 7: Research Augmentation (Agent Task)
**Goal:** Find the latest (2024-2026) consensus on handling the specific genotypes found.

1.  **Optional Template Script (new):** `uv run --script build_research_findings.py [filename]`
    *   Writes `runs/YYYYMMDD/[filename]/research_findings.json` in UTF-8 (no BOM).
    *   **Default behavior:** skips output if no high-priority findings.
    *   Use `--write-empty` to force an empty file.
    *   Fill in `content` + `source` fields after targeted searches.
2.  **Action:** The Agent reviews previous steps, identifies high-priority risks, and performs targeted searches.
3.  **Tool:** use your web search tool
4.  **Output:** Save findings to `runs/YYYYMMDD/[filename]/research_findings.json`.
    *   **Format:** A JSON list of objects: `{"topic": "...", "content": "...", "source": "..."}`.
    *   Empty `content` entries are ignored in the report.
5.  **Example Queries:**
    *   "MTHFR compound heterozygote lifestyle recommendations 2025"
    *   "Prothrombin G20210A lifestyle precautions 2025"
    *   "Lp(a) high risk lifestyle management 2025"
    *   "CHRNA5 AA nicotine dependence implications"
    *   "CYP2C19 poor metabolizer drugs to avoid"
6.  **Rule:** Use local reference tables first; only run web searches when a SNP is clinically actionable or the local table lacks interpretation.

## Step 8: Clinical Trials Search (Conditional)
**Goal:** Only if **critical findings** are detected, search for recruiting clinical trials relevant to those findings.

1.  **Script:** `uv run --script search_trials_for_findings.py [filename]`
    * Optional filters: `--location "City, State"` or `--geo "distance(lat,lon,50mi)"` to narrow recruiting studies.
2.  **Action:** Reads detected **critical findings** (high-severity only), maps them to ClinicalTrials.gov search terms, and queries the **API v2** for `RECRUITING` studies. Writes `runs/YYYYMMDD/[filename]/trials_by_finding.json`. This avoids low-priority queries (e.g., NAT2-only findings) unless explicitly requested.
3.  **Behavior:** If **no critical findings** or **no mapped queries**, no trials file is produced and the report omits the trials section.
4.  **Skill Used:** `clinicaltrials-database`.

## Step 9: Final Report Generation (Markdown & HTML)
**Goal:** Synthesize all findings into polished documents.

1.  **Script:** `uv run --script generate_report.py [filename]`
2.  **Use Data Sources:** `summary.json`, `core_traits.json`, `healthy_aging.json`, `hidden_risks.json`, `expanded_panels.json`, `variant_verification.json`, `trials_by_finding.json`, `research_findings.json`, plus `data/clinical_interpretations.json`.
3.  **Content Must Include:**
    * **Critical Findings** (Actionable items like Prothrombin G20210A, CYP2C9)
    * **Wellness Tables** (Metabolism, Fitness, Functional Health)
    * **Research Augmentation** (Section 8: 2025/2026 Consensus notes from Step 7)
    * **Clinical Trials** (Recruiting studies from Step 8)
    * **Limitations** (Disclaimer and technical caveats)
4.  **Output:** `runs/YYYYMMDD/[filename]/[filename]_Report.md` and `runs/YYYYMMDD/[filename]/[filename]_Report.html`.
5.  **HTML Template:** Uses `report_template.html` in root; sections expand based on JSON availability.

