# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "polars",
# ]
# ///

import re
import sys
from pathlib import Path

import polars as pl

from run_utils import resolve_base_name, run_root, update_summary


_VALID_ALLELES = {"A", "C", "G", "T"}


def _detect_build(input_path: Path) -> dict[str, str | None]:
    build = None
    source_line = None
    pattern = re.compile(r"(grch\s?3[78]|hg1[89]|build\s*3[78](?:\.\d+)?)", re.IGNORECASE)
    try:
        with input_path.open("r", encoding="utf-8", errors="ignore") as handle:
            for _ in range(200):
                line = handle.readline()
                if not line:
                    break
                if not line.startswith("#"):
                    continue
                if "build" not in line.lower() and "grch" not in line.lower() and "hg" not in line.lower():
                    continue
                source_line = line.strip()
                match = pattern.search(source_line)
                if not match:
                    continue
                token = match.group(1).lower().replace(" ", "")
                if token.startswith("build"):
                    if "37" in token:
                        build = "GRCh37"
                    elif "38" in token:
                        build = "GRCh38"
                elif token in {"grch37", "hg19"}:
                    build = "GRCh37"
                elif token in {"grch38", "hg38"}:
                    build = "GRCh38"
                break
    except OSError:
        return {"build_detected": None, "build_source_line": None}
    return {"build_detected": build, "build_source_line": source_line}


def _normalize_chromosome(expr: pl.Expr) -> pl.Expr:
    upper = expr.cast(pl.String).str.strip_chars().str.to_uppercase()
    return (
        pl.when(upper.is_in(["23", "X"]))
        .then(pl.lit("X"))
        .when(upper.is_in(["24", "Y"]))
        .then(pl.lit("Y"))
        .when(upper.is_in(["25", "MT", "M"]))
        .then(pl.lit("MT"))
        .otherwise(upper)
    )


def _infer_sex(x_called: int, x_hetero_rate: float | None, y_called: int) -> tuple[str, str]:
    if x_called == 0:
        return "ambiguous", "No X chromosome calls available."
    if x_hetero_rate is None:
        return "ambiguous", "Unable to compute X heterozygosity."
    if y_called >= 100 and x_hetero_rate < 0.05:
        return "male", "Y calls present with low X heterozygosity."
    if y_called <= 5 and x_hetero_rate > 0.2:
        return "female", "Minimal Y calls with higher X heterozygosity."
    return "ambiguous", "X/Y metrics do not strongly indicate male or female."

def process_dna_file(input_path: str, output_path: str, base_name: str) -> None:
    print(f"Processing {input_path}...")
    
    try:
        # AncestryDNA files are tab-delimited with headers
        # We handle potential type issues by enforcing schema
        # chromosome can be '23', 'X', 'Y', 'MT' -> String
        # position -> Int64
        
        df = pl.read_csv(
            input_path,
            separator="\t",
            comment_prefix="#",
            has_header=True,
            schema_overrides={
                "chromosome": pl.String,
                "position": pl.Int64,
                "allele1": pl.String,
                "allele2": pl.String,
                "rsid": pl.String
            },
            ignore_errors=True
        )
        
        # Normalize column names to lowercase
        df = df.rename({col: col.lower() for col in df.columns})
        df = df.with_columns(
            _normalize_chromosome(pl.col("chromosome")).alias("chr_norm"),
            pl.col("allele1").str.to_uppercase().alias("allele1_u"),
            pl.col("allele2").str.to_uppercase().alias("allele2_u"),
        )

        # Basic validation stats
        total_count = df.height

        missing_mask = (pl.col("allele1_u").is_in(["0", "--"])) | (pl.col("allele2_u").is_in(["0", "--"]))
        valid_mask = pl.col("allele1_u").is_in(list(_VALID_ALLELES)) & pl.col("allele2_u").is_in(list(_VALID_ALLELES))
        called_mask = valid_mask & ~missing_mask
        invalid_mask = ~missing_mask & ~valid_mask

        df = df.with_columns(
            missing_flag=missing_mask,
            invalid_flag=invalid_mask,
            called_flag=called_mask,
            hetero_flag=called_mask & (pl.col("allele1_u") != pl.col("allele2_u")),
        )

        missing_count = df.filter(pl.col("missing_flag")).height
        invalid_count = df.filter(pl.col("invalid_flag")).height
        called_count = df.filter(pl.col("called_flag")).height

        hetero_count = df.filter(pl.col("hetero_flag")).height
        heterozygosity_rate = (hetero_count / called_count) if called_count else 0.0

        ambiguous_mask = (
            (pl.col("allele1_u").is_in(["A", "T"]) & pl.col("allele2_u").is_in(["A", "T"]))
            | (pl.col("allele1_u").is_in(["C", "G"]) & pl.col("allele2_u").is_in(["C", "G"]))
        )
        df = df.with_columns(
            ambiguous_flag=called_mask & ambiguous_mask,
            missing_or_invalid=pl.col("missing_flag") | pl.col("invalid_flag"),
        )
        ambiguous_count = df.filter(pl.col("ambiguous_flag")).height

        print(f"Total SNPs processed: {total_count}")
        print(f"No-calls (0/0): {missing_count}")
        print(f"Invalid allele rows: {invalid_count}")
        print(f"Call Rate:{((called_count) / total_count * 100):.2f}%")

        # Missingness by chromosome
        missing_by_chr = (
            df.group_by("chr_norm")
            .agg(
                pl.len().alias("total"),
                pl.sum("missing_or_invalid").alias("missing"),
                pl.sum("called_flag").alias("called"),
            )
            .sort("chr_norm")
            .to_dicts()
        )

        # Duplicate rsIDs
        dup_df = df.group_by("rsid").len().filter(pl.col("len") > 1)
        duplicate_count = dup_df.height
        duplicate_examples = dup_df.select(pl.col("rsid")).head(5).to_series().to_list()

        # Sex inference
        x_called = df.filter(pl.col("called_flag") & (pl.col("chr_norm") == "X")).height
        x_hetero = df.filter(pl.col("hetero_flag") & (pl.col("chr_norm") == "X")).height
        x_hetero_rate = (x_hetero / x_called) if x_called else None
        y_called = df.filter(pl.col("called_flag") & (pl.col("chr_norm") == "Y")).height
        sex_inference, sex_note = _infer_sex(x_called, x_hetero_rate, y_called)
        
        # Create a combined 'genotype' column for easier querying (e.g., "AG", "CC")
        # We sort alleles alphabetically for consistent matching (e.g., "GA" becomes "AG")
        # except for 0 which we leave or handle.
        
        # Logic: 
        # 1. Replace 0 with ''
        # 2. Concat
        # 3. Sort characters
        
        # For performance in polars, we can just keep allele1/2 and handle sorting at query time, 
        # or do a struct sort. Let's just write the parquet now.
        
        df.write_parquet(output_path)
        print(f"Successfully saved to {output_path}")
        update_summary(
            run_root(base_name),
            {
                "base_name": base_name,
                "input_file": input_path,
                "normalized_parquet": output_path,
                "total_snps": total_count,
                "no_calls": missing_count,
                "invalid_allele_rows": invalid_count,
                "call_rate_percent": round((called_count) / total_count * 100, 2),
                "heterozygosity_rate": round(heterozygosity_rate, 4),
                "ambiguous_snp_count": ambiguous_count,
                "missing_by_chromosome": missing_by_chr,
                "duplicate_rsid_count": duplicate_count,
                "duplicate_rsid_examples": duplicate_examples,
                "sex_inference": sex_inference,
                "sex_inference_note": sex_note,
                "sex_inference_x_called": x_called,
                "sex_inference_x_heterozygosity_rate": (
                    round(x_hetero_rate, 4) if x_hetero_rate is not None else None
                ),
                "sex_inference_y_called": y_called,
                **_detect_build(Path(input_path)),
            },
        )
        
    except Exception as e:
        print(f"CRITICAL ERROR: {e}")
        sys.exit(1)

if __name__ == "__main__":
    base_name = resolve_base_name(sys.argv[1] if len(sys.argv) > 1 else None)
    input_file = f"{base_name}.txt"
    run_dir = run_root(base_name)
    output_file = run_dir / f"{base_name}.normalized.parquet"

    if not input_file or not Path(input_file).exists():
        print(f"Input file not found: {input_file}")
        sys.exit(1)

    process_dna_file(input_file, str(output_file), base_name)
