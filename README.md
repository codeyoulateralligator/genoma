# VCF Explorer

A privacy-first Python tool for exploring VCF and VCF.GZ genome files locally.

VCF Explorer reads variant-call data, summarizes the file structure and genotype statistics, checks a curated set of well-known genetic markers, and generates an educational PDF report with clear tables and explanations.

## Features

* Reads `.vcf` and `.vcf.gz` files
* Detects common genome builds where possible
* Supports single- and multi-sample VCF files
* Prints a structured console summary
* Counts:

  * SNVs
  * insertions
  * deletions
  * heterozygous calls
  * homozygous reference calls
  * homozygous alternate calls
  * missing calls
  * transitions and transversions
* Shows variant counts by chromosome/contig
* Checks curated educational loci, including:

  * CCR5-Δ32
  * APOE markers
  * lactase persistence
  * HFE C282Y and H63D
  * ALDH2
  * MTHFR C677T
  * Factor V Leiden
  * pigmentation and eye-colour markers
* Explains:

  * your called genotype
  * reference and alternate alleles
  * what `0/0`, `0/1`, and `1/1` mean for each marker
  * known limitations and when clinical confirmation is appropriate
* Generates a table-based PDF report and optional JSON output
* Runs locally: your genome file is not uploaded anywhere by the script

## Usage

```bash
python genoma.py Genoom.vcf.gz
```

Optional output paths:

```bash
python genoma.py Genoom.vcf.gz \
  --pdf my_genome_report.pdf \
  --json my_genome_report.json
```

## Requirements

```bash
pip install reportlab
```

The parser otherwise uses Python's standard library.

## Example Output

```text
Factor V Leiden (rs6025; F5)

Your VCF call:
GT=1/1; hom_alt; bases=C/C

Reference / comparison:
REF=T; ALT=C
Reference genotype: T/T (0/0)

Interpretation:
Two Factor V Leiden-associated alleles detected.
This is a potentially important finding to confirm with a clinical-grade test.
```

## Important Limitations

This project is for educational exploration of VCF data only.

* A VCF may not include every genomic position.
* A missing record does not prove that a variant is absent.
* `0/0`, `0/1`, and `1/1` are relative to the VCF's `REF` and `ALT` allele representation.
* Genome build differences, strand orientation, phasing, variant normalization, and call quality can affect interpretation.
* Consumer, research, or filtered VCF files are not automatically clinical-grade.
* Do not use this tool for diagnosis, medication decisions, reproductive decisions, or disease-risk conclusions.
* Confirm medically important findings with a qualified clinician and a clinical laboratory test.

## Privacy

VCF Explorer is designed to run entirely on your own machine. Genetic data is highly sensitive; avoid committing VCF files, generated reports, or JSON exports to public repositories.
