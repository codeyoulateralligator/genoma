#!/usr/bin/env python3
"""
genoma.py - read-only educational VCF explorer with console tables and PDF report.

Examples
--------
    python3 genoma.py Genoom.vcf.gz
    python3 genoma.py Genoom.vcf.gz --pdf my_genetics_report.pdf
    python3 genoma.py Genoom.vcf.gz --no-pdf --json report.json

Dependencies
------------
- Python standard library for parsing and console output.
- reportlab only for PDF output:
      python -m pip install reportlab

Important
---------
This program is educational only. It does not validate a call, measure sequencing
coverage, infer a clinical diagnosis, or make medication decisions. It must not be
used as a substitute for clinical genetic testing or medical advice.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import sys
import textwrap
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Curated educational loci
# ---------------------------------------------------------------------------
# Coordinates are GRCh37/hg19. rsID matching works irrespective of build if the
# rsID is retained by the source VCF. Coordinate fallback is deliberately used
# only when the VCF header confidently identifies GRCh37/hg19.

@dataclass(frozen=True)
class Locus:
    rsid: str
    chrom: str
    pos_b37: int
    label: str
    gene: str
    category: str
    short_description: str
    comparison_note: str
    genotype_0: str
    genotype_1: str
    genotype_2: str
    caution: str


CURATED_LOCI: List[Locus] = [
    Locus(
        "rs333", "3", 46414943, "CCR5-Delta32", "CCR5", "Infectious-disease biology",
        "A 32-base deletion in CCR5. At this locus the VCF reference allele is the non-deletion sequence and the alternate allele is the deletion representation.",
        "Comparison call: 0/0 means two non-deletion alleles at this VCF record.",
        "0 deletion alleles detected. This does not provide special protection from HIV.",
        "1 deletion allele detected. Some studies associate this with partial protection or altered disease course for CCR5-tropic HIV, but it is not protection against HIV and is not a reason to change prevention or testing.",
        "2 deletion alleles detected. Strongly reduces cell-surface CCR5 and is associated with resistance to many CCR5-tropic HIV strains, but not all HIV strains; it is not immunity.",
        "Never use this result to estimate HIV risk or make prevention decisions. A VCF call is not a medical test."
    ),
    Locus(
        "rs429358", "19", 45411941, "APOE marker 1", "APOE", "Neurology / lipid biology",
        "One of two SNPs used to define APOE epsilon alleles together with rs7412.",
        "This single marker has no stand-alone APOE type. Compare it with rs7412 and, if heterozygous, phasing may matter.",
        "Reference allele at this marker. APOE interpretation still requires rs7412.",
        "One alternate allele at this marker. APOE interpretation still requires rs7412 and may require phase.",
        "Two alternate alleles at this marker. APOE interpretation still requires rs7412.",
        "APOE type is a susceptibility marker, not a diagnosis or prediction."
    ),
    Locus(
        "rs7412", "19", 45412079, "APOE marker 2", "APOE", "Neurology / lipid biology",
        "The second SNP used with rs429358 to define APOE epsilon alleles.",
        "This single marker has no stand-alone APOE type. Compare it with rs429358 and, if heterozygous, phasing may matter.",
        "Reference allele at this marker. APOE interpretation still requires rs429358.",
        "One alternate allele at this marker. APOE interpretation still requires rs429358 and may require phase.",
        "Two alternate alleles at this marker. APOE interpretation still requires rs429358.",
        "APOE type is a susceptibility marker, not a diagnosis or prediction."
    ),
    Locus(
        "rs4988235", "2", 136608646, "Lactase-persistence marker", "MCM6 / LCT", "Nutrition trait",
        "A regulatory marker commonly associated with adult lactase persistence in many European populations.",
        "At GRCh37 this is usually G>A; the A allele is the common persistence-associated allele in many Europeans.",
        "No A allele at this marker. In European-derived populations this is more compatible with reduced adult lactase persistence, but symptoms depend on diet, microbiome, dose, and other genetics.",
        "One A allele. Often associated with adult lactase persistence in European-derived populations.",
        "Two A alleles. Often associated with adult lactase persistence in European-derived populations.",
        "This marker does not diagnose lactose intolerance and is less informative in some non-European ancestries."
    ),
    Locus(
        "rs1800562", "6", 26093141, "HFE C282Y", "HFE", "Iron metabolism",
        "A common HFE variant associated with hereditary hemochromatosis risk when relevant combinations are present.",
        "At GRCh37 this is commonly G>A; the A allele corresponds to the C282Y-associated variant.",
        "No C282Y-associated allele detected at this record.",
        "One C282Y-associated allele detected. Usually a carrier state; interpretation depends on H63D, iron tests, sex, family history, and clinical context.",
        "Two C282Y-associated alleles detected. This genotype is associated with increased hemochromatosis risk, but penetrance is incomplete and clinical assessment requires ferritin and transferrin saturation.",
        "Do not diagnose iron overload from genotype alone."
    ),
    Locus(
        "rs1799945", "6", 26091179, "HFE H63D", "HFE", "Iron metabolism",
        "A common HFE variant that can matter in combination with other HFE variants.",
        "At GRCh37 this is commonly C>G; the G allele corresponds to the H63D-associated variant.",
        "No H63D-associated allele detected at this record.",
        "One H63D-associated allele detected. Alone, it is usually not sufficient to cause iron overload.",
        "Two H63D-associated alleles detected. Can be associated with mild iron-parameter changes in some people; clinical relevance is variable.",
        "Interpret jointly with rs1800562 and actual iron studies."
    ),
    Locus(
        "rs671", "12", 112241766, "ALDH2*2 marker", "ALDH2", "Alcohol metabolism",
        "A functional alcohol-metabolism marker. The alternate A allele is commonly called ALDH2*2.",
        "At GRCh37 this is commonly G>A. G/G is the reference comparison call in this VCF representation.",
        "No ALDH2*2-associated allele detected at this record.",
        "One ALDH2*2-associated allele detected. Often associated with reduced acetaldehyde metabolism and flushing after alcohol; alcohol-related health implications depend on behaviour and other factors.",
        "Two ALDH2*2-associated alleles detected. Typically associated with markedly reduced enzyme activity and strong flushing response.",
        "Do not use this as a diagnosis or as permission to consume alcohol."
    ),
    Locus(
        "rs1801133", "1", 11856378, "MTHFR C677T", "MTHFR", "Folate-pathway marker",
        "A common MTHFR marker. The VCF may use the opposite DNA strand, so the GRCh37 G>A representation corresponds to the commonly named C677T variant.",
        "Reference comparison is G/G in this VCF; the alternate A represents the T-equivalent allele on the named transcript strand.",
        "No T-equivalent allele detected at this marker.",
        "One T-equivalent allele detected. Common; usually not clinically decisive on its own.",
        "Two T-equivalent alleles detected. Can modestly affect enzyme activity, but genotype alone generally does not determine folate status or treatment.",
        "Do not use MTHFR genotype alone for treatment decisions."
    ),
    Locus(
        "rs6025", "1", 169519049, "Factor V Leiden", "F5", "Blood clotting",
        "A well-known F5 variant associated with activated-protein-C resistance and venous thrombosis susceptibility.",
        "At GRCh37 this is commonly T>C; C is the Factor V Leiden-associated allele. Reference comparison: T/T (0/0).",
        "No Factor V Leiden-associated allele detected at this record.",
        "One Factor V Leiden-associated allele detected. Risk of venous thrombosis is elevated compared with non-carriers, but most carriers never develop a clot.",
        "Two Factor V Leiden-associated alleles detected. Risk is higher than for one-copy carriers. This is a potentially important finding to clinically confirm.",
        "Confirm an important result with a clinical-grade test; discuss it before estrogen treatment, pregnancy, major surgery, or prolonged immobilisation."
    ),
    Locus(
        "rs4244285", "10", 94781859, "CYP2C19*2 marker", "CYP2C19", "Pharmacogenetics",
        "A no-function CYP2C19 star-allele marker used in pharmacogenetics.",
        "At GRCh37 this is commonly G>A; the A allele contributes to CYP2C19*2. A full phenotype needs a validated multi-variant haplotype panel.",
        "No *2 marker allele detected at this record.",
        "One *2 marker allele detected. May contribute to reduced CYP2C19 activity, but the full phenotype cannot be inferred from this SNP alone.",
        "Two *2 marker alleles detected. Often compatible with a no-function *2/*2 genotype, but a full validated panel is still required.",
        "Never change medication from this report. Use a clinical pharmacogenetic test and prescribing guideline."
    ),
    Locus(
        "rs12248560", "10", 94761900, "CYP2C19*17 marker", "CYP2C19", "Pharmacogenetics",
        "An increased-function CYP2C19 star-allele marker used in pharmacogenetics.",
        "At GRCh37 this is commonly C>T; the T allele contributes to CYP2C19*17. A full phenotype needs a validated multi-variant haplotype panel.",
        "No *17 marker allele detected at this record.",
        "One *17 marker allele detected. May contribute to increased CYP2C19 activity, but a full phenotype cannot be inferred from this SNP alone.",
        "Two *17 marker alleles detected. Often compatible with increased activity, but a full validated panel is still required.",
        "*17 can be in complex combination with other variants; do not infer drug dosing from one or two SNPs."
    ),
    Locus(
        "rs1815739", "11", 66328095, "ACTN3 R577X", "ACTN3", "Muscle-function trait",
        "A common ACTN3 stop-gain marker often discussed in sports genetics.",
        "At GRCh37 this is commonly C>T; T is the stop-gain (X) allele in this marker's standard representation.",
        "Two R-associated alleles at this marker. This does not predict athletic performance.",
        "One R-associated and one X-associated allele. Common and not diagnostic.",
        "Two X-associated alleles. Alpha-actinin-3 is absent in fast fibres; this is common and does not prevent normal activity or athletic success.",
        "Sport performance is highly polygenic and dominated by training, health, opportunity, and environment."
    ),
    Locus(
        "rs12913832", "15", 28365618, "HERC2 / OCA2 eye-colour marker", "HERC2", "Pigmentation trait",
        "A major eye-colour association marker in many Europeans, but not a deterministic predictor.",
        "At GRCh37 this is commonly A>G; G is associated with lower OCA2 expression and, in many Europeans, lighter/blue eye colour.",
        "Two A alleles. In many Europeans this is more often associated with brown/darker eyes, but exceptions are common.",
        "One G allele. Intermediate probabilities; eye colour cannot be reliably predicted from this one marker.",
        "Two G alleles. In many Europeans this is strongly associated with blue/light eyes, though not guaranteed.",
        "This is an ancestry-dependent association, not an identity inference."
    ),
    Locus(
        "rs1426654", "15", 48426484, "SLC24A5 pigmentation marker", "SLC24A5", "Pigmentation trait",
        "A population-associated pigmentation marker with substantial allele-frequency differences across ancestries.",
        "Use the VCF REF/ALT bases shown in the report as the direct comparison; strand conventions can make named alleles look reversed across databases.",
        "Reference/reference at this VCF record. Interpretation is population dependent.",
        "One alternate allele at this VCF record. Interpretation is population dependent.",
        "Two alternate alleles at this VCF record. Interpretation is population dependent.",
        "This marker is not a personal identity inference and does not determine appearance by itself."
    ),
    Locus(
        "rs601338", "19", 49206674, "FUT2 secretor-status marker", "FUT2", "Mucosal biology",
        "One common FUT2 marker related to secretor status, especially in European populations.",
        "The functional interpretation can depend on strand convention and other FUT2 variants; use this marker only as a partial clue.",
        "Reference/reference at this VCF record. This alone does not establish secretor status.",
        "One alternate allele at this VCF record. This alone does not establish secretor status.",
        "Two alternate alleles at this VCF record. This alone does not establish secretor status.",
        "Secretor status should not be inferred confidently from a single marker in every ancestry."
    ),
]

SOURCES = [
    ("MedlinePlus Genetics - Factor V Leiden thrombophilia", "https://medlineplus.gov/genetics/condition/factor-v-leiden-thrombophilia/"),
    ("MedlinePlus Genetics - Direct-to-consumer genetic testing", "https://medlineplus.gov/genetics/understanding/dtcgenetictesting/dtcghrpages/"),
    ("CPIC - CYP2C19 allele definition table", "https://files.cpicpgx.org/data/report/current/allele_definition/CYP2C19_allele_definition_table.xlsx"),
    ("CPIC - CYP2C19 pharmacogenetics guideline resources", "https://cpicpgx.org/guidelines/"),
    ("NCBI dbSNP", "https://www.ncbi.nlm.nih.gov/snp/"),
]


# ---------------------------------------------------------------------------
# VCF parsing
# ---------------------------------------------------------------------------

def open_text_vcf(path: str):
    with open(path, "rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b" or path.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def normalise_chrom(chrom: str) -> str:
    chrom = chrom.strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    return "MT" if chrom.upper() in {"M", "MT"} else chrom


def parse_contig_header(line: str) -> Optional[Tuple[str, int]]:
    if not line.startswith("##contig=<"):
        return None
    inside = line[len("##contig=<"):].rstrip().rstrip(">")
    fields = dict(part.split("=", 1) for part in inside.split(",") if "=" in part)
    try:
        return normalise_chrom(fields["ID"]), int(fields["length"])
    except (KeyError, ValueError):
        return None


def infer_build(contigs: Dict[str, int]) -> str:
    if contigs.get("1") == 249250621 and contigs.get("X") == 155270560:
        return "GRCh37 / hg19"
    if contigs.get("1") == 248956422 and contigs.get("X") == 156040895:
        return "GRCh38 / hg38"
    return "Unknown"


def parse_gt(sample_field: str, format_field: str) -> Tuple[Optional[List[Optional[int]]], bool]:
    if not sample_field or sample_field == "." or not format_field:
        return None, False
    keys = format_field.split(":")
    values = sample_field.split(":")
    try:
        gt = values[keys.index("GT")]
    except (ValueError, IndexError):
        return None, False
    if gt in {".", "./.", ".|."}:
        return None, "|" in gt
    phased = "|" in gt
    parsed: List[Optional[int]] = []
    for value in gt.replace("|", "/").split("/"):
        if value == ".":
            parsed.append(None)
        else:
            try:
                parsed.append(int(value))
            except ValueError:
                return None, phased
    return parsed, phased


def gt_text(alleles: Optional[List[Optional[int]]], phased: bool) -> str:
    if alleles is None:
        return "./."
    separator = "|" if phased else "/"
    return separator.join("." if item is None else str(item) for item in alleles)


def genotype_class(alleles: Optional[List[Optional[int]]]) -> str:
    if alleles is None or any(item is None for item in alleles):
        return "missing"
    if all(item == 0 for item in alleles):
        return "hom_ref"
    if all(item > 0 for item in alleles):
        return "hom_alt"
    return "het"


def alt_copy_count(alleles: Optional[List[Optional[int]]]) -> Optional[int]:
    if alleles is None or any(item is None for item in alleles):
        return None
    return sum(item > 0 for item in alleles)


def allele_strings(ref: str, alt_text: str, alleles: Optional[List[Optional[int]]]) -> List[str]:
    all_alleles = [ref] + ([] if alt_text in {"", "."} else alt_text.split(","))
    if alleles is None:
        return []
    output: List[str] = []
    for item in alleles:
        if item is None:
            output.append(".")
        elif 0 <= item < len(all_alleles):
            output.append(all_alleles[item])
        else:
            output.append(f"<invalid:{item}>")
    return output


def variant_type(ref: str, alts: List[str]) -> str:
    if not alts or alts == ["."]:
        return "unknown"
    if any(item.startswith("<") or item in {"*", "."} for item in alts):
        return "symbolic"
    lengths = {(len(ref), len(item)) for item in alts}
    if all(ref_len == 1 and alt_len == 1 for ref_len, alt_len in lengths):
        return "SNV"
    if all(ref_len == alt_len for ref_len, alt_len in lengths):
        return "MNV"
    if all(alt_len > ref_len for ref_len, alt_len in lengths):
        return "insertion"
    if all(alt_len < ref_len for ref_len, alt_len in lengths):
        return "deletion"
    return "complex"


def is_transition(ref: str, alt: str) -> Optional[bool]:
    if len(ref) != 1 or len(alt) != 1:
        return None
    pair = {ref.upper(), alt.upper()}
    if pair in ({"A", "G"}, {"C", "T"}):
        return True
    if ref.upper() in "ACGT" and alt.upper() in "ACGT":
        return False
    return None


def wrapped(text: str, width: int = 34) -> str:
    return "\n".join(textwrap.wrap(str(text), width=width, break_long_words=False))


def short_allele(allele: str) -> str:
    if len(allele) <= 22:
        return allele
    return f"{allele[:9]}...{allele[-9:]} ({len(allele)} bp)"


def apoe_summary(calls: Dict[str, Dict[str, Any]]) -> Optional[str]:
    """Only give exact APOE type when unambiguous from called alleles."""
    a = calls.get("rs429358")
    b = calls.get("rs7412")
    if not a or not b or not a.get("record") or not b.get("record"):
        return None
    ra, rb = a["record"], b["record"]
    if ra["classification"] == "missing" or rb["classification"] == "missing":
        return "APOE: one or both defining markers are missing."
    # Base pairs in GRCh37 forward reference: rs429358 T/C, rs7412 C/T.
    # haplotype: T/C=e2, T/T=e3, C/C=e4. Exact diplotype requires phase when both are het.
    alleles_a = ra["called_alleles"]
    alleles_b = rb["called_alleles"]
    if len(alleles_a) != 2 or len(alleles_b) != 2:
        return "APOE: non-diploid or incomplete call; not interpreted."
    if ra["phased"] and rb["phased"]:
        haps = list(zip(alleles_a, alleles_b))
        mapping = {("T", "T"): "e2", ("T", "C"): "e3", ("C", "C"): "e4"}
        if all(h in mapping for h in haps):
            return "APOE diplotype (phased): " + "/".join(sorted(mapping[h] for h in haps))
    # Fully homozygous loci give unambiguous diplotypes even unphased.
    if len(set(alleles_a)) == 1 and len(set(alleles_b)) == 1:
        mapping = {("T", "T"): "e2/e2", ("T", "C"): "e3/e3", ("C", "C"): "e4/e4"}
        value = mapping.get((alleles_a[0], alleles_b[0]))
        return f"APOE diplotype: {value}" if value else "APOE marker combination is unusual or not represented by the basic e2/e3/e4 mapping."
    return "APOE: markers are present, but phase is needed for an unambiguous diplotype."


def interpret_locus(locus: Locus, record: Optional[Dict[str, Any]]) -> Dict[str, str]:
    """Create the result, reference comparison and other-genotype descriptions."""
    other = (
        f"0/0: {locus.genotype_0}\n"
        f"0/1: {locus.genotype_1}\n"
        f"1/1: {locus.genotype_2}"
    )
    if not record:
        return {
            "result": "Not present in this VCF. This is not evidence of the reference genotype or of absence of the variant.",
            "reference": "No VCF record to show REF/ALT bases or comparison genotype.",
            "other": other,
        }
    if record["classification"] == "missing":
        return {
            "result": "The VCF record exists, but the genotype is missing or incomplete.",
            "reference": f"Reference allele: {short_allele(record['ref'])}; alternate allele(s): {', '.join(short_allele(x) for x in record['alts'])}.",
            "other": other,
        }

    n_alt = alt_copy_count(record["alleles"])
    if n_alt == 0:
        result = locus.genotype_0
    elif n_alt == 1:
        result = locus.genotype_1
    elif n_alt == 2:
        result = locus.genotype_2
    else:
        result = "Non-standard / multi-allelic genotype; see the exact VCF call."

    called = "/".join(short_allele(x) for x in record["called_alleles"])
    ref = short_allele(record["ref"])
    alt = ", ".join(short_allele(x) for x in record["alts"])
    comparison = (
        f"Your called bases: {called}. VCF reference/base comparison: REF={ref}; ALT={alt}; "
        f"reference genotype is {ref}/{ref} (normally GT 0/0). {locus.comparison_note}"
    )
    return {"result": result, "reference": comparison, "other": other}


def scan_vcf(path: str) -> Dict[str, Any]:
    loci_by_rsid = {locus.rsid: locus for locus in CURATED_LOCI}
    loci_by_coord = {(normalise_chrom(locus.chrom), locus.pos_b37): locus for locus in CURATED_LOCI}

    contigs: Dict[str, int] = {}
    samples: List[str] = []
    sample_stats: Dict[str, Counter] = defaultdict(Counter)
    chrom_records: Counter = Counter()
    variant_types: Counter = Counter()
    filter_counts: Counter = Counter()
    ti = tv = total_records = 0
    hits: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    ccr5_fallback: Optional[Dict[str, Any]] = None

    def make_record(fields: List[str], sample_name: str, sample_field: str) -> Dict[str, Any]:
        chrom = normalise_chrom(fields[0])
        alleles, phased = parse_gt(sample_field, fields[8])
        alts = fields[4].split(",") if fields[4] not in {"", "."} else []
        return {
            "chrom": chrom,
            "pos": int(fields[1]),
            "id": fields[2],
            "ref": fields[3],
            "alt": fields[4],
            "alts": alts,
            "sample": sample_name,
            "alleles": alleles,
            "called_alleles": allele_strings(fields[3], fields[4], alleles),
            "gt": gt_text(alleles, phased),
            "phased": phased,
            "classification": genotype_class(alleles),
        }

    with open_text_vcf(path) as handle:
        for raw in handle:
            if raw.startswith("##"):
                parsed = parse_contig_header(raw)
                if parsed:
                    contigs[parsed[0]] = parsed[1]
                continue
            if raw.startswith("#CHROM"):
                samples = raw.rstrip("\n").split("\t")[9:]
                continue
            if raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            total_records += 1
            chrom = normalise_chrom(fields[0])
            try:
                pos = int(fields[1])
            except ValueError:
                continue
            alts = fields[4].split(",") if fields[4] else []
            chrom_records[chrom] += 1
            variant_types[variant_type(fields[3], alts)] += 1
            filter_counts[fields[6]] += 1
            found_loci = [loci_by_rsid[item] for item in fields[2].split(";") if item in loci_by_rsid]
            ccr5_like = (
                chrom == "3" and 46414920 <= pos <= 46414950 and
                any(len(fields[3]) - len(alt) >= 30 for alt in alts if not alt.startswith("<"))
            )
            if len(fields) >= 10:
                for sample_name, sample_field in zip(samples, fields[9:]):
                    record = make_record(fields, sample_name, sample_field)
                    cls = record["classification"]
                    sample_stats[sample_name][cls] += 1
                    if cls != "missing":
                        for allele_idx in record["alleles"] or []:
                            if allele_idx and 0 < allele_idx <= len(alts):
                                transit = is_transition(fields[3], alts[allele_idx - 1])
                                if transit is True:
                                    ti += 1
                                elif transit is False:
                                    tv += 1
                    for locus in found_loci:
                        hits[locus.rsid].append({"match": "rsID", **record})
                    if ccr5_like:
                        ccr5_fallback = {"match": "CCR5 deletion-size fallback", **record}

    build = infer_build(contigs)

    # Coordinate fallback only after build confirmation. One small second pass avoids accidental cross-build matches.
    missing = {key: locus for key, locus in loci_by_coord.items() if locus.rsid not in hits}
    if build == "GRCh37 / hg19" and missing:
        with open_text_vcf(path) as handle:
            for raw in handle:
                if raw.startswith("#"):
                    continue
                fields = raw.rstrip("\n").split("\t")
                if len(fields) < 10:
                    continue
                try:
                    key = (normalise_chrom(fields[0]), int(fields[1]))
                except ValueError:
                    continue
                locus = missing.get(key)
                if not locus:
                    continue
                for sample_name, sample_field in zip(samples, fields[9:]):
                    record = make_record(fields, sample_name, sample_field)
                    hits[locus.rsid].append({"match": "GRCh37 coordinate", **record})

    if ccr5_fallback is not None and "rs333" not in hits:
        hits["rs333"].append(ccr5_fallback)

    # Select the first record for each locus and sample; standard single-sample VCFs have one.
    selected: Dict[str, Dict[str, Any]] = {}
    for locus in CURATED_LOCI:
        records = hits.get(locus.rsid, [])
        selected[locus.rsid] = {
            "locus": asdict(locus),
            "record": records[0] if records else None,
            "all_records": records,
        }
        selected[locus.rsid].update(interpret_locus(locus, records[0] if records else None))

    report = {
        "file": os.path.abspath(path),
        "file_size_bytes": os.path.getsize(path),
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        "reference_build_inferred": build,
        "samples": samples,
        "contigs_in_header": len(contigs),
        "records": total_records,
        "variant_types": dict(variant_types),
        "filters": dict(filter_counts),
        "sample_genotype_classes": {name: dict(counter) for name, counter in sample_stats.items()},
        "records_by_chromosome": dict(chrom_records),
        "called_snv_transition_count": ti,
        "called_snv_transversion_count": tv,
        "called_snv_titv": round(ti / tv, 4) if tv else None,
        "curated": selected,
        "apoe_summary": apoe_summary(selected),
        "sources": SOURCES,
    }
    return report


# ---------------------------------------------------------------------------
# Console report
# ---------------------------------------------------------------------------

def print_table(headers: List[str], rows: List[List[str]], widths: List[int]) -> None:
    def rule(char: str = "-") -> str:
        return "+" + "+".join(char * (width + 2) for width in widths) + "+"

    def split_cell(value: str, width: int) -> List[str]:
        parts = []
        for line in str(value).splitlines() or [""]:
            parts.extend(textwrap.wrap(line, width=width, break_long_words=False) or [""])
        return parts

    print(rule("="))
    header_lines = [split_cell(header, width) for header, width in zip(headers, widths)]
    for idx in range(max(map(len, header_lines))):
        print("|" + "|".join(f" {lines[idx] if idx < len(lines) else '':<{width}} " for lines, width in zip(header_lines, widths)) + "|")
    print(rule("="))
    for row in rows:
        cell_lines = [split_cell(value, width) for value, width in zip(row, widths)]
        for idx in range(max(map(len, cell_lines))):
            print("|" + "|".join(f" {lines[idx] if idx < len(lines) else '':<{width}} " for lines, width in zip(cell_lines, widths)) + "|")
        print(rule())


def print_console_report(report: Dict[str, Any], top_chromosomes: int) -> None:
    print("=" * 108)
    print("VCF EXPLORER - READ-ONLY EDUCATIONAL REPORT")
    print("=" * 108)
    print(f"File:               {report['file']}")
    print(f"Compressed size:    {report['file_size_bytes']:,} bytes")
    print(f"Reference build:    {report['reference_build_inferred']}")
    print(f"Samples:            {', '.join(report['samples']) or '(none found)'}")
    print(f"Records:            {report['records']:,}")
    print(f"Contigs in header:  {report['contigs_in_header']:,}")

    print("\nVARIANT SUMMARY")
    for name, count in Counter(report["variant_types"]).most_common():
        print(f"  {name:12} {count:>12,}")
    print("\nGENOTYPE CALL SUMMARY")
    for sample_name in report["samples"]:
        stats = report["sample_genotype_classes"].get(sample_name, {})
        missing = stats.get("missing", 0)
        print(f"  {sample_name}: called={report['records'] - missing:,} missing={missing:,} "
              f"hom-ref={stats.get('hom_ref', 0):,} het={stats.get('het', 0):,} hom-alt={stats.get('hom_alt', 0):,}")
    print("\nCALLED SNV Ti/Tv")
    print(f"  transitions={report['called_snv_transition_count']:,} "
          f"transversions={report['called_snv_transversion_count']:,} "
          f"Ti/Tv={report['called_snv_titv'] if report['called_snv_titv'] is not None else 'n/a'}")

    print(f"\nTOP {top_chromosomes} CONTIGS BY RECORD COUNT")
    for chrom, count in Counter(report["records_by_chromosome"]).most_common(top_chromosomes):
        print(f"  {chrom:>8} {count:>12,}")

    print("\nCURATED LOCI - RESULTS AND REFERENCE COMPARISON")
    rows: List[List[str]] = []
    for locus in CURATED_LOCI:
        item = report["curated"][locus.rsid]
        record = item["record"]
        if record:
            call = (f"GT {record['gt']} ({record['classification']})\n"
                    f"bases {'/'.join(short_allele(x) for x in record['called_alleles'])}\n"
                    f"{record['chrom']}:{record['pos']} {short_allele(record['ref'])}>{','.join(short_allele(x) for x in record['alts'])}")
        else:
            call = "Not present in VCF"
        rows.append([
            f"{locus.label}\n{locus.rsid} ({locus.gene})",
            call,
            item["result"],
            item["reference"],
        ])
    print_table(["Locus", "Your VCF call", "Meaning of your call", "Reference / comparison"], rows, [22, 28, 42, 55])

    if report["apoe_summary"]:
        print("\nCOMBINED APOE INTERPRETATION")
        print(f"  {report['apoe_summary']}")

    print("\nOTHER POSSIBLE CALLS")
    other_rows = [[f"{l.label}\n{l.rsid}", report["curated"][l.rsid]["other"]] for l in CURATED_LOCI]
    print_table(["Locus", "How 0/0, 0/1 and 1/1 would be described"], other_rows, [28, 110])

    print("\nIMPORTANT LIMITS")
    for text in [
        "A missing VCF record is not evidence that you have the reference genotype or that a variant is absent.",
        "0/0, 0/1 and 1/1 are relative to this VCF record's REF and ALT alleles. They are not universal labels for normal, healthy, good or bad.",
        "This file has genotype calls but no depth/quality fields in the displayed output; a call cannot be treated as clinically confirmed from this report.",
        "Do not use this output for diagnosis, medication dosing, HIV-risk decisions, or pregnancy/surgical decisions. Confirm anything consequential with an accredited clinical test and a clinician or genetic counsellor.",
    ]:
        print(f"  - {text}")


# ---------------------------------------------------------------------------
# PDF generation
# ---------------------------------------------------------------------------

def make_pdf(report: Dict[str, Any], output_path: str) -> None:
    try:
        from reportlab.lib import colors
        from reportlab.lib.enums import TA_CENTER, TA_LEFT
        from reportlab.lib.pagesizes import A4, landscape
        from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
        from reportlab.lib.units import mm
        from reportlab.platypus import (
            KeepTogether, PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle
        )
    except ImportError as exc:
        raise RuntimeError("PDF generation needs reportlab. Install it with: python -m pip install reportlab") from exc

    page_size = landscape(A4)
    doc = SimpleDocTemplate(
        output_path,
        pagesize=page_size,
        rightMargin=11 * mm,
        leftMargin=11 * mm,
        topMargin=10 * mm,
        bottomMargin=11 * mm,
        title="VCF Explorer Educational Report",
        author="vcf_explorer_v2.py",
    )
    styles = getSampleStyleSheet()
    title = ParagraphStyle("TitleCustom", parent=styles["Title"], fontName="Helvetica-Bold", fontSize=18, leading=22, spaceAfter=6)
    subtitle = ParagraphStyle("Subtitle", parent=styles["Normal"], fontSize=8.5, leading=11, textColor=colors.HexColor("#444444"), spaceAfter=7)
    h1 = ParagraphStyle("H1Custom", parent=styles["Heading1"], fontName="Helvetica-Bold", fontSize=12, leading=15, spaceBefore=8, spaceAfter=5)
    h2 = ParagraphStyle("H2Custom", parent=styles["Heading2"], fontName="Helvetica-Bold", fontSize=9.5, leading=12, spaceBefore=7, spaceAfter=4)
    normal = ParagraphStyle("NormalCustom", parent=styles["Normal"], fontSize=7.4, leading=9.2, spaceAfter=3)
    small = ParagraphStyle("Small", parent=normal, fontSize=6.5, leading=7.8)
    cell = ParagraphStyle("Cell", parent=normal, fontSize=6.3, leading=7.6)
    head = ParagraphStyle("Head", parent=cell, fontName="Helvetica-Bold", textColor=colors.white, alignment=TA_CENTER)
    alert = ParagraphStyle("Alert", parent=normal, fontSize=7.2, leading=9, backColor=colors.HexColor("#FFF4E5"), borderColor=colors.HexColor("#E7A33E"), borderWidth=0.5, borderPadding=5, spaceBefore=4, spaceAfter=6)

    def P(text: Any, style=normal) -> Paragraph:
        # Escape arbitrary VCF text while retaining only the controlled markup that
        # this script emits (<br/>, <b>, </b>). This prevents raw values from
        # breaking ReportLab paragraphs.
        raw = str(text).replace("\n", "<br/>")
        keep = {"<br/>": "__BR__", "<b>": "__BOPEN__", "</b>": "__BCLOSE__"}
        for original, token in keep.items():
            raw = raw.replace(original, token)
        escaped = raw.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
        for original, token in keep.items():
            escaped = escaped.replace(token, original)
        return Paragraph(escaped, style)

    def footer(canvas, doc_obj):
        canvas.saveState()
        canvas.setFont("Helvetica", 6.5)
        canvas.setFillColor(colors.HexColor("#555555"))
        canvas.drawString(11 * mm, 7 * mm, "Educational VCF report - not a clinical test or diagnosis")
        canvas.drawRightString(page_size[0] - 11 * mm, 7 * mm, f"Page {doc_obj.page}")
        canvas.restoreState()

    story: List[Any] = []
    story.append(P("VCF Explorer - Educational Genetic Variant Report", title))
    story.append(P(
        f"Generated {report['generated_utc']} from {os.path.basename(report['file'])}. "
        f"The file contains {report['records']:,} variant records and the inferred reference build is "
        f"{report['reference_build_inferred']}. This report describes a small transparent set of loci; it is not a genome-wide clinical interpretation.",
        subtitle,
    ))
    story.append(P(
        "<b>Important:</b> A VCF call is not a laboratory confirmation. REF/ALT labels are how this particular VCF represents a site; 0/0 is not a universal \"normal\" or \"good\" result. Do not make medical, medication, HIV-prevention, pregnancy, or surgery decisions from this report.",
        alert,
    ))

    story.append(P("1. File and calling summary", h1))
    sample_summary = []
    for sample in report["samples"]:
        stats = report["sample_genotype_classes"].get(sample, {})
        sample_summary.append(f"{sample}: called {report['records'] - stats.get('missing', 0):,}; missing {stats.get('missing', 0):,}; het {stats.get('het', 0):,}; hom-ref {stats.get('hom_ref', 0):,}; hom-alt {stats.get('hom_alt', 0):,}")
    summary_data = [
        [P("Item", head), P("Value", head)],
        [P("Input file", cell), P(report["file"], cell)],
        [P("Compressed size", cell), P(f"{report['file_size_bytes']:,} bytes", cell)],
        [P("Inferred build", cell), P(report["reference_build_inferred"], cell)],
        [P("Sample(s)", cell), P(", ".join(report["samples"]) or "No sample columns found", cell)],
        [P("VCF records", cell), P(f"{report['records']:,}", cell)],
        [P("Variant types", cell), P("; ".join(f"{key}: {value:,}" for key, value in Counter(report["variant_types"]).most_common()), cell)],
        [P("Called SNV Ti/Tv", cell), P(f"{report['called_snv_titv']} ({report['called_snv_transition_count']:,} transitions / {report['called_snv_transversion_count']:,} transversions)", cell)],
        [P("Genotype calls", cell), P("<br/>".join(sample_summary), cell)],
    ]
    table = Table(summary_data, colWidths=[43 * mm, 225 * mm], repeatRows=1)
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#24476A")),
        ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#B7C2CD")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("BACKGROUND", (0, 1), (0, -1), colors.HexColor("#EAF0F5")),
        ("LEFTPADDING", (0, 0), (-1, -1), 4), ("RIGHTPADDING", (0, 0), (-1, -1), 4),
        ("TOPPADDING", (0, 0), (-1, -1), 3), ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ]))
    story.append(table)

    story.append(P("2. Curated loci - your calls, reference comparison, and meanings", h1))
    story.append(P(
        "Reference/comparison reports the actual REF and ALT bases in your VCF and the usual reference genotype (REF/REF). This is a sequence comparison, not a healthy-person reference range. The 'other possible calls' column spells out the standard biallelic 0/0, 0/1, and 1/1 cases.",
        normal,
    ))

    rows: List[List[Any]] = [[P("Locus", head), P("Your call", head), P("Meaning of your call", head), P("Reference / comparison", head), P("Other possible calls", head)]]
    for locus in CURATED_LOCI:
        item = report["curated"][locus.rsid]
        record = item["record"]
        call = "Not present in this VCF"
        if record:
            bases = "/".join(short_allele(x) for x in record["called_alleles"]) if record["called_alleles"] else "not called"
            call = (f"<b>{record['gt']}</b> ({record['classification']})<br/>"
                    f"bases: {bases}<br/>"
                    f"VCF: {record['chrom']}:{record['pos']} {short_allele(record['ref'])}>{', '.join(short_allele(x) for x in record['alts'])}")
        rows.append([
            P(f"<b>{locus.label}</b><br/>{locus.rsid}<br/>{locus.gene}", cell),
            P(call, cell),
            P(item["result"], cell),
            P(item["reference"], cell),
            P(item["other"], cell),
        ])
    results_table = Table(rows, colWidths=[33 * mm, 43 * mm, 65 * mm, 72 * mm, 64 * mm], repeatRows=1, splitByRow=1)
    results_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#24476A")),
        ("GRID", (0, 0), (-1, -1), 0.2, colors.HexColor("#AAB8C5")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F6F9FC")]),
        ("LEFTPADDING", (0, 0), (-1, -1), 3), ("RIGHTPADDING", (0, 0), (-1, -1), 3),
        ("TOPPADDING", (0, 0), (-1, -1), 2.5), ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5),
    ]))
    story.append(results_table)

    story.append(PageBreak())
    story.append(P("3. Detailed explanations", h1))
    for locus in CURATED_LOCI:
        item = report["curated"][locus.rsid]
        record = item["record"]
        status = "Not present in VCF" if not record else f"GT {record['gt']}; {record['classification']}; bases {'/'.join(short_allele(x) for x in record['called_alleles'])}"
        detail_rows = [
            [P("Category", head), P("Explanation", head)],
            [P("Your VCF call", cell), P(status, cell)],
            [P("What this marker is", cell), P(locus.short_description, cell)],
            [P("Interpretation of your call", cell), P(item["result"], cell)],
            [P("Reference / comparison", cell), P(item["reference"], cell)],
            [P("0/0 possibility", cell), P(locus.genotype_0, cell)],
            [P("0/1 possibility", cell), P(locus.genotype_1, cell)],
            [P("1/1 possibility", cell), P(locus.genotype_2, cell)],
            [P("Caution", cell), P(locus.caution, cell)],
        ]
        detail_table = Table(detail_rows, colWidths=[44 * mm, 233 * mm], repeatRows=1)
        detail_table.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#24476A")),
            ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#B7C2CD")),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("BACKGROUND", (0, 1), (0, -1), colors.HexColor("#EAF0F5")),
            ("LEFTPADDING", (0, 0), (-1, -1), 4), ("RIGHTPADDING", (0, 0), (-1, -1), 4),
            ("TOPPADDING", (0, 0), (-1, -1), 3), ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
        ]))
        story.append(P(f"{locus.label} ({locus.rsid}; {locus.gene})", h2))
        story.append(detail_table)
        story.append(Spacer(1, 3))

    if report["apoe_summary"]:
        story.append(P("4. Combined APOE note", h1))
        story.append(P(report["apoe_summary"], normal))
        story.append(P("The exact APOE diplotype is only stated when the two calls make it unambiguous. If either marker is heterozygous and unphased, the two variants can be arranged on chromosomes in more than one way.", normal))

    story.append(P("5. Limitations and sources", h1))
    limits = [
        "This report is limited to the listed loci. It does not screen all genes or all clinically relevant variants.",
        "A VCF can omit a locus entirely; absence of a row is not absence of a variant.",
        "No quality/depth validation, copy-number analysis, structural-variant validation, or clinical annotation was performed.",
        "Allele names can differ by genome build and DNA strand. The report shows the actual bases in this VCF so they can be checked.",
        "Clinical relevance depends on phenotype, family history, ancestry, environmental exposures, co-variants, and laboratory confirmation.",
    ]
    story.append(P("<br/>".join(f"- {item}" for item in limits), normal))
    story.append(P("Background resources used when composing the educational descriptions:", h2))
    source_rows = [[P("Resource", head), P("Link", head)]] + [[P(name, cell), P(url, small)] for name, url in report["sources"]]
    source_table = Table(source_rows, colWidths=[100 * mm, 177 * mm], repeatRows=1)
    source_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#24476A")),
        ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#B7C2CD")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 4), ("RIGHTPADDING", (0, 0), (-1, -1), 4),
        ("TOPPADDING", (0, 0), (-1, -1), 3), ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ]))
    story.append(source_table)

    doc.build(story, onFirstPage=footer, onLaterPages=footer)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def default_pdf_path(vcf_path: str) -> str:
    base = os.path.basename(vcf_path)
    base = re.sub(r"\.vcf(?:\.gz)?$", "", base, flags=re.IGNORECASE)
    return os.path.join(os.path.dirname(os.path.abspath(vcf_path)), f"{base}_genetic_report.pdf")


def main() -> int:
    parser = argparse.ArgumentParser(description="Read-only VCF explorer with console tables and an educational PDF report.")
    parser.add_argument("vcf", help="Path to .vcf or .vcf.gz")
    parser.add_argument("--pdf", nargs="?", const="AUTO", default="AUTO", help="Write a PDF report. Supply a path or omit the value to use the default next to the VCF.")
    parser.add_argument("--no-pdf", action="store_true", help="Do not write a PDF report.")
    parser.add_argument("--json", dest="json_path", help="Optional path to write the structured report as JSON.")
    parser.add_argument("--top-chromosomes", type=int, default=12, help="Number of top contigs to print (default: 12).")
    args = parser.parse_args()

    if not os.path.isfile(args.vcf):
        print(f"ERROR: file not found: {args.vcf}", file=sys.stderr)
        return 2
    if args.no_pdf and args.pdf != "AUTO":
        print("ERROR: --pdf and --no-pdf cannot be used together.", file=sys.stderr)
        return 2

    report = scan_vcf(args.vcf)
    print_console_report(report, args.top_chromosomes)

    if args.json_path:
        with open(args.json_path, "w", encoding="utf-8") as output:
            json.dump(report, output, indent=2, ensure_ascii=False)
        print(f"\nWrote JSON report: {os.path.abspath(args.json_path)}")

    if not args.no_pdf:
        pdf_path = default_pdf_path(args.vcf) if args.pdf == "AUTO" else args.pdf
        try:
            make_pdf(report, pdf_path)
        except RuntimeError as exc:
            print(f"\nPDF NOT WRITTEN: {exc}", file=sys.stderr)
            return 1
        print(f"\nWrote PDF report: {os.path.abspath(pdf_path)}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
