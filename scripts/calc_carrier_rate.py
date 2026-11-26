#!/usr/bin/env python3
"""
General-purpose population carrier rate and allele frequency calculator from
paired sites-only VEP VCF/BCF and multi-sample genotype (GT) VCF/BCF files.

Supported modes (mutually exclusive):
  1. --sites-file : Evaluate a user-supplied TSV/CSV list of target variants.
  2. --region     : Scan one or more genomic regions (comma-separated).
  3. --bed        : Scan regions from a BED file (4th column optionally specifies gene SYMBOL).

Variant filtering (--region / --bed mode):
  --vep-filter    : Comma-separated terms matched against VEP Consequence,
                    IMPACT, and CLIN_SIG (case-insensitive).  When omitted,
                    ALL annotated variants are included.
  --exclude-conflicting : Exclude variants with 'conflicting' in CLIN_SIG.

Author: Shujia Huang & Gemini
Date: 2026-11-25


# 1. 区域模式(统计 PTEN 基因范围内 ClinVar P/LP + VEP 预测的高风险变异)
python calc_carrier_rate.py \
    --vep-vcf cohort.sites.vep.vcf.gz \
    --gt-vcf cohort.gt.bcf \
    --region chr10:87863625-87971930 \
    --gene PTEN \
    --vep-filter pathogenic,likely_pathogenic,HIGH \
    --exclude-conflicting \
    -s subcohort_samples.txt \
    -o pten_summary.tsv \
    --carriers-out pten_carriers.tsv

# 2. 多基因 BED (比如输入整个 ACMG 基因列表 BED)
python calc_carrier_rate.py \
    --vep-vcf cohort.sites.vep.vcf.gz \
    --gt-vcf cohort.gt.bcf \
    --bed acmg_genes.bed \
    --vep-filter pathogenic,likely_pathogenic \
    --exclude-conflicting \
    -o acmg_panel_summary.tsv

# 3. 外部自定义位点表模式 (传入包含 chrom,pos,ref,alt 的 TSV/CSV)
python calc_carrier_rate.py \
    --vep-vcf cohort.sites.vep.vcf.gz \
    --gt-vcf  cohort.gt.bcf \
    --sites-file target_hotspots.tsv \
    -o hotspots_summary.tsv

# 4. 全注释模式（纳入区域内所有有 VEP 注释的变异）
python calc_carrier_rate.py \
    --vep-vcf cohort.sites.vep.vcf.gz \
    --gt-vcf cohort.gt.bcf \
    --region chr10:87863625-87971930 \
    --gene PTEN \
    -o pten_all_variants.tsv \
    --carriers-out pten_all_carriers.tsv

# 5. 按 Consequence 筛选（只看错义 + 蛋白截短变异）
python calc_carrier_rate.py \
    --vep-vcf cohort.sites.vep.vcf.gz \
    --gt-vcf cohort.gt.bcf \
    --bed acmg_genes.bed \
    --vep-filter missense_variant,stop_gained,frameshift_variant,splice_donor_variant,splice_acceptor_variant,start_lost,stop_lost \
    -o acmg_coding.tsv

# 6. 提取全部 ClinVar 收录的变异（无论致病性）
python calc_carrier_rate.py \
    --vep-vcf cohort.sites.vep.vcf.gz \
    --gt-vcf cohort.gt.bcf \
    --bed acmg_genes.bed \
    --vep-filter clinvar \
    -o acmg_clinvar_all.tsv

# 7. ClinVar 全部收录 + 新发 pLoF
python calc_carrier_rate.py \
    --vep-vcf cohort.sites.vep.vcf.gz \
    --gt-vcf cohort.gt.bcf \
    --bed acmg_genes.bed \
    --vep-filter clinvar,HIGH \
    --exclude-conflicting \
    -o acmg_clinvar_plus_plof.tsv
    
一些值得注意的地方：

1. 按照 VEP Consequence 注释信息提取 pLoF 变异，两个方式：
    (1) --vep-filter stop_gained,frameshift_variant,splice_donor_variant,splice_acceptor_variant (这些是gnomAD定义的 pLoF)
    (2) --vep-filter stop_gained,frameshift_variant,splice_donor_variant,splice_acceptor_variant,start_lost,stop_lost,transcript_ablation 
        (这是我的定义，多包括了 tart_lost/stop_lost/transcript_ablation)

2. P/LP & VEP HIGH  

VEP IMPACT=HIGH 指变异对蛋白的预测影响(截短、移码、剪接位点破坏), 由算法自行计算的。

- ClinVar P/LP and HIGH      :  如致病性无义突变、移码缺失
- ClinVar P/LP but not HIGH  :  如致病性错义突变(MODERATE)、致病性同义突变(LOW, 影响剪接)
- HIGH but not ClinVar P/LP  :  如新发 pLoF(未收录)、良性截短(末外显子，不触发 NMD)

"""
import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import pysam

AlleleKey = Tuple[str, int, str, str]  # (chrom_no_chr, pos, ref_upper, alt_upper)
RegionSpec = Tuple[str, int, int, Optional[str]]  # (chrom_with_chr, 1based_start, 1based_end, gene_filter)
CLINVAR_ID_FIELDS = ("CLNVI", "CLNVAR_ID", "CLNVARID", "CLNVID", "CLNDN_ID", "CLNVAR")
BENIGN_TERMS = {"benign", "likely_benign", "benign/likely_benign"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="General-purpose carrier rate calculator from paired VEP + GT VCFs.")
    p.add_argument("--vep-vcf", required=True, type=Path, help="Sites-only VEP-annotated VCF/BCF.")
    p.add_argument("--gt-vcf", required=True, type=Path, help="Multi-sample genotype VCF/BCF.")
    p.add_argument("-o", "--output", required=True, type=Path, help="Output variant-level summary TSV.")
    p.add_argument("-s", "--samples", type=str, default=None, help="Sample list file (1 ID/line) or comma-separated IDs.")
    p.add_argument("--carriers-out", type=Path, default=None, help="Optional sample-level carrier details TSV.")

    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--sites-file", type=Path, help="TSV/CSV file with columns: chrom,pos,ref,alt[,variant_id,hgvsc,hgvsp,phenotype_group].")
    mode.add_argument("--region", type=str, help="Comma-separated genomic regions (e.g. chr10:87863625-87971930).")
    mode.add_argument("--bed", type=Path, help="BED file (0-based half-open; optional 4th column = gene SYMBOL).")

    p.add_argument("--gene", type=str, default=None, help="Optional VEP SYMBOL filter (e.g. PTEN) for --region/--sites-file.")
    p.add_argument("--vep-filter", type=str, default=None,
                   help="Region/BED mode: comma-separated terms matched against VEP Consequence, IMPACT, "
                        "and CLIN_SIG (case-insensitive). Reserved keyword 'clinvar' matches any variant "
                        "with ClinVar annotation. E.g. 'clinvar', 'pathogenic,likely_pathogenic,HIGH', "
                        "'missense_variant,stop_gained'. When omitted, ALL annotated variants are included.")
    p.add_argument("--exclude-conflicting", action="store_true", help="Region/BED mode: exclude variants with 'conflicting' in CLIN_SIG.")
    p.add_argument("--min-gq", type=int, default=0, help="Min GQ threshold (default: 0).")
    p.add_argument("--min-dp", type=int, default=0, help="Min DP threshold (default: 0).")
    p.add_argument("--pass-only", action="store_true", help="Keep only FILTER=PASS records in GT VCF.")
    return p.parse_args()


def clean_chr(c: str) -> str:
    return c[3:] if c.startswith("chr") else c


def fmt_chr(c: str, hdr: pysam.VariantHeader) -> str:
    has_chr = any(x.startswith("chr") for x in hdr.contigs)
    base = clean_chr(c)
    return f"chr{base}" if has_chr else base


def make_key(c: str, pos: int, ref: str, alt: str) -> AlleleKey:
    return (clean_chr(c), pos, ref.upper(), alt.upper())


def parse_region_token(token: str, default_gene: Optional[str]) -> RegionSpec:
    m = re.match(r"^(?:chr)?(.+?):(\d+)[-\u2013:](\d+)$", token.strip())
    if not m:
        sys.exit(f"[ERROR] Cannot parse region '{token}'. Expected 'chr10:start-end'.")
    chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
    if end < start:
        start, end = end, start
    return (f"chr{chrom}", start, end, default_gene)


def load_regions(args: argparse.Namespace) -> List[RegionSpec]:
    if args.region:
        return [parse_region_token(tok, args.gene) 
                for tok in args.region.split(",") if tok.strip()]

    regions: List[RegionSpec] = []
    with args.bed.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            
            parts = line.split("\t")
            if len(parts) < 3:
                continue

            chrom = f"chr{clean_chr(parts)}"
            start_1based = int(parts) + 1
            end_1based = int(parts)
            gene = parts.strip() if len(parts) >= 4 and parts.strip() not in (".", "") else args.gene
            regions.append((chrom, start_1based, end_1based, gene))

    if not regions:
        sys.exit(f"[ERROR] No valid regions parsed from BED file: {args.bed}")
    
    return regions


def load_sites_table(sites_path: Path) -> List[Dict[str, Any]]:
    delim = "," if sites_path.suffix.lower() == ".csv" else "\t"
    targets: List[Dict[str, Any]] = []
    used_ids: Set[str] = set()
    with sites_path.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        if not reader.fieldnames or not {"chrom", "pos", "ref", "alt"}.issubset(set(reader.fieldnames)):
            sys.exit("[ERROR] --sites-file must contain header columns: chrom, pos, ref, alt")
        
        for row in reader:
            chrom = f"chr{clean_chr(row['chrom'].strip())}"
            pos = int(row["pos"].strip())
            ref = row["ref"].strip().upper()
            alt = row["alt"].strip().upper()
            raw_vid = (row.get("variant_id") or "").strip() or f"{clean_chr(chrom)}:{pos}:{ref}:{alt}"
            vid = raw_vid if raw_vid not in used_ids else f"{raw_vid}@{chrom}:{pos}{ref}>{alt}"
            used_ids.add(vid)
            targets.append({
                "variant_id": vid,
                "gene": (row.get("gene") or "NA").strip() or "NA",
                "rsid": (row.get("rsid") or "NA").strip() or "NA",
                "clinvar_id": (row.get("clinvar_id") or "NA").strip() or "NA",
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "hgvsc": (row.get("hgvsc") or "").strip(),
                "hgvsp": (row.get("hgvsp") or "").strip(),
                "exon": (row.get("exon") or "NA").strip() or "NA",
                "variant_tier": "USER_SITE",
                "phenotype_group": (row.get("phenotype_group") or "User Target").strip() or "User Target"
            })
    
    if not targets:
        sys.exit(f"[ERROR] No variants loaded from {sites_path}")
    
    return targets


def load_target_samples(samples_arg: Optional[str], vcf_samples: List[str]) -> List[str]:
    if not samples_arg:
        return vcf_samples

    sample_path = Path(samples_arg)
    if sample_path.is_file():
        requested = [line.strip().split()[0] for line in sample_path.open("r", encoding="utf-8") if line.strip() and not line.startswith("#")]
    else:
        requested = [s.strip() for s in samples_arg.split(",") if s.strip()]

    vcf_set = set(vcf_samples)
    valid = list(dict.fromkeys(s for s in requested if s in vcf_set))
    missing = len(set(requested)) - len(valid)
    if missing > 0:
        print(f"[WARNING] {missing} requested sample(s) not found in --gt-vcf.", file=sys.stderr)
    if not valid:
        sys.exit("[ERROR] None of the specified samples exist in --gt-vcf.")
    
    return valid


def get_csq_fields(hdr: pysam.VariantHeader) -> List[str]:
    if "CSQ" not in hdr.info or "Format:" not in hdr.info["CSQ"].description:
        return []

    return hdr.info["CSQ"].description.split("Format:")[1].strip().strip('"').split("|")


def parse_vep(rec: pysam.VariantRecord, alt_idx: int, fields: List[str], target_gene: Optional[str] = None) -> Dict[str, str]:
    if not fields or "CSQ" not in rec.info:
        return {}
    raw = rec.info["CSQ"]
    entries = [raw] if isinstance(raw, str) else list(raw)
    hits: List[Dict[str, str]] = []
    for e in entries:
        ann = dict(zip(fields, e.split("|")))
        anum = ann.get("ALLELE_NUM", "")
        if anum.isdigit() and int(anum) != alt_idx:
            continue
        if target_gene and ann.get("SYMBOL", "").upper() != target_gene.upper():
            continue
        hits.append(ann)
    for ann in hits:
        if ann.get("MANE_SELECT"):
            return ann
    for ann in hits:
        if ann.get("CANONICAL") == "YES":
            return ann
    return hits if hits else {}


def extract_clinvar_id(ann: Dict[str, str]) -> str:
    for f in CLINVAR_ID_FIELDS:
        v = ann.get(f)
        if v:
            return re.split(r"[&;,]", v).strip()
    return ""


def extract_rsid(ann: Dict[str, str]) -> str:
    ex = ann.get("existing_variation") or ann.get("Existing_variation") or ann.get("RSID") or ""
    for tok in re.split(r"[&;,]", ex):
        tok = tok.strip()
        if tok.lower().startswith("rs") and tok[2:].isdigit():
            return tok
    return ""


def clin_sig_terms(ann: Dict[str, str]) -> List[str]:
    raw = (ann.get("CLIN_SIG") or "").strip()
    return [t.strip() for t in re.split(r"[&;,]", raw) if t.strip()] if raw else []


def classify_region_variant(
    ann: Dict[str, str],
    vep_filter: Optional[Set[str]],
    exclude_conflicting: bool,
) -> Optional[str]:
    """Return tier label if variant passes filters, else None.

    vep_filter: lowercase terms matched against VEP Consequence, IMPACT,
    and CLIN_SIG fields.  The reserved keyword 'clinvar' matches any variant
    with a ClinVar ID or CLIN_SIG annotation.  When None, all annotated
    variants are included.

    Tier labels (based on variant properties, independent of filter):
      TIER1_CLINVAR          - has ClinVar ID or CLIN_SIG annotation
      TIER2_HIGH_IMPACT_PLOF - VEP IMPACT=HIGH and not benign-only
      ALL_ANNOTATED          - other annotated variants
    """
    if not ann:
        return None

    terms = clin_sig_terms(ann)
    terms_lower = [t.lower() for t in terms]
    has_conflicting = any("conflicting" in t for t in terms_lower)
    is_benign_only = bool(terms_lower) and all(t in BENIGN_TERMS for t in terms_lower)

    # Hard exclusion: conflicting ClinVar interpretations
    if exclude_conflicting and has_conflicting:
        return None

    # Inclusion filter: match against Consequence, IMPACT, CLIN_SIG
    if vep_filter:
        csq_terms = {t.strip().lower() for t in re.split(r"[&;,]", ann.get("Consequence", "")) if t.strip()}
        impact = (ann.get("IMPACT") or "").lower()
        matched = (
            bool(csq_terms & vep_filter)
            or impact in vep_filter
            or bool(set(terms_lower) & vep_filter)
            or ("clinvar" in vep_filter and bool(extract_clinvar_id(ann) or terms))
        )
        if not matched:
            return None

    # Tier labeling (priority: TIER1 > TIER2 > ALL_ANNOTATED)
    if extract_clinvar_id(ann) or terms:
        return "TIER1_CLINVAR"
    if ann.get("IMPACT", "").upper() == "HIGH" and not is_benign_only:
        return "TIER2_HIGH_IMPACT_PLOF"
    
    return "ALL_ANNOTATED"


def build_dynamic_target(akey: AlleleKey, ann: Dict[str, str], tier: str, used_ids: Set[str]) -> Dict[str, Any]:
    chrom, pos, ref, alt = akey
    rsid = extract_rsid(ann)
    clinvar_id = extract_clinvar_id(ann)
    sig = "&".join(clin_sig_terms(ann)) or "NA"
    exon = ann.get("EXON", "") or "NA"
    hgvsc = (ann.get("HGVSc", "").split(":")[-1]) or "NA"
    hgvsp = (ann.get("HGVSp", "").replace("%3D", "=").split(":")[-1]) or "NA"
    gene = ann.get("SYMBOL") or "NA"

    base_vid = rsid if rsid else f"{chrom}:{pos}:{ref}:{alt}"
    vid = base_vid if base_vid not in used_ids else f"{base_vid}@{chrom}:{pos}{ref}>{alt}"
    used_ids.add(vid)

    return {
        "variant_id": vid, 
        "gene": gene, 
        "rsid": rsid or "NA", 
        "clinvar_id": clinvar_id or "NA",
        "chrom": f"chr{chrom}", 
        "pos": pos, 
        "ref": ref, 
        "alt": alt,
        "hgvsc": hgvsc, 
        "hgvsp": hgvsp, 
        "exon": (f"Exon {exon}" if exon != "NA" else "NA"),
        "variant_tier": tier, 
        "phenotype_group": f"CLIN_SIG:{sig}",
    }


def match_user_site(akey: AlleleKey, vep: Dict[str, str], targets: List[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    r_chr, r_pos, r_ref, r_alt = akey
    hgvsc = vep.get("HGVSc", "").split(":")[-1]
    hgvsp = vep.get("HGVSp", "").replace("%3D", "=").split(":")[-1]
    for t in targets:
        if r_chr != clean_chr(t["chrom"]):
            continue
        if r_pos == t["pos"] and r_ref == t["ref"] and r_alt == t["alt"]:
            return t
        if (t["hgvsc"] and hgvsc == t["hgvsc"]) or (t["hgvsp"] and hgvsp == t["hgvsp"]):
            return t
    
    return None


def qc_pass(sdata: pysam.VariantRecordSample, min_gq: int, min_dp: int) -> bool:
    gt = sdata.get("GT")
    if not gt or any(a is None for a in gt):
        return False
    if min_gq > 0 and (sdata.get("GQ") is None or sdata.get("GQ") < min_gq):
        return False
    if min_dp > 0 and (sdata.get("DP") is None or sdata.get("DP") < min_dp):
        return False
    return True


def fetch_region(vcf: pysam.VariantFile, chrom: str, start: int, end: int):
    try:
        return vcf.fetch(chrom, max(0, start - 1), end)
    except ValueError:
        return vcf.fetch()


def main() -> None:
    args = parse_args()
    vep_filter = (
        {t.strip().lower() for t in args.vep_filter.split(",") if t.strip()}
        if args.vep_filter else None
    )

    vep_vcf = pysam.VariantFile(str(args.vep_vcf))
    csq_fields = get_csq_fields(vep_vcf.header)
    a2vep: Dict[AlleleKey, Dict[str, str]] = {}
    a2tgt: Dict[AlleleKey, Dict[str, Any]] = {}

    if args.sites_file:
        target_list = load_sites_table(args.sites_file)
        scan_regions: List[RegionSpec] = [
            (t["chrom"], max(1, t["pos"] - 50), t["pos"] + 50, args.gene) for t in target_list
        ]
        for r_chr, r_start, r_end, r_gene in scan_regions:
            for rec in fetch_region(vep_vcf, fmt_chr(r_chr, vep_vcf.header), r_start, r_end):
                if not rec.alts:
                    continue
                for idx, alt in enumerate(rec.alts, start=1):
                    akey = make_key(rec.chrom, rec.pos, rec.ref, alt)
                    ann = parse_vep(rec, idx, csq_fields, target_gene=r_gene)
                    tgt = match_user_site(akey, ann, target_list)
                    if tgt:
                        a2vep[akey] = ann
                        a2tgt[akey] = tgt
    else:
        scan_regions = load_regions(args)
        used_ids: Set[str] = set()
        for r_chr, r_start, r_end, r_gene in scan_regions:
            for rec in fetch_region(vep_vcf, fmt_chr(r_chr, vep_vcf.header), r_start, r_end):
                if not rec.alts:
                    continue
                for idx, alt in enumerate(rec.alts, start=1):
                    akey = make_key(rec.chrom, rec.pos, rec.ref, alt)
                    if akey in a2tgt:
                        continue
                    ann = parse_vep(rec, idx, csq_fields, target_gene=r_gene)
                    tier = classify_region_variant(
                        ann, vep_filter, args.exclude_conflicting,
                    )
                    if tier:
                        a2vep[akey] = ann
                        a2tgt[akey] = build_dynamic_target(akey, ann, tier, used_ids)
        
        target_list = sorted(a2tgt.values(), key=lambda t: (clean_chr(t["chrom"]), t["pos"], t["ref"], t["alt"]))

    vep_vcf.close()

    gt_vcf = pysam.VariantFile(str(args.gt_vcf))
    samples = load_target_samples(args.samples, list(gt_vcf.header.samples))
    if args.samples:
        gt_vcf.subset_samples(samples)
    n_samples = len(samples)

    stats: Dict[str, Dict[str, Any]] = {
        t["variant_id"]: {
            "in_vep": False, "in_gt": False, "pos": f"{t['chrom']}:{t['pos']}",
            "alleles": f"{t['ref']}>{t['alt']}", "called": 0, "het": 0, "hom": 0, "hemi": 0,
            "ac": 0, "an": 0, "csq": "NA", "imp": "NA", "clin": "NA",
            "hgvsc": t["hgvsc"] or "NA", "hgvsp": t["hgvsp"] or "NA", "af": "NA",
            "gene": t["gene"],
        }
        for t in target_list
    }

    for akey, t in a2tgt.items():
        st = stats[t["variant_id"]]
        st["in_vep"] = True
        ann = a2vep.get(akey, {})
        if ann:
            st["gene"] = ann.get("SYMBOL") or st["gene"]
            st["csq"] = ann.get("Consequence") or "NA"
            st["imp"] = ann.get("IMPACT") or "NA"
            st["clin"] = ann.get("CLIN_SIG") or "NA"
            st["hgvsc"] = ann.get("HGVSc") or st["hgvsc"]
            st["hgvsp"] = ann.get("HGVSp", "").replace("%3D", "=") or st["hgvsp"]
            for k in ("gnomADe_AF", "gnomADg_AF", "gnomAD_AF", "AF"):
                if ann.get(k):
                    st["af"] = ann[k]
                    break


    uniq_carriers: Set[str] = set()
    tier1_carriers: Set[str] = set()
    tier2_carriers: Set[str] = set()
    all_annotated_carriers: Set[str] = set()
    seen_records: Set[Tuple[str, int, str, str]] = set()
    carrier_rows: List[Dict[str, str]] = []

    for r_chr, r_start, r_end, _ in scan_regions:
        for rec in fetch_region(gt_vcf, fmt_chr(r_chr, gt_vcf.header), r_start, r_end):
            if args.pass_only and list(rec.filter.keys()) not in ([], ["PASS"]):
                continue
            if not rec.alts:
                continue
            for idx, alt in enumerate(rec.alts, start=1):
                akey = make_key(rec.chrom, rec.pos, rec.ref, alt)
                if akey in seen_records:
                    continue
                tgt = a2tgt.get(akey) or (match_user_site(akey, {}, target_list) if args.sites_file else None)
                if not tgt or tgt["variant_id"] not in stats:
                    continue
                seen_records.add(akey)
                vid = tgt["variant_id"]
                st = stats[vid]
                st["in_gt"] = True
                st["pos"] = f"{rec.chrom}:{rec.pos}"
                st["alleles"] = f"{rec.ref}>{alt}"
                for s in samples:
                    sd = rec.samples[s]
                    if not qc_pass(sd, args.min_gq, args.min_dp):
                        continue
                    gt = sd["GT"]
                    ploidy = len(gt)
                    copies = sum(1 for a in gt if a == idx)
                    st["called"] += 1
                    st["an"] += ploidy
                    st["ac"] += copies
                    if copies > 0:
                        uniq_carriers.add(s)
                        if tgt["variant_tier"] == "TIER1_CLINVAR":
                            tier1_carriers.add(s)
                        elif tgt["variant_tier"] == "TIER2_HIGH_IMPACT_PLOF":
                            tier2_carriers.add(s)
                        elif tgt["variant_tier"] == "ALL_ANNOTATED":
                            all_annotated_carriers.add(s)
                        if ploidy == 1:
                            zyg = "HEMI"
                            st["hemi"] += 1
                        elif copies == ploidy:
                            zyg = "HOM_ALT"
                            st["hom"] += 1
                        else:
                            zyg = "HET"
                            st["het"] += 1
                        carrier_rows.append({
                            "sample_id": s, "gene": st["gene"], "variant_id": vid,
                            "variant_tier": tgt["variant_tier"], "rsid": tgt["rsid"],
                            "clinvar_id": tgt["clinvar_id"], "zygosity": zyg,
                            "gt": "/".join(map(str, gt)), "gq": str(sd.get("GQ", "NA")),
                            "dp": str(sd.get("DP", "NA")), "phenotype_group": tgt["phenotype_group"],
                            "hgvsc": st["hgvsc"], "hgvsp": st["hgvsp"], "clin_sig": st["clin"],
                        })
    gt_vcf.close()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    cols = [
        "gene", "variant_id", "variant_tier", "rsid", "clinvar_id", "exon",
        "phenotype_group", "vcf_pos", "vcf_alleles", "found_in_vep", "found_in_gt",
        "total_samples", "called_samples", "call_rate", "het_carriers",
        "hom_alt_carriers", "hemi_carriers", "total_carriers", "carrier_rate",
        "carrier_rate_pct", "alt_ac", "total_an", "cohort_af", "vep_hgvsc",
        "vep_hgvsp", "vep_consequence", "vep_impact", "vep_clin_sig", "vep_external_af",
    ]
    with args.output.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for t in target_list:
            vid = t["variant_id"]
            st = stats[vid]
            called = st["called"]
            tot_c = st["het"] + st["hom"] + st["hemi"]
            crate = tot_c / called if called else 0.0
            af = st["ac"] / st["an"] if st["an"] else 0.0
            w.writerow({
                "gene": st["gene"], 
                "variant_id": vid, 
                "variant_tier": t["variant_tier"],
                "rsid": t["rsid"], 
                "clinvar_id": t["clinvar_id"], 
                "exon": t["exon"],
                "phenotype_group": t["phenotype_group"], 
                "vcf_pos": st["pos"],
                "vcf_alleles": st["alleles"], 
                "found_in_vep": st["in_vep"],
                "found_in_gt": st["in_gt"], 
                "total_samples": n_samples,
                "called_samples": called, 
                "call_rate": f"{called / n_samples:.4f}",
                "het_carriers": st["het"], 
                "hom_alt_carriers": st["hom"],
                "hemi_carriers": st["hemi"], 
                "total_carriers": tot_c,
                "carrier_rate": f"{crate:.6f}", 
                "carrier_rate_pct": f"{crate * 100:.4f}%",
                "alt_ac": st["ac"], 
                "total_an": st["an"], 
                "cohort_af": f"{af:.6f}",
                "vep_hgvsc": st["hgvsc"], 
                "vep_hgvsp": st["hgvsp"],
                "vep_consequence": st["csq"], 
                "vep_impact": st["imp"],
                "vep_clin_sig": st["clin"], 
                "vep_external_af": st["af"],
            })

    if args.carriers_out:
        args.carriers_out.parent.mkdir(parents=True, exist_ok=True)
        c_cols = [
            "sample_id", "gene", "variant_id", "variant_tier", "rsid", "clinvar_id",
            "zygosity", "gt", "gq", "dp", "phenotype_group", "hgvsc", "hgvsp", "clin_sig",
        ]
        with args.carriers_out.open("w", newline="", encoding="utf-8") as cfh:
            cw = csv.DictWriter(cfh, fieldnames=c_cols, delimiter="\t")
            cw.writeheader()
            cw.writerows(carrier_rows)

    print(f"[DONE] Total unique carriers: {len(uniq_carriers)}/{n_samples} ({len(uniq_carriers)/n_samples*100:.4f}%)")
    if not args.sites_file:
        print(f"       Tier 1 (ClinVar) carriers: {len(tier1_carriers)}/{n_samples} ({len(tier1_carriers)/n_samples*100:.4f}%)")
        print(f"       Tier 2 (HIGH impact pLoF) carriers: {len(tier2_carriers)}/{n_samples} ({len(tier2_carriers)/n_samples*100:.4f}%)")
        print(f"       Other annotated carriers: {len(all_annotated_carriers)}/{n_samples} ({len(all_annotated_carriers)/n_samples*100:.4f}%)")


if __name__ == "__main__":
    main()


