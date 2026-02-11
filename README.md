# Omics Population Genomics Toolkit

**Author:** Dr Tahir Ali  
**Date:** 2025  
**Language:** Python 3  
**License:** MIT  

---

## Overview

A growing collection of Python scripts for population genomics and multi-omics data analysis.  
Designed for VCF manipulation, ancestral allele polarization, and population genetic statistics.

**Current scripts:** 1 | **Planned:** Many more 🔬

---

## Available Scripts

### 1. `polarize_vcf_by_parentOutgroup.py`

**Purpose:**  
Polarize VCF files using a known parent/outgroup individual. Swaps REF/ALT alleles and flips all FORMAT fields (GT, PL, AD, GP) when outgroup is homozygous ALT.

**Features:**
- ✅ Parent/outgroup-based ancestral allele determination  
- ✅ Full FORMAT field flipping (GT, PL, AD, GP)  
- ✅ Handles phased (`|`) and unphased (`/`) genotypes  
- ✅ Gzip support (input/output)  
- ✅ Optional ancestral allele (AA) tag  
- ✅ Optional retention of heterozygous outgroups  

**Usage:**

# ─────────────────────────────────────────────────────────────
# STEP 1: Find your outgroup/parent sample index
# ─────────────────────────────────────────────────────────────
# This returns: [INDEX]    SampleName
# Example:       80        Mil2
# Your -ind value is the number in the first column (0-based)
# ─────────────────────────────────────────────────────────────
zgrep "^#CHROM" input.vcf.gz | head -1 | tr '\t' '\n' | nl -v 0 | grep "Mil2"

# ─────────────────────────────────────────────────────────────
# STEP 2: Run polarisation (use the index from Step 1)
# ─────────────────────────────────────────────────────────────
python polarize_vcf_by_parent.py \
    -vcf input.vcf.gz \
    -out polarized.vcf.gz \
    -ind 80 \          # ← IMPORTANT: Replace 80 with your actual index!
    -add
