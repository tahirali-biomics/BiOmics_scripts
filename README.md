# Population Genomics and Evolutionary Omics Toolkit

**Author:** Dr Tahir Ali  
**Date:** 2024  
**Languages:** Python, Bash, Perl, R  
**License:** MIT  

---

## Overview

A growing collection of multi-language scripts for population genomics and multi-omics data analysis.  
Designed for VCF manipulation, ancestral allele polarisation, population genetic statistics, and downstream visualisation.

**Current scripts:** 1 | **Planned:** Many more

**Supported languages:** 🐍 Python | 🐚 Bash | 🦈 Perl | 📊 R
---

## Available Scripts

### 1. `polarize_vcf_by_parent.py` (Python)

**Purpose:**  
Polarise VCF files using a known parent/outgroup individual. Swaps REF/ALT alleles and flips all FORMAT fields (GT, PL, AD, GP) when outgroup is homozygous ALT.

**Features:**
- ✅ Parent/outgroup-based ancestral allele determination  
- ✅ Full FORMAT field flipping (GT, PL, AD, GP)  
- ✅ Handles phased (`|`) and unphased (`/`) genotypes  
- ✅ Gzip support (input/output)  
- ✅ Optional ancestral allele (AA) tag  
- ✅ Optional retention of heterozygous outgroups  

**Usage:**

```bash
# STEP 1: Find your outgroup/parent sample index (0-based)
zgrep "^#CHROM" input.vcf.gz | head -1 | tr '\t' '\n' | nl -v 0 | grep "Mil2"

# STEP 2: Run polarisation (adjust -ind value from Step 1)
python polarize_vcf_by_parent.py \
    -vcf input.vcf.gz \
    -out polarized.vcf.gz \
    -ind 80 \          # ← Replace 80 with your outgroup index from Step 1
    -add
