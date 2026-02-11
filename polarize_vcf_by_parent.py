#!/usr/bin/python
# -*- coding: UTF-8 -*-

"""
================================================================================
PARENT-BASED REF/ALT SWITCHING & VCF POLARIZATION SCRIPT
================================================================================
Author: Tahir Ali | https://tahirali-biomics.github.io/
Date: 11.02.2026

Adapted from: Kristian Ullrich's polarizeVCFbyOutgroup.py (2021)
             https://github.com/kullrich/bio-scripts/blob/adfef64a803ca2104c60e68b29e241f6b4776405/vcf/polarizeVCFbyOutgroup.py
Original license: MIT

Modifications and improvements:
    • Added comprehensive FORMAT field flipping (GT, PL, AD, GP)
    • Fixed missing genotype handling (./. sites now properly removed)
    • Enhanced documentation for outgroup indexing
    • Added detailed processing statistics
    • Improved error handling and debugging output
    • Added parent-based polarization workflow description

This script is released under the MIT License (maintained from original work).

PURPOSE:
    This script polarizes a VCF file using a specified outgroup/parent individual.
    When the outgroup is homozygous ALT (1/1), it swaps REF and ALT alleles and
    flips all relevant FORMAT fields across ALL samples. This ensures the outgroup
    becomes homozygous REF (ancestral state) after flipping.

WORKFLOW SUMMARY:
    ┌─────────────────────────────────────────────────────────────────────┐
    │  1. Read VCF line by line                                           │
    │  2. Skip header lines (add AA tag if -add specified)                │
    │  3. Extract outgroup genotype (GT field)                            │
    │  4. Apply filtering logic:                                          │
    │     • Missing GT (./.) → REMOVE site                                │
    │     • Homozygous REF (0/0) → KEEP site (no change)                  │
    │     • Heterozygous (0/1,1/0) → REMOVE site (unless -keep)           │
    │     • Homozygous ALT (1/1) → SWITCH alleles & FLIP all samples      │
    │  5. When switching:                                                 │
    │     • Swap REF and ALT columns                                      │
    │     • Flip GT: 0↔1 for all genotypes                                │
    │     • Flip PL: swap 0/0 and 1/1 likelihoods                         │
    │     • Flip AD: swap REF and ALT read counts                         │
    │     • Flip GP: swap 0/0 and 1/1 posterior probabilities             │
    │     • DP, SP, GQ: unchanged (allele-independent)                    │
    │  6. Add Ancestral Allele (AA) tag if -add specified                 │
    │  7. Write processed line to output                                  │
    └─────────────────────────────────────────────────────────────────────┘

KEY FEATURES:
    • Handles phased (|) and unphased (/) genotypes
    • Properly flips all FORMAT fields: GT, PL, AD, GP
    • Preserves DP, SP, GQ (invariant to allele switching)
    • Works with gzipped or plain text VCF
    • Adds ancestral allele info to INFO field (optional)
    • Keeps heterozygous outgroup sites (optional)

FILTERING BEHAVIOR (based on outgroup genotype):
    ┌─────────────┬─────────────────────────┬──────────────────────────────┐
    │ Outgroup GT │      Default (-keep)    │      With -keep             │
    ├─────────────┼─────────────────────────┼──────────────────────────────┤
    │    ./., .   │        REMOVE           │         REMOVE              │
    │    0/0      │         KEEP            │          KEEP               │
    │    0/1,1/0  │        REMOVE           │          KEEP               │
    │    1/1      │    SWITCH & KEEP        │     SWITCH & KEEP           │
    └─────────────┴─────────────────────────┴──────────────────────────────┘

OUTGROUP INDEXING (CRITICAL!):
    The -ind parameter expects a 0-based SAMPLE INDEX (position after FORMAT),
    NOT the absolute column number.
    
    CORRECT: -ind 80 (80th sample after FORMAT column)
    WRONG:   -ind 88 (absolute column position)
    
    TO FIND CORRECT INDEX:
        zgrep "^#CHROM" your.vcf.gz | head -1 | tr '\t' '\n' | nl -v 0 | grep "OutgroupName"

OUTPUT STATISTICS:
    • Parsed:     Total sites processed
    • Removed:    Sites with missing outgroup data
    • Kept:       Homozygous REF and flipped sites
    • Removed:    Heterozygous sites (unless -keep)
    • Switched:   Sites where REF/ALT were swapped
    • Total written: Sites in output VCF

EXAMPLE USAGE:
    # Basic polarization (remove missing/heterozygous, add AA tag)
    python polarizeVCF.py -vcf input.vcf.gz -out polarized.vcf.gz -ind 80 -add
    
    # Keep heterozygous outgroup sites
    python polarizeVCF.py -vcf input.vcf -out polarized.vcf -ind 80 -keep -add
    
    # No ancestral allele annotation
    python polarizeVCF.py -vcf input.vcf.gz -out polarized.vcf.gz -ind 80

NOTE:
    This script is designed for BI-ALLELIC SNPs only. Multi-allelic sites
    and indels should be filtered out before running.
================================================================================
"""

import sys
import argparse
import textwrap
import gzip
import re


def multiple_replace(string, rep_dict):
    """Replace multiple patterns in string using dictionary mapping."""
    pattern = re.compile("|".join([re.escape(k) for k in sorted(rep_dict, key=len, reverse=True)]), flags=re.DOTALL)
    return pattern.sub(lambda x: rep_dict[x.group(0)], string)


def flip_genotype(gt):
    """
    Flip genotype 0↔1 while preserving phasing.
    
    Examples:
        0/0 → 1/1, 0|0 → 1|1
        1/1 → 0/0, 1|1 → 0|0
        0/1 → 1/0, 0|1 → 1|0
        1/0 → 0/1, 1|0 → 0|1
        ./. → ./. (unchanged)
    """
    if gt == '0/0' or gt == '0|0':
        return '1/1' if '/' in gt else '1|1'
    elif gt == '1/1' or gt == '1|1':
        return '0/0' if '/' in gt else '0|0'
    elif gt == '0/1' or gt == '0|1':
        return '1/0' if '/' in gt else '1|0'
    elif gt == '1/0' or gt == '1|0':
        return '0/1' if '/' in gt else '0|1'
    else:
        return gt  # ./., .|., etc.


def flip_pl_values(pl_str):
    """
    Flip PL (Phred-scaled genotype likelihood) values.
    
    FORMAT: PL[0]=0/0, PL[1]=0/1, PL[2]=1/1
    After swap: PL[0]=1/1, PL[1]=0/1, PL[2]=0/0
    """
    if pl_str == '.':
        return pl_str
    pl_values = pl_str.split(',')
    if len(pl_values) == 3:  # For bi-allelic sites
        # Swap PL[0] (0/0) with PL[2] (1/1), keep PL[1] (0/1) the same
        return ','.join([pl_values[2], pl_values[1], pl_values[0]])
    return pl_str


def flip_ad_values(ad_str):
    """
    Flip AD (Allelic depths) values.
    
    FORMAT: AD[0]=REF reads, AD[1]=ALT reads
    After swap: AD[0]=ALT reads, AD[1]=REF reads
    """
    if ad_str == '.':
        return ad_str
    ad_values = ad_str.split(',')
    if len(ad_values) == 2:  # For bi-allelic sites
        # Swap AD[0] (REF) with AD[1] (ALT)
        return ','.join([ad_values[1], ad_values[0]])
    return ad_str


def flip_gp_values(gp_str):
    """
    Flip GP (Genotype posterior probability) values.
    
    FORMAT: GP[0]=0/0, GP[1]=0/1, GP[2]=1/1
    After swap: GP[0]=1/1, GP[1]=0/1, GP[2]=0/0
    """
    if gp_str == '.':
        return gp_str
    gp_values = gp_str.split(',')
    if len(gp_values) == 3:  # For bi-allelic sites
        # Swap GP[0] (0/0) with GP[2] (1/1), keep GP[1] (0/1) the same
        return ','.join([gp_values[2], gp_values[1], gp_values[0]])
    return gp_str


def flip_format_fields(format_str, sample_str):
    """
    Master function to flip all relevant FORMAT fields when alleles are swapped.
    
    Fields flipped:
        - GT: 0↔1 genotype swap
        - PL: 0/0 ↔ 1/1 likelihood swap
        - AD: REF ↔ ALT read count swap
        - GP: 0/0 ↔ 1/1 probability swap
    
    Fields unchanged (allele-independent):
        - DP: Total depth
        - SP: Strand bias p-value
        - GQ: Genotype quality
    """
    if format_str == '.' or sample_str == '.':
        return format_str, sample_str
    
    format_fields = format_str.split(':')
    sample_fields = sample_str.split(':')
    
    # Create a dictionary for easier access
    sample_dict = dict(zip(format_fields, sample_fields))
    
    # Flip each field as needed
    if 'GT' in sample_dict:
        sample_dict['GT'] = flip_genotype(sample_dict['GT'])
    
    if 'PL' in sample_dict:
        sample_dict['PL'] = flip_pl_values(sample_dict['PL'])
    
    if 'AD' in sample_dict:
        sample_dict['AD'] = flip_ad_values(sample_dict['AD'])
    
    if 'GP' in sample_dict:
        sample_dict['GP'] = flip_gp_values(sample_dict['GP'])
    
    # DP, SP, and GQ remain unchanged when swapping alleles
    # (they are invariant to which allele is REF vs ALT)
    
    # Reconstruct the sample string
    new_sample_fields = [sample_dict.get(field, sample_fields[i] if i < len(sample_fields) else '.') 
                         for i, field in enumerate(format_fields)]
    
    return format_str, ':'.join(new_sample_fields)


def parse_lines(fin, fou, ind, keep, add):
    """
    Main VCF parsing and processing function.
    
    Steps for each non-header line:
    1. Extract outgroup genotype
    2. Check for missing data → skip if found
    3. Case 0/0: Keep site (add AA tag if requested)
    4. Case 0/1,1/0: Skip unless -keep specified
    5. Case 1/1: Swap REF/ALT, flip all FORMAT fields for ALL samples
    6. Write processed line
    """
    switchcount = 0
    removecount = 0
    outmissing = 0
    totalcount = 0
    kept_heterozygous = 0
    kept_hom_ref = 0
    
    for line in fin:
        # ------------------- HEADER PROCESSING -------------------
        if line[0] == '#':
            if add:
                if line.split('\t')[0] == '#CHROM':
                    fou.write('##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">\n')
                    fou.write(line)
                else:
                    fou.write(line)
            else:
                fou.write(line)
            continue
        
        # ------------------- DATA LINE PROCESSING -------------------
        totalcount += 1
        linesplit = line.strip().split('\t')
        
        # Check if outgroup column exists
        if len(linesplit) <= 8 + ind:
            outmissing += 1
            continue
        
        # Extract outgroup information
        format_col = linesplit[8]
        sample_col = linesplit[8 + ind]
        
        # Get genotype from GT field
        gt_field = None
        if format_col != '.' and sample_col != '.':
            try:
                format_fields = format_col.split(':')
                if 'GT' in format_fields:
                    gt_index = format_fields.index('GT')
                    sample_fields = sample_col.split(':')
                    if gt_index < len(sample_fields):
                        gt_field = sample_fields[gt_index]
            except:
                pass
        
        # DEBUG: Uncomment to verify outgroup genotypes
        # print(f"DEBUG: Format: {format_col}, Sample: {sample_col}, GT: {gt_field}")
        
        # ------------------- MISSING DATA CHECK -------------------
        if gt_field is None or gt_field == '.' or gt_field.startswith('.') or '..' in gt_field:
            outmissing += 1
            continue  # Skip sites with missing outgroup data
        
        # ------------------- CASE 1: HOMOZYGOUS REF -------------------
        if gt_field == '0/0' or gt_field == '0|0':
            kept_hom_ref += 1
            if add:
                linesplit[7] = 'AA=' + linesplit[3] + ';' + linesplit[7]
                fou.write('\t'.join(linesplit) + '\n')
            else:
                fou.write(line)
        
        # ------------------- CASE 2: HETEROZYGOUS -------------------
        elif gt_field in ['0/1', '1/0', '0|1', '1|0']:
            removecount += 1
            if keep:
                kept_heterozygous += 1
                fou.write(line)
            # If not keep, skip heterozygous sites
        
        # ------------------- CASE 3: HOMOZYGOUS ALT - SWITCH NEEDED -------------------
        elif gt_field == '1/1' or gt_field == '1|1':
            switchcount += 1
            
            # Swap REF and ALT alleles
            REF = linesplit[4]
            ALT = linesplit[3]
            linesplit[3] = REF
            linesplit[4] = ALT
            
            # Flip FORMAT fields for ALL samples (not just outgroup)
            format_col = linesplit[8]
            for i in range(9, len(linesplit)):
                if linesplit[i] != '.':
                    format_col, linesplit[i] = flip_format_fields(format_col, linesplit[i])
            
            if add:
                linesplit[7] = 'AA=' + linesplit[3] + ';' + linesplit[7]
                fou.write('\t'.join(linesplit) + '\n')
            else:
                fou.write('\t'.join(linesplit) + '\n')
        
        # ------------------- CASE 4: OTHER (multi-allelic, etc.) -------------------
        else:
            removecount += 1
            if keep:
                fou.write(line)
            # If not keep, skip these sites
    
    # ------------------- FINAL STATISTICS -------------------
    print('\n' + '='*60)
    print('PROCESSING SUMMARY')
    print('='*60)
    print(f'Total sites parsed:        {totalcount:8d}')
    print(f'Removed (missing outgroup): {outmissing:8d}')
    print(f'Kept (homozygous REF):     {kept_hom_ref:8d}')
    if keep:
        print(f'Kept (heterozygous):       {kept_heterozygous:8d}')
        print(f'Removed (other):           {removecount - kept_heterozygous:8d}')
    else:
        print(f'Removed (heterozygous/other): {removecount:8d}')
    print(f'Switched (REF↔ALT):        {switchcount:8d}')
    print(f'Total sites written:       {kept_hom_ref + switchcount + (kept_heterozygous if keep else 0):8d}')
    print('='*60)


def main():
    parser = argparse.ArgumentParser(
        prog='polarizeVCFbyOutgroup',
        description='''\
        ██████╗ ███████╗███████╗    ███████╗██╗    ██╗██╗████████╗ ██████╗██╗  ██╗
        ██╔══██╗██╔════╝██╔════╝    ██╔════╝██║    ██║██║╚══██╔══╝██╔════╝██║  ██║
        ██████╔╝█████╗  █████╗      ███████╗██║ █╗ ██║██║   ██║   ██║     ███████║
        ██╔═══╝ ██╔══╝  ██╔══╝      ╚════██║██║███╗██║██║   ██║   ██║     ██╔══██║
        ██║     ███████╗██║         ███████║╚███╔███╔╝██║   ██║   ╚██████╗██║  ██║
        ╚═╝     ╚══════╝╚═╝         ╚══════╝ ╚══╝╚══╝ ╚═╝   ╚═╝    ╚═════╝╚═╝  ╚═╝
        
        Parent-based REF/ALT switching and VCF polarization using outgroup genotype.
        
        ==============================================================================
        QUICK START:
        1. Find your outgroup index:
           zgrep "^#CHROM" input.vcf.gz | tr '\t' '\n' | nl -v 0 | grep "Outgroup"
        
        2. Run polarization:
           python %(prog)s -vcf input.vcf.gz -out polarized.vcf.gz -ind 80 -add
        ==============================================================================
        ''',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''\
        EXAMPLES:
          # Basic usage with gzipped files
          python polarizeVCFbyOutgroup.py -vcf variants.vcf.gz -out polarized.vcf.gz -ind 80 -add
          
          # Keep heterozygous outgroup sites
          python polarizeVCFbyOutgroup.py -vcf variants.vcf -out polarized.vcf -ind 80 -keep -add
          
          # No ancestral allele annotation
          python polarizeVCFbyOutgroup.py -vcf variants.vcf.gz -out polarized.vcf.gz -ind 80
        
        CITATION:
          Adapted from Kristian Ullrich (2021)
        '''
    )
    
    parser.add_argument('-vcf', required=True, help='Input VCF file (can be .gz)')
    parser.add_argument('-out', required=True, help='Output VCF file (can be .gz)')
    parser.add_argument('-ind', type=int, required=True,
                       help='0-based SAMPLE INDEX of outgroup/parent (NOT column number!)')
    parser.add_argument('-keep', action='store_true',
                       help='Keep sites where outgroup is heterozygous')
    parser.add_argument('-add', action='store_true',
                       help='Add Ancestral Allele (AA) tag to INFO field')
    
    args = parser.parse_args()
    
    print('\n' + '='*60)
    print('PARENT-BASED REF/ALT SWITCHING')
    print('='*60)
    print(f'Input VCF:  {args.vcf}')
    print(f'Output VCF: {args.out}')
    print(f'Outgroup index: {args.ind} (0-based sample position)')
    print(f'Keep heterozygous: {args.keep}')
    print(f'Add AA tag: {args.add}')
    print('='*60 + '\n')
    
    # Handle gzipped or plain text files
    if args.out.endswith('gz'):
        with gzip.open(args.out, 'wt') as fou:
            if args.vcf.endswith('gz'):
                with gzip.open(args.vcf, 'rt') as fin:
                    parse_lines(fin, fou, args.ind, args.keep, args.add)
            else:
                with open(args.vcf, 'rt') as fin:
                    parse_lines(fin, fou, args.ind, args.keep, args.add)
    else:
        with open(args.out, 'wt') as fou:
            if args.vcf.endswith('gz'):
                with gzip.open(args.vcf, 'rt') as fin:
                    parse_lines(fin, fou, args.ind, args.keep, args.add)
            else:
                with open(args.vcf, 'rt') as fin:
                    parse_lines(fin, fou, args.ind, args.keep, args.add)


if __name__ == '__main__':
    main()
