# LR_Filter
The **LR_Filter** provides information for filtering the false Structural Variations (SVs, VCF format) based on a long-read mapping file (BAM format).
This is developed for ETCHING project.

# Requirements
The **LR_Filer** was written by Python 3 a along with some **Python modules**. The required modules are as follows.

* sys
* argparse
* math
* re
* os
* pickle
* pysam
* numpy
* queue
* threading
* multiprocessing
  
# Usages

''' Python

usage: LR-Filter.v0.4.py [-h] [-i INPUT_VCF] [-t TARGET_TYPE] [-lbt LONG_BAM_T] [-lbn LONG_BAM_N] [-rf REFERENCE_SEQ] [-c CPUS] [-tr TARGET_RANGE] [-o OUTNAME]

LR_Filter is SV filtering tool using long-read BAM file

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_VCF, --input_VCF INPUT_VCF
                        The input VCF file by short-read SV caller, etc ETCHING, Lumpy, Delly ...
  -t TARGET_TYPE, --target_type TARGET_TYPE
                        The target SV type, etc DEL, DUP, INV, or TRA
  -lbt LONG_BAM_T, --long_bam_t LONG_BAM_T
                        A tumor long-read mapping sorted BAM file by Minimap2 or NGMLR
  -lbn LONG_BAM_N, --long_bam_n LONG_BAM_N
                        A normal long-read mapping sorted BAM file by Minimap2 or NGMLR
  -rf REFERENCE_SEQ, --reference_seq REFERENCE_SEQ
                        The reference fasta file, etc hg19 or hg38
  -c CPUS, --cpus CPUS  Set the number of CPUs, this option is for multi-processing
  -tr TARGET_RANGE, --target_range TARGET_RANGE
                        Set the range of SV BP for verifying by long-reads
  -o OUTNAME, --outname OUTNAME
                        Output name

'''
