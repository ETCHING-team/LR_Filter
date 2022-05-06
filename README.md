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

If you want to install required modules, refer to the command below
<pre>
<code>
pip install {"module_name"}
</code>
</pre>

# Usages
python LR_Filter.py. -h

<pre>
<code>
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
  -c CPUS,           --cpus CPUS, default: 1 cpu
                        Set the number of CPUs, this option is for multi-processing
  -tr TARGET_RANGE, --target_range TARGET_RANGE, default: 500bp
                        Set the range of the SV BreakPoints (BPs) for verifying by long-reads
  -o OUTNAME, --outname OUTNAME
                        Output name

</code>
</pre>

# Example
LR_Filter requires long-read BAM file as primary input and two modes are currently available.
1. Somatic SV mode: this mode requires both tumor-normal long-read mapped BAM files
2. General SV mode: in this mode, only a single long-read mapped BAM file is required

## Somatic SV mode
<pre>
<code>

python LR_Filter -i {input_VCF_file} -t {DEL} -lbt {tumor_long_read_sorted_BAM_file} -lbn {normal_long_read_sorted_BAM_file} \
                 -rf {hg19.fa} -c {1} -tr {500} -o {test}

</code>
</pre>

## General SV mode
<pre>
<code>

python LR_Filter -i {input_VCF_file} -t {DEL} -lbt {tumor_long_read_sorted_BAM_file} -rf {hg19.fa} -c {1} -tr {500} -o {test}

</code>
</pre>
