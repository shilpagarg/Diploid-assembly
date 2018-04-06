# Diploid-assembly

The pipeline describes the linear approach to generate haplotigs for every phased block.

Requirements: samtools, bwa, picard, seqtk, whatshap

Customize the paths to initial Illumina, PacBio and reference fasta files
 
Usage:
snakemake -p master
snakemake -p master2
