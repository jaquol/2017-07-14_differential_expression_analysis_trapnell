#!/bin/bash
#$ -N fastq_dump_from_sra_SRR493368_SRR493368
#$ -q long-sl65
#$ -l virtual_free=2G
#$ -l h_rt=48:00:00
#$ -o /users/GR/mb/jquilez/projects/misc/analysis/2017-07-14_differential_expression_analysis_trapnell/job_out/fastq_dump_from_sra_SRR493368_SRR493368_$JOB_ID.out
#$ -e /users/GR/mb/jquilez/projects/misc/analysis/2017-07-14_differential_expression_analysis_trapnell/job_out/fastq_dump_from_sra_SRR493368_SRR493368_$JOB_ID.err
#$ -j y
#$ -M javier.quilez@crg.eu
#$ -m abe
#$ -pe smp 1
/software/mb/bin/fastq-dump SRR493368 --split-files -O /users/GR/mb/jquilez/data/rnaseq/raw/2017-07-14 -DQ '+' --gzip
mv /users/GR/mb/jquilez/data/rnaseq/raw/2017-07-14/SRR493368_1.fastq.gz /users/GR/mb/jquilez/data/rnaseq/raw/2017-07-14/SRR493368_read1.fastq.gz
mv /users/GR/mb/jquilez/data/rnaseq/raw/2017-07-14/SRR493368_2.fastq.gz /users/GR/mb/jquilez/data/rnaseq/raw/2017-07-14/SRR493368_read2.fastq.gz
