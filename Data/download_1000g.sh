#!/bin/bash
# 1000 Genomes phase 3 (build 37) VCFs and sample panel, used as the LD reference of
# Data/prepare_mediation_ldl_cond.R (about 15 GB). Run from the repository root:
#   nohup bash Data/download_1000g.sh > Data/raw/g1000/download.log 2>&1 &
# Interrupted downloads resume (wget -c). Each file is checked with gzip -t at the end.
set -e
cd "$(dirname "$0")"
D=raw/g1000; mkdir -p $D
U=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502
wget -nv -c -O $D/integrated_call_samples_v3.20130502.ALL.panel $U/integrated_call_samples_v3.20130502.ALL.panel
for k in $(seq 1 22); do
  f=ALL.chr$k.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
  wget -nv -c -O $D/$f $U/$f
  wget -nv -c -O $D/$f.tbi $U/$f.tbi
done
for k in $(seq 1 22); do
  gzip -t $D/ALL.chr$k.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz && echo "chr$k OK"
done
echo "DOWNLOAD DONE"
