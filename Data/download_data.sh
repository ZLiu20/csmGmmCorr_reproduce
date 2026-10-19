#!/bin/bash
# Download the public summary statistics used in the applications (about 25 GB in total).
# Run from the repository root on a node with internet access:   bash Data/download_data.sh
#   GLGC (Graham et al. 2021; https://csg.sph.umich.edu/willer/public/glgc-lipids2021/):
#     five trans-ancestry lipid files (pleiotropy, Supplementary Section S5.5) and the LDL files
#     of ten ancestry-by-sex strata (replication; HIS and SAS are held out for validation)
#   UK Biobank, Neale lab round 2 (http://www.nealelab.is/uk-biobank): five traits, both sexes,
#     imputed v3, inverse-rank-normalized (irnt): BMI 21001, SBP 4080, TG 30870, HDL 30760, HbA1c 30750
#     (the three biomarker files are now served only under the name *.varorder.tsv.bgz; these files
#     are byte-identical to the files used in the paper, same md5, checked 2026-10-10)
set -e
cd "$(dirname "$0")"
get() { mkdir -p "$(dirname "$2")"; [ -s "$2" ] && { echo "exists: $2"; return; }; echo "downloading $1"; wget -q --tries=10 --timeout=120 --waitretry=30 -c -O "$2.part" "$1" && mv "$2.part" "$2"; }
G=https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results
for tr in HDL LDL TC logTG nonHDL; do
  f=with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_${tr}_INV_ALL_with_N_1.gz
  get $G/trans_ancestry/$f raw/glgc/pleiotropy/GLGC_trans/$f
done
get $G/trans_ancestry/README.txt raw/glgc/pleiotropy/GLGC_trans/README.txt
for s in AFR_HRC_1KGP3_others EAS_1KGP3 EUR_HRC_1KGP3_others HIS_1KGP3 SAS_HRC_1KGP3_others; do
  for sex in FEMALE MALE; do
    f=LDL_INV_${s}_${sex}.meta.singlevar.results.gz
    get $G/sex_and_ancestry_specific_summary_stats/$f raw/glgc/replication/GLGC_LDL_sex_ancestry/$f
  done
done
B=https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs
for f in 21001_irnt 4080_irnt 30870_irnt 30760_irnt 30750_irnt; do
  v=""; case $f in 30870_irnt|30760_irnt|30750_irnt) v=".varorder";; esac   # biomarker files: only the .varorder name is still served
  get $B/$f.gwas.imputed_v3.both_sexes$v.tsv.bgz raw/ukb/$f.gwas.imputed_v3.both_sexes.tsv.bgz
done
ls -lR raw | grep -E "gz|bgz" | awk '{print $5, $NF}'
echo "DOWNLOAD DONE"
