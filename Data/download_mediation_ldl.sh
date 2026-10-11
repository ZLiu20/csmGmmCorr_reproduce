#!/bin/bash
# Downloads for the three-mediator LDL application (Data/prepare_mediation_ldl.R), about 4 GB:
#   GTEx v10 liver cis-eQTL statistics of all tested variant-gene pairs, eQTL Catalogue dataset
#     QTD000266 (QTD000266.all.tsv.gz, 2.9 GB), and the GENCODE v26 gene annotation of GTEx (TSS)
#   GTEx v8 variant lookup table (build-38 variant IDs with their build-37 positions)
#   GTEx v8 MASHR expression-prediction models (Barbeira et al. 2021; liver model and covariances)
# The LDL-cholesterol statistics are the GLGC EUR female and male files of download_data.sh.
# Run from the repository root on a node with internet access:   bash Data/download_mediation_ldl.sh
# If a download fails on the cluster, download the file on another computer to the same path.
set -e
cd "$(dirname "$0")"
D=raw/mediation_ldl
mkdir -p $D
get() { [ -s "$2" ] && { echo "exists: $2"; return; }; echo "downloading $1"; wget -nv --tries=10 --timeout=120 --waitretry=30 -c -O "$2.part" "$1" && mv "$2.part" "$2"; }
G=https://storage.googleapis.com/adult-gtex
get https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000266/QTD000266.all.tsv.gz $D/QTD000266.all.tsv.gz
get $G/references/v8/reference-tables/gencode.v26.GRCh38.genes.gtf $D/gencode.v26.GRCh38.genes.gtf
get $G/references/v8/reference-tables/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz \
    $D/GTEx_v8_lookup_table.txt.gz
get "https://zenodo.org/records/3518299/files/mashr_eqtl.tar?download=1" $D/mashr_eqtl.tar
# keep only the liver model and its covariances
if [ ! -s $D/mashr_Liver.db ]; then
  tar -xf $D/mashr_eqtl.tar -C $D --wildcards '*mashr_Liver.db' '*mashr_Liver.txt.gz'
  find $D -name 'mashr_Liver.db' -not -path "$D/mashr_Liver.db" -exec mv {} $D/ \;
  find $D -name 'mashr_Liver.txt.gz' -not -path "$D/mashr_Liver.txt.gz" -exec mv {} $D/ \;
fi
ls -l $D
echo "DOWNLOAD DONE"
