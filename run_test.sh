#!/bin/bash

set -e

cd test_data
for f in genomic.gtf rna.fna protein.faa; do
    [ ! -f "$f" ] && gunzip -k "$f.gz"
done
cd ..

./EasyRNASeq \
  -input_dir test_data \
  -reference test_data/rna.fna \
  -gtf test_data/genomic.gtf \
  -meta test_data/metadata.txt \
  -ref_group Control \
  -out_dir test_result \
  -threads 4 


# ./EasyRNASeq \
#   -input_dir test_data \
#   -reference test_data/rna.fna \
#   -gtf test_data/genomic.gtf \
#   -meta test_data/metadata.txt \
#   -ref_group Control \
#   -out_dir test_result \
#   -threads 4 \
#   -run_enrichment \
#   -pep test_data/protein.faa \
#   -eggnog_db /path/to/eggnog_db \
#   -run_wgcna


echo "Results saved in ./test_result"
