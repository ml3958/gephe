#!/bin/bash

# Author: Menghan Liu
# Date: September 29, 2021

# select phenotype-conferring proteins,
# proteins-->POG,
# POG annotation

echo "  Select proteins"  [$(date --rfc-3339=seconds)]
python $gephe_dir/pog/select_proteins.py \
  $DIR_ALIGNMENT $DIR_ASSOCIATION $DIR_POG \
  -p $PREFIX_PROTEIN -percentile ${TOP_PERCENT}

echo "  Phylogenetic profile file for selected proteins"  [$(date --rfc-3339=seconds)]
python $gephe_dir/pog/generate_pp_for_selectedproteins.py \
  $DIR_ASSOCIATION $DIR_POG -p $PREFIX_PROTEIN


# echo "  .abc file for selected proteins"  [$(date --rfc-3339=seconds)]   ---- Retired in gephe v5
cut -d, -f1 ${PREFIX_PROTEIN}.csv | seqtk subseq ${DIR_ALIGNMENT_MASTER}/input.fasta > ${PREFIX_PROTEIN}.protseq.faa
run_diamond_abc(){
  f=$1
  if [ ! -f $DIR_POG/${f}.diamond.out ] || [ ! -s $DIR_POG//${f}.diamond.out  ]
    then
      echo $f started $(date)
          diamond blastp \
          -q ${DIR_POG}/${f}.faa \
          --out ${DIR_POG}/${f}.diamond.out \
          --db ${DIR_POG}/${f}.faa \
          -e ${ALIGNMENT_EVALUE} -k ${ALIGNMENT_MAX} \
          --query-cover ${ALIGNMENT_QUERY_COVERAGE} \
          --subject-cover ${ALIGNMENT_SUBJECT_COVERAGE} \
          --outfmt 6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore \
	        -b8 -c1
          echo $f finished $(date)
    else
      echo $f already exists
  fi
}
export -f run_diamond
run_diamond_abc ${PREFIX_PROTEIN}.protseq.faa

cat ${PREFIX_PROTEIN}.protseq.faa.diamond.out| cut -f1,2,13 > ${PREFIX_PROTEIN}.abc



echo "  mxcload"  [$(date --rfc-3339=seconds)]
if [ ! -f $DIR_POG/${PREFIX_PROTEIN}.mci ]
then
  mcxload \
    -abc $DIR_POG/${PREFIX_PROTEIN}.abc \
    -o $DIR_POG/${PREFIX_PROTEIN}.mci \
    -write-tab $DIR_POG/${PREFIX_PROTEIN}.dic \
    --stream-neg-log10 -stream-tf 'ceil(200)' \
    --write-binary -ri max
fi
