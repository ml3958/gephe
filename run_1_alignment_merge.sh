#!/bin/bash
test specific things

# DIR_INPUT=/mnt/data1/menghanliu/gephe_jgi/0_data/ # JGI data specific

echo "  build reference db"[$(date --rfc-3339=seconds)]
rm -rf $DIR_ALIGNMENT/input.fasta
rm -rf  $DIR_ALIGNMENT/input_ipr.txt
for i in $(cut -f1 $METADATA)
  do
    if [ -f $DIR_INPUT/${i}.faa ]
      then
        cat $DIR_INPUT/${i}.faa  >> $DIR_ALIGNMENT/input.fasta
      else
        echo ${i}.faa does not exist
    fi
  done
diamond makedb --in $DIR_ALIGNMENT/input.fasta --db $DIR_ALIGNMENT/input_database

# -----------------------------
echo " merge .faa files"[$(date --rfc-3339=seconds)] # added as part of V4
python $gephe_dir/alignment/diamond_to_pickle.py $METADATA_POS $DIR_INPUT $DIR_FAA_MERGE
# -----------------------------

echo "  Aligning..."[$(date --rfc-3339=seconds)]
run_diamond(){
  f=$1
  if [ ! -f $DIR_ALIGNMENT/${f}.diamond.out ] || [ ! -s $DIR_ALIGNMENT/${f}.diamond.out  ]
    then
      echo $f started $(date)
          diamond blastp \
          -q $DIR_INPUT/${f}.faa \
          --db $DIR_ALIGNMENT/input_database \
          --out $DIR_ALIGNMENT/${f}.diamond.out \
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
parallel -j ${ALIGNMENT_NJOBS} run_diamond  ::: `sed '1d' $METADATA_POS |cut -f1` ::: ${ALIGNMENT_EVALUE} ::: ${ALIGNMENT_MAX} ::: $DIR_ALIGNMENT


# 08/02/2022 [This is added in V4 because I always notice empty alignment output]
echo " Scan for empty files, and realign"[$(date --rfc-3339=seconds)]
parallel -j ${ALIGNMENT_NJOBS} run_diamond  ::: `sed '1d' $METADATA_POS |cut -f1` ::: ${ALIGNMENT_EVALUE} ::: ${ALIGNMENT_MAX} ::: $DIR_ALIGNMENT

# -----------------------------
echo " divide .merged diamond output"[$(date --rfc-3339=seconds)] # added as part of V4
python $gephe_dir/alignment/divide_diamondout.py $DIR_ALIGNMENT_MERGE $DIR_ALIGNMENT
# -----------------------------


echo "  Pickle diamond output"[$(date --rfc-3339=seconds)]
python $gephe_dir/alignment/diamond_to_pickle.py ${DIR_ALIGNMENT}
