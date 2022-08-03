#!/bin/bash

# DIR_INPUT=/mnt/data1/menghanliu/gephe_jgi/0_data/ # JGI data specific

echo "  Copy input sequence"[$(date --rfc-3339=seconds)]
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
          -e ${ALIGNMENT_EVALUE} -k ${ALIGNMENT_MAX} --query-cover ${ALIGNMENT_COVERAGE}
          echo $f finished $(date)
    else
      echo $f already exists
  fi
}
export -f run_diamond
parallel -j ${ALIGNMENT_NJOBS} run_diamond  ::: `sed '1d' $METADATA_POS |cut -f1` ::: ${ALIGNMENT_EVALUE} ::: ${ALIGNMENT_MAX} ::: $DIR_ALIGNMENT


echo "  Pickle diamond output"[$(date --rfc-3339=seconds)]
python $gephe_dir/alignment/diamond_to_pickle.py ${DIR_ALIGNMENT}
