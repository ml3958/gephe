#!/bin/bash

input_dir=/mnt/data1/menghanliu/gephe_jgi/0_data/ # JGI data specific

echo "  Copy input sequence"[$(date --rfc-3339=seconds)]
rm -rf $DIR_ALIGNMENT/input.fasta
rm -rf  $DIR_ALIGNMENT/input_ipr.txt
for i in $(cut -f1 $METADATA)
  do
    if [ -f $input_dir/${i}.genes.faa_corrected ]
      then
        cat $input_dir/${i}.genes.faa_corrected  >> $DIR_ALIGNMENT/input.fasta
      else
        echo ${i}.genes.faa_corrected does not exist
    fi
  done
diamond makedb --in $DIR_ALIGNMENT/input.fasta --db $DIR_ALIGNMENT/input_database


echo "  Aligning..."[$(date --rfc-3339=seconds)]
run_diamond(){
  f=$1
  input_dir=/mnt/data1/menghanliu/gephe_jgi/0_data/
  DIR_ALIGNMENT=$4
  if [ ! -f $DIR_ALIGNMENT/${f}.diamond.out ] || [ ! -s $DIR_ALIGNMENT/${f}.diamond.out  ]
    then
      echo $f started $(date)
          diamond blastp \
          -q $input_dir/${f}.genes.faa_corrected \
          --db $DIR_ALIGNMENT/input_database \
          --out $DIR_ALIGNMENT/${f}.diamond.out \
          -e ${ALIGNMENT_EVALUE} -k ${ALIGNMENT_MAX} --query-cover ${ALIGNMENT_COVERAGE}
          echo $f finished $(date)
    else
      echo $f already exists
  fi
}
export -f run_diamond
parallel -j ${ALIGNMENT_NJOBS} run_diamond  ::: `sed '1d' $METADATA_POS |cut -f1` ::: $diamond_e ::: $diamond_k ::: $DIR_ALIGNMENT


echo "  Pickle diamond output"[$(date --rfc-3339=seconds)]
python $gephe_dir/alignment/diamond_to_pickle.py ${DIR_ALIGNMENT}
