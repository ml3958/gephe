#!/bin/bash

# DIR_INPUT=/mnt/data1/menghanliu/gephe_jgi/0_data/ # JGI data specific

echo "  1.1 cp input .faa files"[$(date --rfc-3339=seconds)]
for i in $(cut -f1 $METADATA_POS)
  do
    if [ ! -f $DIR_FAA/${i}.faa ]
      then
        if [ ! -f $DIR_INPUT/${i}.faa ]
          then
            cp $DIR_INPUT/${i}.faa $DIR_FAA/
        else
          echo ${i}.faa does not exist
        fi
    fi
done


# -----------------------------
echo " 1.2 merge .faa files"[$(date --rfc-3339=seconds)] # added as part of V4
python $gephe_dir/alignment/merge_faa.py ${DIR_FAA} ${DIR_FAA_MERGE} ${N_FAA_TO_MERGE}
# -----------------------------


echo "  1.3 build reference db"[$(date --rfc-3339=seconds)]
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


echo "  1.4 Aligning..."[$(date --rfc-3339=seconds)]
run_diamond(){
  f=$1
  if [ ! -f $DIR_ALIGNMENT_MERGE/${f}.diamond.out ] || [ ! -s $DIR_ALIGNMENT_MERGE/${f}.diamond.out  ]
    then
      echo $f started $(date)
          diamond blastp \
          -q ${DIR_FAA_MERGE}/${f}.faa \
          --out ${DIR_ALIGNMENT_MERGE}/${f}.diamond.out \
          --db ${DIR_ALIGNMENT}/input_database \
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
# parallel -j ${ALIGNMENT_NJOBS} run_diamond  ::: `ls ${DIR_FAA_MERGE}| grep faa| sed 's/.faa//g'`
parallel -j 2 run_diamond  ::: `ls ${DIR_FAA_MERGE}| grep faa| sed 's/.faa//g'`

# 08/02/2022 [This is added in V4 because I always notice empty alignment output]
echo " 1.5 Scan for empty files, and realign"[$(date --rfc-3339=seconds)]
parallel -j 2 run_diamond  ::: `ls ${DIR_FAA_MERGE}| grep faa| sed 's/.faa//g'`

# -----------------------------
echo " 1.6 Divide .merged diamond output"[$(date --rfc-3339=seconds)] # added as part of V4
parallel "python $gephe_dir/alignment/divide_diamondout.py $DIR_ALIGNMENT_MERGE/{.} $DIR_ALIGNMENT" ::: `ls $DIR_ALIGNMENT_MERGE`
# -----------------------------


echo "  1.7 Pickle diamond output"[$(date --rfc-3339=seconds)]
python $gephe_dir/alignment/diamond_to_pickle.py ${DIR_ALIGNMENT}
