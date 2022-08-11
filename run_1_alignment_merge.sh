#!/bin/bash

# DIR_INPUT=/mnt/data1/menghanliu/gephe_jgi/0_data/ # JGI data specific

echo "  1.1 cp input .faa files"[$(date --rfc-3339=seconds)]
if [ -d $DIR_FAA ]
then
    echo $DIR_FAA exists, skip coping [$(date --rfc-3339=seconds)]
else
    mkdir $DIR_FAA
    for i in $(cut -f1 $METADATA_POS)
      do
        if [ ! -f $DIR_FAA/${i}.faa ]
          then
            echo coping to $DIR_FAA/${i}.faa....
            if [ -f $DIR_INPUT/${i}.faa ]
              then
                cp $DIR_INPUT/${i}.faa $DIR_FAA/
              else
                echo $DIR_INPUT/${i}.faa does not exist
            fi
        fi
    done
fi


# -----------------------------
echo " 1.2 merge .faa files"[$(date --rfc-3339=seconds)] # added as part of V4
if [ -d ${DIR_FAA_MERGE} ]
then
  echo $DIR_FAA_MERGE exists, skip faa merging [$(date --rfc-3339=seconds)]
else
  mkdir $DIR_FAA_MERGE
  python $gephe_dir/alignment/merge_faa.py ${DIR_FAA} ${DIR_FAA_MERGE} ${N_FAA_TO_MERGE}
fi
# -----------------------------


echo "  1.3 build reference db"[$(date --rfc-3339=seconds)]
if [ ! -f ${DIR_ALIGNMENT_MASTER}/input.fasta ] || [ ! -s ${DIR_ALIGNMENT_MASTER}/input.fasta ]
then
  # for i in $(cut -f1 $METADATA)
  #   do
  #     if [ -f $DIR_INPUT/${i}.faa ]
  #       then
  #         cat $DIR_INPUT/${i}.faa  >> ${DIR_ALIGNMENT_MASTER}/input.fasta
  #       else
  #         echo ${i}.faa does not exist
  #     fi
  #   done
  echo coping input to ${DIR_ALIGNMENT_MASTER}/input.fasta[$(date --rfc-3339=seconds)]
else
  echo ${DIR_ALIGNMENT_MASTER}/input.fasta exists, skip coping....[$(date --rfc-3339=seconds)]
fi
# rm -rf ${DIR_ALIGNMENT_MASTER}/input.fasta
# rm -rf  ${DIR_ALIGNMENT_MASTER}/input_ipr.txt
if [ ! -f ${DIR_ALIGNMENT_MASTER}/input_database.dmnd ] || [ ! -s ${DIR_ALIGNMENT_MASTER}/input_database.dmndinput.fasta ]
then
  diamond makedb --in ${DIR_ALIGNMENT_MASTER}/input.fasta --db ${DIR_ALIGNMENT_MASTER}/input_database
else
  echo ${DIR_ALIGNMENT_MASTER}/input_database.dmnd exists, skip coping....[$(date --rfc-3339=seconds)]
fi


echo "  1.4 Aligning..."[$(date --rfc-3339=seconds)]
run_diamond(){
  f=$1
  if [ ! -f $DIR_ALIGNMENT_MERGE/${f}.diamond.out ] || [ ! -s $DIR_ALIGNMENT_MERGE/${f}.diamond.out  ]
    then
      echo $f started $(date)
          diamond blastp \
          -q ${DIR_FAA_MERGE}/${f}.faa \
          --out ${DIR_ALIGNMENT_MERGE}/${f}.diamond.out \
          --db ${DIR_ALIGNMENT_MASTER}/input_database \
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
parallel -j 1 run_diamond  ::: `ls ${DIR_FAA_MERGE}| grep faa| sed 's/.faa//g'`

# 08/02/2022 [This is added in V4 because I always notice empty alignment output]
echo " 1.5 Scan for empty files, and realign"[$(date --rfc-3339=seconds)]
parallel -j 1 run_diamond  ::: `ls ${DIR_FAA_MERGE}| grep faa| sed 's/.faa//g'`

# -----------------------------
echo " 1.6 Divide merged .diamond.out"[$(date --rfc-3339=seconds)] # added as part of V4
if [ -d ${DIR_ALIGNMENT} ]
then
  echo ${DIR_ALIGNMENT} exists, skip dividing merged*.diamond.out [$(date --rfc-3339=seconds)]
else
  mkdir ${DIR_ALIGNMENT}
  parallel "python $gephe_dir/alignment/divide_diamondout.py $DIR_ALIGNMENT_MERGE/{} $DIR_ALIGNMENT" ::: `ls $DIR_ALIGNMENT_MERGE`
fi
# -----------------------------


echo "  1.7 Pickle diamond output"[$(date --rfc-3339=seconds)]
parapllel -j ${ALIGNMENT_NJOBS} \
  "[ ! -f {}.pickle ] && python $gephe_dir/alignment/diamond_to_pickle.py {} || echo {}.pickle exists, skipping...." \
  ::: `ls ${DIR_ALIGNMENT}/*diamond.out`
