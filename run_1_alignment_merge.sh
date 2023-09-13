#!/bin/bash



# -----------------------------
echo "  1.1 build reference db"[$(date --rfc-3339=seconds)]
# -----------------------------
if [ ! -f ${DIR_ALIGNMENT_MASTER}/input.fasta ] || [ ! -s ${DIR_ALIGNMENT_MASTER}/input.fasta ]
then
  for i in $(cut -f2 $METADATA)
    do
      if [ -f ${i} ]
        then
          cat ${i}  >> ${DIR_ALIGNMENT_MASTER}/input.fasta
        else
          echo $(basename ${i}) does not exist
      fi
    done
  echo coping input to ${DIR_ALIGNMENT_MASTER}/input.fasta[$(date --rfc-3339=seconds)]
else
  echo ${DIR_ALIGNMENT_MASTER}/input.fasta exists, skip coping....[$(date --rfc-3339=seconds)]
fi
# rm -rf ${DIR_ALIGNMENT_MASTER}/input.fasta
# rm -rf  ${DIR_ALIGNMENT_MASTER}/input_ipr.txt
if [ ! -f ${DIR_ALIGNMENT_MASTER}/input.dmnd ] || [ ! -s ${DIR_ALIGNMENT_MASTER}/input.dmnd ]
then
  diamond makedb --in ${DIR_ALIGNMENT_MASTER}/input.fasta --db ${DIR_ALIGNMENT_MASTER}/input
else
  echo ${DIR_ALIGNMENT_MASTER}/input.dmnd exists, skip coping....[$(date --rfc-3339=seconds)]
fi

# -----------------------------
echo "  1.2 cp input .faa files"[$(date --rfc-3339=seconds)]
# -----------------------------
# determine whether to overwrite exisiting
if [ -d "$DIR_FAA" ]
then
    # read -p "Directory $DIR_FAA already exists. Do you want to overwrite it? (y/n) " answer
    echo "Directory $DIR_FAA already exists. Do you want to overwrite it? (y/n)" > /dev/tty
    # Wait for input for up to 1 minute or automatically respond "yes"
    read -t 60 answer < /dev/tty || answer="yes"



    if [[ "$answer" =~ [yY](es)* ]]
    # echo "Enter your name: " > /dev/tty
    # read name < /dev/tty
    then
        echo "Overwriting $DIR_FAA..."
        rm -rf $DIR_FAA
    else
        echo "Skipping copying"
    fi
else
    echo "Creating directory $DIR_FAA..."
fi

# actual copying
if [ ! -d "$DIR_FAA" ]
then
  mkdir -p $DIR_FAA
  while read -r line
  do
      file=$(echo "$line" | cut -f2)
      if [ ! -f "$file" ]
      then
          echo "Error: file $file not found"
          # exit 1
      fi
      echo "Copying $(basename "$file") to $DIR_FAA..."
      cp "$file" "$DIR_FAA/"
  done < "$METADATA_POS"
fi

# -----------------------------
echo " 1.3 merge .faa files"[$(date --rfc-3339=seconds)] # added as part of V4
# -----------------------------

# determine whether to overwrite exisiting
if [ -d "$DIR_FAA_MERGE" ]
then
    echo "Directory $DIR_FAA_MERGE already exists. Do you want to overwrite it? (y/n)" > /dev/tty

    # Wait for input for up to 1 minute or automatically respond "yes"
    read -t 60 answer < /dev/tty || answer="yes"

    if [[ "$answer" =~ [yY](es)* ]]
    then
        echo "Overwriting $DIR_FAA_MERGE..."
        rm -rf "$DIR_FAA_MERGE"
    else
        echo "Skipping copying"
        # exit 0
    fi
else
    echo "Creating directory $DIR_FAA_MERGE..."
fi


# actual merging
if [ ! -d "$DIR_FAA_MERGE" ]
then
    mkdir -p "$DIR_FAA_MERGE"
    python "$gephe_dir/alignment/merge_faa.py" -n_faa_to_merge "$N_FAA_TO_MERGE" "$DIR_FAA/" "$DIR_FAA_MERGE/"
fi



# -----------------------------
echo "  1.4 Aligning..."[$(date --rfc-3339=seconds)]
# -----------------------------
run_diamond(){
  f=$1
  if [ ! -f $DIR_ALIGNMENT_MERGE/${f}.diamond.out ] || [ ! -s $DIR_ALIGNMENT_MERGE/${f}.diamond.out  ]
    then
      echo $f started $(date)
          diamond blastp \
          -q ${DIR_FAA_MERGE}/${f}.faa \
          --out ${DIR_ALIGNMENT_MERGE}/${f}.diamond.out \
          --db ${DIR_ALIGNMENT_MASTER}/input \
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
mkdir -p ${DIR_ALIGNMENT_MERGE}
parallel -j 1 run_diamond  ::: `ls ${DIR_FAA_MERGE}| grep faa$| sed 's/.faa//g'`

# -----------------------------
echo " 1.5 Scan for empty files, and realign"[$(date --rfc-3339=seconds)]  # 08/02/2022 [This is added in V4 because I always notice empty alignment output]
# -----------------------------
parallel -j 1 run_diamond  ::: `ls ${DIR_FAA_MERGE}| grep faa| sed 's/.faa//g'`

# -----------------------------
echo " 1.6 Divide merged .diamond.out"[$(date --rfc-3339=seconds)] # added as part of V4  # 08/02/2022 [This is added in V4 because I always notice empty alignment output]
# -----------------------------

# determine whether to overwrite exisiting
if [ -d "$DIR_ALIGNMENT" ]
then
    echo "Directory $DIR_ALIGNMENT already exists. Do you want to overwrite it? (y/n)" > /dev/tty

    # Wait for input for up to 1 minute or automatically respond "yes"
    read -t 60 answer < /dev/tty || answer="yes"

    if [[ "$answer" =~ [yY](es)* ]]
    then
        echo "Overwriting $DIR_ALIGNMENT..."
        rm -rf "$DIR_ALIGNMENT"
    else
        echo ${DIR_ALIGNMENT} unchanged [$(date --rfc-3339=seconds)]
    fi
else
    echo "Creating directory $DIR_ALIGNMENT..."
fi

if [ -d ${DIR_ALIGNMENT} ]
then
  echo ${DIR_ALIGNMENT} exists, skip dividing merged*.diamond.out [$(date --rfc-3339=seconds)]
else
  mkdir -p ${DIR_ALIGNMENT}
  parallel "python $gephe_dir/alignment/divide_diamondout.py $DIR_ALIGNMENT_MERGE/{} $DIR_ALIGNMENT" ::: `ls $DIR_ALIGNMENT_MERGE`
fi


# -----------------------------
echo "  1.7 Pickle diamond output"[$(date --rfc-3339=seconds)]
# -----------------------------
parallel -j ${ALIGNMENT_NJOBS} \
  "[ ! -f {}.pickle ] && python $gephe_dir/alignment/diamond_to_pickle.py {} || echo {}.pickle exists, skipping...." \
  ::: `ls ${DIR_ALIGNMENT}/*diamond.out`
