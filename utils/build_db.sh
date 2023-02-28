#!/bin/bash

# DIR_FAA=/mnt/data1/menghanliu/gephe_jgi/0_data/ # JGI data specific
DIR_FAA=${1}     # a directory of proteome files
METADATA=${2}    # a metadata with first column being the target genomes/proteomes of interests
PREFIX=${3}
DIR_OUT=${4} # a directory to store intermediate files


# -----------------------------
echo "  1.1 build reference db"[$(date --rfc-3339=seconds)]
if [ ! -f ${DIR_OUT}/${PREFIX}.faa ] || [ ! -s ${DIR_OUT}/${PREFIX}.faa ]
then
  for i in $(cut -f1 $METADATA)
    do
      if [ -f $DIR_FAA/${i}.faa ]
        then
          cat $DIR_FAA/${i}.faa  >> ${DIR_OUT}/${PREFIX}.faa
        else
          echo ${i}.faa does not exist
      fi
    done
  echo coping input to ${DIR_OUT}/${PREFIX}.faa[$(date --rfc-3339=seconds)]
else
  echo ${DIR_OUT}/${PREFIX}.faa exists, skip coping....[$(date --rfc-3339=seconds)]
fi
# rm -rf ${DIR_OUT}/${PREFIX}.faa
# rm -rf ${DIR_OUT}/input_ipr.txt
if [ ! -f ${PREFIX}.dmnd ] || [ ! -s ${PREFIX}.dmnd ]
then
  diamond makedb --in ${DIR_OUT}/${PREFIX}.faa --db ${DIR_OUT}/${PREFIX}
else
  echo ${DIR_OUT}/${PREFIX}.dmnd exists, skip coping....[$(date --rfc-3339=seconds)]
fi
