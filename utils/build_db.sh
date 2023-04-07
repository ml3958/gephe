#!/bin/bash

# DIR_FAA=/mnt/data1/menghanliu/gephe_jgi/0_data/ # JGI data specific
INPUT_FILE_PATH=${1}     # a directory of proteome files
DIR_OUT=${2}
PREFIX=${3}
# the INPUT_FILE_PATH needs to contain a column names input_path


mkdir -p ${DIR_OUT}
# -----------------------------
echo "  1.1 build reference db"[$(date --rfc-3339=seconds)]
# parallel "echo {}" ::: $(awk 'NR==1{for(i=1;i<=NF;i++) if($i=="input_path") col=i} {print $col}' input_2881_5species_per_GTDB_2428gem.txt))

for path in $(awk 'NR==1{for(i=1;i<=NF;i++) if($i=="input_path") col=i} {print $col}' ${INPUT_FILE_PATH})
do
  if [ -f $path]; then
    echo coping $path
    cat $path >> ${DIR_OUT}/${PREFIX}.faa
  else
    echo $path not exists
  fi
done



if [ ! -f ${PREFIX}.dmnd ] || [ ! -s ${PREFIX}.dmnd ]
then
  diamond makedb --in ${DIR_OUT}/${PREFIX}.faa --db ${DIR_OUT}/${PREFIX}
else
  echo ${DIR_OUT}/${PREFIX}.dmnd exists, skip coping....[$(date --rfc-3339=seconds)]
fi


# # -----------------------------
# echo "  1.1 build reference db"[$(date --rfc-3339=seconds)]
# if [ ! -f ${DIR_OUT}/${PREFIX}.faa ] || [ ! -s ${DIR_OUT}/${PREFIX}.faa ]
# then
#   for i in $(cut -f1 $METADATA)
#     do
#       if [ -f $DIR_FAA/${i}.faa ]
#         then
#           cat $DIR_FAA/${i}.faa  >> ${DIR_OUT}/${PREFIX}.faa
#         else
#           echo ${i}.faa does not exist
#       fi
#     done
#   echo coping input to ${DIR_OUT}/${PREFIX}.faa[$(date --rfc-3339=seconds)]
# else
#   echo ${DIR_OUT}/${PREFIX}.faa exists, skip coping....[$(date --rfc-3339=seconds)]
# fi
# # rm -rf ${DIR_OUT}/${PREFIX}.faa
# # rm -rf ${DIR_OUT}/input_ipr.txt
# if [ ! -f ${PREFIX}.dmnd ] || [ ! -s ${PREFIX}.dmnd ]
# then
#   diamond makedb --in ${DIR_OUT}/${PREFIX}.faa --db ${DIR_OUT}/${PREFIX}
# else
#   echo ${DIR_OUT}/${PREFIX}.dmnd exists, skip coping....[$(date --rfc-3339=seconds)]
# fi
