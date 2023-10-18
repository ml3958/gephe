
PREFIX=$1
FAA_DB=$2
N_PROTEIN=$3
MCL_I=$4


PROTEIN_FILE=${PREFIX}.faa
DIAMOND_OUT=${PREFIX}.faa.diamond.out
ABC_FILE=${PREFIX}.faa.abc
MCI_FILE=${PREFIX}.faa.mci
DIC_FILE=${PREFIX}.faa.dic
MCL_FILE=${PREFIX}.faa_I${MCL_I}.mcl

seqtk sample -s $RANDOM ${FAA_DB} ${N_PROTEIN} > ${PROTEIN_FILE}

diamond blastp \
          -q ${PROTEIN_FILE} \
          --db ${PROTEIN_FILE} \
          --out ${DIAMOND_OUT}| \
          -e 1E-10 -k 100000 \
          --query-cover 60 \
          --subject-cover 50 \
          --outfmt 6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore \
	  -b8 -c1

cat ${DIAMOND_OUT}| cut -f1,2,13 > ${ABC_FILE}


echo "  mxcload"  [$(date --rfc-3339=seconds)]
mcxload \
  -abc ${ABC_FILE} \
  -o ${MCI_FILE} \
  -write-tab ${DIC_FILE} \
  --stream-neg-log10 -stream-tf 'ceil(200)' \
  --write-binary -ri max


for MCL_I in {1.4,2.5,3}
do
  mcl \
    ${DIC_FILE}  \
    -o ${MCL_FILE} \
    -I ${MCL_I} -te 10
done
