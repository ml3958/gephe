#!/bin/bash

# DIR_FAA=/mnt/data1/menghanliu/gephe_jgi/0_data/ # JGI data specific
QUERY_PROTEIN=${1}  # a .faa FILE
DB=${2}     # a directory of proteome files
METADATA=${3}    # a metadata with first column being the target genomes/proteomes of interests
PREFIX=${4}
DIR_OUT=${5} # a directory to store intermediate files
ALIGNMENT_MAX=${6}  # max hit per protein, usually the # of genomes x 5

# DIR_FAA_MERGE=$DIR_OUT/faa_merge/
# DIR_ALIGNMENT_MERGE=$DIR_OUT/diamond_merge/
# DIR_ALIGNMENT=$DIR_OUT/diamond/
# PREFIX=$(basename ${QUERY_PROTEIN} .faa) # needs to be predefined
ALIGNMENT_EVALUE='1e-10'
ALIGNMENT_QUERY_COVERAGE='66'
ALIGNMENT_SUBJECT_COVERAGE='50'


# -----------------------------
echo "  2 Aligning..."[$(date --rfc-3339=seconds)]
diamond blastp \
    -q ${QUERY_PROTEIN} \
    --out ${DIR_OUT}/${PREFIX}.diamond.out \
    --db ${DB} \
    -e ${ALIGNMENT_EVALUE} -k ${ALIGNMENT_MAX} \
    --query-cover ${ALIGNMENT_QUERY_COVERAGE} \
    --subject-cover ${ALIGNMENT_SUBJECT_COVERAGE} \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore \
    -b8 -c1


# # -----------------------------
# echo "  1.4 Pickle diamond output"[$(date --rfc-3339=seconds)]
# [ ! -f ${PREFIX}.diamond.out.pickle ] && python $gephe_dir/alignment/diamond_to_pickle.py ${PREFIX}.diamond.out || echo ${PREFIX}.diamond.out.pickle exists, skipping....

# -----------------------------
echo "  1.3 Convert diamond.out to phylogenetic profiles"[$(date --rfc-3339=seconds)]
diamond_to_pp.py ${DIR_OUT}${PREFIX}.diamond.out ${METADATA} ${DIR_OUT}
diamond_to_pp_count.py ${DIR_OUT}${PREFIX}.diamond.out ${METADATA} ${DIR_OUT}
