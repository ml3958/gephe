#!/bin/bash

# Author: Menghan Liu
# Date: September 29, 2021

echo "  mcl"[$(date --rfc-3339=seconds)]
if [ ! -f $DIR_POG/${PREFIX_POG}.mcl ]
then
  mcl \
    $DIR_POG/${PREFIX_PROTEIN}.mci \
    -o $DIR_POG/${PREFIX_POG}.mcl \
    -I ${MCL_I} -te 10
fi

#
echo "  POG annotate"[$(date --rfc-3339=seconds)]
python $gephe_dir/pog/annotate_pog.py \
  $DIR_ALIGNMENT $DIR_POG $PREFIX_PROTEIN $PREFIX_POG -f ${REP_PROTEIN_FRAC}

echo "  POG phylogenetic profile file"[$(date --rfc-3339=seconds)]
python $gephe_dir/pog/generate_pog_pp.py \
  $DIR_POG \
  -prefix_protein $PREFIX_PROTEIN \
  -prefix_pog $PREFIX_POG \
  -c ${PHYLOGENETIC_PROFILE_BINARY_CUTOFF}
