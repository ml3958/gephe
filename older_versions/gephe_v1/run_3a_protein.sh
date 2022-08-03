#!/bin/bash

# Author: Menghan Liu
# Date: September 29, 2021

# select phenotype-conferring proteins,
# proteins-->POG,
# POG annotation



echo "  Select proteins"  [$(date --rfc-3339=seconds)]
python $gephe_dir/pog/select_proteins.py \
  $DIR_ALIGNMENT $DIR_ASSOCIATION $DIR_POG \
  -p $PREFIX_PROTEIN -percentile ${TOP_PERCENT}

echo "  Phylogenetic profile file for selected proteins"  [$(date --rfc-3339=seconds)]
python $gephe_dir/pog/generate_pp_for_selectedproteins.py \
  $DIR_ASSOCIATION $DIR_POG -p $PREFIX_PROTEIN


echo "  .abc file for selected proteins"  [$(date --rfc-3339=seconds)]
python $gephe_dir/pog/generate_abcfile_for_selectedproteins.py \
  $DIR_ALIGNMENT $DIR_ASSOCIATION $DIR_POG -p $PREFIX_PROTEIN


echo "  mxcload"  [$(date --rfc-3339=seconds)]
mcxload \
  -abc $DIR_POG/${PREFIX_PROTEIN}.abc \
  -o $DIR_POG/${PREFIX_PROTEIN}.mci \
  -write-tab $DIR_POG/${PREFIX_PROTEIN}.dic \
  --stream-neg-log10 -stream-tf 'ceil(200)' \
  --write-binary -ri max
