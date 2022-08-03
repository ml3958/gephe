#!/bin/bash

echo "  Association each genome" [$(date --rfc-3339=seconds)]
parallel -j ${ALIGNMENT_NJOBS} \
        "python $gephe_dir/association/association.py \
          $METADATA $PHENOTYPE_COLNAME $DIR_ALIGNMENT/{1}.diamond.out.pickle $DIR_ASSOCIATION/ " \
        ::: `cut -f1 ${METADATA_POS} | sed '1d' `


echo "  Summarize association results"[$(date --rfc-3339=seconds)]
python $gephe_dir/association/association_summary.py $DIR_ASSOCIATION $METADATA $PHENOTYPE_COLNAME
