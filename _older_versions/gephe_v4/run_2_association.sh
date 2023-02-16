#!/bin/bash
echo "  Association each genome" [$(date --rfc-3339=seconds)]
parallel -j ${ALIGNMENT_NJOBS} \
        "python $gephe_dir/association/association.py \
          $METADATA $PHENOTYPE_COLNAME {} $DIR_ASSOCIATION/ " \
        ::: `ls $DIR_ALIGNMENT/*.diamond.out.pickle`

echo "  Summarize association results"[$(date --rfc-3339=seconds)]
python $gephe_dir/association/association_summary.py $DIR_ASSOCIATION $METADATA $PHENOTYPE_COLNAME ${ALIGNMENT_MAX}
