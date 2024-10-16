#!/bin/bash


Rscript --vanilla $gephe_dir/module/cocluster_freq.R \
    --dir_pog ${DIR_POG} --dir_module ${DIR_MODULE} \
    --prefix_pog_pp ${PREFIX_POG_PP}

Rscript --vanilla $gephe_dir/module/module_plot.R \
    -m $METADATA -c $PHENOTYPE_COLNAME \
    --dir_pog ${DIR_POG} --dir_module ${DIR_MODULE} \
    --prefix_protein ${PREFIX_PROTEIN} --prefix_pog ${PREFIX_POG} --prefix_pog_pp ${PREFIX_POG_PP} --prefix_module ${PREFIX_MODULE} \
    -n ${MODULE_N}


# mv ${DIR_POG}/*pdf ${DIR_MODULE}/
# mv ${DIR_POG}/*cocluster* ${DIR_MODULE}/
