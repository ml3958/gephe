---------
Author: Menghan Liu (ml4640@columbia.edu)
Saeed Tavazoie lab, Department of Biology, Columbia University
---------


---------
Usage:
---------


Cross-species microbial genotype phenotype association for functional module discovery

Author: Menghan Liu
Tavazoie lab @ Columibia University

SYNTAX:
  gephe COMMAND [Options] <metadata> <metadata_pos> <phenotype_colname> <alignment_dir> <result_dir>

COMMAND:
  all         Run complete pipeline [align -> associate -> pog -> module]
  align       Align proteins
  associate   Associate proteins with the phenotype
  pog         Select phenotype-conferring proteins and convert to POGs
  module      Assign POGs to modules
  help        Show this message

OPTIONS:
  general:
    -h, --help 		show help message

  COMMAND all:
    --skip-align
    --skip-align-associate

  COMMAND align:
    -q, --align_query_c   alignment query coverage cutoff <default: 66>
    -s, --align_subject_c alignment query coverage cutoff <default: 60>
    -e, --align_e   alignment evalue cutoff   <default: 1e-10>
    -k, --align_k   alignment max hits        <default: 10000>
    -n, --align_n   # of parallel jobs        <default: 50>


  COMMAND pog:
    -t, --top_perc  top percentage (%) as phenotype-conferring proteins <default: 0.2>
    -i, --mcl_i     MCL I value                                         <default: 1.4>
    -p, --pp_mode   phylogenetic profile "binary|proportion|evalue"     <default: binary>
    -b, --pp_bc     cutoff for binary phylogenetic profiles             <default: 0.5>
    -f, --pog_frac  fraction of proteins as representatice sequence     <default: 1>

  COMMAND module:
    -m, --module_n 	# of modules  <default: 20>
