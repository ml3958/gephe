# [help function]
function usage(){

cat << EOF

Cross-species microbial genotype phenotype association for functional module discovery

Author: Menghan Liu
Tavazoie lab @ Columibia University

SYNTAX:
  gephe COMMAND [Options] <metadata> <metadata_pos> <phenotype_colname> <input_dir> <alignment_dir> <result_dir>

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
    -k, --align_k   alignment max hits        <default: 500>
    -n, --align_n   # of parallel jobs        <default: 50>


  COMMAND pog:
    -t, --top_perc  top percentage (%) as phenotype-conferring proteins <default: 0.2>
    -i, --mcl_i     MCL I value                                         <default: 1.4>
    -p, --pp_mode   phylogenetic profile "binary|proportion|evalue"     <default: binary>
    -b, --pp_bc     cutoff for binary phylogenetic profiles             <default: 0.5>
    -f, --pog_frac  fraction of proteins as representatice sequence     <default: 1>

  COMMAND module:
    -m, --module_n 	# of modules  <default: 20>
EOF
}


# [optional variables with default values]
ALIGNMENT_MAX=1000
ALIGNMENT_EVALUE='1e-10'
ALIGNMENT_QUERY_COVERAGE='66'
ALIGNMENT_SUBJECT_COVERAGE='60'
ALIGNMENT_NJOBS=50
TOP_PERCENT=0.2  # percentile of top proteins
MCL_I=1.4
PHYLOGENETIC_PROFILE_MODE='binary' #  `binary`/`evalue`/`count`
PHYLOGENETIC_PROFILE_BINARY_CUTOFF=0.5 #  `binary`/`evalue`/`count`
MODULE_N=12
PHENOTYPE_COLNAME='t'
REP_PROTEIN_FRAC=1

# [replace default with user-supplied value]
while [ $# -gt 0 ]
do
    unset OPTIND
    unset OPTARG
    while getopts "hq:s:e:k:n:t:i:p:b:m:f:j:" opt; do
        case $opt in
            q ) ALIGNMENT_QUERY_COVERAGE=$OPTARG ;;
            s ) ALIGNMENT_SUBJECT_COVERAGE=$OPTARG ;;
            e ) ALIGNMENT_EVALUE=$OPTARG ;;
            k ) ALIGNMENT_MAX=$OPTARG ;;
            n ) ALIGNMENT_NJOBS=$OPTARG ;;
            t ) TOP_PERCENT=$OPTARG ;;
            i ) MCL_I=$OPTARG ;;
            p ) PHYLOGENETIC_PROFILE_MODE=$OPTARG ;;
            b ) PHYLOGENETIC_PROFILE_BINARY_CUTOFF=$OPTARG ;;
            m ) MODULE_N=$OPTARG ;;
            f ) REP_PROTEIN_FRAC=$OPTARG ;;
	    j ) ALIGNMENT_NJOBS=$OPTARG ;;
            # s ) PHENOTYPE_COLNAME=$OPTARG ;;
            h ) usage
            exit 0;;
            *) usage
            exit 1;;
        esac
      done
   shift $((OPTIND-1))
   ARGS="${ARGS} $1"     # store positional values into a string "COMMANDS <metadata> <metadata_pos> <align_dir> <result_dir>"
   shift
done

# retrive each of the positional arguments
IFS=' ' read -r -a ARGS <<< "${ARGS}" # space delimited
export COMMAND=${ARGS[0]}
export METADATA=${ARGS[1]}
export METADATA_POS=${ARGS[2]}
export PHENOTYPE_COLNAME=${ARGS[3]}
export DIR_INPUT=${ARGS[4]}
export DIR_ALIGNMENT=${ARGS[5]}
export DIR=${ARGS[6]}


# export
export ALIGNMENT_MAX=${ALIGNMENT_MAX}
export ALIGNMENT_EVALUE=${ALIGNMENT_EVALUE}
export ALIGNMENT_QUERY_COVERAGE=${ALIGNMENT_QUERY_COVERAGE}
export ALIGNMENT_SUBJECT_COVERAGE=${ALIGNMENT_SUBJECT_COVERAGE}
export ALIGNMENT_NJOBS=${ALIGNMENT_NJOBS}
export PHENOTYPE_COLNAME=${PHENOTYPE_COLNAME}
export TOP_PERCENT=${TOP_PERCENT}  # percent of top proteins
export MCL_I=${MCL_I}
export PHYLOGENETIC_PROFILE_MODE=${PHYLOGENETIC_PROFILE_MODE} #  `binary`/`evalue`/`count`
export PHYLOGENETIC_PROFILE_BINARY_CUTOFF=${PHYLOGENETIC_PROFILE_BINARY_CUTOFF}
export MODULE_N=${MODULE_N}
export REP_PROTEIN_FRAC=${REP_PROTEIN_FRAC}
export ALIGNMENT_NJOBS=${ALIGNMENT_NJOBS}

# [prefix]
export gephe_dir=/home/menghan/Tools/gephe/
export PREFIX_PROTEIN=top${TOP_PERCENT}
export PREFIX_POG=${PREFIX_PROTEIN}_I${MCL_I}
export PREFIX_POG_PP=${PREFIX_POG}_binary${PHYLOGENETIC_PROFILE_BINARY_CUTOFF}
export PREFIX_MODULE=${PREFIX_POG_PP}_module${MODULE_N}
export DIR_FAA_MERGE=${DIR_ALIGNMENT}/faa_merge/
export DIR_ALIGNMENT_MERGE=${DIR_ALIGNMENT}/align_merge/
export DIR_ASSOCIATION=$DIR/association/
export DIR_POG=$DIR/pog/
export DIR_MODULE=$DIR/module/
export DIR_LOG=$DIR/logs_${PREFIX_MODULE}/


# # #[examine]
# echo "<$COMMAND> <$METADATA> <$METADATA_POS> <$DIR_ALIGNMENT> <$DIR>"
# echo  $ALIGNMENT_QUERY_COVERAGE
# echo  $ALIGNMENT_EVALUE
# echo  $ALIGNMENT_MAX
# echo  $ALIGNMENT_NJOBS
# echo  $TOP_PERCENT
# echo  $MCL_I
# echo  $PHYLOGENETIC_PROFILE_MODE
# echo  $PHYLOGENETIC_PROFILE_BINARY_CUTOFF
# echo  $MODULE_N
# # exit 1


function run_align(){
  echo --------------- --------------- --------------- ---------------
  echo STEP1: Alignment STARTED... [$(date --rfc-3339=seconds)]
  bash $gephe_dir/run_1_alignment.sh > $DIR_LOG/alignment.out 2> $DIR_LOG/alignment.err
  echo STEP1: Alignment FINISHED   [$(date --rfc-3339=seconds)]
  echo --------------- --------------- --------------- ---------------

}

function run_associate(){
  echo --------------- --------------- --------------- ---------------
  echo STEP2: Association STARTED... [$(date --rfc-3339=seconds)]
  bash $gephe_dir/run_2_association.sh \
   $METADATA $METADATA_POS $DIR_ALIGNMENT $DIR_ASSOCIATION > $DIR_LOG/association.out 2> $DIR_LOG/association.err
  echo STEP2: Association FINISHED   [$(date --rfc-3339=seconds)]
  echo --------------- --------------- --------------- ---------------
}

function run_pog(){
  echo --------------- --------------- --------------- ---------------
  echo STEP3: POG STARTED... [$(date --rfc-3339=seconds)]
  bash $gephe_dir/run_3a_protein.sh > $DIR_LOG/pog.out 2> $DIR_LOG/pog.err
  bash $gephe_dir/run_3b_pog.sh > $DIR_LOG/pog_pp.out 2> $DIR_LOG/pog_pp.err
  echo STEP3: POG FINISHED   [$(date --rfc-3339=seconds)]
  echo --------------- --------------- --------------- ---------------

}

function run_module(){
  echo --------------- --------------- --------------- ---------------
  echo STEP3: MODULE STARTED... [$(date --rfc-3339=seconds)]
  bash $gephe_dir/run_4_module.sh > $DIR_LOG/module.out 2> $DIR_LOG/module.err
  echo STEP3: MODULE FINISHED... [$(date --rfc-3339=seconds)]
  echo --------------- --------------- --------------- ---------------
}

# [run]
# [[ -z $COMMAND ]] && echo "! ERROR: command empty, exiting... "; usage; exit 1
# [[ -z "$METADATA" ]] && echo "! ERROR: <metadata> empty, exiting... "; usage; exit 1
# [[ -z "$METADATA_POS" ]] && echo "! ERROR: <metadata_pos> empty, exiting... "; usage; exit 1
# [[ -z "$DIR_ALIGNMENT" ]] && echo "! ERROR: <dir_align> empty, exiting... "; usage; exit 1
# [[ -z "$DIR" ]] && echo "! ERROR: <dir> empty, exiting... "; usage; exit 1

mkdir -p  $DIR_ALIGNMENT $DIR $DIR_ASSOCIATION $DIR_ASSOCIATION/summary/ $DIR_POG $DIR_MODULE $DIR_LOG

if [ $COMMAND == "all" ]; then
    run_align
    run_associate
    run_pog
    run_module
  elif [ $COMMAND == "align" ]; then
    run_align
  elif [ $COMMAND == "associate" ]; then
    run_associate
  elif [ $COMMAND == "pog" ]; then
    run_pog
  elif [ $COMMAND == "module" ]; then
    run_module
  elif [ $COMMAND == "help" ]; then
    usage
fi
