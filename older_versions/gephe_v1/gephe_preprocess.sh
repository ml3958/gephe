function usage(){

cat << EOF

Cross-species microbial genotype phenotype association for functional module discovery

Author: Menghan Liu
Tavazoie lab @ Columibia University

SYNTAX:
  gephe_preprocess [Options] <input_dir> <output_dir> <genome_prefix>

EOF
}

while [ $# -gt 0 ]
do
    unset OPTARG
    while getopts "h" opt; do
      case $opt in
            h ) usage
            exit 0;;
          esac
    done
    ARGS="${ARGS} $1"     # store positional values into a string "COMMANDS <metadata> <metadata_pos> <align_dir> <result_dir>"
    shift
done

# Positional arguments
IFS=' ' read -r -a ARGS <<< "${ARGS}" # space delimited
export gephe_dir=/home/menghan/Tools/gephe/
export DIR_INPUT=${ARGS[0]}
export DIR_OUTPUT=${ARGS[1]}
export PREFIX=${ARGS[2]}


# Function
function run_preprocess(){
  echo --------------- --------------- --------------- ---------------
  echo Preprocess STARTED... [$(date --rfc-3339=seconds)]
  python $gephe_dir/preprocess/preprocess.py $DIR_INPUT $DIR_OUTPUT $PREFIX
  echo PreProcess FINISHED   [$(date --rfc-3339=seconds)]
  echo --------------- --------------- --------------- ---------------

}

# Main script
if [[ -z $DIR_INPUT || -z $DIR_OUTPUT || -z $PREFIX ]]; then
  echo "Error: Missing arguments!"
  usage
  exit 1
fi

mkdir -p $DIR_OUTPUT
run_preprocess
