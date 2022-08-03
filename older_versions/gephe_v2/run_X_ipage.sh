RUNFILE=$1
EXPTYPE=$2
BINNR=$3
DB_PREFIX=$4 # database prefix
PAGEDIR=~/Tools/iPAGE/

if [ $EXPTYPE=='continuous' ]; then
	perl $PAGEDIR/page.pl --independence=0 --expfile=$RUNFILE  --minr=-1 --exptype=$EXPTYPE --species=$DB_PREFIX --ebins=$BINNR --max_p 0.05 > $RUNFILE.out 2> $RUNFILE.err
elif [ $EXPTYPE=='discrete' ]; then
	perl $PAGEDIR/page.pl --independence=0 --expfile=$RUNFILE  --minr=-1 --exptype=$EXPTYPE --species=$DB_PREFIX --ebins=$BINNR --max_p 0.05 > $RUNFILE.out 2> $RUNFILE.err
fi
