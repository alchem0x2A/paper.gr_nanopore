#Generic script for submitting comsol job
#usage: sub_job.sh [COMSOL_IN_NAME] [optional:MEM]
INFILE=$1
OUTFILE=$(echo $INFILE | cut -f 1 -d ".")
OUTFILE+="_out.mph"
OUTFILE=$SCRATCH/$OUTFILE


#Memory alloc
MEM=4096
if [ $# -ge 2 ]; then
	MEM=$2
fi
echo "Rmem=$MEM"

NP=4
if [ $# -ge 3 ]; then
        NP=$3
fi

JOB_STRING="comsol batch -np $NP -inputfile $INFILE -outputfile $OUTFILE"

#Now the submission line
SUB_STRING="bsub -n $NP -W 24:00 -R \"rusage[mem=$MEM]\" \"$JOB_STRING\""
echo "Final check:" $SUB_STRING

while true; do
	read -p  "Ready for sumission?[y/n]" OPTION
	case $OPTION in
		[Yy]* ) eval $SUB_STRING; break;;
		[Nn]* ) break;;
		* ) echo "y or n?";;
	esac
done
