##!/bin/bash
## extraction of total, clipped, and genome-specific aligned read pair counts
## INPUT: folder names and reference genomes
## OUTPUT: table with read pair counts (text file)
## Author: Doris Chen 
## version 240318


function usage()
{  echo
   echo "Get read counts from clipping and specific aligments"
   echo
   echo "Usage: $0 -i SAMPLE_ID_INDEX(e.g. 1) -c CLIPPED_LOG_DIR -m MAPPED_DIR -1 GENOME1(CBS138|ZP591.22) -2 GENOME2(R64|ASMv1) -o OUTFILEBASE"
	  echo
			echo "  -i ... index of unique sample id in file name after separation by underscore"
			echo "  -c ... directory of clipping (cutadapt) *.out log files"
			echo "  -m ... directory of mapped files"
			echo "  -1 ... calibration genome"
			echo "  -2 ... experimental genome"
			echo "  -o ... prefix of output file name"
		    echo "  -h ... help"
			echo
			exit
}


## OPTIONS
while getopts ":i:c:m:1:2:o:h" OPTION; do
	 case $OPTION in
		  i ) ID_INDEX="$OPTARG" ;;
		  c ) CLIP_DIR="$OPTARG" ;;
		  m ) MAP_DIR="$OPTARG" ;;
		  1 ) GENOME1="$OPTARG" ;;
		  2 ) GENOME2="$OPTARG" ;;
		  o ) OUTFILEBASE="$OPTARG" ;;
		  h ) usage ;;
		  * ) echo "Unrecognized argument. Use '-h' for usage information."
		  exit 1 ;;
	 esac
done

case "$GENOME1" in
	CBS138 ) SPECIES1=ngla ;;
  ZP591.22 ) SPECIES1=skud ;;
		 * ) echo "Unknown genome "$GENOME1" !!"
			 exit 1 ;;
esac

case "$GENOME2" in
	   R64 ) SPECIES2=scer ;;
     ASMv1 ) SPECIES2=scer ;;	  
		 * ) echo "Unknown genome "$GENOME2" !!"
			 exit 1 ;;
esac

SCRIPTSUFFIX=240318
OUTFILE=${OUTFILEBASE}_${GENOME1}_${GENOME2}_counts${SCRIPTSUFFIX}.txt


## SCRIPT
# print header 
echo -e "Sample id\tTotal\tClipped\t${GENOME2} specific\t${GENOME1} specific\tCommon" | tee $OUTFILE

# go through files in input directory
for INPUTPATH in $(ls ${CLIP_DIR}/*.out); do
   INPUTFILE=$(echo $INPUTPATH | awk -F/ '{print $NF}')
   SAMPLE_ID=$(echo $INPUTFILE | cut -d'_' -f${ID_INDEX})  # extract sample id
   MAPLOGFILE=*${SAMPLE_ID}*.log   
  	
   TOTAL=$(grep 'Total read pairs processed' $INPUTPATH | xargs | cut -d' ' -f5 | tr -dc '0-9') # extract numbers only
   CLIPPED=$(grep 'Pairs written' $INPUTPATH | xargs | cut -d' ' -f5 | tr -dc '0-9')  # remove all whitespaces except one space between words
   SPEC1COUNT=$(grep "_${GENOME1}_woCommon-${GENOME2}.sam aligned" ${MAP_DIR}/log_files/$MAPLOGFILE | cut -d':' -f2) 
   SPEC2COUNT=$(grep "_${GENOME1}_unaligned_${GENOME2}.sam aligned" ${MAP_DIR}/log_files/$MAPLOGFILE | cut -d':' -f2)
   COMMONCOUNT=$(grep "_${GENOME1}_${GENOME2}.sam aligned" ${MAP_DIR}/log_files/$MAPLOGFILE | cut -d':' -f2)
   
   echo -en "$SAMPLE_ID\t$TOTAL\t$CLIPPED\t$SPEC2COUNT\t$SPEC1COUNT\t$COMMONCOUNT\n" | tee -a $OUTFILE
done
 

exit 0
