#!/bin/bash
## INPUT: sample table (with Sample Id, Index1 Id, Index2 Id)
## INPUT: adapter table (with Id, Index (barcode sequence), Sequence (adapter sequence), IndexRC (index rev. compl.), SequenceRC (adapter rev. compl.)
## OUTPUT: cutadapt command lines
## Author: Doris Chen (2016)
## version 240712




function usage()
{  echo
   echo "Create cutadapt commands"
   echo
   echo "Usage: $0 -d INPUT_DIR -i SAMPLE_TABLE -a ADAPTER_TABLE -r READ_LENGTH -s OUT_FILE_SUFFIX(e.g. 160405) -o MIN_CLIP_LENGTH -e ERROR_RATE(e.g. 0.1) -m MIN_CLIPPED_READ_LENGTH -n MAX_N_COUNT -t NR_OF_THREADS"
	  echo
			echo "  -d ... INPUT_DIR with fastq.gz files, also used as output directory"
			echo "  -i ... SAMPLE_TABLE with sample name (prefix before '(_1).fastq.gz', for identification of fastq files), adapter index id 1, adapter index id 2"
			echo "  -a ... ADAPTER_TABLE with adapter id, index, sequence, index rev. compl., sequence rev. compl."
			echo "  -r ... read length" 
			echo "  -s ... suffix added to output file name"
			echo "  -o ... minimal clipping length"
			echo "  -e ... error rate (0-1)"
			echo "  -m ... minimal clipped read length"
			echo "  -n ... maximal N count tolerated (otherwise read will be discarded)"
			echo "  -t ... number of threads/CPU cores to be used"
			echo "  -h ... help"
			echo
			exit
}


## OPTIONS
while getopts ":d:i:a:r:s:o:e:m:n:t:h" OPTION; do
	 case $OPTION in
	      d ) INPUTDIR="$OPTARG" ;;
		  i ) SAMPLETABLE="$OPTARG" ;;
		  a ) ADAPTERTABLE="$OPTARG" ;;
		  r ) READLENGTH="$OPTARG" ;;
		  s ) SUFFIX="$OPTARG" ;;
		  o ) O="$OPTARG" ;;
		  e ) E="$OPTARG" ;;
		  m ) M="$OPTARG" ;;
		  n ) N="$OPTARG" ;;
		  t ) THREADS="$OPTARG" ;;
		  h ) usage ;;
		  * ) echo "Unrecognized argument. Use '-h' for usage information."
		  exit 1 ;;
	 esac
done


## DEFAULTS

## SCRIPT
EP$(echo "$E" | awk '{print $1 * 100}')  # error percentage   

while read LINE; do SAMPLE=$(echo $LINE | cut -d' ' -f1); if (ls $INPUTDIR | grep -q $SAMPLE); then INDEX1=$(echo $LINE | cut -d' ' -f2); SEQ1=$(grep -w $INDEX1 $ADAPTERTABLE | cut -f3); SUBSEQ1=${SEQ1:0:((READLENGTH+1))}; INDEX2=$(echo $LINE | cut -d' ' -f3); SEQ2=$(grep -w $INDEX2 $ADAPTERTABLE | cut -f5); SUBSEQ2=${SEQ2:0:((READLENGTH+1))};\
    CMD="cutadapt --cores $THREADS -a $SUBSEQ1 -A $SUBSEQ2 -O $O -e $E -m $M --max-n $N -o ${INPUTDIR}/${SAMPLE}_R1_clipped-aO${O}e${EP}m${M}N${N}.fastq.gz -p ${INPUTDIR}/${SAMPLE}_R2_clipped-aO${O}e${EP}m${M}N${N}.fastq.gz ${INPUTDIR}/${SAMPLE}_1.fastq.gz ${INPUTDIR}/${SAMPLE}_2.fastq.gz | tee ${INPUTDIR}/${SAMPLE}_cutadapt-aA${INDEXTYPE}-O${O}-e${EP}-m${M}-N${N}_${SUFFIX}.out"; \
	echo $CMD; fi; done < $SAMPLETABLE	 
   
  
exit 0
