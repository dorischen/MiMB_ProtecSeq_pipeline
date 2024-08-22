##!/bin/bash
## INPUT: (clipped) fastq(.gz) files containing reads from two species (experimental and calibration reference genomes)
## OUTPUT: species-specific and common alignments 
## Author: Doris Chen (April, 2016)
## version 240821


## FUNCTIONS
source get_specific_alignments_pe_config.sh  # containing paths to tools

function ngm_align_paired()
{  REFGENOME=$1
   SPECIES=$2
   INFILE1=$3
   INFILE2=$4
   PREFIX=$5
   OUTSUFFIX=$6
   LOGFILE=$7
   THREADS=$8
   ID=$9
   HITS=${10}  # only 1 allowed with ngm
   MINLENGTH=${11}
   MAXINSERT=${12}
   KMERLENGTH=${13}
   CONFIG=${14}
   CSUFFIX=${15}
  
   LONGLOGFILE=${INPUTDIR}/${PREFIX}${CSUFFIX}_${SPECIES}_${OUTSUFFIX}.log
   
   # align
   if [ "$CONFIG" = "global" ]; then
     ${ALIGNPATH}ngm -r $REFGENOME -1 $INFILE1 -2 $INFILE2 -o ${INPUTDIR}/${PREFIX}${CSUFFIX}_${SPECIES}_withUnal.sam --end-to-end -k $KMERLENGTH --kmer-skip 1 -n $HITS --strata -t $THREADS -i $ID -R 1 -X$MAXINSERT 2>&1 | tee $LONGLOGFILE
   elif [ "$CONFIG" = "local" ]; then
     ${ALIGNPATH}ngm -r $REFGENOME -1 $INFILE1 -2 $INFILE2 -o ${INPUTDIR}/${PREFIX}${CSUFFIX}_${SPECIES}_withUnal.sam --local -k $KMERLENGTH --kmer-skip 1 -n $HITS --strata -t $THREADS -i $ID -R $MINLENGTH -X$MAXINSERT 2>&1 | tee $LONGLOGFILE
   fi
   
   # create short log file
   sed '/Time:/d' $LONGLOGFILE | tee -a $LOGFILE;     
   rm $LONGLOGFILE
}

function hisat_align_paired()
{  REFGENOME=$1
   SPECIES=$2
   INFILE1=$3
   INFILE2=$4
   PREFIX=$5
   OUTSUFFIX=$6
   LOGFILE=$7
   THREADS=$8
   HITS=$9  
   MAXINSERT=${10}
   CONFIG=${11}
   CSUFFIX=${12}
   
   LONGLOGFILE=${INPUTDIR}/${PREFIX}${CSUFFIX}_${SPECIES}_${OUTSUFFIX}.log
   
   # align
   if [ "$CONFIG" = "global" ]; then
     ${ALIGNPATH}hisat2 -x $REFGENOME -q -1 $INFILE1 -2 $INFILE2 --no-spliced-alignment --end-to-end --very-sensitive -X $MAXINSERT -k $HITS --no-discordant --no-mixed --threads $THREADS --met-file ${INPUTDIR}/${PREFIX}${CSUFFIX}_${SPECIES}_metrics.txt -S ${INPUTDIR}/${PREFIX}${CSUFFIX}_${SPECIES}_withUnal.sam 2>&1 | tee $LONGLOGFILE
   elif [ "$CONFIG" = "local" ]; then
     ${ALIGNPATH}hisat2 -x $REFGENOME -q -1 $INFILE1 -2 $INFILE2 --no-spliced-alignment --very-sensitive -X $MAXINSERT -k $HITS --no-discordant --no-mixed --threads $THREADS --met-file ${INPUTDIR}/${PREFIX}${CONFIG}_${SPECIES}_metrics.txt -S ${INPUTDIR}/${PREFIX}${CSUFFIX}_${SPECIES}_withUnal.sam 2>&1 | tee $LONGLOGFILE
   fi    
   
   # create short log file
   sed '/Time:/d' $LONGLOGFILE | tee -a $LOGFILE;     
   rm $LONGLOGFILE
}

function extract_aligned_and_count()
{  OUTFILEBASE=$1
   INFILE=$2
   LOGFILE=$3
   OUTFILE=${OUTFILEBASE}.sam   
  
   # extract and count
   awk '$1 ~ /^@/ || $2==83 || $2==99 || $2==147 || $2==163' $INFILE > $OUTFILE   # include only properly paired
   echo "$OUTFILE aligned pairs: $(($(${SAMTOOLSPATH}samtools view -c ${OUTFILE})/2))" | tee -a $LOGFILE
}

function extract_unaligned_and_count()
{  OUTFILEBASE=$1
   INFILE=$2
   LOGFILE=$3
   OUTFILE=${OUTFILEBASE}_unaligned.sam
   
   # extract
   awk '$1 ~ /^@/ || !($2==83 || $2==99 || $2==147 || $2==163 || $2==339 || $2==355 || $2==403 || $2==419)' $INFILE > $OUTFILE  # don't include properly paired
   echo "$OUTFILE unaligned pairs: $(($(${SAMTOOLSPATH}samtools view -c ${OUTFILE})/2))" | tee -a $LOGFILE
}

function extract_ids_and_remove()
{  OUTFILEBASE=$1
   COMMONFILE=$2
   ALIGNFILE=$3
   SPECIES2=$4
   LOGFILE=$5
   MEM=$6
   OUTFILESUFFIX=${OUTFILEBASE}_woCommon-${SPECIES2}
   OUTFILE=${OUTFILESUFFIX}.sam
   
   # extract reads ids from common
   awk '!($1 ~ /^@/)' $COMMONFILE | cut -f1 | sort | uniq > ${OUTFILEBASE}_reads.txt
   java -Xmx${MEM}000m -jar ${JAVASAMPATH}SAMManipulator.jar -sam $ALIGNFILE -ids ${OUTFILEBASE}_reads.txt -option delete -out $OUTFILESUFFIX
   
   echo "$OUTFILE aligned pairs: $(($(${SAMTOOLSPATH}samtools view -c ${OUTFILE})/2))" | tee -a $LOGFILE
   
   rm ${OUTFILEBASE}_reads.txt  
}

function convert_sam_to_fastq_seqtk()  # better for unaligned pairs
{  INFILEBASE=$1
   
   # extract ids from sam file
   IDFILE=${INFILEBASE}_ids.txt
   awk '!($1 ~ /^@/)' ${INFILEBASE}.sam | cut -f1 | sort | uniq > ${IDFILE}
   
   # extract fastq based on ids
   sed -e 's/$/\/1/' $IDFILE > ${IDFILE}.1
   ${SEQTKPATH}seqtk subseq $INPUTPATH1 ${IDFILE}.1 > ${INFILEBASE}_R1.fastq
   rm ${IDFILE}.1
   gzip ${INFILEBASE}_R1.fastq
 
   sed -e 's/$/\/2/' $IDFILE > ${IDFILE}.2
   ${SEQTKPATH}seqtk subseq $INPUTPATH2 ${IDFILE}.2 > ${INFILEBASE}_R2.fastq
   rm ${IDFILE}.2 
   gzip ${INFILEBASE}_R2.fastq  

   wc -l ${IDFILE}
   rm ${IDFILE}  
}

function convert_sam_to_fasta()
{  INFILEBASE=$1
   cat ${INFILEBASE}.sam | grep -v ^@ | awk '{print ">"$1"\n"$10"\n"}' > ${INFILEBASE}.fa 
}

function usage()
{  echo
   echo "Get species-specific alignments"
   echo
   echo "Usage: $0 -i INPUT_PATH_1 -j INPUT_PATH_2 -o OUTPUT_BASE -s OUTPUT_SUFFIX -1 SPECIES_1(e.g. ZP591.22|CBS138) -2 SPECIES_2(e.g. ASMv1|R64) -t NR_OF_THREADS -m RAM_MEMORY_GB(e.g. 64) -n NGM|HISAT -e MIN_PERC_IDENTITY(e.g. 0.95) -c [global|local] -l MIN_PERC_READLENGTH(e.g. 0.9) -r [TRUE|FALSE] [-k KMER_LENGTH(10-14)] [-z MAX_INSERT_SIZE(default:1000)] [-x MAX_HITS(>1 only for HISAT2] [-f JOB_NR_1(0-10) -g JOB_NR_2(0-10)]"
	  echo
			echo "  -i ... read 1 fastq(.gz) file path; file name should start with SAMPLE_NAME and contain 'R1'"
			echo "  -j ... read 2 fastq(.gz) file path; file name should start with SAMPLE_NAME and contain 'R2'"
			echo "  -o ... base of output filenames"
			echo "  -s ... output suffix, e.g. 240122; will be appended to log file"
			echo "  -1 ... first genome to align to, calibration genome"
			echo "  -2 ... second, specific genome to align to, 'experimental' genome"
			echo "  -t ... number of threads"
			echo "  -m ... GB max. RAM"  # for SAMManipulator.jar
			echo "  -n ... alignment tool, NextGenMap ('NGM') or HISAT2 ('HISAT')"
			echo "  -e ... minimum percent identity for alignment (max.=1)"
			echo "  -c ... either local or global alignment"
			echo "  -l ... minimum percent read length required (max.=1); in case of global alignment no difference"
			echo "  -r ... if TRUE, all unnecessary files removed"
			echo "  -k ... k-mer length (integer 10 - 14) in case of NextGenMap alignment tool, the smaller the more sensitive and slower is the alignment"
			echo "  -z ... maximum insert size, default=1000nt"
			echo "  -x ... maximum number of hits in genome (for multiple-hit read pears), only valid for HISAT2 (NGM only allows for uniquely aligned read pairs)"
			echo "  -f ... job nr. to start with; 0-2..align to first genome and count; 3-6..extract unaligned and align to second genome and count, exp. specific alignment; 7-9..extract aligned and align to second genome and count; 10 .. remove common reads, cal. specific alignment"
			echo "  -g ... job nr. to end with"
			echo "  -h ... help"
			echo
			exit
}


## OPTIONS
while getopts ":i:j:o:s:1:2:t:m:n:e:c:l:r:k:z:x:f:g:h" OPTION; do
	 case $OPTION in
		  i ) INPUTPATH1="$OPTARG" ;;
		  j ) INPUTPATH2="$OPTARG" ;;
		  o ) OUTBASE="$OPTARG" ;;
		  s ) OUTSUFFIX="$OPTARG" ;;
		  1 ) SPECIES1="$OPTARG" ;;
		  2 ) SPECIES2="$OPTARG" ;;
		  t ) THREADS="$OPTARG" ;;
		  m ) MEM="$OPTARG" ;;
		  n ) TOOL="$OPTARG" ;;
		  e ) MINID="$OPTARG" ;;
		  c ) CONFIG="$OPTARG" ;;
		  l ) MINLENGTH="$OPTARG" ;;
		  r ) CLEAN="$OPTARG" ;;
		  k ) KMERLENGTH="$OPTARG" ;;
		  z ) MAXINSERT="$OPTARG" ;;
		  x ) MAXHITS="$OPTARG" ;;
		  f ) FROM="$OPTARG" ;;
		  g ) TO="$OPTARG" ;;
		  h ) usage ;;
		  * ) echo "Unrecognized argument. Use '-h' for usage information."
		  exit 1 ;;
	 esac
done


## DEFAULTS
SCRIPTSUFFIX=240821

if [ "$MEM" = "" ]; then
  echo "-m RAM_MEMORY_GB missing !!"  # otherwise later cryptic error message when running Java script
  exit 1
fi

INPUTDIR=$(echo $INPUTPATH1 | rev | cut -d'/' -f2- | rev)
INPUTFILE1=$(echo $INPUTPATH1 | awk -F/ '{print $NF}')
INPUTFILE2=$(echo $INPUTPATH2 | awk -F/ '{print $NF}')
EMPTYSUFFIX=

if [ "$MAXINSERT" = "" ]; then
  MAXINSERT=1000
fi
if [ "$FROM" = "" ] || [ "$TO" = "" ]; then
  FROM=0
  TO=10
fi

if [ "$TOOL" = "HISAT" ]; then
  ALIGNSUFFIX=hisat
  if [ "$MAXHITS" = "" ]; then
    MAXHITS=100
  fi
  CONFIGSUFFIX=-${CONFIG}n${MAXHITS}
elif [ "$TOOL" = "NGM" ]; then
  ALIGNSUFFIX=ngm
  MAXHITS=1
  CONFIGSUFFIX=-${CONFIG}R${MINLENGTH}i${MINID}k${KMERLENGTH}
fi

# create log file
LOGFILE=${INPUTDIR}/specAlign-${SCRIPTSUFFIX}_${OUTBASE}_${ALIGNSUFFIX}${CONFIGSUFFIX}_${SPECIES1}-${SPECIES2}_${OUTSUFFIX}.log
LOG_FIRSTLINE="$0 -i $INPUTPATH1 -j $INPUTPATH2 -o $OUTBASE -s $OUTSUFFIX -1 $SPECIES1 -2 $SPECIES2 -t $THREADS -m $MEM -n $TOOL -e $MINID -c $CONFIG -l $MINLENGTH -r $CLEAN -k $KMERLENGTH -z $MAXINSERT -f $FROM -g $TO"
if [ -z "$FROM" ]; then
  echo $LOG_FIRSTLINE | tee $LOGFILE
else
  echo $LOG_FIRSTLINE | tee -a $LOGFILE
fi

# configuration paths
if [ -n "$ALIGNPATH" ] && [[ "$ALIGNPATH" != */ ]]; then ALIGNPATH=${ALIGNPATH}/; fi
if [ -n "$SAMTOOLSPATH" ] && [[ "$SAMTOOLSPATH" != */ ]]; then SAMTOOLSPATH=${SAMTOOLSPATH}/; fi
if [ -n "$SEQTKPATH" ] && [[ "$SEQTKPATH" != */ ]]; then SEQTKPATH=${SEQTKPATH}/; fi
if [ -n "$GENOMEPATH" ] && [[ "$GENOMEPATH" != */ ]]; then GENOMEPATH=${GENOMEPATH}/; fi
echo "ALIGNPATH = "$ALIGNPATH | tee -a $LOGFILE
echo "SAMTOOLSPATH = "$SAMTOOLSPATH | tee -a $LOGFILE
echo "SEQTKPATH = "$SEQTKPATH | tee -a $LOGFILE
echo "JAVASAMPATH = "$JAVASAMPATH | tee -a $LOGFILE
echo "GENOMEPATH = "$GENOMEPATH | tee -a $LOGFILE

# assign reference genome index file
if [ "$TOOL" = "NGM" ]; then
	case "$SPECIES1" in
		CBS138 ) GENOME1=C_glabrata_CBS138.fa ;;
	  ZP591.22 ) GENOME1=Skud-ZP591_GCA947243785.1_221101.fa ;;
			 * ) echo "Unknown species "$SPECIES1" !!"
			 exit 1 ;;
	esac
    case "$SPECIES2" in
		 ASMv1 ) GENOME2=SK1_ASMv1_chroms.fa ;;	  
		   R64 ) GENOME2=S288C-R64-1-1.fa ;;
			 * ) echo "Unknown species "$SPECIES2" !!"
			 exit 1 ;;
	esac
fi
if [ "$TOOL" = "HISAT" ]; then
	case "$SPECIES1" in  
		CBS138 ) GENOME1=C_glabrata_CBS138 ;;
	  ZP591.22 ) GENOME1=Skud-ZP591_GCA947243785.1_221101 ;;
			 * ) echo "Unknown species "$SPECIES1" !!"
				 exit 1 ;;
	esac
    case "$SPECIES2" in
		 ASMv1 ) GENOME2=SK1_ASMv1 ;;	  
		   R64 ) GENOME2=S288C-R64-1-1 ;;
			 * ) echo "Unknown species "$SPECIES2" !!"
				 exit 1 ;;
	esac
fi

# save time and date
T="$(date +%s)"
echo $(date) | tee -a $LOGFILE


## COMMANDS
PREFIX0=${OUTBASE}_${ALIGNSUFFIX}
PREFIX1=${OUTBASE}_${ALIGNSUFFIX}${CONFIGSUFFIX}_${SPECIES1}
PREFIX2=${INPUTDIR}/${PREFIX1}_unaligned_${SPECIES2}
PREFIX3=${INPUTDIR}/${PREFIX1}_${SPECIES2}

# align to first genome
if [ "$TOOL" = "NGM" ]; then
  CMD[0]="ngm_align_paired ${GENOMEPATH}$GENOME1 $SPECIES1 $INPUTPATH1 $INPUTPATH2 $PREFIX0 $OUTSUFFIX $LOGFILE $THREADS $MINID $MAXHITS $MINLENGTH $MAXINSERT $KMERLENGTH $CONFIG $CONFIGSUFFIX"
elif [ "$TOOL" = "HISAT" ]; then
  CMD[0]="hisat_align_paired ${GENOMEPATH}$GENOME1 $SPECIES1 $INPUTPATH1 $INPUTPATH2 $PREFIX0 $OUTSUFFIX $LOGFILE $THREADS $MAXHITS $MAXINSERT $CONFIG $CONFIGSUFFIX"
fi

# extract aligned sequences and count
CMD[1]="extract_aligned_and_count ${INPUTDIR}/${PREFIX1} ${INPUTDIR}/${PREFIX1}_withUnal.sam $LOGFILE"
# extract unaligned sequences
CMD[2]="extract_unaligned_and_count ${INPUTDIR}/${PREFIX1} ${INPUTDIR}/${PREFIX1}_withUnal.sam $LOGFILE"

# convert sam to fastq (unaligned)
CMD[3]="convert_sam_to_fastq_seqtk ${INPUTDIR}/${PREFIX1}_unaligned"
EXT=fastq.gz
# align to second genome
if [ "$TOOL" = "NGM" ]; then
  CMD[4]="ngm_align_paired ${GENOMEPATH}$GENOME2 $SPECIES2 ${INPUTDIR}/${PREFIX1}_unaligned_R1.fastq.gz ${INPUTDIR}/${PREFIX1}_unaligned_R2.fastq.gz ${PREFIX1}_unaligned $OUTSUFFIX $LOGFILE $THREADS $MINID $MAXHITS $MINLENGTH $MAXINSERT $KMERLENGTH $CONFIG $EMPTYSUFFIX"
elif [ "$TOOL" = "HISAT" ]; then
  CMD[4]="hisat_align_paired ${GENOMEPATH}$GENOME2 $SPECIES2 ${INPUTDIR}/${PREFIX1}_unaligned_R1.fastq.gz ${INPUTDIR}/${PREFIX1}_unaligned_R2.fastq.gz ${PREFIX1}_unaligned $OUTSUFFIX $LOGFILE $THREADS $MAXHITS $MAXINSERT $CONFIG $EMPTYSUFFIX"
fi

# extract aligned sequences and count
CMD[5]="extract_aligned_and_count $PREFIX2 ${PREFIX2}_withUnal.sam $LOGFILE"
# extract unaligned sequences
CMD[6]="extract_unaligned_and_count $PREFIX2 ${PREFIX2}_withUnal.sam $LOGFILE"

# convert sam to fastq (aligned)
CMD[7]="convert_sam_to_fastq_seqtk ${INPUTDIR}/${PREFIX1} $PAIRED"

# align to second genome
if [ "$TOOL" = "NGM" ]; then
  CMD[8]="ngm_align_paired $CONFIGFILE ${GENOMEPATH}$GENOME2 $SPECIES2 ${INPUTDIR}/${PREFIX1}_R1.fastq.gz ${INPUTDIR}/${PREFIX1}_R2.fastq.gz $PREFIX1 $OUTSUFFIX $LOGFILE $THREADS $MINID $MAXHITS $MINLENGTH $MAXINSERT $KMERLENGTH $CONFIG $EMPTYSUFFIX"
elif [ "$TOOL" = "HISAT" ]; then
  CMD[8]="hisat_align_paired ${GENOMEPATH}$GENOME2 $SPECIES2 ${INPUTDIR}/${PREFIX1}_R1.fastq.gz ${INPUTDIR}/${PREFIX1}_R2.fastq.gz $PREFIX1 $OUTSUFFIX $LOGFILE $THREADS $MAXHITS $MAXINSERT $CONFIG $EMPTYSUFFIX"
fi

# extract aligned sequences and count
CMD[9]="extract_aligned_and_count $PREFIX3 ${PREFIX3}_withUnal.sam $LOGFILE"

# extract ids of aligned reads to remove them from first alignment (to get specific alignments for first genome)
CMD[10]="extract_ids_and_remove $PREFIX1 ${PREFIX3}.sam ${INPUTDIR}/${PREFIX1}.sam $SPECIES2 $LOGFILE $MEM"


# files to check before starting next command
CHECKFILE[0]=$INPUTPATH
CHECKFILE[1]=${INPUTDIR}/${PREFIX1}_withUnal.sam
CHECKFILE[2]=${INPUTDIR}/${PREFIX1}_withUnal.sam
CHECKFILE[3]=${INPUTDIR}/${PREFIX1}_unaligned.sam
CHECKFILE[4]=${INPUTDIR}/${PREFIX1}_unaligned_R2.fastq.gz
CHECKFILE[5]=${PREFIX2}_withUnal.sam
CHECKFILE[6]=${PREFIX2}_withUnal.sam
CHECKFILE[7]=${INPUTDIR}/${PREFIX1}.sam
CHECKFILE[8]=${INPUTDIR}/${PREFIX1}_R2.fastq.gz
CHECKFILE[9]=${PREFIX3}_withUnal.sam
CHECKFILE[10]=${PREFIX3}.sam


## SCRIPT
for i in $(seq $FROM $TO); do   
  if [ -r ${CHECKFILE[$i]} ]; then
    T="$(date +%s)"; \
	echo $i | tee -a $LOGFILE; \
    echo ${CMD[$i]} | tee -a $LOGFILE; \
    ${CMD[$i]} || exit 1
    T="$(($(date +%s)-T))"; \
    printf "Elapsed time: %02d:%02d:%02d\n" "$((T/3600%24))" "$((T/60%60))" "$((T%60))" | tee -a $LOGFILE; \
  else
    echo "Error: ${CHECKFILE[$i]} not existing!"; \
    exit 1
  fi
done

# remove unneccessary files (optional)
if [ "$CLEAN" = "TRUE" ]; then
  rm ${INPUTDIR}/${PREFIX1}_withUnal.sam
  rm ${INPUTDIR}/${PREFIX1}_unaligned.sam
  rm ${PREFIX2}_withUnal.sam
  rm ${PREFIX3}_withUnal.sam
  rm ${CHECKFILE[4]}
  rm ${CHECKFILE[8]}
  rm ${INPUTDIR}/${PREFIX1}_unaligned_R1.fastq.gz
  rm ${INPUTDIR}/${PREFIX1}_R1.fastq.gz
fi

exit 0