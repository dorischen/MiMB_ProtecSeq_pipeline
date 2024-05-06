#! /bin/bash 
## Run in folder with cutadapt-clipped fastq.gz files for automated sorting of clipped files into subfolder /clipped and cutadapt output files (*.out) into subfolder /clipped/log_files_cutadapt.
## Author: Doris Chen
## version 190816


files=$(ls *_clipped*.fastq.gz 2> /dev/null | wc -l)
if [ "$files" != "0" ]; then
  mkdir -p clipped
  mv *_clipped*.fastq.gz clipped  
  echo "*_clipped*.fastq.gz moved to clipped"
fi

files=$(ls *_cutadapt*.out 2> /dev/null | wc -l)
if [ "$files" != "0" ]; then
  mkdir -p clipped/log_files_cutadapt
  mv *_cutadapt*.out clipped/log_files_cutadapt  
  echo "*_cutadapt*.out moved to clipped/log_files_cutadapt"
fi

files=$(ls cutadapt_* 2> /dev/null | wc -l)
if [ "$files" != "0" ]; then
  rm cutadapt_*
  echo "cutadapt_* log files removed"
fi



exit 0
