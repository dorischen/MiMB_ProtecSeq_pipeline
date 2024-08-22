## Input for convert_sam.pe_commandLine.R, v.240416


# reference genome
genome <- "ASMv1" 
genomeDataPath <- paste0(workingDirINC,"/genomes/config_file_genomes.R")   # with chromosome ids (saved in output files), names (as in fasta and alignment files), and lengths

# sample file (with columns in the following order: sample id), sample name, color, rank; sample id should be a unique identifier in the input file name)
workingDirINS <- workingDirINC
inputfileS <- "MiMB_dDSB_sampleInfo_forWiggle.txt"  

# output
workingDirOUT <- workingDirIND  # output directory, parent folder for "wiggle_files"
sepStrands <- FALSE    # FALSE | TRUE .. one wiggle file for each strand (only possible if input file contains column strand)

# normalization (optional)
withCal <- TRUE   # FALSE | TRUE .. signals (depths) will be multiplied by calibration factor provided in factorColumn, output file with suffix '_cal'
withRpm <- FALSE  # FALSE | TRUE .. normalisation according aligned reads, depths multiplied by factor provided in factorColumn, output file with suffix '_rpm'
workingDirINCF <- workingDirINC
inputfileCF <- "MiMB_dDSB_ZP591.22_ASMv1_calFactors240724-rp10M.txt"  # (with column "Cal.Factor")
idColumnCF <- "Sample.id"
factorColumn <- "Cal.Factor"  # e.g. "Cal.Factor" or "RP10M"
digitsAC <- 3  # digits after comma in depth values of wiggle file