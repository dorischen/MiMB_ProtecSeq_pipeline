## Input for convert_sam.pe_commandLine.R, v.240503

# reference genome
genome <- "ZP591.22" 
genomeDataPath <- paste0(workingDirINC,"/genomes/config_file_genomes.R")   # with chromosome ids (saved in output files), names (as in fasta and alignment files), and lengths

# alignment format
inputFormat <- "hisat"  # "ngm" | "hisat"; alignment tool

# output
workingDirOUT <- workingDirINS  # output directory, parent folder for subfolders (/dpp_files, /filledDpp_files, /fd_files)
outputOptions <- c("fd")   # "filled", "fd" and/or "5prime"

saveTableParsed <- FALSE   # if TRUE, saving of alignment table
saveLengthTable <- FALSE   # if TRUE, saving of fragment length distribution table
saveLengthPlot <- TRUE   # if TRUE, pdf file of legnth distribution created
saveTableOutward <- FALSE  # if TRUE, outward oriented read pair alignments saved in file *_outward.txt

# further parameters
minScore <- -3  # alignment tool specific !!
maxRowCountPC <- 3000000   # for splitting up of large files and faster processing
maxLengthDist <- 400  # max. length in length distribution plot; if 0 .. derived from data; default: 500


