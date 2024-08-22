## Input for convert_sam.pe_commandLine.R, v.240727

# reference genome
genome <- "ASMv1" 
genomeDataPath <- paste0(workingDirINC,"/genomes/config_file_genomes.R")   # add absolute path; with chromosome ids (saved in output files), names (as in fasta and alignment files), and lengths

# alignment format
inputFormat <- "ngm"  # "ngm" | "hisat"; alignment tool

# output
workingDirOUT <- workingDirINS  # output directory, parent folder for subfolders (/dpp_files, /filledDpp_files, /fd_files)
outputOptions <- c("filled", "fd", "5prime")   # "filled", "fd" and/or "5prime"

saveTableParsed <- FALSE   # if TRUE, saving of alignment table
saveLengthTable <- TRUE   # if TRUE, saving of fragment length distribution table
saveLengthPlot <- TRUE   # if TRUE, pdf file of legnth distribution created
saveTableOutward <- FALSE  # if TRUE, outward oriented read pair alignments saved in file *_outward.txt

# further parameters
minScore <- 0  # applied if > 0
maxRowCountPC <- 3000000   # for splitting up of large files and faster processing
maxLengthDist <- 400  # max. length in length distribution plot; if 0 .. derived from data; default: 500
