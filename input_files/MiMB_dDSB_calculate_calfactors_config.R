## Input for calculate_calibration_factors_commandLine.R, v.240724

## read count table *_counts*.txt (Sample id, Total, Clipped, EXP_GENOME specific, CAL_GENOME specific)
inputfileR <- "MiMB_dDSB_ZP591.22_ASMv1_counts240318.txt" 
idColumnR <- "Sample.id"  # replace space by '.' !!

# sample table with IP/WCE information and sample names
workingDirINS <- paste0(workingDirINR,"/input_files")
inputfileS <- "MiMB_dDSB_sample_table.txt"
idColumnS <- "SRA.id"   # replace space by '.' !!
nameColumnS <- "Sample.name"  # replace space by '.' !!
typeColumnS <- "Type"  # replace space by '.' !!
parentColumnS <- "Parent.id"  # replace space by '.' !!; has to match idColumnS and idColumnR !!

genomeExp <- "ASMv1"
genomeCal <- "ZP591.22"
addPerc <- TRUE    # add column with % of total read pairs
rpmCount <- 10000000   # count to which experimental genome specific read counts are normalized to
saveTables <- TRUE