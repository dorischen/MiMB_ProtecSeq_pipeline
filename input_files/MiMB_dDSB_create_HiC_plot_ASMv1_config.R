## Input for create_HiCstyle_plot_commandLine.R, v.240424

# input
filePatternF <- "_fd.txt"

# output directory
workingDirOUT <- paste0(workingDirINF,"/HiCplots")

# plot
selRegion <- c("7", "832000", "836000") # chrom, start, end
plotHeight0 <- 0 # inch; if 0, derived from data
plotWidth0 <- 0  # inch; if 0, derived from data
maxY0 <- 1000  # max. Y value (dDSB fragment length) in plot; if '0', derived from the data
symbolSize <- 0.6

# colors
rainbow <- TRUE    # if FALSE, then grey
fixedColCount <- 512  # number of colors in gradient; 0 for automatic count

# normalization (optional)
withCal <- FALSE   # FALSE | TRUE .. signals (depths) will be multiplied by calibration factor provided in factorColumn, output file with suffix '_cal'
withRpm <- FALSE  # FALSE | TRUE .. normalisation according aligned reads, depths multiplied by factor provided in factorColumn, output file with suffix '_rpm'
workingDirINCF <- workingDirINC
inputfileCF <- "MiMB_dDSB_ZP591.22_ASMv1_calFactors240724-rp10M.txt"  # (with column "Cal.Factor")
idColumnCF <- "Sample.id"
factorColumn <- "Cal.Factor"  # e.g. "Cal.Factor" or "RP10M"

# additional DSB track (vertical lines)
plusDSB <- FALSE  # FALSE or TRUE
workingDirIND <- workingDirINC
inputfileD <- "Spo11.oligos_t4_ASMv1_5prime_dpp.txt"
sampleNameD <- "Mohibullah (2017),\nwt t4"
sampleNameDF <- "Moh1" # for output filename
scaleFactorD <- 0.6  # scaling factor for signal depths
lineWidthD <- 2
selColorD <- "grey40" 
alphaFactorD <- 0.6  # 0 - 1; the higher, the lower the transparence
plotOrder <- "DSBtoBack"  # "DSBtoFront" | "DSBtoBack"
