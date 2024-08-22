## Input for create_arc_plot, v.240801


# input file (fd format: chrom, start, length, depth)
workingDirIN <- "Z:/forMiMB/test"
workingDirINF <- paste0(workingDirIN,"/mapping_results/scer_specific/fd_files")
filePatternF <- "_fd.txt"
delimF <- "_"  # delimiter in input file name for finding sample id
idIndexF <- 2  # index within input file name where sample id can be found after splitting at delimiter delimC

# reference file, displayed as bars (optional; dpp format: chrom, position, depth), e.g. Spo11 oligo 5prime ends
withRef <- FALSE  # TRUE | FALSE
plotOrder <- "refToFront" # "refToFront" .. bars will be plotted on top of arcs | "refToBack" .. arcs will be plotted on top of bars
workingDirINR <- workingDirINC
inputfileR <- "Spo11.oligos_t4_ASMv1_5prime_dpp.txt"
sampleNameR <- "Mohibullah (2017), wt t4"
sampleNameRF <- "Moh.wt.t4" # will be added to output file name
selColorR <- "#E67300"  # (orange, recommendable for grey-shaded arcs)
optionR <- "sumStrands"  # "sepStrands" .. Watson and Crick strands shown separately above or below x-axis) | "sumStrands" .. both strands shown above x-axis
lineWidthR <- 3    # line width of bars
dppScaleFactorR <- 0.25  # depth scaling factor; if dppScaleFactorR=1, max. depth of reference file = max. y
selColorRT <- add_alpha(selColorR, alpha=1*255)  # alpha .. between 0 and 255; the higher the less transparent
#selColorRT <- add_alpha(selColorR, alpha=0.7*255)  # recommended for plotOrder="refToFront"

# color table (optional); with sample id (same as in input file name) and color (R color name or #rrggbb)
withColorTable <- TRUE  # if FALSE vColors needed !
workingDirINS <- workingDirINC
inputfileS <- "MiMB_dDSB_sample_table.txt"
idColumnS <- "SRA.id"  # column with sample ids for identification of input samples and corresponding colors
colorColumnS <- "Color"
vColors <- c()   # needed if colorTable=FALSE; default: "black"

# calibration factor table (optional)
withCal <- TRUE
workingDirINCF <- workingDirINC
inputfileCF <- "MiMB_dDSB_ZP591.22_ASMv1_calFactors240724-rp10M.txt"
idColumnCF <- "Sample.id"   # column with sample ids for identification of input samples and corresponding calibration factors
factorColumnCF <- "Cal.Factor"  # column containing calibration factor

# plot settings
plotOption <- "pdf"   # if not "pdf", plots will be printed to the R Console or RStudio Plots pane
workingDirOUT <- paste0(workingDirINF, "/arcPlots")
plotHeight <- 35 

# plot contrast settings
fixedPlotWidth <- 76  # plot width in inch; recommendable in case comparison between multiple samples needed; if 0 -> plot width calculated according data
circleLwd <- 0.1  # line width for arcs; default: 0.1
fixedMaxDepth <- 400  # set to 0 if max. depth should be extracted from the (region-specific) data; will be rounded

# chromosomal regions to plot (chrom, start, end)
selRegionList <- list(c("7" ,"833500", "834000"))  # multiple regions (vectors) possible
vAlphaOffset <- c(0.01)   # between 0 and 1; so that depth=1 is visible (the higher, the darker); can be set for each region separately (if the same for all regions, just indicate one value)


