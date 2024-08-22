## Hi-C-type heatmap representation of dDSB fragments
## INPUT: fd file (chrom, start, length, depth)
## OUTPUT: Hi-C-style 2D heatmap plot (pdf file)
## Author: Doris Chen, Franz Klein
## version 240424


## PACKAGES
## get library path
libLoc <- .libPaths()[1]   # [1] in case of multiple library locations
repository <- "https://cloud.r-project.org/"

## install packages if not yet installed and load
if(!require("unix", quietly=TRUE))
  install.packages("unix", lib=libLoc, repo=repository, dependencies=TRUE)
library(unix, lib.loc=libLoc)
rlimit_as(1e15)  # memory limit increases to ~96GB

if(!require("data.table", quietly=TRUE))
  install.packages("data.table", lib=libLoc, repo=repository, dependencies=TRUE) 
library(data.table, lib.loc=libLoc)


## SUB-FUNCTIONS
process_dpp_c <- function(dppTable, option="", cal=FALSE, calFactor=1)
{  setkey(dppTable, chrom, position)
   if(dppTable[position<0, .N]>0)
     warning("Negative positions !")

   if(option=="sum")
   {  # reduce to 3 columns
      if(dim(dppTable)[2]==3)
      {  vDepthSum <- dppTable[,3,with=FALSE]
      } else
      {  vDepthSum <- dppTable[,3,with=FALSE] + dppTable[,4,with=FALSE]
      }
      resultTable <- dppTable[,1:2,with=FALSE]
      resultTable[,depth:=vDepthSum]
      if(cal)
        resultTable[,depth:=depth*calFactor]
   } else
   {  resultTable <- dppTable
      if(cal)
        resultTable[,norm_depth_plus:=norm_depth_plus*calFactor, norm_depth_minus:=norm_depth_minus*calFactor]
   }

   return(resultTable)
}

process_fd_c <- function(fdTable, createIndex=FALSE, withCal=FALSE, calFactor=1)  # addition column "end", optionally addition of "index", optionally multiplying depth with calFactor
{  resultTable <- data.table::copy(fdTable)
   resultTable[,end := start + length - 1L]
   if(createIndex)
   {  resultTable[,index := seq(1, nrow(resultTable))]
      setkey(resultTable, index)
   }
   if(withCal)
     resultTable[,depth:=depth*calFactor]
   return(resultTable)
}

add_alpha <- function(vColors, alpha=255)  # alpha between 0 and 1!!
{  if(missing(vColors))
     stop("Please provide a vector of colours.")
   
   vColorsT <- as.character(apply(sapply(vColors, col2rgb), 2, function(x) rgb(x[1], x[2], x[3], maxColorValue=255, alpha=alpha)))
   return(vColorsT)
}



## INPUT 
## COMMAND-LINE INPUT
argCount <- 5
args <- commandArgs(trailingOnly=TRUE)
print("Parameters:")
print(args)
if(length(args)<argCount)
{  print("Rscript create_HiCstyle_plot_commandLine.R CONFIG_FOLDER CONFIG_FILE INPUT_FOLDER INPUT_FILE(fd file,*_fd.txt) SAMPLE_ID")
}
workingDirINC <- file.path(args[1])
inputfileC <- args[2]
workingDirINF <- file.path(args[3])
inputfileF <- args[4]
sampleId <- args[5]

# ## INTERACTIVE INPUT
# workingDirINC <- "Z:/forMiMB/input_files"
# inputfileC <- "MiMB_dDSB_create_HiC_plot_ASMv1_config.R"
# workingDirINF <- "Z:/forMiMB/mapping_results/scer_specific/fd_files"
# inputfileF <- "sp101.IP_SRR14093065_Spo11-MYC18_wt_t4_1358_ASMv1_fd.txt"  
# sampleId <- "SRR14093065"

source(paste0(workingDirINC,"/",inputfileC))  



## DEFAULTS
# creation of output directory (if not existing yet)
if(!file.exists(workingDirOUT))
{  dir.create(workingDirOUT, recursive=TRUE)
   print(paste(workingDirOUT,"created"))
}

# output file suffix
scriptSuffix <- "240424"
fileSuffix <- paste0("_HiC",scriptSuffix)
if(maxY0>0)
  fileSuffix <- paste0(fileSuffix,"-maxY",maxY0)
if(!rainbow)
  fileSuffix <- paste0(fileSuffix,"-grey")
if(plusDSB)
  fileSuffix <- paste0(fileSuffix,"-",sampleNameDF,"DSB")
if(withCal)
{  fileSuffix <- paste0(fileSuffix,"_cal")
} else
if(withRpm)
  fileSuffix <- paste0(fileSuffix,"_rpm")

# calculation of data-derived plot width
if(plotWidth0==0)
  plotWidth <- 8 + (as.numeric(selRegion[3]) - as.numeric(selRegion[2])) * 7/6000



## FUNCTION
# load DSB data (optional)
if(plusDSB)
{  # get reference files
   setwd(workingDirIND)
   dppTableSource <- fread(inputfileD, header=TRUE, sep="\t", stringsAsFactors=FALSE)
   dppTable <- process_dpp_c(dppTableSource, option="sum")  # chrom - position - depth
   print(dppTable)
}

# get calibration or normalisation factor (optional)
if(withCal | withRpm)
{  setwd(workingDirINCF)
   calTable <- fread(inputfileCF, header=TRUE, sep="\t", stringsAsFactors=FALSE)
   setnames(calTable, gsub(" ", ".", names(calTable)))
   calTable 
   
   normFactor <- calTable[base::get(idColumnCF)==sampleId, base::get(factorColumn)]
   if(length(normFactor)>0)
   {  print(paste("Normalisation factor:",normFactor))
   } else
     stop(paste0("Sample '",sampleId, "' not found in column '",idColumnCF,"' of table ",workingDirINCF,"/",inputfileCF))
}
  
# create output file name base
fileNameBase <- gsub("_fd.txt", "", inputfileF)
currChrom <- selRegion[1]
currStart <- as.numeric(selRegion[2])
currEnd <- as.numeric(selRegion[3])
print(paste0("chrom", currChrom, ":", currStart, "-", currEnd))
fileNameBaseP <- paste0(fileNameBase,fileSuffix,"_chr",currChrom,".",currStart,"-",currEnd)
print(fileNameBaseP)

# get fd table (chrom - start - length - depth - end)
setwd(workingDirINF)
fdTable0 <- fread(inputfileF, header=TRUE, sep="\t", stringsAsFactors=FALSE)
fdTable <- process_fd_c(fdTable0, withCal=(withCal|withRpm), calFactor=normFactor)
print(fdTable)

# get fd in region
fdTableSel <- fdTable[chrom==currChrom & start>=currStart & end<=currEnd,]
fdTableSel <- fdTableSel[order(fdTableSel$depth),]  # for plotting order (strongest fragments last)
print(fdTableSel)

# get dpp in region
if(plusDSB)
{  dppTableSel <- dppTable[chrom==currChrom & position>=currStart & position<=currEnd,]
   dppTableSel[,depth:=depth*scaleFactorD]
   
   # color opacity 
   selColorDT <- add_alpha(selColorD, alpha=alphaFactorD*255)  # the lower alpha value, the higher the transparence
}

# set max. depth
maxDepth <- ceiling(max(fdTableSel$depth))

# create plot title (with output file name base, max. depth and configurtion file path)
normString <- ""
if(withCal)
{  normString <- paste0(", calibration factor=",round(normFactor,3))
} else
if(withRpm)
{  normString <- paste0(", normalisation factor=",round(normFactor,3))
}
plotTitle <- paste0(fileNameBaseP,"\n(max. depth=",maxDepth,normString,")\n",paste0(workingDirINC,"/",inputfileC))
print(plotTitle)

# create color palette
if(rainbow)
{  makeramp <- colorRampPalette(c(rgb(0,0,1,0.2),rgb(0,1,0,0.4), rgb(1,1,0,0.6),rgb(1,0.7,0,0.7),rgb(1,0.3,0,0.8), rgb(1,0.15,0,0.9), rgb(1,0,0,1)), alpha=TRUE, bias=4)
} else
  makeramp <- colorRampPalette(c(rgb(0,0,0,0.2),rgb(0,0,0,0.4), rgb(0,0,0,0.6),rgb(0,0,0,0.8), rgb(0,0,0,1)), alpha=TRUE, bias=1)
if(fixedColCount==0) 
{  fixedColCount <- maxDepth
}
colorRange <- makeramp(fixedColCount)
colCountLog <- ceiling(log2(fixedColCount))
selColCounts <- 2^seq(0,colCountLog)  # for legend

# generate plot parameters
factorX <- 10^floor(log10(max(fdTableSel$start)-min(fdTableSel$end)))
minX <- currStart
maxX <- currEnd
factorY <- 10^floor(log10(max(fdTableSel$length)-min(fdTableSel$length))) * 2
if(maxY0==0)
{  maxY <- ceiling(max(fdTableSel$length)/factorY)*factorY
} else
  maxY <- maxY0
   
if(plotHeight0==0)
  plotHeight <- max(3.5, plotWidth * maxY/(as.numeric(selRegion[3])-as.numeric(selRegion[2])))

# create plot and save in pdf file (old pdf with same name is overwritten)
oldpar <- par()  # save old plot settings

setwd(workingDirOUT)
outfile <- paste0(fileNameBaseP,".pdf")
if(file.exists(outfile))
{  file.remove(outfile)
   print(paste(outfile,"removed"))
}
pdf(file=outfile, bg="transparent", height=plotHeight, width=plotWidth) 
par(mar=c(2,6,5,2), cex=1, xaxs="i", yaxs="i", bty="n", mgp=c(0,0.8,0))  

# preparation of plot
labelSize <- 1.7
axisFontSize <- 1.4
plot(0, 0, xlim=c(minX,maxX), ylim=c(0,maxY), pch=".", main=plotTitle, font.main=1, cex.main=0.8, cex.lab=labelSize, xaxt="n", yaxt="n", xlab="", ylab="", asp=0.5)
mtext(text=paste("Position on chromosome",currChrom,"(kb)"), side=1, cex=labelSize, line=1)
axis(side=1, at=seq(minX, maxX, factorX), labels=seq(minX, maxX, factorX)/factorX, tcl=-0.5, cex.axis=axisFontSize, lwd=2, pos=0) # thick ticks
#axis(side=1, pos=yCross-0.02*yLength, at=seq(minX, maxX, factorX/10), cex.axis=1, lwd=0) # labels of thick ticks
axis(side=1, at=seq(minX, maxX, factorX/10), labels=FALSE, lwd=1, pos=0) # thin ticks
mtext(text="dDSB fragment\nlength (bp)", side=2, line=2.5, cex=labelSize)
axis(side=2, at=seq(0, maxY, factorY), labels=TRUE, cex.axis=axisFontSize, lwd=2, pos=minX) 

# DSBs (optional)
if(plusDSB & plotOrder=="DSBtoBack")
{  lines(x=dppTableSel$position, y=dppTableSel$depth, type="h", lwd=lineWidthD, col=selColorDT, asp=1)
   legend(x="topleft", inset=c(0,0), sampleNameD, cex=labelSize, adj=c(0,0.2), pt.cex=labelSize*1.5, pch="-", col=selColorD, text.col=selColorD, bty="n")  # legend at top left
}

# dDSBs
points((fdTableSel$start+fdTableSel$length/2), fdTableSel$length, pch=18, col=colorRange[fdTableSel$depth], cex=symbolSize) 

# DSBs (optional)
if(plusDSB & plotOrder=="DSBtoFront")
{  lines(x=dppTableSel$position, y=dppTableSel$depth, type="h", lwd=lineWidthD, col=selColorDT, asp=1)
   legend(x="topleft", inset=c(0,0), sampleNameD, cex=labelSize, adj=c(0,0.2), pt.cex=labelSize*1.5, pch="-", col=selColorD, text.col=selColorD, bty="n")  # legend at top left
}

# legend
legend(x="topright", horiz=TRUE, legend=as.character(selColCounts), col=colorRange[selColCounts], pch=15, bty="n", pt.cex=labelSize, cex=axisFontSize)

# saving of pdf file
dev.off()
print(paste(outfile,"saved in",workingDirOUT))


par(oldpar)  # reset par settings







