## Arc representation of dDSB fragments, shade correlating with depth
## INPUT: fd files (chrom, start, length, depth), color table (optional), opacity setting (alpha factor), DSB reference (dpp) file (optional)
## OUTPUT: arc plots (pdf)
## Author: Doris Chen, Franz Klein
## version 240801


## PACKAGES
## get library path
libLoc <- .libPaths()[1]   # [1] in case of multiple library locations
repository <- "https://cloud.r-project.org/"

## install packages if not yet installed and load
if(!require("data.table", quietly=TRUE))
  install.packages("data.table", lib=libLoc, repo=repository, dependencies=TRUE) 
library(data.table, lib.loc=libLoc)


## SUB-FUNCTIONS
process_dpp <- function(dppTable, option="", cal=FALSE, calFactor=1)  # optional summing of strand-specific depths and reduction to 3 columns; optional calibration
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

process_fd <- function(fdTable, createIndex=FALSE, withCal=FALSE, calFactor=1)  # addition of column "end", optional addition of "index", optional calibration
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

round_number <- function(range, option, value, corrFactor=0) # option="max|min"
{  result <- 0

   factor <- 10^(floor(log10(range))-corrFactor) # [change?!]
   if(option=="min")
   {  result <- floor(value/factor) * factor
   } else
   if(option=="maxR")
   {  result <- round(value/factor) * factor
   }
   if(option=="maxC")
   {  result <- ceiling(value/factor) * factor
   }
   
   return(result)
}

# modified from http://www.magesblog.com/2013/04/how-to-change-alpha-value-of-colours-in.html
add_alpha <- function(vColors, alpha=255)  # alpha between 0 and 255; the lower alpha the more transparent
{  if(missing(vColors))
     stop("Please provide a vector of colours.")
   
   vColorsT <- as.character(apply(sapply(vColors, col2rgb), 2, function(x) rgb(x[1], x[2], x[3], maxColorValue=255, alpha=alpha)))
   return(vColorsT)
}

create_legend_text <- function(depthMaxL)
{  if(depthMaxL<=10)
   {  stepL <- 1
      startDepth <- 1
   } else
   if(depthMaxL<20)
   {  stepL <- 2
      startDepth <- 1
   } else
   {  stepL <- depthMaxL/10
      startDepth <- 0
   }
   legText <- seq(startDepth, depthMaxL, stepL)
   if(stepL>2)
     legText[1] <- 1  # replace 0

   return(legText)
}

plot_arcs <- function(fdTableSel, fixedMaxDepth, depthMax0, alphaOffset, selColor, circleLwd, maxY, minX, maxX, plotTitle, plotOrder, plotOption, workingDirOUT, outfilePrefix, fixedPlotWidth, labelPV, vLengthsP, vLengthsV)
{  oldpar <- par()

   # assign color to depth
   if(fixedMaxDepth>0)
   {  depthMax <- fixedMaxDepth
   } else
     depthMax <- depthMax0
   
   vDepthProp0 <- alphaOffset + fdTableSel$depth/depthMax  # alpha offset add to avoid too low alpha (leading to invisible lines)
   vDepthProp <- ifelse(vDepthProp0>1, 1, vDepthProp0)  # max. is 1
   fdTableSel[,color:=add_alpha(selColor, ifelse(vDepthProp>1, 255, 255*vDepthProp))]  # all depths larger than depthMax will be plotted with most opaque color
   print(fdTableSel)

   # create plot (pdf file)
   if(plotOption=="pdf")
   {  setwd(workingDirOUT)
      if(fixedPlotWidth>0)
      {  plotWidth <- fixedPlotWidth
      } else
      {  plotWidth <- plotHeight*xLength*1.2/maxY
      }
      outfile <- paste0(outfilePrefix,".pdf")
      if(file.exists(outfile))  # existing pdf files with the same name will be removed (if accessible)
      {  file.remove(outfile)
         print(paste(outfile,"removed"))
      }
      pdf(file=outfile, bg="transparent", height=plotHeight, width=plotWidth) 
   }
   
   yLength <- maxY + abs(minY)
   par(yaxs="i", bty="n", tcl=-0.7, mar=c(0,0,8,0), mai=c(0,0,0,0), lend="butt", pin=c(plotHeight, plotHeight*yLength/xLength))  # -> y-axis starts exactly at 0
   par(plt=c(0.04,0.96,0.1,0.96))
   print(par("plt"))
   
   vRadius <- (fdTableSel$length-1)/2 # set radii of arcs
   
   if(plotOrder=="refToFront")
   {  symbols(fdTableSel$start+vRadius, y=rep(0, nrow(fdTableSel)), circles=vRadius, inches=FALSE, xlim=c(minX, maxX), ylim=c(minY,maxY), xaxt="n", yaxt="n", main=plotTitle, font.main=1, col.main=selColor, xlab="", ylab="", fg=fdTableSel$color, lwd=circleLwd, cex.main=4.6)
   } else
     symbols(0, 0, circles=0.001, inches=FALSE, xlim=c(minX, maxX), ylim=c(minY,maxY), xaxt="n", yaxt="n", main=plotTitle, font.main=1, col.main=selColor, xlab="", ylab="", fg="transparent", lwd=circleLwd, cex.main=4.6) # "empty" plot with axes only

   # reference (vertical lines)
   if(withRef)
   {  if(optionR=="sumStrands")
      {  lines(x=dppTableRSelNorm$position, y=dppTableRSelNorm$depth, type="h", lwd=lineWidthR, col=selColorRT, asp=1)
      } else
      {  lines(x=dppTableRSelNorm$position, y=dppTableRSelNorm$norm_depth_plus, type="h", lwd=lineWidthR, col=selColorRT, asp=1)
         lines(x=dppTableRSelNorm$position, y=-dppTableRSelNorm$norm_depth_minus, type="h", lwd=lineWidthR, col=selColorRT, asp=1)
      }
      legend(x="topleft", inset=c(0.04,0.02), sampleNameR, cex=4.3*1.5, adj=c(0,0.2), pt.cex=9, pch="-", col=selColorR, text.col=selColorR, bty="n")  # legend at top left
   }
 
   if(plotOrder=="refToBack")
   {  par(new=TRUE)  # plot on top of existing plot
      symbols(fdTableSel$start+vRadius, y=rep(0, nrow(fdTableSel)), circles=vRadius, inches=FALSE, xlim=c(minX, maxX), ylim=c(minY,maxY), xaxt="n", yaxt="n", xlab="", ylab="", fg=fdTableSel$color, bg="transparent", lwd=circleLwd)
   }
   
   # white rectangle for covering lower half of circle
   rect(minX-xLength*0.05, -maxY, maxX+xLength*0.05, 0, col="white", border="white")
   
   # x axis
   yCross <- -(maxY/100)
   factorX <- 10^floor(log10(maxX-minX))
   axis(side=1, pos=yCross, at=seq(minX, maxX, factorX/10*2), lwd=3, labels=FALSE, tcl=-3) # thick ticks
   axis(side=1, pos=yCross-0.04*yLength, at=seq(minX, maxX, factorX*2), cex.axis=3.5*2, lwd=0) # labels of thick ticks
   
   # legend with depth color scale
   depthMaxL <- round_number(depthMax, "maxR", depthMax)  # for legend (see below) and outfile name
   legText <- create_legend_text(depthMaxL)
   vAlphaL0 <- alphaOffset + legText/depthMaxL  # add offset to make lightest color visible
   vAlphaL <- ifelse(vAlphaL0>1, 1, vAlphaL0)   # restrict max. alpha to 1
   if(depthMax0>depthMaxL)  # in case max. depth is larger than fixedMaxDepth or last (rounded) number in legend
     legText[length(legText)] <- paste0(depthMaxL," - ",depthMax0)   # range indicated
   legCols <- add_alpha(selColor, 255*vAlphaL)
   legend(x="topright", inset=c(0,0), legend=legText, cex=4.3*1.5, pt.cex=6*1.5, pch=15, col=legCols, bty="n", title="Depth")  # legend at top right
   
   # save pdf file
   dev.off()
   print(paste(outfile,"saved in",workingDirOUT))

   par(oldpar)  # reset par settings
}



## INPUT (user-defined parameters)
workingDirINC <- "Z:/forMiMB/test/input_files"
inputfileC <- "MiMB_dDSB_create_arc_plot_config.R"

source(paste0(workingDirINC,"/",inputfileC))  



## DEFAULTS
if(!file.exists(workingDirOUT))
{  dir.create(workingDirOUT, recursive=TRUE)
   print(paste(workingDirOUT,"created"))
}

fileSuffix <- paste0("lwd",circleLwd)
if(fixedMaxDepth>0)
  fileSuffix <- paste0(fileSuffix,"-mD",fixedMaxDepth)
if(withRef)
{  fileSuffix <- paste0(fileSuffix,"-",gsub(" ","",sampleNameRF))
   if(optionR=="sepStrands")
     fileSuffix <- paste0(fileSuffix,"Sep")
}
if(withCal)
  fileSuffix <- paste0(fileSuffix,"-cal")
if(plotOrder!="")
  fileSuffix <- paste0(fileSuffix,"-",plotOrder)



## FUNCTION
# get input file list
fileListF <- list.files(path=workingDirINF, pattern=filePatternF)
print(fileListF)

# get sample ids
vIds <- unlist(lapply(strsplit(fileListF, delimF), "[[", idIndexF))
print(vIds)

# get sample colors
if(withColorTable)
{  setwd(workingDirINS)
   sampleTable <- fread(inputfileS, header=TRUE, sep="\t", stringsAsFactors=FALSE)
   setnames(sampleTable, new=gsub(" ", ".", colnames(sampleTable)))
   print(sampleTable)
   
   vColors <- sampleTable[match(vIds,base::get(idColumnS)), base::get(colorColumnS)]
}
if(length(vColors)>0)
{  print(vColors)
} else
  vColors <- rep("black", length(vIds))

# get calibration factors (optional)
if(withCal)
{  setwd(workingDirINCF)
   calTable <- fread(inputfileCF, header=TRUE, sep="\t", stringsAsFactors=FALSE)
   setnames(calTable, gsub(" ", ".", names(calTable)))
   print(calTable) 
   
   vCals <- calTable[match(vIds,base::get(idColumnCF)), base::get(factorColumnCF)]
   if(length(vCals)>0)
   {  print(vCals)
   } else
     stop(print("No calibration factors found. Please check idIndexF, delimF, and the calibration table settings in your configuration file."))
}

# get reference file, displayed as bars (optional)
if(withRef)
{  # get reference files
   setwd(workingDirINR)
   dppTableR0 <- fread(inputfileR, header=TRUE, sep="\t", stringsAsFactors=FALSE)
   if(optionR=="sumStrands")  
   {  dppTableR <- process_dpp(dppTableR0, option="sum")  # chrom - position - depth
   } else
   {  dppTableR <- dppTableR0
   }
   print(head(dppTableR))
   print(paste(nrow(dppTableR),"rows found."))
}

# go through input fd files
for(f in seq_along(fileListF))
{  inputfileF <- fileListF[f]
   fileNameBase <- gsub("_fd.txt", "", inputfileF)
   print(fileNameBase) # for output file name
   
   # sample color
   selColor <- vColors[f]
   print(selColor)
   
   # sample calibration factor
   if(withCal)
   {  calFactor <- vCals[f]
      print(paste("Calibration factor:", calFactor))
   }
   
   # get fd table (chrom - start - length - depth)
   setwd(workingDirINF)
   fdTable0 <- fread(inputfileF, header=TRUE, sep="\t", stringsAsFactors=FALSE)
   fdTable <- process_fd(fdTable0, withCal=withCal, calFactor=calFactor)
   print(fdTable)
   
   # go through chromosomal regions
   for(r in seq_along(selRegionList))
   {  # get region
      selRegion <- selRegionList[[r]]
      currChrom <- selRegion[1]
      currStart <- as.numeric(selRegion[2])
      currEnd <- as.numeric(selRegion[3])
      print(paste0("Selected region: chrom ", currChrom, ", ", currStart, " - ", currEnd))
      
      # get fd subtable
      fdTableSel <- fdTable[chrom==currChrom & end>=currStart & start<=currEnd,]
      fdTableSel <- fdTableSel[order(fdTableSel$depth),]  # sort by depth (ascending)
      fragCount <- nrow(fdTableSel)
      
      if(fragCount>0)
      {  print(paste(fragCount, "fragments found."))
      
         # max. depth
         depthMax0 <- ceiling(max(fdTableSel$depth))
         
         # get reference dpp subtable
         dppRCount <- 0
         if(withRef)
         {  dppTableRSel <- dppTableR[chrom==currChrom & position>=currStart & position<=currEnd,]
            dppRCount <- nrow(dppTableRSel)
            print(paste(dppRCount,"rows from",sampleNameRF,"found."))
         }
            
         # prepare plot
         # get opacity setting
         alphaOffset <- vAlphaOffset[r]
         print(paste("Alpha offset:",alphaOffset))
         
         # set plot size
         xLength <- currEnd - currStart
         minX <- round_number(xLength, "min", currStart, corrFactor=1)
         maxX <- round_number(xLength, "maxC", currEnd, corrFactor=1)
         maxY <- max(ceiling((fdTableSel$length-1)/2))
         minY <- -maxY * 0.5
         
         # scale reference signals
         if(dppRCount>0)
         {  if(optionR=="sepStrands")
            {  minY <- -maxY * dppScaleFactor
               maxDR <- max(dppTableRSel$norm_depth_plus, dppTableRSel$norm_depth_minus)
               scaleFactorR <- dppScaleFactorR * maxY / maxDR 
               dppTableRSelNorm <- data.table::copy(dppTableRSel)
               dppTableRSelNorm[,`:=`(norm_depth_plus=norm_depth_plus*scaleFactorR, norm_depth_minus=norm_depth_minus*scaleFactorR)]
            } else 
            {  maxDR <- max(dppTableRSel$depth)
               scaleFactorR <- dppScaleFactorR * maxY / maxDR
               dppTableRSelNorm <- data.table::copy(dppTableRSel)
               dppTableRSelNorm[,depth:=depth*scaleFactorR]
            } 
         }
         
         # set output filename and plot title
         if(selColor!="black")
         {  selColorString <- paste0("-",gsub("#","",selColor))
         } else
           selColorString <- ""
         outfilePrefix <- paste0(fileNameBase,"_chr",currChrom,"-",currStart,"-",currEnd,"_",fileSuffix,"-a",alphaOffset,selColorString)
         plotTitle <- paste0(fileNameBase,"\nchrom. ",unique(fdTableSel$chrom)," (",fragCount," fragments, max. depth=",depthMax0,")")
     
         # create arc plot
         plot_arcs(fdTableSel, fixedMaxDepth, depthMax0, alphaOffset, selColor, circleLwd, maxY, minX, maxX, plotTitle, plotOrder, plotOption, workingDirOUT, outfilePrefix, fixedPlotWidth)
     
     } else
       stop("No signals found in this region !")
   }
}





