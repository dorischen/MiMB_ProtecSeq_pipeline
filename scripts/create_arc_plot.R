# draw fragments as arcs with depth sum curves 
# version 240425

library(stringr)
library(plyr)
library(shape)
library(Hmisc)
library(data.table)


## FUNCTIONS
process_dpp_c <- function(dppTable, option="", cal=FALSE, calFactor)
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

process_fd_c <- function(fdTable, createIndex=FALSE, withCal=FALSE, calFactor)
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
 

# modified from http://www.magesblog.com/2013/04/how-to-change-alpha-value-of-colours-in.html
add_alpha <- function(vColors, alpha=255)  # alpha between 0 and 1!!
{  if(missing(vColors))
     stop("Please provide a vector of colours.")
   
   vColorsT <- as.character(apply(sapply(vColors, col2rgb), 2, function(x) rgb(x[1], x[2], x[3], maxColorValue=255, alpha=alpha)))
   return(vColorsT)
}

# transparent colors _ note: always pass alpha on the 0-255 scale
make_transparent <- function(inputColor, alpha=200)
{  rgbColor <- col2rgb(inputColor)
   alpha[alpha<0] <- 0
   alpha[alpha>255] <- 255
   
   transpColor <- apply(rgbColor, 2, function(x) {rgb(red=x[1], green=x[2], blue=x[3], alpha=alpha, maxColorValue=255)})
   return(transpColor)
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

plot_gene_region_simple <- function(geneTableW, geneTableC, posLIndex, posRIndex, symbolIndex, chrom, start, end, withTc, tcTableW, tcTableC, offsetG)
{  # get genes in selected region
   geneTableWsub <- geneTableW[which(geneTableW[,"chrom"]==chrom & geneTableW[,posRIndex]>=start & geneTableW[,posLIndex]<=end),]
   wCount <- nrow(geneTableWsub)
   geneTableCsub <- geneTableC[which(geneTableC[,"chrom"]==chrom & geneTableC[,posRIndex]>=start & geneTableC[,posLIndex]<=end),]
   cCount <- nrow(geneTableCsub)
   
   # draw genes
   if(wCount>0)
   {  rect(geneTableWsub[,posLIndex], offsetG+5, geneTableWsub[,posRIndex], offsetG+23, col="snow2")
      text(geneTableWsub[,posLIndex]+5, offsetG+15, offset=0, geneTableWsub[,symbolIndex], cex=4, adj=0)
      arrows(geneTableWsub[,posRIndex], offsetG+14, geneTableWsub[,"tc_midpoint"], offsetG+14, angle=90, code=2, length=0.06, lwd=2)
      if(withTc)
      {  tcTableWC <- tcTableW[which(tcTableW[,"chrom"]==chrom),]
         tcOffset <- 30
         lines(x=tcTableWC[,"posL"], y=tcTableWC[,"posR"]+tcOffset, type="l")
      }
   }
   
   if(cCount>0)
   {  rect(geneTableCsub[,posLIndex], offsetG-5, geneTableCsub[,posRIndex], offsetG-23, border="sienna4", col="moccasin")
      text(geneTableCsub[,posRIndex]-5, offsetG-13, offset=0, geneTableCsub[,symbolIndex], cex=4, adj=1, col="sienna4")
      arrows(geneTableCsub[,posLIndex], offsetG-14, geneTableCsub[,"tc_midpoint"], offsetG-14, angle=90, code=2, col="sienna4", length=0.06, lwd=2)
      if(withTc)
      {  tcTableCC <- tcTableC[which(tcTableC[,"chrom"]==chrom),]
         tcOffset <- 30
         lines(x=tcTableCC[,"posL"], y=-tcTableCC[,"posR"]-tcOffset, col="sienna4")
      }
   }
}

#plot_half_circles(fdTableSel, fixedMaxDepth, depthMax0, alphaOffset, selColor, maxY, minX, maxX, plotTitle=paste0(fileNameBase,"\nchrom. ",unique(fdTableSel$chrom)," (",fragCount," fragments, max. depth = ",depthMax0,")"), plotOrder, plotOption, workingDirOUT, outfilePrefix, fixedPlotWidth)
plot_half_circles <- function(fdTableSel, fixedMaxDepth, depthMax0, alphaOffset, selColor, circleLwd, maxY, minX, maxX, plotTitle, plotOrder, plotOption, workingDirOUT, outfilePrefix, fixedPlotWidth, labelPV, vLengthsP, vLengthsV)
{  # assign color to depth
   if(fixedMaxDepth>0)
   {  depthMax <- fixedMaxDepth
   } else
   if(depthMaxQ>0)
   {  depthMax <- quantile(fdTableSel$depth, depthMaxQ)
   } else
   {  depthMax <- depthMax0
   }
   overallMaxDepth <- max(depthMax, overallMaxDepth)
   vDepthProp0 <- alphaOffset + fdTableSel$depth/depthMax  # add to avoid too low alpha
   vDepthProp <- ifelse(vDepthProp0>1, 1, vDepthProp0)  # max. is 1
   fdTableSel[,color:=add_alpha(selColor, ifelse(vDepthProp>1, 255, 255*vDepthProp))]
   if(labelPV)
   {  fdTableSel[which(length %in% vLengthsP),color:=add_alpha(colorP, ifelse(vDepthProp>1, 255, 255*vDepthProp))]
   }
   print(head(fdTableSel))

   # create plot
   depthMaxL <- round_number(depthMax, "maxR", depthMax)  # for legend (see below) and outfile name
   
   if(plotOption=="pdf")
   {  setwd(workingDirOUT)
      if(fixedPlotWidth>0)
      {  plotWidth <- fixedPlotWidth
      } else
      {  plotWidth <- plotHeight*xLength*1.2/maxY
      }
      outfile <- paste0(outfilePrefix,"-mD",depthMaxL,"-pW",plotWidth,".pdf")
      if(file.exists(outfile))
      {  file.remove(outfile)
         print(paste(outfile,"removed"))
      }
      pdf(file=outfile, bg="transparent", height=plotHeight, width=plotWidth) # increasing width (compared to plot width) looks nicer but results in cutoff of y-axis
   }
   
   yLength <- maxY + abs(minY)
   par(yaxs="i", bty="n", tcl=-0.7, mar=c(0,0,8,0), mai=c(0,0,0,0), lend="butt", pin=c(plotHeight, plotHeight*yLength/xLength))  # -> y-axis starts exactly at 0
   par(plt=c(0.04,0.96,0.1,0.96))
   print(par("plt"))
   
   vRadius <- (fdTableSel$length-1)/2
   
   if(!woCircles & plotOrder=="DSBtoFront")
   {  symbols(fdTableSel$start+vRadius, y=rep(0, nrow(fdTableSel)), circles=vRadius, inches=FALSE, xlim=c(minX, maxX), ylim=c(minY,maxY), xaxt="n", yaxt="n", main=plotTitle, font.main=1, col.main=selColor, xlab="", ylab="", fg=fdTableSel$color, lwd=circleLwd, cex.main=4.6/2)
   } else
     symbols(0, 0, circles=0.001, inches=FALSE, xlim=c(minX, maxX), ylim=c(minY,maxY), xaxt="n", yaxt="n", main=plotTitle, font.main=1, col.main=selColor, xlab="", ylab="", fg="transparent", lwd=circleLwd, cex.main=4.6/2) # "empty" plot

   # dpp
   if(dpp1)
   {  if(optionD1=="sumStrands")
      {  lines(x=dppTable1SelNorm$position, y=dppTable1SelNorm$depth, type="h", lwd=lineWidthD1, col=selColorD1T, asp=1)
      } else
      {  lines(x=dppTable1SelNorm$position, y=dppTable1SelNorm$norm_depth_plus, type="h", lwd=lineWidthD1, col=selColorD1T, asp=1)
         lines(x=dppTable1SelNorm$position, y=-dppTable1SelNorm$norm_depth_minus, type="h", lwd=lineWidthD1, col=selColorD1T, asp=1)
      }
      legend(x="topleft", inset=c(0.04,0.02), sampleNameD1, cex=4.3*1.5, adj=c(0,0.2), pt.cex=9, pch="-", col=selColorD1, text.col=selColorD1, bty="n")  # legend at top left
   }
 
   if(withD2)
   {  if(optionD2=="sumStrands")
      {  lines(x=dppTable2SelNorm$position, y=dppTable2SelNorm$depth, type="h", lwd=lineWidthD2, col=selColorD2T)
      } else
      {  lines(x=dppTable2SelNorm$position, y=dppTable2SelNorm$norm_depth_plus, type="h", lwd=lineWidthD2, col=selColorD2T)
         lines(x=dppTable2SelNorm$position, y=-dppTable2SelNorm$norm_depth_minus, type="h", lwd=lineWidthD2, col=selColorD2T)
      }
      legend(x="topleft", inset=c(0.04,0.04), sampleNameD2, cex=4.3*1.5, adj=c(0,0.2), pt.cex=9*1.5, pch="-", col=selColorD2, text.col=selColorD2, bty="n")  # legend at top right
   }
   
   if(!woCircles & plotOrder=="DSBtoBack")
   {  par(new=TRUE)  # plot on top of existing plot
      symbols(fdTableSel$start+vRadius, y=rep(0, nrow(fdTableSel)), circles=vRadius, inches=FALSE, xlim=c(minX, maxX), ylim=c(minY,maxY), xaxt="n", yaxt="n", xlab="", ylab="", fg=fdTableSel$color, bg="transparent", lwd=circleLwd)
   }
   
   # white rectangle for covering lower half of circle
   rect(minX-xLength*0.05, -maxY, maxX+xLength*0.05, 0, col="white", border="white")
   
   # x axis
   #yCross <- -(maxY/990)
   yCross <- -(maxY/100)
   factorX <- 10^floor(log10(maxX-minX))
   axis(side=1, pos=yCross, at=seq(minX, maxX, factorX/10*2), lwd=3, labels=FALSE, tcl=-3) # thick ticks
   axis(side=1, pos=yCross-0.04*yLength, at=seq(minX, maxX, factorX*2), cex.axis=3.5*2, lwd=0) # labels of thick ticks
   #axis(side=1, pos=yCross, at=seq(minX, maxX, factorX/100), labels=FALSE, lwd=1) # thin ticks
   #mtext(side=1, line=-29, text="Position (nt)", cex=4)  # .. very difficult to position
   
   #if(withGenes)
   #   plot_gene_region_simple(geneTableW, geneTableC, 4, 5, 7, currChrom, currStart, currEnd, withTc=FALSE, tcTableW, tcTableC, offsetG=-50)
   # 
   
   if(!woCircles)
   {  # legend
      legText <- create_legend_text(depthMaxL)
      vAlphaL0 <- alphaOffset + legText/depthMaxL  # add offset to make lightest color visible
      vAlphaL <- ifelse(vAlphaL0>1, 1, vAlphaL0)   # restrict max. alpha to 1
      if(depthMax0>depthMaxL)  # in case fixed max. depth set and real max. depth larger
        legText[length(legText)] <- paste0(depthMaxL," - ",depthMax0)
      legCols <- add_alpha(selColor, 255*vAlphaL)
      legend(x="topright", inset=c(0,0), legend=legText, cex=4.3*1.5, pt.cex=6*1.5, pch=15, col=legCols, bty="n", title="Depth")  # legend at top right
   }
 
   if(plotOption!="")
   {  dev.off()
      print(paste(outfile,"saved in",workingDirOUT))
   }
}



## USER INPUT
#workingDirINF <- "/Volumes/vbc/CHROM/KLEIN/Elisa VBC/160323_160526_fd Spo11 IP"
#workingDirINF <- "D:/forFK/solexa/ddsb_fd_files"
#workingDirINF <- "E:/forFK/solexa/ddsb_fd_files_v3/transient"
#filePatternF <- "_fd.txt"
# workingDirINF <- "E:/forFK/solexa/ddsb_mapping_results_v3/R64_pe/scer_specific/fd_files"
# filePatternF <- "_fd.txt"
#workingDirINF <- "E:/forFK/solexa/sequences_180920-R6722_Cel/sub/fd_files"
#filePatternF <- "_fd.txt"
#workingDirINF <- "E:/forFK/solexa/sequences_180920-R6722_Cel/mapping_results/cel_specific/sorted/fd_files"
#filePatternF <- "_fd.txt"
#workingDirINF <- "Z:/forFK/ddsb_fd_v3/ASMv1/sub210628"
#filePatternF <- "_fd.txt"
#workingDirINF <- "Z:/forTransfer/R9525/mapping_results_ASMv1.Li.compl/scer_specific/fd_files"
filePatternF <- "_fd.txt"
#workingDirINF <- "Z:/forFK/ddsb_fd_v3/ASMv1/sub200116"
workingDirINF <- "Z:/forMiMB/mapping_results/scer_specific/fd_files"

# OPTIONAL: DSB dpp file (visible as bars)
withD1 <- FALSE  # TRUE or FALSE
# workingDirIND1 <- "/Volumes/vbc/CHROM/KLEIN/Elisa VBC/Keeney DSBs"  
# workingDirIND1 <- "D:/forPS/bernd/scer/keeney/mapping_results_sk1-mvo1/dpp_files"
# inputfileD1 <- "keeney_ngm-local-i96-r75_MVO1_pmlmin60-mm4_5prime_dpp.txt"
# sampleNameD1 <- "Keeney DSB"
# selColorD1 <- "red3"  # recommendable for 'black'
# dppScaleFactor <- 0.1  # scaling factor for dpp1 (e.g. Keeney DSB)
# workingDirIND1 <- "E:/forFK/keeney/Spo11_tel1d/mapping_results/clipped3/R1/dpp_files"
# inputfileD1 <- "keeneyMoh1_SRR3942949_Spo11-PrA_wt_t4_clipped3_MVO1-R1_5prime_dpp.txt"
# workingDirIND1 <- "E:/forFK/keeney/Spo11_tel1d/mapping_results/clipped3/ASMv1/n100/dpp_files/average"
# inputfileD1 <- "keeneyMoh1_SRR3942949_Spo11-PrA_wt_t4_ASMv1-n100_5prime_dpp.txt"
# sampleNameD1 <- "Mohibullah (2017), wt t4"
# sampleNameD1F <- "keeneyMoh1" # for file name
# selColorD1 <- "#E67300"  # recommendable for 'black'
# #selColorD1 <- "black"  # recommendable for 'black'
# optionD1 <- "sepStrands"  # "sepStrands" (Watson and Crick strands separate) or "sumStrands"
# lineWidthD1 <- 8 # default is 5; if plotOrder="DSBtoBack"
# #lineWidthD1 <- 6 # default is 5
# dppScaleFactor <- 0.25  # scaling factor for dpp1 (e.g. Keeney DSB); sumStrands:0.25, sepStrands:0.4-0.5
# # selColorD1T <- add_alpha(selColorD1, alpha=0.8*255)  # if plotOrder="DSBtoBack"
# selColorD1T <- add_alpha(selColorD1, alpha=0.7*255)  # if plotOrder="DSBtoFront"
# workingDirIND1 <- "E:/forFK/keeney/Spo11_tel1d/mapping_results/clipped3/R64/dpp_files"
# inputfileD1 <- "keeneyMoh1_SRR3942949_Spo11-PrA_wt_t4_clipped3_R64_5prime_dpp.txt"
# sampleNameD1 <- "Keeney-Mohibullah DSB (wt, t4)"
# sampleNameD1F <- "keeneyMoh1" # for file name
# selColorD1 <- "#E67300"  # recommendable for 'black'
# optionD1 <- "sepStrands"  # "sepStrands" (Watson and Crick strands separate) or "sumStrands"
# lineWidthD1 <- 10  # default is 5
# dppScaleFactor <- 0.3  # scaling factor for dpp1 (e.g. Keeney DSB)
# selColorD1T <- add_alpha(selColorD1, alpha=0.8*255)
# workingDirIND1 <- "E:/forFK/keeney/mapping_results_v3/spo11/scer_specific/dpp_files"
# inputfileD1 <- "Keeney_ngm050-localR0.75i0.96n100_ZP591_unaligned_MVO1_pml-mm4_5prime_dpp.txt"
# sampleNameD1 <- "Pan et al. (2011)"
# sampleNameD1F <- "keeneyPan-n100" # for file name
# selColorD1 <- "#FF7F00"  # recommendable for 'black'
# optionD1 <- "sepStrands"  # "sepStrands" (Watson and Crick strands separate) or "sumStrands"
# lineWidthD1 <- 6  # default is 5
# dppScaleFactor <- 0.25  # scaling factor for dpp1 (e.g. Keeney DSB)
# selColorD1T <- add_alpha(selColorD1, alpha=0.7*255)
# workingDirIND1 <- "E:/forFK/keeney/mapping_results"
# inputfileD1 <- "Keeney_S288C_R64-1-1_ngm3-config454-ws9_pmls35min52_5prime_sam_dpp.tab"
# sampleNameD1 <- "Keeney-Pan DSB (wt, t4)"
# sampleNameD1F <- "keeneyPan" # for file name
# selColorD1 <- "#E67300"  # recommendable for 'black'
# optionD1 <- "sepStrands"  # "sepStrands" (Watson and Crick strands separate) or "sumStrands"
# lineWidthD1 <- 10  # default is 5
# dppScaleFactor <- 0.3  # scaling factor for dpp1 (e.g. Keeney DSB)
# selColorD1T <- add_alpha(selColorD1, alpha=0.8*255)
# workingDirIND1 <- "E:/forPS/bernd/scer/mapping_results/R64-1-1/dpp_files/5prime"
# inputfileD1 <- "8399_15_S288C-R64_ngm3_pmls37min60_5prime_sam_dpp.tab"
# sampleNameD1 <- "Bernd, 8399"
# sampleNameD1F <- "bernd8399" # for file name
# selColorD1 <- "#E67300"  # recommendable for 'black'
# optionD1 <- "sepStrands"  # "sepStrands" (Watson and Crick strands separate) or "sumStrands"
# lineWidthD1 <- 10  # default is 5
# dppScaleFactor <- 0.3  # scaling factor for dpp1 (e.g. Keeney DSB)
# selColorD1T <- add_alpha(selColorD1, alpha=0.8*255)
workingDirIND1 <- "Z:/literature_data/Mohibullah/ASMv1/n100/dpp_files/average"
inputfileD1 <- "keeneyMoh1_SRR3942949_Spo11-PrA_wt_t4_ASMv1-n100_5prime_dpp.txt"
sampleNameD1 <- "Mohibullah (2017), wt t4"
sampleNameD1F <- "keeneyMoh1" # for file name
selColorD1 <- "#E67300"  # recommendable for 'black'
#selColorD1 <- "black"  # recommendable for 'black'
optionD1 <- "sumStrands"  # "sepStrands" (Watson and Crick strands separate) or "sumStrands"
lineWidthD1 <- 3 # default is 5; if plotOrder="DSBtoBack"
dppScaleFactor <- 0.25  # scaling factor for dpp1 (e.g. Keeney DSB); sumStrands:0.25, sepStrands:0.4-0.5
#selColorD1T <- add_alpha(selColorD1, alpha=0.7*255)  # if plotOrder="DSBtoFront"
selColorD1T <- add_alpha(selColorD1, alpha=1*255)  # if plotOrder="DSBtoBack"

plotOrder <- "DSBtoFront" # "DSBtoFront"|"DSBtoBack"to

# OPTIONAL: dDSB dpp folder (visible as bars)
withD2 <- FALSE
workingDirIND2 <- "D:/forFK/solexa/sequences_160526_R3828/mapping_results/MVO1/scer_specific/dpp_files"  
filePatternD2 <- "_5prime_dpp.txt"
sampleNameD2 <- "5prime Ends"  # for result file name and plot legend
optionD2 <- "sumStrands"   # "sepStrands" .. W and C separated OR "sumStrands" .. W and C summed up
lineWidthD2 <- 5  # default is 5
vColorsD2 <- c("tomato4", "tomato4", "black", "tomato4", "red4", "blue4", "blue4", "blue4", "red4", "red4", "black", "black")  
woCircles <- FALSE  # TRUE or FALSE; if TRUE, half circles will not be plotted

# selected lengths file
selectLengths <- FALSE
# workingDirINL <- "E:/forFK/solexa/ddsb_fd_files_v3/lengths/periodicity/peak_lists_max400_p_nonOv/tables"
# filePatternL <- "_table.txt"
workingDirINL <- "E:/forFK/solexa/ddsb_fd_v3/ASMv1/lengthDist/periodicity/tables"
filePatternL <- "sp036_39957_SPO11-MYC18_rec114-8A_rad50S_t4_5456_ASMv1_lengths180926-0-500-freq-cal-calB_countsAll_pv180626-bw3res1-thresh0.0001-maxPeak400-maxX400.txt"
selLenCutoff <- 0
selLenOptions <- "bw3"
selLenWindow <- 1
selLenFormat <- "sameMaxD"

labelPV <- FALSE
workingDirINP <- "E:/forFK/solexa/ddsb_fd_v3/ASMv1/lengthDist/periodicity/tables"
inputfileP <- "sp036_39957_SPO11-MYC18_rec114-8A_rad50S_t4_5456_ASMv1_lengths180926-0-500-freq-cal-calB_countsAll_pv180626-bw3res1-thresh0.0001-maxPeak400-maxX400.txt"
window <- 1
minPeak <- 20 # for peak distance calculation, affects max. periodic length
peakThresh <- 400  # max. peak length; for calculation of peak periodicity
periodThresh1 <- 8  # for calculation of periodicity (incl.)
periodThresh2 <- 15 # for calculation of periodicity (incl.)
colorP <- "darkred"
colorV <- "palegreen"

# selected colors in order of list of files (alphabetically); s.a. http://www.bxhorn.com/Downloads/RColors1_bxhorn.pdf or http://research.stowers-institute.org/efg/R/Color/Chart/ColorChart.pdf
#vColors <- c("tomato4", "tomato4", "black", "tomato4", "red4", "blue4", "blue4", "blue4", "red4", "red4", "black", "black")  
#vColors <- rep("black",12)  
#vColors <- rep("navy",20)  
#vColors <- rep("black",3)  

# color table
# colorTable <- FALSE  # if FALSE vColors and vSamples needed
# vColors <- rep("black",9)
# workingDirINS <- "E:/forFK/solexa/samples"
#inputfileS <- "samples_SP_170509.txt"
#inputfileS <- "samples_SP_170821.txt"
#inputfileS <- "samples_SP_180205.txt"
#inputfileS <- "samples_SP_210216.txt"
colorTable <- TRUE  # if FALSE vColors and vSamples needed
workingDirINS <- "Z:/forMiMB/input_files"
inputfileS <- "MiMB_dDSB_sample_table.txt"
idColumnS <- "SRA.id"
idIndexS <- 2

# calibration factors (for total counts)
withCal <- TRUE
#workingDirINC <- "E:/forFK/solexa/samples"
#inputfileC <- "samples_SP_170216_calFactors170216.txt"
#inputfileC <- "samples_SP_calFactors_180508.txt"
#workingDirINC <- "Z:/forTransfer/R9525"
#inputfileC <- "200401-R9525_ZP591_ASMv1_counts190828_calFactors200323-rp10M.txt"
# workingDirINC <- "Z:/forFK/samples"
# inputfileC <- "samples_dDSB_ASMv1_counts_210324_calFactors190813-rp10M.txt"
workingDirINC <- "Z:/forMiMB/input_files"
inputfileC <- "MiMB_dDSB_ZP591.22_ASMv1_calFactors240221-rp10M.txt"
idColumnC <- "Sample.id"
idIndexC <- 2

# fixedMaxDepth <- 0  # set to 0 if max. depth should be extracted from data (region); will be rounded
# alphaOffset <- 0.08 # between 0 and 1; so that depth=1 color is visible (the higher, the darker)
# plotHeight <- 70 # min. recommended plot height=20 (max. recommended length around 20kb); 50
# fixedMaxDepth <- 5000  # set to 0 if max. depth should be extracted from data (region); will be rounded
# alphaOffset <- 0 # recommendable for 'black'
# plotHeight <- 50 # recommended height for length approx. 800nt
# dppScaleFactor <- 0.7
#fixedMaxDepth <- 400  # set to 0 if max. depth should be extracted from data (region); will be rounded
#alphaOffset <- 0.1 # between 0 and 1; so that depth=1 color is visible (the higher, the darker)
#plotHeight <- 50 # min. recommended plot height=20 (max. recommended length around 20kb); 50
plotHeight <- 35 
fixedPlotWidth <- 76  # inch; recommendable in case comparison between multiple samples needed
#workingDirOUT <- paste0("E:/forFK/solexa/halfCirclePlots_mD",fixedMaxDepth,"_withCol")
workingDirOUT <- paste0(workingDirINF, "/halfCirclePlots_mD",fixedMaxDepth)
if(selectLengths)
  workingDirOUT <- paste0(workingDirOUT, "/selLen-",selLenOptions,"-",selLenFormat)
if(labelPV)
  workingDirOUT <- paste0(workingDirOUT, "/labPV")

plotOption <- "pdf"
depthMaxQ <- 0  # between 0 and 1; for setting max. depth for darkest color (highest opacity); set to 0 if real maximum should be taken
withGenes <- FALSE  # not working yet

# chrom, start, end; multiple regions possible
#selRegionList <- list(c("4", "826300", "827100"))
#selRegionList <- list(c("4", "967500", "968500"))
#selRegionList <- list(c("12", "432000", "453000"))
#selRegionList <- list(c("16", "665950", "666300"))
#selRegionList <- list(c("1", "174600", "175600"), c("4", "967500", "968400"), c("6", "109300", "110300"), c("8", "408800", "409800"))
#vAlphaOffset <- c(0.04, 0.06, 0, 0.02) # between 0 and 1; so that depth=1 color is visible (the higher, the darker)
#selRegionList <- list(c("3", "236500", "237500"))
#selRegionList <- list(c("4", "826400", "827100"),c("4", "967500", "968500"))
#selRegionList <- list(c("4", "826400", "827100"),c("3", "236800", "237400"))
#vAlphaOffset <- c(0.01,0.01) # for color and hotspots chr4 and chr3
#vAlphaOffset <- c(0,0) # between 0 and 1; so that depth=1 color is visible (the higher, the darker); one per region!!
#selRegionList <- list(c("3", "216300", "217700"),c("12", "793500", "794100"),c("16", "388500", "389400"))
# selRegionList <- list(c("3", "216500", "217300"),c("4", "823850", "824500"),c("12", "793500", "794100"))  # ASMv1
# vAlphaOffset <- c(0.01,0.01,0.03) # between 0 and 1; so that depth=1 color is visible (the higher, the darker); one per region!! ; ASMv1
#selRegionList <- list(c("1","15067100","15067700"), c("2","8289300","8289900"), c("3","7444400","7445000"),c("4","11073000","11073500"),c("5","5285200","5285600"))
#vAlphaOffset <- c(0.01,0.05,0.05,0.05,0.05) 
# selRegionList <- list(c("3", "214700", "215200"),c("3", "215600", "216000"),c("3", "218900", "219200"))  # ASMv1
# vAlphaOffset <- c(0.01,0.01,0.03) # between 0 and 1; so that depth=1 color is visible (the higher, the darker); one per region!! ; ASMv1
# circleLwd <- 3 # 0.1 default
# selRegionList <- list(c("STE50_FUS1_S2210" ,"3000", "7000"),c("STE50_FUS1_S4921", "3000", "7000"), c("STE50_FUS1_S2210" ,"7000", "10500"), c("STE50_FUS1_S4921", "7000", "10500"))  # ASMv1
# vAlphaOffset <- c(0.01,0.01,0.01,0.01) # between 0 and 1; so that depth=1 color is visible (the higher, the darker); one per region!! ; ASMv1
# circleLwd <- 3 # 0.1 default
# selRegionList <- list(c("4" ,"690400", "691400"))  # ASMv1
# vAlphaOffset <- c(0) # between 0 and 1; so that depth=1 color is visible (the higher, the darker); one per region!! ; ASMv1
# circleLwd <- 0.1 # 0.1 default

# fixedMaxDepth <- 400  # set to 0 if max. depth should be extracted from data (region); will be rounded
# selRegionList <- list(c("4" ,"823700", "824600"))  # ASMv1
# vAlphaOffset <- c(0) # between 0 and 1; so that depth=1 color is visible (the higher, the darker); one per region!! ; ASMv1
# circleLwd <- 0.1 # 0.1 default
#fixedMaxDepth <- 400  # set to 0 if max. depth should be extracted from data (region); will be rounded
#selRegionList <- list(c("5" ,"480500", "482000"))  # ASMv1
#vAlphaOffset <- c(0.01) # between 0 and 1; so that depth=1 color is visible (the higher, the darker); one per region!! ; ASMv1
# fixedMaxDepth <- 400  # set to 0 if max. depth should be extracted from data (region); will be rounded
# selRegionList <- list(c("4" ,"1052400", "1053000"))  # ASMv1
# vAlphaOffset <- c(0.02) # between 0 and 1; so that depth=1 color is visible (the higher, the darker); one per region!! ; ASMv1
# circleLwd <- 0.1 # 0.1 default
# fixedMaxDepth <- 400  # set to 0 if max. depth should be extracted from data (region); will be rounded
# selRegionList <- list(c("3" ,"216800", "217100"))  # ASMv1
# vAlphaOffset <- c(0.01) # between 0 and 1; so that depth=1 color is visible (the higher, the darker); one per region!! ; ASMv1
# circleLwd <- 0.1 # 0.1 default
fixedMaxDepth <- 400  # set to 0 if max. depth should be extracted from data (region); will be rounded
selRegionList <- list(c("7" ,"833500", "834000"))  # ASMv1
vAlphaOffset <- c(0.01) # between 0 and 1; so that depth=1 color is visible (the higher, the darker); one per region!! ; ASMv1
circleLwd <- 0.1 # 0.1 default



plusAnnotation <- FALSE
genome <- "ASMv1"
#genome <- "Cel235.r93"



## DEFAULTS
if(!file.exists(workingDirOUT))
{  dir.create(workingDirOUT, recursive=TRUE)
   print(paste(workingDirOUT,"created"))
}

scriptSuffix <- "240425"
fileSuffix <- paste0("-lwd",circleLwd)
if(depthMaxQ>0)
  fileSuffix <- paste0(fileSuffix,"-dmq",depthMaxQ*100)
if(withD1)
{  fileSuffix <- paste0(fileSuffix,"-",str_replace(sampleNameD1F," ",""))
   if(optionD1=="sepStrands")
     fileSuffix <- paste0(fileSuffix,"Sep")
}
if(withD2)
{  fileSuffix <- paste0(fileSuffix,"-",str_replace(sampleNameD2," ",""))
   if(optionD2=="sepStrands")
     fileSuffix <- paste0(fileSuffix,"Sep")
}
if(woCircles)
  fileSuffix <- paste0(fileSuffix,"-woCircles")
if(withGenes)
  fileSuffix <- paste0(fileSuffix,"-genes")
if(withCal)
  fileSuffix <- paste0(fileSuffix,"-cal")
if(plotOrder!="")
  fileSuffix <- paste0(fileSuffix,"-",plotOrder)

titleSuffix <- ""
if(selLenCutoff>0)
  titleSuffix <- paste0(", <=",selLenCutoff,"nt")

if(plusAnnotation)
{  fileSuffix <- paste0(fileSuffix,"-plusAnnot")
   cutOption <- "depth"  # "depth"|"rpb"
   exprCutQ <- 0.1  # quantile cutoff (either depth or rpb)
   maxGap <- 3000 # max. gap (or overlap) between convergent transcripts
   
   if(genome=="sk1-mvo1")
   {  # transcript files
      workingDirINT <- "E:/genomes/Scer/sources/sk1_revised_gbrar/MVO1/ngm_out/dpp_files"
      inputFileT0 <- "SRR387845_mRNA-traditional-I_ngm050_MVO1-local-i0_pmlmin60-mm4_dpp.txt"
      
      saveTableMidpoint <- TRUE
      workingDirOUTM <- "E:/genomes/Scer/transcripts"
      outfileM <- "sk1_gbrar_tc_meio_convTc-150827_midpoints_MVO1.txt"
      
      # annotation file (bed format)
      workingDirINA <- "E:/genomes/Scer/sources/sgd/version_20150113"
      inputFileA <- "scerR64-2-1_genes-160925_outerEnds-161004_bwa_MVO1_b.bed"   # from process_sgd_transcript_annotation_160927.R; extraction of R64.2 sequences, alignment to sk1-mvo1
      inputFileTrans <- "saccharomyces_cerevisiae_160925_outerEnds-160930-transGenes.txt"  # ids only
   
      # recombination sites (bed file)
      workingDirINRec <- "E:/genomes/Scer/sources/steinmetz/recombination_sites_2013"
      inputFileRec <- "steinmetzRecSitesR63_bwa_MVO1_sorted_b.bed"
   }

   # colors
   selColorsT <- c("black", "red3")
   selColorTrans <- "turquoise4"
   selColorMid <- "grey30"
   selColorRec <- "purple4"
   
   selColorsTT <- make_transparent(selColorsT, alpha=200)
   selColorTransT <- make_transparent(selColorTrans, alpha=200)
   selColorMidT <- make_transparent(selColorMid, alpha=150)
   selColorRecT <- make_transparent(selColorRec, alpha=100)
}




## FUNCTION
oldpar <- par()

fileListF <- list.files(path=workingDirINF, pattern=filePatternF)
print(fileListF)

if(colorTable)
{  setwd(workingDirINS)
   sampleTable <- fread(inputfileS, header=TRUE, sep="\t", stringsAsFactors=FALSE)
   colNamesS <- str_replace_all(colnames(sampleTable), " ", ".")
   colnames(sampleTable) <- colNamesS
   sampleTable 
   
   vIds <- unlist(lapply(str_split(fileListF, "_"), "[[", idIndexS))
   
   vColors <- sampleTable[match(vIds,base::get(idColumnS)), Color]
}
print(vColors)

if(withCal)
{  setwd(workingDirINC)
   calTable <- fread(inputfileC, header=TRUE, sep="\t", stringsAsFactors=FALSE)
   calTable 
   
   vSampleIds <- unlist(lapply(str_split(fileListF, "_"), "[[", idIndexC))
   vCals <- calTable[match(vSampleIds,base::get(idColumnC)), Cal.Factor]
   vCals[is.na(vCals)] <- 1  # in case of missing OR fill in with 1
} else
  vCals <- rep(1, length(fileListF))
vCals

if(withD1)
{  # get reference files
   setwd(workingDirIND1)
   dppTable1Source <- fread(inputfileD1, header=TRUE, sep="\t", stringsAsFactors=FALSE)
   if(optionD1=="sumStrands")
   {  dppTable1 <- process_dpp_c(dppTable1Source, option="sum")  # chrom - position - depth
   } else
   {  dppTable1 <- dppTable1Source
   }
   print(head(dppTable1))
   print(paste(nrow(dppTable1),"rows found."))
}

if(withD2)
{  fileListD2 <- list.files(path=workingDirIND2, pattern=filePatternD2)
   print(fileListD2)
}

# get length selection file (optional)
if(selectLengths)
{  fileListL <- list.files(path=workingDirINL, pattern=filePatternL)
   fileListL <- rep(fileListL,length(fileListF))
   print(fileListL)
   vSampleIds <- unlist(lapply(str_split(fileListF, "_"), "[[", 1))
}

if(labelPV)
{  # load peak valley table
   pvTable0 <- load_table(workingDirINP, inputfileP, plusHeader=TRUE)
   pvTable <- pvTable0[Value>=minPeak & Value<=peakThresh & (peakDiff==0 | is.na(peakDiff) | (peakDiff>=periodThresh1 & peakDiff<=periodThresh2))]
   pvTable
   
   vLengthsP0 <- pvTable[Type=="peak",Value]
   vLengthsP <- sort(c(vLengthsP0, vLengthsP0+1, vLengthsP0-1))
   
   vLengthsV0 <- pvTable[Type=="valley",Value]
   vLengthsV <- sort(c(vLengthsV0, vLengthsV0+1, vLengthsV0-1))
}

overallMaxDepth <- 0
for(f in seq_along(fileListF))
{  inputfileF <- fileListF[f]
   fileNameBase <- str_replace(inputfileF, "_fd.txt", "")
   print(fileNameBase)
   
   selColor <- vColors[f]
   print(selColor)
   calFactor <- 1
   if(withCal)
   {  calFactor <- vCals[f]
      print(paste("Calibration factor:", calFactor))
   }
   
   # get fd table (chrom - start - length - depth - end)
   setwd(workingDirINF)
   fdTable0 <- fread(inputfileF, header=TRUE, sep="\t", stringsAsFactors=FALSE)
   fdTable <- process_fd_c(fdTable0, withCal=withCal, calFactor=calFactor)
   print(head(fdTable))
   print(paste(nrow(fdTable),"rows found."))
   
   if(withD2)
   {  inputfileD2 <- fileListD2[f]
      setwd(workingDirIND2)
      dppTable2Source <- fread(inputfileD2, header=TRUE, sep="\t", stringsAsFactors=FALSE)
      if(optionD2=="sumStrands")
      {  dppTable2 <- process_dpp_c(dppTable2Source, option="sum", cal=withCal, calFactor)  # chrom - position - depth
      } else
      {  dppTable2 <- process_dpp_c(dppTable2Source, cal=withCal, calFactor)
      }
      print(head(dppTable2))
      print(paste(nrow(dppTable2),"rows found."))
      
      selColorD2 <- vColorsD2[f]
      selColorD2T <- add_alpha(selColorD2, alpha=0.75*255)
   }
   
   if(selectLengths)
   {  #inputfileL <- grep(vSampleIds[f], fileListL, value=TRUE)
      inputfileL <- fileListL[f]
      setwd(workingDirINL)
      selLenTable0 <- fread(inputfileL, header=TRUE, sep="\t", stringsAsFactors=FALSE)
      
      if(selLenCutoff>0)
        selLenTable0 <- selLenTable0[Value<=selLenCutoff,]
   }
  
   # plot regions
   for(r in seq_along(selRegionList))
   {  selRegion <- selRegionList[[r]]
      currChrom <- selRegion[1]
      currStart <- as.numeric(selRegion[2])
      currEnd <- as.numeric(selRegion[3])
      print(paste0("chrom", currChrom, ":", currStart, "-", currEnd))
      
      alphaOffset <- vAlphaOffset[r]
      
      # get fd in region
      fdTableSel <- fdTable[chrom==currChrom & end>=currStart & start<=currEnd,]
      fdTableSel <- fdTableSel[order(fdTableSel$depth),]
    
      depthMax0 <- ceiling(max(fdTableSel$depth))
      
      # get dpp in region
      dpp1Count <- 0
      dpp2Count <- 0
      dpp1 <- FALSE
      if(withD1)
      {  dppTable1Sel <- dppTable1[chrom==currChrom & position>=currStart & position<=currEnd,]
         dpp1Count <- nrow(dppTable1Sel)
         print(paste(dpp1Count,sampleNameD1F,"found"))
      }
      if(withD2)
      {  dppTable2Sel <- dppTable2[chrom==currChrom & position>=currStart & position<=currEnd,]
         dpp2Count <- nrow(dppTable2Sel)
         print(paste(dpp2Count," rows found"))
      }
         
      # prepare plot
      xLength <- currEnd - currStart
      minX <- round_number(xLength, "min", currStart, corrFactor=1)
      maxX <- round_number(xLength, "maxC", currEnd, corrFactor=1)
      maxY <- max(ceiling((fdTableSel$length-1)/2))
      #dppScaleFactor <- 0.7  # fraction of maxY that will be assigned to maxD1
      if(optionD1=="sepStrands" | optionD2=="sepStrands")
      {  minY <- -maxY * dppScaleFactor
      } else
        minY <- -maxY * 0.5
      if(dpp1Count>0)
      {  dpp1 <- TRUE
         if(optionD1=="sumStrands")
         {  maxD1 <- max(dppTable1Sel$depth)
            scaleFactorD1 <- dppScaleFactor * maxY / maxD1 
            #scaleFactorD1 <- dppScaleFactor
            dppTable1SelNorm <- data.table::copy(dppTable1Sel)
            dppTable1SelNorm$depth <- dppTable1Sel$depth * scaleFactorD1
         } else
         {  maxD1 <- max(dppTable1Sel$norm_depth_plus, dppTable1Sel$norm_depth_minus)
            scaleFactorD1 <- dppScaleFactor * maxY / maxD1 
            #scaleFactorD1 <- dppScaleFactor
            dppTable1SelNorm <- data.table::copy(dppTable1Sel)
            dppTable1SelNorm$norm_depth_plus <- dppTable1Sel$norm_depth_plus * scaleFactorD1
            dppTable1SelNorm$norm_depth_minus <- dppTable1Sel$norm_depth_minus * scaleFactorD1
         }
      }
      if(dpp2Count>0)
      {  if(optionD2=="sumStrands")
         {  maxD2 <- max(dppTable2Sel$depth)
            scaleFactorD2 <- dppScaleFactor * maxY / maxD2
            dppTable2SelNorm <- data.table::copy(dppTable2Sel)
            dppTable2SelNorm$depth <- dppTable2Sel$depth * scaleFactorD2
         } else
         {  maxD2 <- max(dppTable2Sel$norm_depth_plus, dppTable2Sel$norm_depth_minus)
            scaleFactorD2 <- dppScaleFactor * maxY / maxD2
            dppTable2SelNorm <- data.table::copy(dppTable2Sel)
            dppTable2SelNorm$norm_depth_plus <- dppTable2Sel$norm_depth_plus * scaleFactorD2
            dppTable2SelNorm$norm_depth_minus <- dppTable2Sel$norm_depth_minus * scaleFactorD2
         }
      }
      
      # subselect certain lengths
      if(selectLengths)
      {  for(selLenType in c("peak","valley"))
         {  print(selLenType)
            selLenTable <- selLenTable0[Type==selLenType,]
            vSelLengths <- unique(c(selLenTable$Value, selLenTable$Value+selLenWindow, selLenTable$Value-selLenWindow))
            fdTableSelL <- fdTableSel[length %in% vSelLengths,]
            print(fdTableSelL)
      
            if(nrow(fdTableSelL)>0)
            {  if(selLenFormat!="sameMaxD")
                  depthMax0 <- ceiling(max(fdTableSelL$depth))
               
               fragCount <- nrow(fdTableSelL)
               print(paste(fragCount, "fragments found."))
            
               fileSuffixL <- paste0(fileSuffix,"_selLen-",selLenOptions,"-w",selLenWindow,"-",selLenCutoff,"-",selLenFormat,"-",selLenType)
               outfilePrefix <- paste0(fileNameBase,"_chr",currChrom,"-",currStart,"-",currEnd,"_fd",scriptSuffix,fileSuffixL,"-aO",alphaOffset,"-",selColor)
               plot_half_circles(fdTableSelL, fixedMaxDepth, depthMax0, alphaOffset, selColor, circleLwd, maxY, minX, maxX, plotTitle=paste0(fileNameBase,"\nchrom. ",unique(fdTableSelL$chrom)," (",fragCount," fragments, max. depth = ",max(fdTableSelL$depth),"), ", selLenType,"s (+/-",selLenWindow, "nt) only",titleSuffix), plotOrder, plotOption, workingDirOUT, outfilePrefix, fixedPlotWidth, labelPV, vLengthsP, vLengthsV)
            }
         }
      } else
      {  if(nrow(fdTableSel)>0)
         {  fragCount <- nrow(fdTableSel)
            print(paste(fragCount, "fragments found."))
         
            outfilePrefix <- paste0(fileNameBase,"_chr",currChrom,"-",currStart,"-",currEnd,"_fd",scriptSuffix,fileSuffix,"-aO",alphaOffset,"-",selColor)
            plot_half_circles(fdTableSel, fixedMaxDepth, depthMax0, alphaOffset, selColor, circleLwd, maxY, minX, maxX, plotTitle=paste0(fileNameBase,"\nchrom. ",unique(fdTableSel$chrom)," (",fragCount," fragments, max. depth = ",depthMax0,")"), plotOrder, plotOption, workingDirOUT, outfilePrefix, fixedPlotWidth, labelPV, vLengthsP, vLengthsV)
         }
      }
   }  
}

par(oldpar)  # reset par settings



