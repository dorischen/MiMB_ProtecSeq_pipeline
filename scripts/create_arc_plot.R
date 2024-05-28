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
# selColorD1 <- "#F