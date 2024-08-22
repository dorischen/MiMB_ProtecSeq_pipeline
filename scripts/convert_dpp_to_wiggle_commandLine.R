## Conversion of dpp tables to wiggle files (tracks for genome browsers like IGV)
## INPUT: (5-prime or filled) dpp file, sample table with colors and sample names, (optional) normalisation factor table
## OUTPUT: wiggle files (1-based; optionally strand-separated and/or calibrated/normalized)
## Author: Doris Chen 
## version 240416

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
normalize_dpp <- function(dppTable, normFactor, digitsAC)
{  # dpp table with strand-specific depths
   if(dim(dppTable)[2]>3) 
   {  resultTable <- dppTable[,`:=`(norm_depth_plus=round(norm_depth_plus*normFactor,digitsAC), norm_depth_minus=round(norm_depth_minus*normFactor,digitsAC))]
   } else 
     resultTable <- dppTable[,depth:=round(depth*normFactor,digitsAC)]
   
   setkey(resultTable, chrom, position)
   print(resultTable)
   
   return(resultTable)
}

write_wiggle_file <- function(wigTable, id, desc, colorString, colorRGB, rank, outfilePath, suffix="")
{  # write first track lines
   trackLine <- paste0("track type=wiggle_0 name=",id,suffix," description=",desc,suffix," visibility=full ",colorString,"=",colorRGB," priority=",rank," autoScale=on gridDefault=on graphType=bar","\n")
   cat(trackLine, file=outfilePath)
   
   # go through chromosomes
   vChroms <- unique(wigTable$chrom)
   for(chr in vChroms)
   {  chrTableD <- wigTable[chrom==chr,.SD,.SDcols=c(2,3)]
      cat(paste0("variableStep chrom=",chromNames[chr],"\n"), file=outfilePath, append=TRUE)
      write.table(chrTableD, file=outfilePath, append=TRUE, quote=FALSE, sep="\t", dec=".", row.names=FALSE, col.names=FALSE)
   }
   
   print(paste(outfilePath,"created."))
}

create_wiggle <- function(dppTable, sepStrands, id, desc, colorRGB, rank, workingDirOUT, fileNameBase)
{  colorString <- "color"
   if(sepStrands)
   {  for(strand in c(1,2))
      {  if(strand==1)
         {  wigTable <- dppTable[norm_depth_plus>0,.(chrom,position,norm_depth_plus)]
         } else
         {  wigTable <- dppTable[norm_depth_minus>0,.(chrom,position,-norm_depth_minus)]
            colorString <- "altColor"
            rank <- rank + 1
         }
    
         suffix <- paste0("_",strand)
         outfilePath <- paste0(workingDirOUT,"/",fileNameBase,suffix,".wig")
         write_wiggle_file(wigTable, id, desc, colorString, colorRGB, rank, outfilePath, suffix)
      }
   } else
   {  suffix <- "" 
      outfilePath <- paste0(workingDirOUT,"/",fileNameBase,".wig")
      write_wiggle_file(dppTable, id, desc, colorString, colorRGB, rank, outfilePath, suffix)
   }
}
   
move_file <- function(from, to)
{  todir <- dirname(to)
   
   # check if directory is existing, if not then create
   if (!isTRUE(file.info(todir)$isdir)) 
     dir.create(todir, recursive=TRUE)

   file.rename(from=from,  to=to)
}

move_files_to_folder <- function(filePattern, workingDirOUT, subfolder)
{  fileList <- list.files(path=workingDirOUT, pattern=filePattern)
   setwd(workingDirOUT)
   for(file in fileList)
   {  print(file)
      move_file(from=file, to=paste0(workingDirOUT,"/",subfolder,"/",file))
   }
}


### INPUT
## COMMAND-LINE INPUT
argCount <- 5
args <- commandArgs(trailingOnly=TRUE)
print("Parameters:")
print(args)
if(length(args)<argCount)
{  print("Rscript convert_to_wiggle_commandLine.R CONFIG_FOLDER CONFIG_FILE INPUT_FOLDER INPUT_FILE(5-prime or filled dpp) SAMPLE_ID")
}
workingDirINC <- file.path(args[1])
inputfileC <- args[2]
workingDirIND <- file.path(args[3])
inputfileD <- args[4]
sampleId <- args[5]

## INTERACTIVE INPUT
# workingDirINC <- "Z:/forMiMB/input_files"
# inputfileC <- "MiMB_dDSB_convert_dpp_to_wiggle_ASMv1_dpp5_config.R"
# workingDirIND <- "Z:/forMiMB/mapping_results/scer_specific/dpp_files"
# inputfileD <- "SRR14093066_ngm-globalR1i0.96k11_ZP591.22_unaligned_ASMv1_conv240228_5prime_dpp.txt"  # sam file (..1-based!!; BAM files are 0-based)
# sampleId <- "SRR14093066"

source(paste0(workingDirINC,"/",inputfileC))  
source(genomeDataPath)


## DEFAULTS
if(!file.exists(workingDirOUT))
{  dir.create(workingDirOUT, recursive=TRUE)
   print(paste(workingDirOUT,"created"))
}

scriptSuffix <- "240416"
fileSuffix <- ""
if(withCal)
{  fileSuffix <- paste0(fileSuffix,"_cal")
} else
if(withRpm)
  fileSuffix <- paste0(fileSuffix,"_rpm")

options(scipen=9999) # avoid scientific notation

sampleTableHeader <- c("Sample.id", "Sample.name", "Color", "Rank")


## FUNCTION
# get file name base
fileNameBase <- gsub(".txt", "", inputfileD)
print(paste("File name base:",fileNameBase))

# load dpp table
setwd(workingDirIND)
dppTable0 <- fread(inputfileD, header=TRUE, sep="\t", stringsAsFactors=FALSE)
dppTable <- dppTable0[chrom %in% chroms,]  # restrict to chromosomes indicated in genome configuration file
dppTableW <- dppTable[,position:=position+1]
print(dppTableW)

# get sample information 
setwd(workingDirINS)
sampleTable <- fread(inputfileS, header=TRUE, sep="\t", stringsAsFactors=FALSE)
setnames(sampleTable, new=sampleTableHeader)  
sampleInfo <- sampleTable[Sample.id==sampleId,]
if(nrow(sampleInfo)>0)
{  print(sampleInfo)
} else
  stop(paste0("Sample '",sampleId, "' not found in first column of table ",workingDirINS,"/",inputfileS))

# get calibration or normalisation factor (optional)
if(withCal | withRpm)
{  setwd(workingDirINCF)
   calTable <- fread(inputfileCF, header=TRUE, sep="\t", stringsAsFactors=FALSE)
   setnames(calTable, gsub(" ", ".", names(calTable)))
   calTable 
   
   normFactor <- calTable[base::get(idColumnCF)==sampleId, base::get(factorColumn)]
   if(length(normFactor)>0)
   {  print(paste("Normalisation factor:",normFactor))
      dppTable <- normalize_dpp(dppTable, normFactor, digitsAC)
   } else
     stop(paste0("Sample '",sampleId, "' not found in column '",idColumnCF,"' of table ",workingDirINCF,"/",inputfileCF))
}
  
# create wiggle file
selColorRGB <- paste0(as.vector(col2rgb(sampleInfo$Color)[,1]), collapse=",")
create_wiggle(dppTableW, sepStrands, id=sampleInfo$Sample.id, desc=sampleInfo$Sample.name, colorRGB=selColorRGB, rank=sampleInfo$Rank, workingDirOUT, paste0(fileNameBase,fileSuffix))

# move wiggle files to subfolder
move_files_to_folder("\\.wig", workingDirOUT, "/wiggle_files")







