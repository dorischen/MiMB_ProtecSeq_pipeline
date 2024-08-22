## calculation of calibration factors from read (pair) counts
## INPUT: read (pair) count table
## OUTPUT: (optional) read (pair) count table with percentages and sample names; table with counts, RPM, OR and calibration factors
## v. 240724

## get library path
libLoc <- .libPaths()[1]   # [1] in case of multiple library locations
repository <- "https://cloud.r-project.org/"

## install packages if not yet installed and load
if(!require("data.table", quietly=TRUE))
  install.packages("data.table", lib=libLoc, repo=repository, dependencies=TRUE) 
library(data.table, lib.loc=libLoc)


## SUB-FUNCTIONS
save_table <- function(inputTable, workingDirOUT, outfile, withColNames=TRUE, withRowNames=FALSE)
{  setwd(workingDirOUT)
   write.table(inputTable, file=outfile, append=FALSE, sep="\t", dec=".", row.names=withRowNames, col.names=withColNames)
   print(paste(outfile,"saved in ",workingDirOUT))
}


### INPUT
## COMMAND-LINE INPUT
argCount <- 4
args <- commandArgs(trailingOnly=TRUE)
print("Parameters:")
print(args)
if(length(args)<argCount)
{  print("Rscript calculate_calibration_factors_commandLine.R CONFIG_FOLDER CONFIG_FILE READCOUNTTABLE_FOLDER OUTPUT_FOLDER")
}
workingDirINC <- file.path(args[1])  # directory with config file
inputfileC <- args[2]   # config file (with input parameters)
workingDirINR <- file.path(args[3]) # directory with read count table (Sample id, Total, Clipped, EXP_GENOME specific, CAL_GENOME specific)
workingDirOUT <- file.path(args[4])

## INTERACTIVE INPUT
# workingDirINC <- "Z:/forMiMB/input_files"
# inputfileC <- "MiMB_dDSB_calculate_calfactors_config.R"
# workingDirINR <- "Z:/forMiMB/"
# workingDirOUT <- workingDirINR

source(paste0(workingDirINC,"/",inputfileC))  


## DEFAULTS
scriptSuffix <- "240724"
fileSuffix <- paste0("_calFactors",scriptSuffix)
rpmString <- as.character(rpmCount/1000000)
if(rpmString!="1")
  fileSuffix <- paste0(fileSuffix,"-rp",rpmString,"M")

if(!file.exists(workingDirOUT))
{  dir.create(workingDirOUT)
   print(paste(workingDirOUT,"created"))
}


## SCRIPT
# get output filename base
fileNameBase0 <- unlist(strsplit(inputfileR, ".txt"))[1]  # for count table with sample names and percentages (optional)
print(fileNameBase0)
fileNameBase <- paste0(unlist(strsplit(inputfileR, "_counts"))[1],fileSuffix)   # for calibration factor table
print(fileNameBase)

# load alignment counts
setwd(workingDirINR)
countTable0 <- fread(inputfileR, header=TRUE, sep="\t", stringsAsFactors=FALSE)
colNamesR <- gsub(" ", ".", colnames(countTable0))   # replace spaces by '.' in column names
colnames(countTable0) <- colNamesR
countTable <- countTable0[,(idColumnR):=as.character(base::get(idColumnR))]

# add column "Unaligned" if not existing
if(length(grep("Unaligned",colnames(countTable)))==0)
{  countTable[,Unaligned:=Clipped-sum(base::get(paste0(genomeExp,".specific")), base::get(paste0(genomeCal,".specific")), Common), by=eval(idColumnR)]
   print("Column 'Unaligned' added.")
}

# add percentages (optional)
if(addPerc)
{  countTable[,`:=`(Clipped.percOfTotal=round(100*Clipped/Total,1), genomeExp.spec.percOfTotal=round(100*base::get(paste0(genomeExp,".specific"))/Total,1), genomeCal.spec.percOfTotal=round(100*base::get(paste0(genomeCal,".specific"))/Total,1), Common.percOfTotal=round(100*Common/Total,2), Unaligned.percOfTotal=round(100*Unaligned/Total,1))]
   setnames(countTable, old=c("genomeExp.spec.percOfTotal", "genomeCal.spec.percOfTotal"), new=c(paste0(genomeExp,".spec.percOfTotal"), paste0(genomeCal,".spec.percOfTotal")))
}

# load sample table with IP/WCE information and sample names
setwd(workingDirINS)
sampleTable0 <- fread(inputfileS, header=TRUE, sep="\t", stringsAsFactors=FALSE)
colNamesS <- gsub(" ", ".", colnames(sampleTable0))
colnames(sampleTable0) <- colNamesS
sampleTable0

sampleTable <- data.table::copy(sampleTable0[,.(base::get(idColumnS), base::get(nameColumnS), base::get(typeColumnS), base::get(parentColumnS))])
setnames(sampleTable, c("Sample.id", "Sample", "Type", "Parent.id"))
sampleTable[,`:=`(Sample.id=as.character(Sample.id), Parent.id=as.character(Parent.id))]
sampleTable

# add sample names
countTableN <- merge.data.table(countTable, sampleTable[,.(Sample.id, Sample)], by.x=idColumnR, by.y="Sample.id")  # table with percentages and sample names
print(countTableN)
if(saveTables)
  save_table(countTableN, workingDirOUT, outfile=paste0(fileNameBase0,"_ext.txt"))

# create table with chIP and WCE next to each other
calTable0 <- merge.data.table(sampleTable, countTable[,c(1,4,5),with=FALSE], by.x="Sample.id", by.y=idColumnR)  # for calibration factor calculation
calTableI <- calTable0[Type=="IP", -c("Parent.id"), with=FALSE]
calTableW <- calTable0[Type=="WCE"]
calTable <- merge.data.table(calTableI[,-c("Type"), with=FALSE], calTableW[,-c("Type", "Sample"), with=FALSE], by.x="Sample.id", by.y="Parent.id", suffixes=c(".IP",".WCE"))
setnames(calTable, old=names(calTable), new=gsub(".specific",".spec",names(calTable)))

# calculate occupancy ratio (OR)
calTable[,OR:=(base::get(paste0(genomeCal,".spec.WCE"))/1 * base::get(paste0(genomeExp,".spec.IP"))/1) / as.numeric(base::get(paste0(genomeExp,".spec.WCE"))/1 * base::get(paste0(genomeCal,".spec.IP"))/1)]  # Nc * Fx / Nx * Fc; OR = Wc IPx / Wx IPc 

# calculate RPM and calibration factor
calTable[,(paste0("RP",rpmString,"M")):=rpmCount/base::get(paste0(genomeExp,".spec.IP"))/1]
calTable[,Cal.Factor:=OR*base::get(paste0("RP",rpmString,"M"))]
print(calTable)

# save table
if(saveTables)
  save_table(calTable, workingDirOUT, outfile=paste0(fileNameBase,".txt"))

