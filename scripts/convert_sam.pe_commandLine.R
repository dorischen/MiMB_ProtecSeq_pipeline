## Conversion of sam alignments (of NextGenMap or HISAT2) to dpp, fd or filledDpp tables
## INPUT: sam alignment file 
## OUTPUT: dpp, fd or filledDpp tables
## Author: Doris Chen 
## version 240727

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
convert_to_numeric <- function(inputTable, vColumns)
{  resultTable <- data.table::copy(inputTable)
   for(currCol in vColumns)
   {  resultTable[,match(currCol,colnames(resultTable))] <- as.numeric(unlist(resultTable[,match(currCol,colnames(resultTable)), with=FALSE]))
   }
   return(resultTable)
} 

extract_deletion_length <- function(cigarString)
{  # split CIGAR string into numbers and letters
   vCigar <- unlist(strsplit(cigarString, "(?=[MIDS])(?<=[0-9])|(?=[0-9])(?<=[MIDS])", perl=TRUE))
   vNumbers <- as.numeric(vCigar[which(vCigar=="D")-1])
   return(sum(vNumbers))
}

extract_insertion_length <- function(cigarString)
{  # split CIGAR string into numbers and letters
   vCigar <- unlist(strsplit(cigarString, "(?=[MIDS])(?<=[0-9])|(?=[0-9])(?<=[MIDS])", perl=TRUE))
   vNumbers <- as.numeric(vCigar[which(vCigar=="I")-1])
   return(sum(vNumbers))
}

extract_qlength_fast <- function(cigarString) 
{  sumLengths <- 0
   tempNum <- ""
   # Loop through each character in the CIGAR string
   for (i in seq(1,nchar(cigarString))) 
   {  char <- substr(cigarString, i, i)
      if(grepl("[0-9]", char))  # if number, appended to temporary variable
      {  tempNum <- paste0(tempNum, char)
      } else
      if(grepl("[SD]", char))  # ignore soft clipping and deletion lengths
      {  tempNum <- ""
      } else
      if(grepl("[MI]", char)) # add (mis)matches and insertions
      {  # convert the number and add to sum
         sumLengths <- sumLengths + as.numeric(tempNum)
         tempNum <- "" # Reset the temporary variable for numbers
      } 
   }
   return(sumLengths)
}

replace_na_in_table <- function(inputTable, colNames, replaceString, numeric=FALSE)
{  for(colName in colNames)
   {  if(numeric)
      {  inputTable[is.na(base::get(colName)), (colName):=as.numeric(replaceString)]
      } else
      {  inputTable[is.na(base::get(colName)), (colName):=replaceString]
      }
   }
   return(inputTable)
}

save_table <- function(inputTable, workingDirOUT, outfile)
{  setwd(workingDirOUT)
   write.table(inputTable, file=outfile, append=FALSE, sep="\t", dec=".", row.names=FALSE, col.names=TRUE, quote=FALSE)
   print(paste(outfile,"saved in ",workingDirOUT))
}

#convert_sam_to_table_pe(workingDirINS, inputfileS, inputFormat, seqFormat, outputOptions, minScore)
convert_sam_to_table_pe <- function(workingDirIN, inputfile, inputFormat, seqFormat, outputOptions, minScore=0)  # add option for HISAT!!
{  print("Loading of SAM alignment file...")  
   setwd(workingDirIN)
   samTable0 <- data.table(read.csv(inputfile, header=FALSE, sep="\t", stringsAsFactors=FALSE, comment.char="@", col.names=samHeaderOri))  # number of columns depending on alignment tool
   samTable1 <- samTable0[,names(samTable0)[dropVector]:=NULL] # remove unnecessary columns
    
   if(inputFormat=="ngm")
   {  setnames(samTable1, samHeader)
   } else
   if(inputFormat=="hisat")
   {  subTableZS <- samTable1[grep("ZS",V13),]  # extract rows with additional ZS column
      setnames(subTableZS, samHeader0)
      subTableZS[,c("ZS"):=NULL]
      subTable <- samTable1[!grep("ZS",V13),.SD,.SDcols=-(ncol(samTable1))]  # extract rest
      setnames(subTable, samHeader0[-11])
      samTable1 <- rbindlist(list(subTable, subTableZS), use.names=TRUE)  # remove ZS column and merge subtables
   }
   
   # remove tags
   for(tagName in tagTable$name)
   {  samTable1[,eval(tagName):=gsub(tagTable[name==tagName,tag], "", base::get(tagName))]
   }
   samTable1 <- convert_to_numeric(samTable1, vNumeric)
   
   # set strand according flag of first read
   samTable1[,`:=`(strand=ifelse(FLAG %in% vFlagsPlus,"+","-"))]
   
   # filter by score
   samTable2 <- samTable1[SCORE>=minScore,]
   print(paste0("Removal of ",nrow(samTable1[SCORE<minScore,])," rows with score < ",minScore,"."))
   
   # get read lengths
   samTable2[,`:=`(index=seq(1,nrow(samTable2)), del=0, ins=0)]
   if(inputFormat=="hisat")
   {  samTable2[,QLENGTH:=extract_qlength_fast(CIGAR), by=index]
      samTable2 <- samTable2[,`:=`(IDENT=(QLENGTH-EDITDIST)/QLENGTH)]  # select columns
   }
   samTable2[,`:=`(del=0, ins=0)]
   samTable2[grep("D|I",CIGAR),`:=`(del=extract_deletion_length(CIGAR), ins=extract_insertion_length(CIGAR)), by=index]  # for reference alignment lengths 
   samTable2[,`:=`(rlength=(QLENGTH+del-ins))]  # read length on reference; QLENGTH and TLEN report sequence read and insert lengths, respectively
   setkey(samTable2, index)
   
   # read pair information added
   samTable <- samTable2[,.(index=index, read=QNAME, chrom=match(RNAME, chromNames), strand=strand, rstart=POS-1, rend=POS-1+rlength, perc_identity=IDENT, score=SCORE, hit_count=(HITS+1), seq_changes=MD)] # 0-based
   samTableW <- samTable[strand=="+",] 
   samTableC <- samTable[strand=="-",] 
   setkey(samTableW, read)  # for merging
   setkey(samTableC, read)
   
   # remove "outward" pairs
   tol <- -1   # tolerance for extruding reads
   pairTable0 <- samTableW[samTableC]  
   setnames(pairTable0, old=c("rstart", "rend", "i.rstart", "i.rend"), new=c("rstart.1", "rend.1", "rstart.2", "rend.2"))
   pairTable <- pairTable0[((rend.2-rend.1) >= tol) & ((rstart.2-rstart.2) >= tol)]  # accept only read pairs with read2 overlapping or downstream of read1, with tolerance of tol nt
   print(paste0((nrow(pairTable0)-nrow(pairTable))," alignments discarded due to wrong location of pairs."))
   if(saveTableOutward)
     save_table(pairTable0[((rend.2-rend.1) < tol) & ((rstart.2-rstart.1) < tol)], workingDirOUT, outfile=paste0(fileNameBase,"_outward.txt"))
   
   # calculate lengths
   length_table <- pairTable[,.(read, start=rstart.1, end=rend.2, ref_length=rend.2-rstart.1)]  # paired lengths (insert sizes)
   setkey(length_table, read)
   
   # concatenate strand tables, add reference length information
   samTableWC <- rbindlist(list(samTableW, samTableC))
   setkey(samTableWC, read)
   parsedTable0 <- samTableWC[length_table]
   print(paste(length(unique(parsedTable0$read)),"read pairs found at",nrow(unique(parsedTable0[strand=="+",.(chrom, start, end)])),"genomic sites."))  # parsed table contains 2 rows per read pair with pair (insert) start and ends
   parsedTable0[,count:=1/hit_count]   # actually only relevant if alignments with multi-mapping read pairs
   
   # for 5prime ends only start of reads relevant
   parsedTableF <- data.table()
   parsedTable5p <- data.table()
   if(("filled" %in% outputOptions | "fd" %in% outputOptions))
   {  parsedTableF <- parsedTable0[order(chrom, start, ref_length, read)]
   } 
   if("5prime" %in% outputOptions)
   {  parsedTable05p <- parsedTable0[,.(index, read, chrom, strand, start=rstart, end=rend, perc_identity, score, count, seq_changes, ref_length)]
      parsedTable5p <- parsedTable05p[order(chrom, start, ref_length, read)]
   }
  
   return(list(parsedTableF, parsedTable5p))
}

convert_to_dpp <- function(expandedTable, outputOption, createFile=FALSE, workingDirOUT, outfileD) 
{   dppTable0 <- expandedTable[,.(depth=sum(count)), by=.(chrom, strand, position)]

    if(outputOption=="filled")
    {  dppTable <- dppTable0[,.(chrom, position, depth)]
    } else
    {  dppTableW <- dppTable0[strand=="+",.(chrom, position, norm_depth_plus=depth)]
       dppTableC <- dppTable0[strand=="-",.(chrom, position, norm_depth_minus=depth)]
       dppTable <- merge.data.table(dppTableW, dppTableC, by=c("chrom", "position"), all=TRUE)
       dppTable <- replace_na_in_table(dppTable, c("norm_depth_plus", "norm_depth_minus"), "0", numeric=TRUE)
    }
    setkey(dppTable, chrom, position)
    print(dppTable)
      
    if(createFile)
    {  save_table(dppTable, workingDirOUT, outfile=outfileD)  # create file
    } else
       write.table(dppTable, file=outfileD, append=TRUE, sep="\t", dec=".", row.names=FALSE, col.names=FALSE, quote=FALSE)  # append to file
}

convert_and_save_dpp <- function(parsedTable, maxRowCountPC, outputOption, workingDirOUT, fileNameBase, outputSuffix)
{  createFile <- TRUE  # new table (and file) created 
   outfile <- paste0(fileNameBase,"_",outputSuffix)
 
   for(chr in unique(parsedTable$chrom))  # to avoid memory problems
   {  print(chr)
      parsedTableChr <- parsedTable[chrom==chr,]
      setkey(parsedTableChr, start, end)
      rowCountPC <- nrow(parsedTableChr)
           
      if(rowCountPC>maxRowCountPC)
      {  parsedTableChr[,rowIndex:=.I]  # get gaps which split overlapping read pairs (should be processed together)
         setkey(parsedTableChr, rowIndex)
         parsedTableChr[,bin:=maxRowCountPC*ceiling(rowIndex/maxRowCountPC)]  # bins by row number
         
         # get bin borders (ends)
         binTable <- parsedTableChr[,tail(end,1), by=bin]
         setnames(binTable, c("bin","last"))  
         
         for(lastEnd in binTable$last)  # go through each bin
         {  #print(lastEnd)
            parsedTableChrSub <- parsedTableChr[start<=lastEnd,]
             
            # remove rows from table
            parsedTableChr <- parsedTableChr[rowIndex>max(parsedTableChrSub$rowIndex),]  # remaining rows
            
            # expand positions, convert to dpp 
            if(outputOption=="5prime")
            {  expTableChrSub <- parsedTableChrSub[,.(position=ifelse(strand=="+", start, end)), by=.(index, chrom, strand, count)]
            } else
              expTableChrSub <- parsedTableChrSub[,.(position=seq(start,end)), by=.(index, chrom, strand, count)]
            
            # save in file
            convert_to_dpp(expTableChrSub, outputOption, createFile, workingDirOUT, outfile)  # conversion to dpp and saving of data to file
            createFile <- FALSE
         }
      } else
      {  # expand positions, convert to dpp   
         if(outputOption=="5prime")
         {  expTableChr <- parsedTableChr[,.(position=ifelse(strand=="+", start, end)), by=.(index, chrom, strand, count)]
         } else
           expTableChr <- parsedTableChr[,.(position=seq(start,end)), by=.(index, chrom, strand, count)]
         
         # save in file
         convert_to_dpp(expTableChr, outputOption, createFile, workingDirOUT, outfile) # conversion to dpp and saving of data to file
         createFile <- FALSE
      }
   }
   
   print(paste0(outfile," was saved to ", workingDirOUT))
  
}

create_and_save_pdf <- function(plotTable, plotType, workingDirOUT, fileNameBase, plotWidth, plotHeight)
{  oldpar <- par()  # save old plot settings
   
   setwd(workingDirOUT)
   outfile <- paste0(fileNameBase,".pdf")
   if(file.exists(outfile))
   {  file.remove(outfile)
      print(paste(outfile,"removed"))
   }
   pdf(file=outfile, bg="transparent", width=plotWidth, height=plotHeight) # transparent background
   par(mar=c(5,6,2,2))
 
   plot(plotTable, type=plotType, lwd=2, cex.axis=2, cex.lab=2.5)
 
   dev.off()
   print(paste(outfile,"saved in",workingDirOUT))
   
   par(oldpar)  # reset par settings
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
argCount <- 4
args <- commandArgs(trailingOnly=TRUE)
print("Parameters:")
print(args)
if(length(args)<argCount)
{  print("Rscript convert_sam.pe_commandLine.R CONFIG_FOLDER CONFIG_FILE SAM_FOLDER SAM_FILE")
}
workingDirINC <- file.path(args[1])
inputfileC <- args[2]
workingDirINS <- file.path(args[3])
inputfileS <- args[4]

## INTERACTIVE INPUT
# workingDirINC <- "Z:/forMiMB/input_files"
# inputfileC <- "MiMB_dDSB_convert_sam.pe_ASMv1_config.R"
# workingDirINS <- "Z:/forMiMB/mapping_results/scer_specific/sub2/"
# inputfileS <- "SRR14093079_ngm-globalR1i0.96k11_ZP591.22_unaligned_ASMv1.sam"  # sam file (..1-based!!; BAM files are 0-based)

source(paste0(workingDirINC,"/",inputfileC))  
source(genomeDataPath)


## DEFAULTS
if(!file.exists(workingDirOUT))
{  dir.create(workingDirOUT, recursive=TRUE)
   print(paste(workingDirOUT,"created"))
}

scriptSuffix <- "240727"
fileSuffix <- paste0("_conv",scriptSuffix)
if(minScore>0)
  fileSuffix <- paste0(fileSuffix,"-minS",minScore)

# alignment tool specific parameters
if(inputFormat=="ngm")
{  samHeader <- c("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SCORE", "NM", "HITS", "IDENT", "QLENGTH", "MD")
   samHeaderOri <- paste0("V",seq(1,19))
   dropVector <- c(10,11,16,17) # columns in sam file which are ignored: 10..SEQ; 11..* or quality; 16..X0, nr. of equal scoring hits, same as NH; 17..XE, nr. of seeds
   tagTable <- data.table(name=c("SCORE", "NM", "HITS", "IDENT", "QLENGTH", "MD"), tag=c("AS:i:", "NM:i:", "NH:i:", "XI:f:", "XR:i:", "MD:Z:"))  
   vNumeric <- c("SCORE", "NM", "HITS", "IDENT", "QLENGTH")  
} else
if(inputFormat=="hisat")  
{  samHeaderOri <- paste0("V",seq(1,22))
   samHeader0 <- c("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SCORE", "ZS", "XN", "NM", "GAPOPEN", "GAPEXT", "EDITDIST", "MD", "YS", "PAIRCLASS", "HITS")
   dropVector <- c(10,11) # # columns in sam file which are ignored: 10..SEQ; 11..SEQ2 (rest removed after dealing with optional ZS column)
   tagTable <- data.table(name=c("SCORE", "EDITDIST", "HITS", "PAIRCLASS", "MD", "NM"), tag=c("AS:i:", "NM:i:", "NH:i:","YT:Z:", "MD:Z:", "XM:i:"))
   vNumeric <- c("SCORE", "NM", "EDITDIST", "HITS")  
}

# sam flags for forward and reverse paired-end alignments
vFlagsPlus <- c(99,163,355,419)
vFlagsMinus <- c(83,147,339,403)

# set default for max. length in length distribution plot
if(is.na(maxLengthDist) | maxLengthDist==0)
  maxLengthDist <- 1000


## FUNCTION
# create log file
fileNameBase <- paste0(gsub(".sam", "", inputfileS),fileSuffix)
print(paste("File name base:",fileNameBase))

setwd(workingDirOUT)
logfile <- paste0(fileNameBase,".log")
sink(file=logfile, append=FALSE, type="output", split=TRUE)
   
# load sam file and convert to table with one read pair per row
parsedTables <- convert_sam_to_table_pe(workingDirINS, inputfileS, inputFormat, seqFormat, outputOptions, minScore) # 0-based table
parsedTableF <- parsedTables[[1]]
parsedTable5p <- parsedTables[[2]]

# save table
if(saveTableParsed)
{  if(("filled" %in% outputOptions | "fd" %in% outputOptions))
     save_table(parsedTableF, workingDirOUT, outfile=paste0(fileNameBase,"_table.txt"))
   if("5prime" %in% outputOptions)
     save_table(parsedTable5p, workingDirOUT, outfile=paste0(fileNameBase,"_5prime_table.txt"))
}

# get insert sizes
if(nrow(parsedTableF)>0)
{  parsedTable <- data.table::copy(parsedTableF)
} else
  parsedTable <- data.table::copy(parsedTable5p)
lengthCountTable <- data.table(table(parsedTable[strand=="+",ref_length]))
setnames(lengthCountTable, c("Length", "Count"))
lengthCountTable[,Length:=as.numeric(Length)]
if(saveLengthTable)
  save_table(lengthCountTable, workingDirOUT, outfile=paste0(fileNameBase,"_lengthTable.txt"))
if(saveLengthPlot)
{  create_and_save_pdf(lengthCountTable[Length<=maxLengthDist,], plotType="l", workingDirOUT, paste0(fileNameBase,"_lengthDist"), plotWidth=8, plotHeight=4.9)
}

# convert to (filled) dpp table(s)
print("Conversion to dpp table...")  
     
if(("filled" %in% outputOptions) | ("fd" %in% outputOptions))
{  parsedTable <- parsedTableF[strand=="+",]  # count each read pair only once (already contains whole insert information)
   
   if("filled" %in% outputOptions)
     convert_and_save_dpp(parsedTable, maxRowCountPC, outputOption="filled", workingDirOUT, fileNameBase, "filledDpp.txt")

   if("fd" %in% outputOptions)
   {  # convert to fd table (0-based)
      fdTable0 <- parsedTable[strand=="+",.(chrom, start, length=ref_length, count)]  # count each read pair only once (already contains whole insert information)
      fdTable <- fdTable0[,.(depth=sum(count)), by=.(chrom, start, length)]
      setkey(fdTable, chrom, start, length)
      print(fdTable)
   
      save_table(fdTable, workingDirOUT, outfile=paste0(fileNameBase,"_fd.txt"))
   }
}
if("5prime" %in% outputOptions)
{  parsedTable <- parsedTable5p  # count each read pair only once (already contains whole insert information)
   convert_and_save_dpp(parsedTable, maxRowCountPC, outputOption="5prime", workingDirOUT, fileNameBase, "5prime_dpp.txt")
}



print("finished.")
print("Session info:")
sessionInfo()

sink(file=NULL)
print(paste0("Log saved in ",workingDirOUT,"/",logfile))


# move files to subfolders
move_files_to_folder(paste0(fileNameBase,".*\\.log"), workingDirOUT, "/log_files")

if("filled" %in% outputOptions)
  move_files_to_folder(paste0(fileNameBase,".*_filledDpp.txt"), workingDirOUT, "/filledDpp_files")

if("5prime" %in% outputOptions)
  move_files_to_folder(paste0(fileNameBase,".*_5prime_dpp.txt"), workingDirOUT, "/dpp_files")

if("fd" %in% outputOptions)
  move_files_to_folder(paste0(fileNameBase,".*_fd.txt"), workingDirOUT, "/fd_files")

if(saveLengthTable)
  move_files_to_folder(paste0(fileNameBase,".*_lengthTable.txt"), workingDirOUT, "/lengths/tables")

if(saveLengthPlot)
  move_files_to_folder(paste0(fileNameBase,".*_lengthDist.pdf"), workingDirOUT, "/lengths")

if(saveTableOutward)
  move_files_to_folder(paste0(fileNameBase,".*outward.txt"), workingDirOUT, "/out_tables")

if(saveTableParsed)
  move_files_to_folder(paste0(fileNameBase,".*_table.txt"), workingDirOUT, "/parsed_tables")


