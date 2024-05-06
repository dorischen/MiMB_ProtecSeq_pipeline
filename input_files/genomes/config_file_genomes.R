## configuration file for R scripts using genome information
## chroms .. sorted vector with 1 number for each chromosome, used as chromosome id in tables
## chromNames .. sorted vector with chromosome names as indicated in the fasta and alignment files
## chromLengths .. lengths of each chromosome in base pairs (in same order as above vectors)

if(genome=="R64")
{  chroms <- seq(1,16)   
   chromNames <- c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI", "mitochondrion")
   chromLengths <- c(230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 562643, 439888, 745751, 666816, 1078177, 924431, 784333, 1091291, 948066)  # chrom. lengths 1 - 16  
} else
if(genome=="ASMv1")
{  chroms <- seq(1,17)
   chromNames <- c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI", "mitochondrion")
   chromLengths <- c(228861, 829469, 340914, 1486921, 589812, 299318, 1080440, 542723, 449612, 753937, 690901, 1054145, 923535, 791982, 1053869, 946846, 84638)
 } else
if(genome=="CBS138")
{  chroms <- seq(1,14)
   chromNames <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14")
   chromLengths <- c(491328, 502101, 558804, 651701, 687738, 927101, 992211, 1050361, 1100349, 1195129, 1302831, 1455689, 1402899, 20063)
} else
if(genome=="ZP591.22")
{  chroms <- seq(1,16)
   chromNames <- paste0("OX3658",seq(80.1,95.1))
   chromLengths <- c(204266, 807096, 335459, 1464898, 540048, 287366, 1097743, 519556, 450862, 709397, 663429, 1029505, 916890, 784778, 1039997, 908134)
} else
  stop(paste("Unknown genome",genome))