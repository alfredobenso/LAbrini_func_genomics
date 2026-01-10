

require(stringr)
require(xlsx)
require(stringr)
require(GenomicRanges)
require(DescTools)
require(seqinr) 
require(bamsignals) # bam files to coverage


#-------------------------------------------------------------------------------
# Check the read coverage of the genomic regionscorresponding to operons 
#-------------------------------------------------------------------------------

# Load gene (operon) annotations
f = read.xlsx("Operon.description.gene.annot.xlsx", sheetIndex=1)

# Load BAM files
syngas.bam.fs = dir("./02.Bam/", pattern="^F58.*total.bam$")
fruct.bam.fs = dir("./02.Bam/", pattern="^FF1.*total.bam$")

# Load the FASTA file of the genome
fastafile<- read.fasta(file = "D:/A Re - Work/co-polito/PoliTo - Research/2023/RBP/study/analysis-main/LT_LAbrini_ncbi_dataset/ncbi_dataset/data/GCA_040166795.1/GCA_040166795.1_LT_LAbrini_genomic.fna", 
                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)

# Genome-wide genomic range
df = cbind.data.frame(start = 1, 
                      end = nchar(fastafile$CP110420.1),
                      strand = "*", 
                      seqid="NC_999999999.1")
cur.genomic.range = makeGRangesFromDataFrame(df, seqnames.field="seqid", start.field="start",
                                             end.field = "end", strand.field  = "strand", ignore.strand = TRUE)

#
#Count the operons where the read coverage along the genomic region occupied by
#the operon is consistent. To make the decision, the script computes the median 
#read coverage over sliding windows of 200 nt and computes the median read 
#coverage over the total region. We conclude that the read coverage is consistent
#if the median of the read coverages of the sliding windows does not deviate from
#the median read coverage over the total region by more than 10%. 
#If the genomic region of the operon is shorther or equal to 200 nt, then the 
#check is automatically passed. 
#
#Fructose samples
rpm = matrix(NA, nrow=3, ncol=nchar(fastafile$CP110420.1))
for(s in 1 : length(fruct.bam.fs)){
  readCovGenome = bamCoverage(bampath = paste0("./02.Bam/",fruct.bam.fs[s]), gr = cur.genomic.range , mapqual = 0,
                              paired.end = "extend", tlenFilter = NULL,
                              filteredFlag = -1, verbose = TRUE)[1]
  scaling.factor = sum(readCovGenome)/1000000
  rpm[s,] = readCovGenome/scaling.factor
}
RPM = colSums(rpm)
# Check if an operon has consistent read coverage over its region
lgcUnif = array(FALSE, length(opGRs))
for(i in 1 : length(opGRs)){
  
  if(f$GenesInOperon[i] == 1){
    lgcUnif[i] = TRUE
  }
  
  if(f$GenesInOperon[i] > 1){
    cur.rpm = RPM[OpFrom[i]:OpTo[i]]
    cur.median = median(RPM[OpFrom[i]:OpTo[i]])
    
    slidingWindow = 200
    
    if(length(cur.rpm) <= 200){
      lgcUnif[i] = TRUE
    }
    
    if(length(cur.rpm) > 200){
      
      rollmedian <- rollapply(cur.rpm, width = slidingWindow, FUN = median, align = "right", fill = NA)
      medianRollmedian = median(rollmedian[is.na(rollmedian)==FALSE])
      if( (medianRollmedian >= .90 * cur.median) & (medianRollmedian <= 1.10 * cur.median) ){
        lgcUnif[i] = TRUE
      }
    }
  }
  
}
names(lgcUnif) = paste0("Operon-", OpFrom)
lgcUnifFruct = lgcUnif

#
#Syngas samples
rpm = matrix(NA, nrow=3, ncol=nchar(fastafile$CP110420.1))
for(s in 1 : length(syngas.bam.fs)){
  readCovGenome = bamCoverage(bampath = paste0("./02.Bam/",syngas.bam.fs[s]), gr = cur.genomic.range , mapqual = 0,
                              paired.end = "extend", tlenFilter = NULL,
                              filteredFlag = -1, verbose = TRUE)[1]
  scaling.factor = sum(readCovGenome)/1000000
  rpm[s,] = readCovGenome/scaling.factor
}
RPM = colSums(rpm)
# Check if if an operon has consistent read coverage over its region for each operon
lgcUnif = array(FALSE, length(opGRs))
for(i in 1 : length(opGRs)){
  
  if(f$GenesInOperon[i] == 1){
    lgcUnif[i] = TRUE
  }
  
  if(f$GenesInOperon[i] > 1){
    cur.rpm = RPM[OpFrom[i]:OpTo[i]]
    cur.median = median(RPM[OpFrom[i]:OpTo[i]])
    
    slidingWindow = 200
    
    if(length(cur.rpm) <= 200){
      lgcUnif[i] = TRUE
    }
    
    if(length(cur.rpm) > 200){
      rollmedian <- rollapply(cur.rpm, width = slidingWindow, FUN = median, align = "right", fill = NA)
      medianRollmedian = median(rollmedian[is.na(rollmedian)==FALSE])
      if( (medianRollmedian >= .90 * cur.median) & (medianRollmedian <= 1.10 * cur.median) ){
        lgcUnif[i] = TRUE
      }
    }
  }
  
  show(i) 
}
names(lgcUnif) = paste0("Operon-", OpFrom)
lgcUnifSyngas = lgcUnif

op.unif.read.cov = union(names(lgcUnifFruct[lgcUnifFruct==TRUE]), names(lgcUnifSyngas[lgcUnifSyngas==TRUE]))
save(op.unif.read.cov, file="op.unif.read.cov.rda")


