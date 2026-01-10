
require(stringr)
require(xlsx)
require(seqinr) 
require(stringr)
require(colorRamp2)
require(GenomicRanges)
require(bamsignals) 
require(zoo) 


#-------------------------------------------------------------------
# Retrieve NCBI annotation for Clostridium autoethanogenum
#-------------------------------------------------------------------

ncbif = read.delim("../data/GCA_040166795.1/genomic.gff", header=FALSE, comment.char = "#", sep="\t")

# Subset the rows of the NCBI files that contains the annotations of the protein coding sequences
annotf = ncbif[ncbif[,3] %in% "CDS",]
# CDS IDs
cdsIds = as.vector(sapply(annotf[,9], function(x) strsplit(x,";")[[1]][1]))
cdsIds = gsub("ID=cds-", "", cdsIds)
# Locus tags
locusTags = as.vector(sapply(annotf[,9], function(x) strsplit(x,";")[[1]][2]))
locusTags = gsub("Parent=gene-", "", locusTags)
# CDS start
cdsStart = annotf[,4]
# CDS end
cdsSEnd = annotf[,5]
# CDS length
cdsLength = cdsSEnd-cdsStart


#------------------------------------------------------------------------
# Load operons displaying uniform read coverage along their region
#------------------------------------------------------------------------

load("../results/op.unif.read.cov.rda")


#-------------------------------------------------------------------------
# Obtain transcript leves from RNA-seq data 
#-------------------------------------------------------------------------

# Load gene (operon) annotations
f = read.xlsx("../data/Operon.description.gene.annot.xlsx", sheetIndex=1)

# Load transcript raw read counts in LAbrini fructose and syngas chemostats
ftranscript = read.delim("../data/Gene_count_transcriptome.txt")
ftranscript = ftranscript[grepl("ERCC",ftranscript$Gene)==FALSE,]

# Get RPKM values
m = ftranscript[,2:7]
rpm = matrix(NA, nrow=dim(m)[1], ncol=dim(m)[2], dimnames=list(ftranscript$Gene, colnames(m)))
scalingFactor = colSums(m)/1000000
for(i in 1 : dim(m)[2]){
  rpm[,i] = m[,i]/scalingFactor[i]
}
rpkm = matrix(NA, nrow=dim(rpm)[1], ncol=dim(rpm)[2], dimnames=list(ftranscript$Gene, colnames(rpm)))
for(i in 1 : dim(rpm)[2]){
  cur.len = cdsLength[match(rownames(rpm)[i], locusTags)]/1000
  rpkm[,i] = rpm[,i]/cur.len
}


# Compute median and median absolute deviation (MAD) of transcript levels by growth substrate
operon.fruct.median = c(); operon.fruct.mad = c()
operon.syngas.median = c(); operon.syngas.mad = c()
for(i in 1 : dim(f)[1]){
  
  genesInOperon = str_split(f$GenesInOperon_CP110420[i],",")[[1]]
  genesInOperon = gsub("LABRINI","LAbrini",genesInOperon)
  
  # Retrieve transcript levels by growth condition
  fruct.vec = c(); syngas.vec = c()
  for(j in 1 : length(genesInOperon)){
    
    tmpf = rpkm[match(genesInOperon[j], rownames(rpkm)),]
    
    cur.fruct.vec = c(as.numeric(tmpf[1]), as.numeric(tmpf[3]), as.numeric(tmpf[5]))
    fruct.vec = c(fruct.vec, cur.fruct.vec)
    
    cur.syngas.vec = c(as.numeric(tmpf[2]), as.numeric(tmpf[4]), as.numeric(tmpf[6])) 
    syngas.vec = c(syngas.vec, cur.syngas.vec)
    
  }
  
  # Compute the median of transcript levels in an operon in syngas and fructose
  operon.fruct.median[paste0("Operon-",f$OperonID[i])] = median(fruct.vec)
  operon.syngas.median[paste0("Operon-",f$OperonID[i])] = median(syngas.vec)
  
  # Compute the MAD of transcript levels in an operon in syngas and fructose
  operon.fruct.mad[paste0("Operon-",f$OperonID[i])] = median(abs(fruct.vec - median(fruct.vec)))
  operon.syngas.mad[paste0("Operon-",f$OperonID[i])] = median(abs(syngas.vec - median(syngas.vec)))
  
} 

# Discretize operon median read coverages in syngas and fructose
step.quantile.exp = .20
bin.operon.fruct.median = quantile(operon.fruct.median, probs=seq(0,1,step.quantile.exp))
bin.operon.fruct.median[1] = bin.operon.fruct.median[1] - 1 

bin.operon.syngas.median = quantile(operon.syngas.median, probs=seq(0,1,step.quantile.exp))
bin.operon.syngas.median[1] = bin.operon.syngas.median[1] - 1 

# Label quantiles on fructose
op.fruct.median.quant.end.point = levels(cut(operon.fruct.median, bin.operon.fruct.median))
op.fruct.median.quant.label = seq(length(op.fruct.median.quant.end.point),1,by=-1)

# Label quantiles on syngas
op.syngas.median.quant.end.point = levels(cut(operon.syngas.median, bin.operon.syngas.median))
op.syngas.median.quant.label = seq(length(op.syngas.median.quant.end.point),1,by=-1)

expQuByOperon = op.fruct.median.quant.label[match(cut(operon.fruct.median, bin.operon.fruct.median), 
                                                  op.fruct.median.quant.end.point)]
names(expQuByOperon) = names(operon.fruct.median)


#-----------------------------------------------------------------------
# Load processed data from Cappable-seq and Term-seq 
#-----------------------------------------------------------------------

load("../data/cap.term.data.rda")


#-----------------------------------------------------------------------
# Filter out TISs whose positions coincide with CDS starts 
#-----------------------------------------------------------------------

lgc.tss = array(TRUE, dim(tbseq.cap.term$cap_site)[1])
for(i in 1 : dim(tbseq.cap.term$cap_site)[1]){
  
  tss.strand = tbseq.cap.term$cap_site$Strand[i]
  if(is.na(match(tbseq.cap.term$cap_site$Strand[i], "+"))==FALSE){
    cur.delta = tbseq.cap.term$cap_site$OperonFrom[i] - tbseq.cap.term$cap_site$TIS_pos[i]
    if(cur.delta==0){lgc.tss[i] = FALSE}
  }
  if(is.na(match(tbseq.cap.term$cap_site$Strand[i], "-"))==FALSE){
    cur.delta = tbseq.cap.term$cap_site$TIS_pos[i] - tbseq.cap.term$cap_site$OperonTo[i]
    if(cur.delta == 0){lgc.tss[i] = FALSE}
  }
  
}
cap.site = tbseq.cap.term$cap_site[lgc.tss==TRUE,]

#-----------------------------------------------------------------------
# Assemble TISs by gene (operon)
#-----------------------------------------------------------------------

# Define gene (operon) identifiers by their "OperonFrom" positions
operonIDs = unique(cap.site$OperonFrom)

# Retrieve TISs per operon 
tssByOperon = list()
for(i in 1 : length(operonIDs)){
  tmp = cap.site[is.na(match(cap.site$OperonFrom, operonIDs[i]))==FALSE,]
  tssByOperon[[paste0("Operon-",operonIDs[i])]] = tmp$TIS_ID
}
tssByOperonCounts = as.numeric(lapply(tssByOperon, function(x) length(x)))


#--------------------------------------------------------------------------
# Site filtering by site coverage and distance from CDS start  
#--------------------------------------------------------------------------

step.quantile = .20
thrs.quantile = 5

#
# Discretize the distance between TISs and CDS start positions.
# Discretize TISs read coverage gauged on fructose and syngas. 
#

# Collect distances between TISs and the start codons of the downstream CDSs
delta = c()
# Collect the coverage of TISs per operon normalized by the average coverage of all sites from all samples
readcov = c()
# Collect the coverage of TISs per operon normalized by the average coverage of all sites from fructose samples
fruct.readcov = c()
# Collect the coverage of TISs per operon normalized by the average coverage of all sites from syngas samples
syngas.readcov = c()
# Compute the median read coverage of TISs per operon on fructose and syngas
fruct.median = c()
syngas.median = c()
# Compute the MAD of read coverage of TISs per operon on fructose and syngas
fruct.mad = c()
syngas.mad = c()
for(i in 1 : dim(cap.site)[1]){
  
  cur_delta=NA
  if(is.na(match(cap.site$Strand[i], "+"))==FALSE){
    cur.delta = abs(cap.site$OperonFrom[i] - cap.site$TIS_pos[i])
  }
  if(is.na(match(cap.site$Strand[i], "-"))==FALSE){
    cur.delta = abs(cap.site$TIS_pos[i] - cap.site$OperonTo[i])
  }
  delta = c(delta, cur.delta)
  
  #TIS coverage normalized by the average coverage of all sites from all samples
  readcov = c(readcov, cap.site$All_L2RelCov[i])
  #TIS coverage normalized by the average coverage of all sites from fructose samples
  fruct.readcov = c(fruct.readcov, cap.site$B_L2RelCov[i])
  #TIS coverage normalized by the average coverage of all sites from syngas samples
  syngas.readcov = c(syngas.readcov, cap.site$BR_L2RelCov[i])
  
  fruct.read.cov = c(cap.site$B1_L2RelCov[i], cap.site$B3_L2RelCov[i], cap.site$B4_reads[i])
  syngas.read.cov = c(cap.site$BR1_L2RelCov[i], cap.site$BR3_L2RelCov[i], cap.site$BR4_reads[i])
  
  fruct.median = c(fruct.median, median(fruct.read.cov))
  syngas.median = c(syngas.median, median(syngas.read.cov))
  
  fruct.mad = c(fruct.mad, median(abs(fruct.read.cov - median(fruct.read.cov))))
  syngas.mad = c(syngas.mad, median(abs(syngas.read.cov - median(syngas.read.cov))))
}

# Create quantiles of Dist(TIS, CDS start): Q(Dist)
bin.delta = quantile(delta, probs=seq(0,1,step.quantile))
bin.delta[1] = bin.delta[1] - 1 

# Label distance quantiles
delta.quant.end.point = levels(cut(delta,bin.delta))
delta.quant.label = 1:length(delta.quant.end.point)

# Create TIS coverage quantiles: Q(SiteCov)
bin.read.cov = quantile(readcov, probs=seq(0,1,step.quantile))
bin.read.cov[1] = bin.read.cov[1] - 1 

# Label TIS coverage quantiles
readcov.quant.end.point = levels(cut(readcov,bin.read.cov))
readcov.quant.label = seq(length(readcov.quant.end.point),1,by=-1)

# Create TIS coverage quantiles on fructose samples: Q_fructose(SiteCov)
fruct.bin.read.cov = quantile(fruct.readcov, probs=seq(0,1,step.quantile))
fruct.bin.read.cov[1] = fruct.bin.read.cov[1] - 1 

# Label TIS coverage quantiles on fructose samples 
fruct.readcov.quant.end.point = levels(cut(fruct.readcov,fruct.bin.read.cov))
fruct.readcov.quant.label = seq(length(fruct.readcov.quant.end.point),1,by=-1)

# Create TIS coverage quantiles on syngas samples: Q_syngas(SiteCov)
syngas.bin.read.cov = quantile(syngas.readcov, probs=seq(0,1,step.quantile))
syngas.bin.read.cov[1] = syngas.bin.read.cov[1] - 1 

# Label TIS coverage quantiles on syngas samples
syngas.readcov.quant.end.point = levels(cut(syngas.readcov,syngas.bin.read.cov))
syngas.readcov.quant.label = seq(length(syngas.readcov.quant.end.point),1,by=-1)


#
# For each TIS of a gene (operon), Q(SiteCov) and Q(Dist) are summed in the composite 
# quantile (CQ). For each gene (operon) sites displaying the lowest CQ are selected. 
# If the CQ of a site is lower than the central composite quantile, the site is accepted. 
#
fruct.tss.operon.dist = list(); syngas.tss.operon.dist = list(); tss.operon.dist = list()
fruct.tss.operon.site = list(); syngas.tss.operon.site = list(); tss.operon.site = list()
fruct.tss.operon.qu = list(); syngas.tss.operon.qu = list(); tss.operon.qu = list()
for(i in 1 : length(tssByOperon)){ # we consider TSSs outside and inside the operon
  
  if(length(tssByOperon[[i]]) >= 1){
    
    # Vector of TIS distances per operon
    delta.vec = c()
    # Vector of TIS coverages over all samples
    read.cov.vec = c()
    # Vector of TIS coverages on fructose samples 
    fruct.read.cov.vec = c()
    # Vector of TIS coverages on syngas samples 
    syngas.read.cov.vec = c()
    
    syngas.median.vec = c(); syngas.mad.vec = c()
    fruct.median.vec = c(); fruct.mad.vec = c()
    for(j in 1 : length(tssByOperon[[i]])){
      
      cur.delta = "NA"
      cur.read.cov = "NA"
      idx = match(tssByOperon[[i]][j], cap.site$TIS_ID)
      
      if(is.na(match(cap.site$Strand[idx], "+"))==FALSE){
        cur.delta = abs(cap.site$OperonFrom[idx] - cap.site$TIS_pos[idx])
        delta.vec = c(delta.vec, cur.delta)
      }
      
      if(is.na(match(cap.site$Strand[idx], "-"))==FALSE){
        cur.delta = abs(cap.site$TIS_pos[idx] - cap.site$OperonTo[idx])
        delta.vec = c(delta.vec, cur.delta)
      }
      
      # Site coverage over all samples
      cur.read.cov = cap.site$All_L2RelCov[idx]
      read.cov.vec = c(read.cov.vec, cur.read.cov)
      
      # Site coverage over samples grown on fructose
      cur.read.cov = cap.site$B_L2RelCov[idx]
      fruct.read.cov.vec = c(fruct.read.cov.vec, cur.read.cov)
      
      # Site coverage over samples grown on syngas
      cur.read.cov = cap.site$BR_L2RelCov[idx]
      syngas.read.cov.vec = c(syngas.read.cov.vec, cur.read.cov)
      
      fruct.read.cov = c(cap.site$B1_L2RelCov[idx], cap.site$B3_L2RelCov[idx], 
                         cap.site$B4_reads[idx])
      syngas.read.cov = c(cap.site$BR1_L2RelCov[idx], cap.site$BR3_L2RelCov[idx], 
                          cap.site$BR4_reads[idx])
      
      fruct.median.vec = c(fruct.median.vec, median(fruct.read.cov))
      syngas.median.vec = c(syngas.median.vec, median(syngas.read.cov))
      
      fruct.mad.vec = c(fruct.mad.vec, median(abs(fruct.read.cov - median(fruct.read.cov))))
      syngas.mad.vec = c(syngas.mad.vec, median(abs(syngas.read.cov - median(syngas.read.cov))))
      
    }
    
    # Map coverage of TISs per operon to Q(SiteCov)
    cur.readcov.quant.label = readcov.quant.label[match(cut(read.cov.vec,bin.read.cov), 
                                                        readcov.quant.end.point)]
    # Map coverage of TISs per operon on fructose to Q_fructose(SiteCov)
    cur.fruct.readcov.quant.label = fruct.readcov.quant.label[match(cut(fruct.read.cov.vec,fruct.bin.read.cov), 
                                                                    fruct.readcov.quant.end.point)]
    # Map coverage on of TISs per operon on syngas to Q_syngas(SiteCov)
    cur.syngas.readcov.quant.label = syngas.readcov.quant.label[match(cut(syngas.read.cov.vec,syngas.bin.read.cov), 
                                                                      syngas.readcov.quant.end.point)]
    # Map distance of TISs per operon to Q(Dist)
    cur.delta.quant.label = delta.quant.label[match(cut(delta.vec,bin.delta), 
                                                    delta.quant.end.point)]
    
    # Select TISs belonging to the lowest CQ (fructose)
    cur.fruct.label = rbind.data.frame(cur.fruct.readcov.quant.label, cur.delta.quant.label)
    tmp = colSums(cur.fruct.label)
    fruct.tss.operon.dist[[names(tssByOperon)[i]]] = delta.vec[is.na(match(tmp, min(tmp)))==FALSE]
    fruct.tss.operon.site[[names(tssByOperon)[i]]] = tssByOperon[[i]][is.na(match(tmp, min(tmp)))==FALSE]
    fruct.tss.operon.qu[[names(tssByOperon)[i]]] = tmp[is.na(match(tmp, min(tmp)))==FALSE]
    
    # Select TISs belonging to the lowest CQ (syngas)
    cur.syngas.label = rbind.data.frame(cur.syngas.readcov.quant.label, cur.delta.quant.label)
    tmp = colSums(cur.syngas.label)
    syngas.tss.operon.dist[[names(tssByOperon)[i]]] = delta.vec[is.na(match(tmp, min(tmp)))==FALSE]
    syngas.tss.operon.site[[names(tssByOperon)[i]]] = tssByOperon[[i]][is.na(match(tmp, min(tmp)))==FALSE]
    syngas.tss.operon.qu[[names(tssByOperon)[i]]] = tmp[is.na(match(tmp, min(tmp)))==FALSE]
    
    # Select TISs belonging to the lowest CQ 
    cur.label = rbind.data.frame(cur.readcov.quant.label, cur.delta.quant.label)
    tmp = colSums(cur.label)
    tss.operon.dist[[names(tssByOperon)[i]]] = delta.vec[is.na(match(tmp, min(tmp)))==FALSE]
    tss.operon.site[[names(tssByOperon)[i]]] = tssByOperon[[i]][is.na(match(tmp, min(tmp)))==FALSE]
    tss.operon.qu[[names(tssByOperon)[i]]] = tmp[is.na(match(tmp, min(tmp)))==FALSE]
    
  }
}

#
# Filter out the operons whose selected TISs display CQs lower than the central CQ. 
#
# Q(SiteCov) is based on fructose samples
distCQ.fruct = list()
siteCQ.fruct = list()
for(i in 1 : length(fruct.tss.operon.dist)){
  
  #check if there is at least a TIS whose CQ fulfills the criterion
  tmp.qu = fruct.tss.operon.qu[[i]]
  lgc = array(FALSE, length(tmp.qu))
  lgc[tmp.qu <= thrs.quantile] = TRUE
  
  tmp.dist = fruct.tss.operon.dist[[i]]
  tmp.site = fruct.tss.operon.site[[i]]
  
  if( length(lgc[lgc==TRUE]) > 0 ){
    distCQ.fruct[[names(fruct.tss.operon.dist)[i]]] = tmp.dist[lgc==TRUE]
    siteCQ.fruct[[names(fruct.tss.operon.site)[i]]] = tmp.site[lgc==TRUE]
  }
  
}

# Q(SiteCov) is based on syngas samples
distCQ.syngas = list()
siteCQ.syngas = list()
for(i in 1 : length(syngas.tss.operon.dist)){
  
  #check if there is at least a TIS whose CQ fulfills the criterion
  tmp.qu = syngas.tss.operon.qu[[i]]
  lgc = array(FALSE, length(tmp.qu))
  lgc[tmp.qu <= thrs.quantile] = TRUE
  
  tmp.dist = syngas.tss.operon.dist[[i]]
  tmp.site = syngas.tss.operon.site[[i]]
  
  if( length(lgc[lgc==TRUE]) > 0 ){
    distCQ.syngas[[names(syngas.tss.operon.dist)[i]]] = tmp.dist[lgc==TRUE]
    siteCQ.syngas[[names(syngas.tss.operon.site)[i]]] = tmp.site[lgc==TRUE]
  }
  
}

# Q(SiteCov) is based on all samples
distCQ = list()
siteCQ = list()
for(i in 1 : length(tss.operon.dist)){
  
  #check if there is at least a TIS whose CQ fulfills the criterion
  tmp.qu = tss.operon.qu[[i]]
  lgc = array(FALSE, length(tmp.qu))
  lgc[tmp.qu <= thrs.quantile] = TRUE
  
  tmp.dist = tss.operon.dist[[i]]
  tmp.site = tss.operon.site[[i]]
  
  if( length(lgc[lgc==TRUE]) > 0 ){
    distCQ[[names(tss.operon.dist)[i]]] = tmp.dist[lgc==TRUE]
    siteCQ[[names(tss.operon.site)[i]]] = tmp.site[lgc==TRUE]
  }
  
}


#-------------------------------------------------------------------------------
#
# Site filtering by the consistency of 5’UTR read-coverage profile with TIS usage.
# Site filtering by read-coverage profile over the genomic region spanned by gene (operon). 
#
#-------------------------------------------------------------------------------

# Load the FASTA file of the genome
fastafile<- read.fasta(file = "../data/GCA_040166795.1/GCA_040166795.1_LT_LAbrini_genomic.fna", 
                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)

#
# Check the validity of the filtering condition in fructose growth condition
#
# Load BAM files
fruct.bam.fs = dir("../data/Transcriptome_raw_data/", pattern="^FF1.*total.bam$")

# List containing for each sample a vector of logical values about whether an operon has at least a site meeting the criteria
lgcRCLstFruct = list()
# List containing for each sample the sites meeting the criteria for each operon
siteRCLstFruct = list()
# List containing for each sample the distances from the CDS start of the sites meeting the criteria for each operon
distRCLstFruct = list()
for(s in 1 : length(fruct.bam.fs)){
  
  lgc.by.operon = array(FALSE, length(distCQ.fruct))
  names(lgc.by.operon) = names(distCQ.fruct)
  
  site.by.operon <- vector(mode = "list", length = length(distCQ.fruct))
  names(site.by.operon) <- names(distCQ.fruct)
  
  dist.by.operon <- vector(mode = "list", length = length(distCQ.fruct))
  names(dist.by.operon) <- names(distCQ.fruct)
  
  outf = paste0("rc.tss.", fruct.bam.fs[s],".pdf")
  pdf(outf)
  
  for(i in 1 : length(distCQ.fruct)){
    
    #
    # Select the genomic range corresponding to the Region Of Interest (ROI) and
    # to the control region which is defined as the region upstream of the TSS of 100 nt
    #
    cur.op.strand = cap.site$Strand[match(names(distCQ.fruct)[i], 
                                          paste0("Operon-",cap.site$OperonFrom))]
    cur.op.from = cap.site$OperonFrom[match(names(distCQ.fruct)[i], 
                                            paste0("Operon-",cap.site$OperonFrom))]
    cur.op.to = cap.site$OperonTo[match(names(distCQ.fruct)[i], 
                                        paste0("Operon-",cap.site$OperonFrom))]
    
    cur.seq = ""; cur.label = ""; cur.genomic.range = ""; cur.genomic.range.ctrl = ""
    
    lgc.by.site = c()
    
    for(j in 1 : length(distCQ.fruct[[i]])){
      
      if(is.na(match(cur.op.strand, "+"))==FALSE){
        
        seqfrag = getFrag(object=fastafile, 
                          begin = cur.op.from - as.numeric(distCQ.fruct[[i]][j]), 
                          end = cur.op.from
        )
        cur.seq = getSequence(seqfrag, as.string = TRUE)[[1]][[1]][1]
        cur.label = names(distCQ.fruct)[i]
        
        # Region Of Interest (ROI)
        df = cbind.data.frame(start = cur.op.from - as.numeric(distCQ.fruct[[i]][j]), 
                              end = cur.op.from,
                              strand = cur.op.strand, 
                              seqid="NC_999999999.1"
        )
        
        cur.genomic.range = makeGRangesFromDataFrame(df,
                                                     seqnames.field="seqid", 
                                                     start.field="start",
                                                     end.field="end",
                                                     strand.field="strand",
                                                     ignore.strand=FALSE
        )
        
        # Control region
        df.ctrl = cbind.data.frame(start = cur.op.from - as.numeric(distCQ.fruct[[i]][j])-100, 
                                   end = cur.op.from - as.numeric(distCQ.fruct[[i]][j])-1,
                                   strand = cur.op.strand, 
                                   seqid="NC_999999999.1"
        )
        
        cur.genomic.range.ctrl = makeGRangesFromDataFrame(df.ctrl,seqnames.field="seqid", start.field="start",
                                                          end.field="end",strand.field="strand", ignore.strand=FALSE
        )
        
      }
      
      if(is.na(match(cur.op.strand, "-"))==FALSE){
        
        seqfrag = getFrag(object=fastafile, 
                          begin = cur.op.to, 
                          end = cur.op.to + as.numeric(distCQ.fruct[[i]][j]))
        cur.seq = getSequence(seqfrag, as.string = TRUE)[[1]][[1]][1]
        cur.label = names(distCQ.fruct[i])
        
        # Region Of Interest (ROI)
        df = cbind.data.frame(start = cur.op.to, 
                              end = cur.op.to + as.numeric(distCQ.fruct[[i]][j]),
                              strand = cur.op.strand, 
                              seqid="NC_999999999.1"
        )
        
        cur.genomic.range = makeGRangesFromDataFrame(df,
                                                     seqnames.field="seqid", 
                                                     start.field="start",
                                                     end.field="end",
                                                     strand.field="strand", 
                                                     ignore.strand=FALSE
        )
        
        # Control region
        df.ctrl = cbind.data.frame(start = cur.op.to + as.numeric(distCQ.fruct[[i]][j]) + 1, 
                                   end = cur.op.to + as.numeric(distCQ.fruct[[i]][j])+100,
                                   strand = cur.op.strand, 
                                   seqid="NC_999999999.1"
        )
        
        cur.genomic.range.ctrl = makeGRangesFromDataFrame(df.ctrl,
                                                          seqnames.field="seqid", 
                                                          start.field="start",
                                                          end.field="end",
                                                          strand.field="strand", 
                                                          ignore.strand=FALSE
        )
      }
      
      
      # Retrieve the read coverage of the ROI and of the control region 
      readCov = bamCoverage(bampath = paste0("../data/Transcriptome_raw_data/",fruct.bam.fs[s]), gr = cur.genomic.range, mapqual = 0, 
                            paired.end = "extend", tlenFilter = NULL,
                            filteredFlag = -1, verbose = TRUE)[1]
      
      readCovCtrl = bamCoverage(bampath = paste0("../data/Transcriptome_raw_data/",fruct.bam.fs[s]), gr = cur.genomic.range.ctrl, mapqual = 0, 
                                paired.end = "extend", tlenFilter = NULL,
                                filteredFlag = -1, verbose = TRUE)[1]
      
      
      #
      # Site filtering by the consistency of 5’UTR read-coverage profile with TIS usage: 
      #
      # - The ROI features a stretch of length at least equal to 50% of the ROI length 
      #  with a read coverage higher than 50% of that of the contiguous control positions 
      #  This condition corresponds to lgcSignal.
      #
      # - The ROI features a stretch of length greater or equal to 50% of the ROI length
      #  with a median read coverage higher than 5. 
      #  This condition corresponds to lgcWin
      #
      
      #
      # Compute the value of the lgcSignal variable 
      #
      lgcSignal = FALSE
      medianCtrl = median(readCovCtrl)
      # Calculate the moving average with a window size of width equal to slidingWindow
      slidingWindow = round(0.5 * length(readCov))
      rollresult <- rollapply(readCov, width = slidingWindow, FUN = median, align = "right", fill = NA)
      # Verify the validity of the condition 
      for(w in 1 : length(rollresult)){
        
        if( (is.na(rollresult[w])==FALSE) & (medianCtrl > 0) ){
          if( (rollresult[w]/medianCtrl) >= 1.5 ){
            lgcSignal = TRUE
            break
          }
        }
        if( (is.na(rollresult[w])==FALSE) & (rollresult[w] > 0) & (medianCtrl == 0) ){
          lgcSignal = TRUE
          break  
        }
      }
      
      #
      # Compute the value of the lgcWin variable 
      #
      
      # Calculate the moving average with a window size of width equal to slidingWindow
      slidingWindow = round(0.5 * length(readCov))
      lgcWin = FALSE
      rollresult <- rollapply(readCov, width = slidingWindow, FUN = median, align = "right", fill = NA)
      # Verify the validity of the condition 
      if(length(rollresult[(is.na(rollresult) == FALSE) & (rollresult > 5)]) >= 1 ){lgcWin=TRUE}
      
      # 
      # Combine the logical variables
      #
      lgc = FALSE
      if( (lgcSignal == TRUE) & (lgcWin == TRUE)  ){
        lgc = TRUE
      }
      if( (lgcSignal == FALSE) | (lgcWin == FALSE) ){
        lgc = FALSE
      }
      lgc.by.site[j] = lgc
      
      #
      # Plot section
      #
      if(is.na(match(cur.op.strand, "+"))==FALSE){
        
        readCovPlot = c(readCovCtrl, readCov)
        tss.idx = length(readCovCtrl+1)
        plot(readCovPlot, ylim=c(0, max(readCovPlot)), ylab="", xlab="", 
             main=paste0(names(distCQ.fruct)[i], "(", cur.op.strand, ")"), 
             type="l", axes=FALSE)
        
        #set the labels on the x-axis
        x.labels = seq(-length(readCovCtrl), length(readCov), by=20)
        #set the ticks on the x-axis
        x.ticks = seq(0, length(readCovPlot), by=20)
        axis(1, at = x.ticks, labels = x.labels)
        
        #set the ticks and labels on the y-axis
        stepY = ceiling(max(readCovPlot)/10)
        axis(2, at = seq(0, max(readCovPlot), by = stepY), labels = seq(0, max(readCovPlot), by = stepY))
        
        segments(x0 = tss.idx, y0 = 0, x1 = tss.idx, y1 = max(readCovPlot), col="red")
        
        mtext(side=3, line=-5, paste0(siteCQ.fruct[[i]][j], ": ", lgc.by.site[j]))
        
      }
      
      if(is.na(match(cur.op.strand, "-"))==FALSE){
        
        readCovPlot = c(readCov, readCovCtrl)
        tss.idx = length(readCov)
        plot(readCovPlot, ylim=c(0, max(readCovPlot)), ylab="", xlab="", 
             main=paste0(names(distCQ.fruct)[i], "(", cur.op.strand, ")"), 
             type="l", axes=FALSE)
        
        #set the labels on the x-axis
        x.labels = seq(-length(readCov), length(readCovCtrl), by=20)
        #set the ticks on the x-axis
        x.ticks = seq(0, length(readCovPlot), by=20)
        axis(1, at = x.ticks, labels = x.labels, las=2)
        
        #set the ticks and labels on the y-axis
        stepY = ceiling(max(readCovPlot)/10)
        axis(2, at = seq(0, max(readCovPlot), by = stepY), labels = seq(0, max(readCovPlot), by = stepY))
        
        segments(x0 = tss.idx, y0 = 0, x1 = tss.idx, y1 = max(readCovPlot), col="red")
        
        mtext(side=3, line=-5, paste0(siteCQ.fruct[[i]][j], ": ", lgc.by.site[j]))
        
      }
      
    }
    
    
    #
    # Site filtering by read-coverage profile over the genomic region spanned by gene (operon). 
    # The gene (operon) has consistent read coverage over the corresponding 
    # genomic region admitting deviations only in the presence of TSSs that are
    # present within the CDS and retained by the analysis
    # This condition corresponds to lgc.read.cov.by.operon.
    #
    
    lgc.constant.RC = FALSE
    if(is.na(match(names(distCQ.fruct)[i], op.unif.read.cov))==FALSE){lgc.constant.RC = TRUE}
    
    # Check if the operon has sites within CDS 
    retained.tss = siteCQ.fruct[[i]][lgc.by.site=TRUE]
    lgcInternal = FALSE
    if(is.na(match("cds", retained.tss))==FALSE){lgcInternal=TRUE}
    
    # Verify the validity of the condition
    lgc.read.cov.by.operon = FALSE
    if(lgc.constant.RC==TRUE){lgc.read.cov.by.operon = TRUE}
    if( (lgc.constant.RC==FALSE) & (lgcInternal==TRUE) ){ lgc.read.cov.by.operon = TRUE}
    
    
    #
    # Establish if an operon is assigned a site or not based on the previous filters 
    #
    
    if( (length(lgc.by.site[lgc.by.site==TRUE]) > 0) & (lgc.read.cov.by.operon == TRUE) ){
      lgc.by.operon[i] = TRUE
      site.by.operon[[i]] = siteCQ.fruct[[i]][lgc.by.site==TRUE]
      dist.by.operon[[i]] = distCQ.fruct[[i]][lgc.by.site==TRUE]
    }
    
  }
  dev.off()
  
  lgcRCLstFruct[[s]] = lgc.by.operon
  siteRCLstFruct[[s]] = site.by.operon
  distRCLstFruct[[s]] = dist.by.operon
  
}


#
# Check the validity of the filtering condition in syngas growth condition
#
syngas.bam.fs = dir("../data/Transcriptome_raw_data/", pattern="^F58.*total.bam$")

# List containing for each sample a vector of logical values about whether an operon has at least a site meeting the criteria
lgcRCLstSyngas = list()
# List containing for each sample the sites meeting the criteria for each operon
siteRCLstSyngas = list()
# List containing for each sample the distances from the CDS start of the sites meeting the criteria for each operon
distRCLstSyngas = list()
for(s in 1 : length(syngas.bam.fs)){
  
  lgc.by.operon=array(FALSE, length(distCQ.syngas))
  names(lgc.by.operon) = names(distCQ.syngas)
  
  site.by.operon <- vector(mode = "list", length = length(distCQ.syngas))
  names(site.by.operon) <- names(distCQ.syngas)
  
  dist.by.operon <- vector(mode = "list", length = length(distCQ.syngas))
  names(dist.by.operon) <- names(distCQ.syngas)
  
  outf = paste0("rc.tss.", syngas.bam.fs[s],".pdf")
  pdf(outf)
  
  for(i in 1 : length(distCQ.syngas)){
    
    #
    # Select the genomic range corresponding to the Region Of Interest (ROI) and
    # to the control region which is defined as the region upstream of the TSS of 100 nt
    #
    cur.op.strand = cap.site$Strand[match(names(distCQ.syngas)[i], 
                                          paste0("Operon-",cap.site$OperonFrom))]
    cur.op.from = cap.site$OperonFrom[match(names(distCQ.syngas)[i], 
                                            paste0("Operon-",cap.site$OperonFrom))]
    cur.op.to = cap.site$OperonTo[match(names(distCQ.syngas)[i], 
                                        paste0("Operon-",cap.site$OperonFrom))]
    
    cur.seq = ""; cur.label = ""; cur.genomic.range = ""; cur.genomic.range.ctrl = ""
    
    lgc.by.site = c()
    
    for(j in 1 : length(distCQ.syngas[[i]])){
      
      if(is.na(match(cur.op.strand, "+"))==FALSE){
        seqfrag = getFrag(object=fastafile, 
                          begin = cur.op.from - as.numeric(distCQ.syngas[[i]][j]), 
                          end = cur.op.from
        )
        cur.seq = getSequence(seqfrag, as.string = TRUE)[[1]][[1]][1]
        cur.label = names(distCQ.syngas)[i]
        
        # Region Of Interest (ROI)
        df = cbind.data.frame(start = cur.op.from - as.numeric(distCQ.syngas[[i]][j]), 
                              end = cur.op.from,
                              strand = cur.op.strand, 
                              seqid="NC_999999999.1"
        )
        
        cur.genomic.range = makeGRangesFromDataFrame(df,seqnames.field="seqid", 
                                                     start.field="start",
                                                     end.field="end",
                                                     strand.field="strand", 
                                                     ignore.strand=FALSE
        )
        
        # Control region
        df.ctrl = cbind.data.frame(start = cur.op.from - as.numeric(distCQ.syngas[[i]][j]) - 100, 
                                   end = cur.op.from - as.numeric(distCQ.syngas[[i]][j]) - 1,
                                   strand = cur.op.strand, 
                                   seqid="NC_999999999.1"
        )
        
        cur.genomic.range.ctrl = makeGRangesFromDataFrame(df.ctrl,
                                                          seqnames.field="seqid", 
                                                          start.field="start",
                                                          end.field="end",
                                                          strand.field="strand", 
                                                          ignore.strand=FALSE
        )
        
      }
      
      if(is.na(match(cur.op.strand, "-"))==FALSE){
        
        seqfrag = getFrag(object=fastafile, 
                          begin = cur.op.to, 
                          end = cur.op.to + as.numeric(distCQ.syngas[[i]][j])
        )
        cur.seq = getSequence(seqfrag, as.string = TRUE)[[1]][[1]][1]
        cur.label = names(distCQ.syngas[i])
        
        # Region of interest (ROI)
        df = cbind.data.frame(start = cur.op.to, 
                              end = cur.op.to + as.numeric(distCQ.syngas[[i]][j]),
                              strand = cur.op.strand, 
                              seqid="NC_999999999.1"
        )
        
        cur.genomic.range = makeGRangesFromDataFrame(df,
                                                     seqnames.field="seqid", 
                                                     start.field="start",
                                                     end.field="end",
                                                     strand.field="strand", 
                                                     ignore.strand=FALSE
        )
        # Control region 
        df.ctrl = cbind.data.frame(start = cur.op.to + as.numeric(distCQ.syngas[[i]][j]) + 1, 
                                   end = cur.op.to + as.numeric(distCQ.syngas[[i]][j]) + 100,
                                   strand = cur.op.strand, 
                                   seqid="NC_999999999.1"
        )
        
        cur.genomic.range.ctrl = makeGRangesFromDataFrame(df.ctrl,
                                                          seqnames.field="seqid", 
                                                          start.field="start",
                                                          end.field="end",
                                                          strand.field="strand", 
                                                          ignore.strand=FALSE
        )
      }
      
      # Retrieve the read coverage of the ROI and of the control region 
      readCov = bamCoverage(bampath = paste0("../data/Transcriptome_raw_data/",syngas.bam.fs[s]), gr = cur.genomic.range, mapqual = 0, 
                            paired.end = "extend", tlenFilter = NULL,
                            filteredFlag = -1, verbose = TRUE)[1]
      
      readCovCtrl = bamCoverage(bampath = paste0("../data/Transcriptome_raw_data/",syngas.bam.fs[s]), gr = cur.genomic.range.ctrl, mapqual = 0, 
                                paired.end = "extend", tlenFilter = NULL,
                                filteredFlag = -1, verbose = TRUE)[1]
      
      
      #
      # Site filtering by the consistency of 5’UTR read-coverage profile with TIS usage: 
      #
      # - The ROI features a stretch of length at least equal to 50% of the ROI length 
      #  with a read coverage higher than 50% of that of the contiguous control positions 
      #  This condition corresponds to lgcSignal.
      #
      # - The ROI features a median read coverage of at least 5 over at least 50% of the ROI.
      #  This condition corresponds to lgcWin
      #
      
      #
      # Compute the value of the lgcSignal variable 
      #
      lgcSignal = FALSE
      medianCtrl = median(readCovCtrl)
      slidingWindow = round(0.5 * length(readCov))
      rollresult <- rollapply(readCov, width = slidingWindow, FUN = median, align = "right", fill = NA)
      # Verify the validity of the condition 
      for(w in 1 : length(rollresult)){
        
        if( (is.na(rollresult[w])==FALSE) & (medianCtrl > 0) ){
          if( (rollresult[w]/medianCtrl) >= 1.5 ){
            lgcSignal = TRUE
            break
          }
        }
        if( (is.na(rollresult[w])==FALSE) & (rollresult[w] > 0) & (medianCtrl == 0) ){
          lgcSignal = TRUE
          break
        }
      }
      
      #
      # Compute the value of the lgcWin variable 
      #
      
      # Calculate the moving average with a window size of width equal to slidingWindow
      slidingWindow = round(0.5 * length(readCov))
      lgcWin = FALSE
      rollresult <- rollapply(readCov, width = slidingWindow, FUN = median, align = "right", fill = NA)
      # Verify the validity of the condition
      if(length(rollresult[(is.na(rollresult) == FALSE) & (rollresult > 5)]) >= 1 ){lgcWin=TRUE}
      
      # 
      # Combine the logical variables
      #
      lgc = FALSE
      if( (lgcSignal == TRUE) & (lgcWin == TRUE) ){
        lgc = TRUE
      }
      if( (lgcSignal == FALSE) | (lgcWin == FALSE) ){
        lgc = FALSE
      }
      lgc.by.site[j] = lgc
      
      
      #
      # Plot section
      #
      if(is.na(match(cur.op.strand, "+"))==FALSE){
        
        readCovPlot = c(readCovCtrl, readCov)
        tss.idx = length(readCovCtrl+1)
        plot(readCovPlot, ylim=c(0, max(readCovPlot)), ylab="", xlab="", 
             main=paste0(names(distCQ.syngas)[i], "(", cur.op.strand, ")"), type="l", axes=FALSE)
        
        #set the labels on the x-axis
        x.labels = seq(-length(readCovCtrl), length(readCov), by=20)
        #set the ticks on the x-axis
        x.ticks = seq(0, length(readCovPlot), by=20)
        axis(1, at = x.ticks, labels = x.labels)
        
        #set the ticks and labels on the y-axis
        stepY = ceiling(max(readCovPlot)/10)
        axis(2, at = seq(0, max(readCovPlot), by = stepY), labels = seq(0, max(readCovPlot), by = stepY))
        
        segments(x0 = tss.idx, y0 = 0, x1 = tss.idx, y1 = max(readCovPlot), col="red")
        
        mtext(side=3, line=-5, paste0(siteCQ.syngas[[i]][j], ": ", lgc.by.site[j]))
        
      }
      
      if(is.na(match(cur.op.strand, "-"))==FALSE){
        
        readCovPlot = c(readCov, readCovCtrl)
        tss.idx = length(readCov)
        plot(readCovPlot, ylim=c(0, max(readCovPlot)), ylab="", xlab="", 
             main=paste0(names(distCQ.syngas)[i], "(", cur.op.strand, ")"), type="l", axes=FALSE)
        
        #set the labels on the x-axis
        x.labels = seq(-length(readCov), length(readCovCtrl), by=20)
        #set the ticks on the x-axis
        x.ticks = seq(0, length(readCovPlot), by=20)
        axis(1, at = x.ticks, labels = x.labels, las=2)
        
        #set the ticks and labels on the y-axis
        stepY = ceiling(max(readCovPlot)/10)
        axis(2, at = seq(0, max(readCovPlot), by = stepY), labels = seq(0, max(readCovPlot), by = stepY))
        
        segments(x0 = tss.idx, y0 = 0, x1 = tss.idx, y1 = max(readCovPlot), col="red")
        
        mtext(side=3, line=-5, paste0(siteCQ.syngas[[i]][j], ": ", lgc.by.site[j]))
        
        
      }
      
    }
    
    
    #
    # Site filtering by read-coverage profile over the genomic region spanned by gene (operon). 
    # The gene (operon) has consistent read coverage over the corresponding 
    # genomic region admitting deviations only in the presence of TSSs that are
    # present within the CDS and retained by the analysis
    # This condition corresponds to lgc.read.cov.by.operon.
    # 
    
    lgc.constant.RC = FALSE
    if(is.na(match(names(distCQ.syngas)[i], op.unif.read.cov))==FALSE){lgc.constant.RC = TRUE}
    #Check if the operon has sites within the CDS
    retained.tss = siteCQ.syngas[[i]][lgc.by.site=TRUE]
    lgcInternal = FALSE
    if(is.na(match("cds", retained.tss))==FALSE){lgcInternal=TRUE}
    #Define the logical variable about consistent read coverage in the light of TSSs
    lgc.read.cov.by.operon = FALSE
    if(lgc.constant.RC==TRUE){lgc.read.cov.by.operon = TRUE}
    if( (lgc.constant.RC==FALSE) & (lgcInternal==TRUE) ){ lgc.read.cov.by.operon = TRUE}
    
    
    #
    # Establish if a gene (operon) can be assigned a site or not based on the previous filters 
    #
    
    if( (length(lgc.by.site[lgc.by.site==TRUE]) > 0) & (lgc.read.cov.by.operon == TRUE) ){
      lgc.by.operon[i] = TRUE
      site.by.operon[[i]] = siteCQ.syngas[[i]][lgc.by.site==TRUE]
      dist.by.operon[[i]] = distCQ.syngas[[i]][lgc.by.site==TRUE]
    }
    
  }
  dev.off()
  
  lgcRCLstSyngas[[s]] = lgc.by.operon
  siteRCLstSyngas[[s]] = site.by.operon
  distRCLstSyngas[[s]] = dist.by.operon
}

#
# save the results of the check of TSSs according to 5' UTR read coverage
#
tssRC = list(ops.syngas = ops.syngas, lgcRCLstSyngas = lgcRCLstSyngas, 
             distRCLstSyngas = distRCLstSyngas, siteRCLstSyngas = siteRCLstSyngas, 
             ops.fruct = ops.fruct, lgcRCLstFruct = lgcRCLstFruct, 
             distRCLstFruct = distRCLstFruct, siteRCLstFruct = siteRCLstFruct
)

save(ttsRC, file = "../results/TIS-T3E/site.by.5utr.read.coverage.rda")


#--------------------------------------------------------------------------------------
# Site clustering by physical proximity. From sites closer than 5 nt from each other, 
# the TIS closest to the start codon of the downstream CDS is the primary site. 
#--------------------------------------------------------------------------------------

load("../results/TIS-T3E/site.by.5utr.read.coverage.rda")

#
# Compile the final set of sites per gene (operon) by considering representative sites per site clusters in fructose
#
final.tss.by.operon = list()
final.tss.short.dist.by.operon = c()
for(i in 1 : length(tssRC$ops.fruct)){
  
  final.site.by.operon = c()
  
  # Unique sites per operon over all samples (including sites that are close to each other)
  site.by.operon = c()
  for(j in 1 : length(tssRC$siteRCLstFruct)){
    idx = match(tssRC$ops.fruct[i], names(tssRC$siteRCLstFruct[[j]]))
    site.by.operon = c(site.by.operon, tssRC$siteRCLstFruct[[j]][idx])
  }
  site.by.operon = unique(unlist(site.by.operon))
  
  
  #
  # Consider the scenario when the operon has multiple TSSs
  #
  
  if(length(site.by.operon)>1){
    site.pos.vec = c()
    for(j in 1 : length(site.by.operon)){
      site.pos.vec[j] = tbseq.cap.term$cap_site$TIS_pos[match(site.by.operon[j],
                                                              tbseq.cap.term$cap_site$TIS_ID)]
    }
    
    inter.dist = c()
    for(j in 2 : length(site.pos.vec)){
      inter.dist = c( inter.dist, abs(site.pos.vec[j] - site.pos.vec[j-1]) )
    }
    
    # Identify the indices of the sites that are distant by more than 5 nt
    retained.idx.by.operon = c()
    #
    # Consider the case where an operon displays more than two sites
    #
    if(length(inter.dist)>1){
      for(j in 1 : length(inter.dist)){
        if( (j!=1) & (j!= length(inter.dist)) & inter.dist[j]>5 ){
          retained.idx.by.operon = c(retained.idx.by.operon, c(j,j+1))
        }
        if( (j==1) & (inter.dist[j]>5) ){
          retained.idx.by.operon = c(retained.idx.by.operon, c(j, j+1))
        }
        if( (inter.dist[j]>5) && j==(length(inter.dist)) && (inter.dist[j-1] > 5) ){
          retained.idx.by.operon = c(retained.idx.by.operon, c(j,j+1))
        }
        if( (inter.dist[j]>5) && (j==length(inter.dist)) && (inter.dist[j-1] <= 5) ){
          retained.idx.by.operon = c(retained.idx.by.operon, (j+1))
        }
        if( length(inter.dist[inter.dist<=5]) == length(inter.dist) ){
          retained.idx.by.operon = c(retained.idx.by.operon, length(site.pos.vec))
        }
      }
    }
    #
    # Consider the case where an operon displays two sites
    #
    if(length(inter.dist)==1){
      if(inter.dist>5){
        retained.idx.by.operon = c(retained.idx.by.operon, 1:2)
      }
      if(inter.dist<=5){
        retained.idx.by.operon = c(retained.idx.by.operon, 2)
      }
    }
    
    # Collect sites corresponding to distinct positions
    final.site.by.operon = site.by.operon[unique(retained.idx.by.operon)]
    final.tss.by.operon[[tssRC$ops.fruct[i]]] = final.site.by.operon
  }
  
  #
  # Consider the scenario where the operon displays a single TSS
  #
  if(length(site.by.operon)==1){
    final.site.by.operon = site.by.operon
    final.tss.by.operon[[tssRC$ops.fruct[i]]] = final.site.by.operon
  }
  
  #
  # Compute the shortest distance of the final sites per operon
  #
  final.dist.pos.vec = c()
  for(j in 1 : length(final.site.by.operon)){
    idx = match(final.site.by.operon[j], siteCQ.fruct[[tssRC$ops.fruct[i]]])
    final.dist.pos.vec[j] = distCQ.fruct[[tssRC$ops.fruct[i]]][idx]
  }
  final.tss.short.dist.by.operon[tssRC$ops.fruct[i]] = min(final.dist.pos.vec)
  
}
final.tss.by.operon.fruct = final.tss.by.operon
final.tss.short.dist.by.operon.fruct = final.tss.short.dist.by.operon


#
# Compile the final set of sites per gene (operon) by considering representative sites per site clusters in syngas
#
final.tss.by.operon = list()
final.tss.short.dist.by.operon = c()
for(i in 1 : length(tssRC$ops.syngas)){
  
  final.site.by.operon = c()
  
  # Unique sites per operon over all samples (including sites that are close to each other)
  site.by.operon = c()
  for(j in 1 : length(tssRC$siteRCLstSyngas)){
    idx = match(tssRC$ops.syngas[i], names(tssRC$siteRCLstSyngas[[j]]))
    site.by.operon = c(site.by.operon, tssRC$siteRCLstSyngas[[j]][idx])
  }
  site.by.operon = unique(unlist(site.by.operon))
  
  #
  # Consider the scenario where an operon has multiple TSSs
  #
  if(length(site.by.operon)>1){
    site.pos.vec = c()
    for(j in 1 : length(site.by.operon)){
      site.pos.vec[j] = tbseq.cap.term$cap_site$TIS_pos[match(site.by.operon[j],
                                                              tbseq.cap.term$cap_site$TIS_ID)]
    }
    
    inter.dist = c()
    for(j in 2 : length(site.pos.vec)){
      inter.dist = c( inter.dist, abs(site.pos.vec[j] - site.pos.vec[j-1]) )
    }
    
    # Identify the indices of the sites that are distant by more than 5 nt
    retained.idx.by.operon = c()
    #
    # Conisder the case where an operon displays more than two sites
    #
    if(length(inter.dist)>1){
      for(j in 1 : length(inter.dist)){
        if( (j!=1) & (j!= length(inter.dist)) & inter.dist[j]>5 ){
          retained.idx.by.operon = c(retained.idx.by.operon, c(j,j+1))
        }
        if( (j==1) & (inter.dist[j]>5) ){
          retained.idx.by.operon = c(retained.idx.by.operon, c(j, j+1))
        }
        if( (inter.dist[j]>5) && j==(length(inter.dist)) && (inter.dist[j-1] > 5) ){
          retained.idx.by.operon = c(retained.idx.by.operon, c(j,j+1))
        }
        if( (inter.dist[j]>5) && (j==length(inter.dist)) && (inter.dist[j-1] <= 5) ){
          retained.idx.by.operon = c(retained.idx.by.operon, (j+1))
        }
        if( length(inter.dist[inter.dist<=5]) == length(inter.dist) ){
          retained.idx.by.operon = c(retained.idx.by.operon, length(site.pos.vec))
        }
      }
    }
    
    #
    # Consider the case where an operon displays two sites
    #
    if(length(inter.dist)==1){
      if(inter.dist>5){
        retained.idx.by.operon = c(retained.idx.by.operon, 1:2)
      }
      if(inter.dist<=5){
        retained.idx.by.operon = c(retained.idx.by.operon, 2)
      }
    }
    # Collect sites corresponding to distinct positions
    final.site.by.operon = site.by.operon[unique(retained.idx.by.operon)]
    final.tss.by.operon[[tssRC$ops.syngas[i]]] = final.site.by.operon
  }
  
  #
  # Consider the scenario where the operon displays a single TSS
  #
  if(length(site.by.operon)==1){
    final.site.by.operon = site.by.operon
    final.tss.by.operon[[tssRC$ops.syngas[i]]] = final.site.by.operon
  }
  
  #
  # Compute the shortest distance of the final sites per operon
  #
  final.dist.pos.vec = c()
  for(j in 1 : length(final.site.by.operon)){
    
    idx = match(final.site.by.operon[j], siteCQ.syngas[[tssRC$ops.syngas[i]]])
    final.dist.pos.vec[j] = distCQ.syngas[[tssRC$ops.syngas[i]]][idx]
    
  }
  final.tss.short.dist.by.operon[tssRC$ops.syngas[i]] = min(final.dist.pos.vec)
  
}
final.tss.by.operon.syngas = final.tss.by.operon
final.tss.short.dist.by.operon.syngas = final.tss.short.dist.by.operon


#-------------------------------------------------------------------------------
#
# Plot the results of TIS filtering
#
#-------------------------------------------------------------------------------

# Load gene (operon) annotations
f = read.xlsx("../data/OperonDescriptions.xlsx", sheetIndex = 1)
gene.set.by.operon = list()
for(i in 1 : dim(f)[1]){
  geneTagTmp <- strsplit(f$GeneTags[i], "([^0-9])")[[1]]
  GeneTags = c()
  for(j in 1 : length(geneTagTmp)){
    if(nchar(geneTagTmp[j])>0){
      GeneTags = c(GeneTags, geneTagTmp[j])
    }
  }
  LABRINI = paste0("LAbrini_", GeneTags)
  gene.set.by.operon[[paste0("Operon-",f$OperonID[i])]] = LABRINI
  
}
operon.size = as.numeric(unlist(lapply(gene.set.by.operon, function(x) length(x))))
names(operon.size) = names(gene.set.by.operon)


jpeg("../results/TIS-T3E/TIS.jpeg", res=600, width=4000, height=2500)
par(mfrow=c(1,2))
par(mar=c(7.5,6,1,1))

#
# Display the 5' UTR length in syngas and fructose
#
# Fructose
final.tss.short.dist.by.gene.fruct = c()
for( i in 1 : length(final.tss.short.dist.by.operon.fruct)){
  tmp = gene.set.by.operon[match(names(final.tss.short.dist.by.operon.fruct)[i], names(gene.set.by.operon))][[1]]
  for(j in 1 : length(tmp)){
    final.tss.short.dist.by.gene.fruct[tmp[j]] = as.numeric(final.tss.short.dist.by.operon.fruct[i])
  }
}

counts = c()
bins = seq(0,2000,by=100)
for(i in 1 : (length(bins)-1)){
  counts[i] = length(final.tss.short.dist.by.gene.fruct[final.tss.short.dist.by.gene.fruct > bins[i] & 
                                                          final.tss.short.dist.by.gene.fruct < bins[i+1]])
}
counts.fruct = counts
# Syngas
final.tss.short.dist.by.gene.syngas = c()
for( i in 1 : length(final.tss.short.dist.by.operon.syngas)){
  tmp = gene.set.by.operon[match(names(final.tss.short.dist.by.operon.syngas)[i], names(gene.set.by.operon))][[1]]
  for(j in 1 : length(tmp)){
    final.tss.short.dist.by.gene.syngas[tmp[j]] = as.numeric(final.tss.short.dist.by.operon.syngas[i])
  }
}

counts = c()
bins = seq(0,2000,by=100)
for(i in 1 : (length(bins)-1)){
  counts[i] = length(final.tss.short.dist.by.operon.syngas[final.tss.short.dist.by.operon.syngas > bins[i] & 
                                                             final.tss.short.dist.by.operon.syngas < bins[i+1]])
}
counts.syngas = counts

# Plotting 
mpl = t(as.matrix(cbind.data.frame(counts.fruct, counts.syngas)))

no.bins.plot = 5
b = barplot(mpl[,1:no.bins.plot], plot=FALSE)
barplot(mpl[,1:no.bins.plot], ylab="Frequency", beside=TRUE, col=c("#e08214", "#8073ac"),las=2, 
        border=NA, cex.axis=1.25, cex.names=1.65)
legend("topright", legend=c("Fructose","Syngas"), fill=c("#e08214", "#8073ac"), bty="n", 
       border="white", ncol=1, cex=1.15)

binlabels = c()
for(i in 1 : (length(bins)-1)){
  binlabels[i] = paste0("[", bins[i], ",", bins[i+1], "]")
}
axis(1, at=seq(2,length(b)*3,by=3), labels = binlabels[1:no.bins.plot], las=2, cex=.85)
mtext(side=1, line=6, text="TIS distance from start codon")

#
# Display the number of genes by the number of selected TIS sites per gene (fructose and syngas)
#
final.tss.by.gene.syngas = list()
for(i in  1 : length(final.tss.by.operon.syngas)){
  tmp = gene.set.by.operon[match(names(final.tss.by.operon.syngas)[i], names(gene.set.by.operon))][[1]]
  for(j in 1 : length(tmp)){
    final.tss.by.gene.syngas[[tmp[j]]] = final.tss.by.operon.syngas[[i]]
  }
}

final.tss.by.gene.fruct = list()
for(i in  1 : length(final.tss.by.operon.fruct)){
  tmp = gene.set.by.operon[match(names(final.tss.by.operon.fruct)[i], names(gene.set.by.operon))][[1]]
  for(j in 1 : length(tmp)){
    final.tss.by.gene.fruct[[tmp[j]]] = final.tss.by.operon.fruct[[i]]
  }
}

final.gene = union(names(final.tss.by.gene.syngas), names(final.tss.by.gene.fruct))
final.tss.by.gene = list()
for(i in 1 : length(final.gene)){
  final.tss.by.gene[[final.gene[i]]] = unique(c(final.tss.by.gene.syngas[[final.gene[i]]], final.tss.by.gene.fruct[[final.gene[i]]]))
  if(length(final.tss.by.gene[[final.gene[i]]])==5){show(final.gene[i])}
}

gene.freq.by.tss.no.fruct = table(unlist(lapply(final.tss.by.gene.fruct, function(x) length(x))))
gene.freq.by.tss.no.syngas = table(unlist(lapply(final.tss.by.gene.syngas, function(x) length(x))))
gene.freq.by.tss.no = table(unlist(lapply(final.tss.by.gene, function(x) length(x))))
tss.no = unique(c(names(gene.freq.by.tss.no.fruct), names(gene.freq.by.tss.no.syngas), names(gene.freq.by.tss.no)))

m = matrix(0, nrow=length(tss.no), ncol=3, dimnames=list(tss.no, c("Fructose","Syngas","Either subs.")))
for(i in 1 : length(tss.no)){
  m[i,1] = as.numeric(gene.freq.by.tss.no.fruct[match(tss.no[i], names(gene.freq.by.tss.no.fruct))])
  m[i,2] = as.numeric(gene.freq.by.tss.no.syngas[match(tss.no[i], names(gene.freq.by.tss.no.syngas))])
  m[i,3] = as.numeric(gene.freq.by.tss.no[match(tss.no[i], names(gene.freq.by.tss.no))])
}

barplot(t(m[,1:2]),beside=TRUE, col=c("#e08214", "#8073ac"), ylab="Frequency", las=1, cex.axis=1.25, cex.names=1.45)
legend("topright", legend=c("Fructose","Syngas"), fill=c("#e08214", "#8073ac"), bty="n", 
       border="white", ncol=1, cex=1.15)
mtext(side=1, line=3, text="TIS no. by gene")

dev.off()


#-------------------------------------------------------------------
# Print the final TISs
#-------------------------------------------------------------------

# Load gene (operon) annotations
foperon = read.xlsx("../data/Operon.description.gene.annot.xlsx", sheetIndex=1)
# Obtain gene sets per operon
gene.set.by.operon = list()
for(i in 1 : dim(foperon)[1]){
  LABRINI = foperon$GenesInOperon_CP110420[i]
  gene.set.by.operon[[paste0("Operon-",f$OperonID[i])]] = strsplit(LABRINI, ",")[[1]]
}


# Collect the final TISs for each operon on either growth condition 
ops = union(names(final.tss.by.operon.syngas), names(final.tss.by.operon.fruct))

tssByOperon = list()
tssByGene = list()
for(i in 1 : length(ops)){
  cur.tss = c()
  
  idx = match(ops[i], names(final.tss.by.operon.fruct))
  cur.tss = c(cur.tss, unlist(final.tss.by.operon.fruct[[idx]]))
  
  idx = match(ops[i], names(final.tss.by.operon.syngas))
  cur.tss = c(cur.tss, unlist(final.tss.by.operon.syngas[[idx]]))
  
  # Assign TISs to operons
  tssByOperon[[ops[i]]] = unique(cur.tss)
  
  # Assign TISs to genes
  cur.gene.by.operon = gene.set.by.operon[[ops[i]]]
  for( j in 1 : length(cur.gene.by.operon) ){
    tssByGene[[cur.gene.by.operon[j]]] = unique(cur.tss)
  }
}

# Print TISs by operon
out.name = names(tssByOperon)
out.signal = as.vector(
  unlist(lapply(tssByOperon, function(x) paste(x, collapse=",")))
)
out = cbind(out.name, out.signal)
colnames(out) = c("Operon","TIS")
write.xlsx(out, file = "../results/TIS-T3E/TIS.xlsx", row.names=FALSE, sheetName = "TIS_by_operon")

# Print TISs by gene
out.name = names(tssByGene)
out.signal = as.vector(
  unlist(lapply(tssByGene, function(x) paste(x, collapse=",")))
)
out = cbind(out.name, out.signal)
colnames(out) = c("Gene","TIS")
write.xlsx(out, file = "../results/TIS-T3E/TIS.xlsx", row.names=FALSE, sheetName = "TIS_by_gene", append=TRUE) 

