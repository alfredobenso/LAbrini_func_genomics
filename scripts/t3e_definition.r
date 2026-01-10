

require(stringr)
require(xlsx)
require(seqinr) 
require(stringr)
require(colorRamp2)
require(GenomicRanges)
require(bamsignals) 
require(zoo) 


#-------------------------------------------------------------------
# Retrieve NCBI annotations for Clostridium autoethanogenum
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
# Filter out T3Es whose positions coincide with CDS stops
#-----------------------------------------------------------------------

lgc.tts = array(TRUE, dim(tbseq.cap.term$term_site)[1])
for(i in 1 : dim(tbseq.cap.term$term_site)[1]){
  
  tts.strand = tbseq.cap.term$term_site$Strand[i]
  if(is.na(match(tbseq.cap.term$term_site$Strand[i], "+"))==FALSE){
    cur.delta = tbseq.cap.term$term_site$TTS_pos[i] - tbseq.cap.term$term_site$OperonTo[i]
    if(cur.delta == 0){lgc.tts[i] = FALSE}
  }
  if(is.na(match(tbseq.cap.term$term_site$Strand[i], "-"))==FALSE){
    cur.delta = tbseq.cap.term$term_site$OperonFrom[i] - tbseq.cap.term$term_site$TTS_pos[i]
    if(cur.delta == 0){lgc.tts[i] = FALSE}
  }
}
term.site = tbseq.cap.term$term_site[lgc.tts==TRUE,]

#-----------------------------------------------------------------------
# Assemble T3Es by gene (operon)
#-----------------------------------------------------------------------

# Define gene (operon) identifiers by their "OperonFrom" positions
operonIDs = unique(term.site$OperonFrom)

# Retrieve the TTS per operon 
ttsByOperon = list()
for(i in 1 : length(operonIDs)){
  tmp = term.site[is.na(match(term.site$OperonFrom, operonIDs[i]))==FALSE,]
  ttsByOperon[[paste0("Operon-",operonIDs[i])]] = tmp$TTS_ID
}
ttsByOperonCounts = as.numeric(lapply(ttsByOperon, function(x) length(x)))


#--------------------------------------------------------------------------
# Site filtering by site coverage and distance from CDS stop  
#--------------------------------------------------------------------------

step.quantile = .20
thrs.quantile = 5

#
# Discretize the distance between T3Es and CDS stop positions.
# Discretize T3Es read coverage gauged on fructose and syngas. 
#

# Collect distances between T3Es and the start codons of the downstream CDSs
delta = c()
# Collect the coverage of T3Es per operon normalized by the average coverage of all sites from all samples
readcov = c()
# Collect the coverage of T3Es per operon normalized by the average coverage of all sites from fructose samples
fruct.readcov = c()
# Collect the coverage of T3Es per operon normalized by the average coverage of all sites from syngas samples
syngas.readcov = c()
# Compute the median read coverage of T3Es per operon on fructose and syngas
fruct.median = c()
syngas.median = c()
# Compute the MAD of read coverage of T3Es per operon on fructose and syngas
fruct.mad = c()
syngas.mad = c()
for(i in 1 : dim(term.site)[1]){
  
  cur_delta=NA
  if(is.na(match(term.site$Strand[i], "+"))==FALSE){
    cur.delta = abs(term.site$TTS_pos[i] - term.site$OperonTo[i])
  }
  if(is.na(match(term.site$Strand[i], "-"))==FALSE){
    cur.delta = abs(term.site$OperonFrom[i] - term.site$TTS_pos[i])
  }
  delta = c(delta, cur.delta)
  
  # T3E coverage normalized by the average coverage of all sites from all samples
  readcov = c(readcov, term.site$All_L2RelCov[i])
  # T3E coverage normalized by the average coverage of all sites from fructose samples
  fruct.readcov = c(fruct.readcov, term.site$B_L2RelCov[i])
  # T3E coverage normalized by the average coverage of all sites from fructose samples
  syngas.readcov = c(syngas.readcov, term.site$BR_L2RelCov[i])
  
  fruct.read.cov = c(term.site$B1_L2RelCov[i], term.site$B3_L2RelCov[i], 
                     term.site$B4_reads[i])
  syngas.read.cov = c(term.site$BR1_L2RelCov[i], term.site$BR3_L2RelCov[i], 
                      term.site$BR4_reads[i])
  
  fruct.median = c(fruct.median, median(fruct.read.cov))
  syngas.median = c(syngas.median, median(syngas.read.cov))
  
  fruct.mad = c(fruct.mad, median(abs(fruct.read.cov - median(fruct.read.cov))))
  syngas.mad = c(syngas.mad, median(abs(syngas.read.cov - median(syngas.read.cov))))
  
}

# Create quantiles of Dist(T3E, CDS stop): Q(Dist)
bin.delta = quantile(delta, probs=seq(0,1,step.quantile))
bin.delta[1] = bin.delta[1] - 1 

# Label distance quantiles
delta.quant.end.point = levels(cut(delta,bin.delta))
delta.quant.label = 1:length(delta.quant.end.point)

# Create T3E coverage quantiles: Q(SiteCov)
bin.read.cov = quantile(readcov, probs=seq(0,1,step.quantile))
bin.read.cov[1] = bin.read.cov[1] - 1 

# Label T3E coverage quantiles
readcov.quant.end.point = levels(cut(readcov,bin.read.cov))
readcov.quant.label = seq(length(readcov.quant.end.point),1,by=-1)

# Create T3E coverage quantiles on fructose samples: Q_fructose(SiteCov)
fruct.bin.read.cov = quantile(fruct.readcov, probs=seq(0,1,step.quantile))
fruct.bin.read.cov[1] = fruct.bin.read.cov[1] - 1 

# T3E coverage quantiles on fructose samples
fruct.readcov.quant.end.point = levels(cut(fruct.readcov,fruct.bin.read.cov))
fruct.readcov.quant.label = seq(length(fruct.readcov.quant.end.point),1,by=-1)

# Create T3E coverage quantiles on syngas samples: Q_syngas(SiteCov)
syngas.bin.read.cov = quantile(syngas.readcov, probs=seq(0,1,step.quantile))
syngas.bin.read.cov[1] = syngas.bin.read.cov[1] - 1 

# Label T3E coverage quantiles on syngas samples
syngas.readcov.quant.end.point = levels(cut(syngas.readcov,syngas.bin.read.cov))
syngas.readcov.quant.label = seq(length(syngas.readcov.quant.end.point),1,by=-1)


#
# For each T3E of a gene (operon), Q(SiteCov) and Q(Dist) are summed in the composite 
# quantile (CQ). For each gene (operon) sites displaying the lowest CQ are selected. 
# If the CQ of a site is lower than the central composite quantile, the site is accepted. 
#
fruct.tts.operon.dist = c(); syngas.tts.operon.dist = c(); tts.operon.dist = c()
fruct.tts.operon.site = list(); syngas.tts.operon.site = list(); tts.operon.site = list()
fruct.tts.operon.qu = list(); syngas.tts.operon.qu = list(); tts.operon.qu = list()
for(i in 1 : length(ttsByOperon)){
  
  if(length(ttsByOperon[[i]]) >= 1){
    
    # Vector of T3E distances per operon
    delta.vec = c()
    # Vector of T3E coverages
    read.cov.vec = c()
    # Vector of T3E coverages on samples grown on fructose
    fruct.read.cov.vec = c()
    # Vector of T3E coverages on samples grown on syngas
    syngas.read.cov.vec = c()
    
    syngas.median.vec = c(); syngas.mad.vec = c()
    fruct.median.vec = c(); fruct.mad.vec = c()
    
    for(j in 1 : length(ttsByOperon[[i]])){
      
      cur.delta = "NA"
      cur.read.cov = "NA"
      idx = match(ttsByOperon[[i]][j], term.site$TTS_ID)
      
      if(is.na(match(term.site$Strand[idx], "+"))==FALSE){
        cur.delta =  abs(term.site$TTS_pos[idx] - term.site$OperonTo[idx])
        delta.vec = c(delta.vec, cur.delta)
      }
      
      if(is.na(match(term.site$Strand[idx], "-"))==FALSE){
        cur.delta =  abs(term.site$OperonFrom[idx] - term.site$TTS_pos[idx])
        delta.vec = c(delta.vec, cur.delta)
      }
      
      # Site coverage over all samples
      cur.read.cov = term.site$All_L2RelCov[idx]
      read.cov.vec = c(read.cov.vec, cur.read.cov)
      
      # Site coverage over samples grown on fructose
      cur.read.cov = term.site$B_L2RelCov[idx]
      fruct.read.cov.vec = c(fruct.read.cov.vec, cur.read.cov)
      
      # Site coverage over samples grown on syngas
      cur.read.cov = term.site$BR_L2RelCov[idx]
      syngas.read.cov.vec = c(syngas.read.cov.vec, cur.read.cov)
      
      fruct.read.cov = c(term.site$B1_L2RelCov[idx], term.site$B3_L2RelCov[idx], 
                         term.site$B4_reads[idx])
      syngas.read.cov = c(term.site$BR1_L2RelCov[idx], term.site$BR3_L2RelCov[idx], 
                          term.site$BR4_reads[idx])
      
      fruct.median.vec = c(fruct.median.vec, median(fruct.read.cov))
      syngas.median.vec = c(syngas.median.vec, median(syngas.read.cov))
      
      fruct.mad.vec = c(fruct.mad.vec, median(abs(fruct.read.cov - median(fruct.read.cov))))
      syngas.mad.vec = c(syngas.mad.vec, median(abs(syngas.read.cov - median(syngas.read.cov))))
      
    }
    
    # Map coverage of T3Es per operon to Q(SiteCov)
    cur.readcov.quant.label = readcov.quant.label[match(cut(read.cov.vec,bin.read.cov), 
                                                        readcov.quant.end.point)]
    
    # Map coverage of T3Es per operon on fructose samples to Q_fructose(SiteCov)
    cur.fruct.readcov.quant.label = fruct.readcov.quant.label[match(cut(fruct.read.cov.vec,fruct.bin.read.cov), 
                                                                    fruct.readcov.quant.end.point)]
    
    # Map coverage of T3Es per operon on syngas samples to Q_syngas(SiteCov)
    cur.syngas.readcov.quant.label = syngas.readcov.quant.label[match(cut(syngas.read.cov.vec,syngas.bin.read.cov), 
                                                                      syngas.readcov.quant.end.point)]
    
    # Map the distance of T3Es per operon to Q(Dist)
    cur.delta.quant.label = delta.quant.label[match(cut(delta.vec,bin.delta), 
                                                    delta.quant.end.point)]
    
    # Select the T3Es belonging to the lowest CQ (fructose)
    cur.fruct.label = rbind.data.frame(cur.fruct.readcov.quant.label, cur.delta.quant.label)
    tmp = colSums(cur.fruct.label)
    fruct.tts.operon.dist[[names(ttsByOperon)[i]]] = delta.vec[is.na(match(tmp, min(tmp)))==FALSE]
    fruct.tts.operon.site[[names(ttsByOperon)[i]]] = ttsByOperon[[i]][is.na(match(tmp, min(tmp)))==FALSE]
    fruct.tts.operon.qu[[names(ttsByOperon)[i]]] = tmp[is.na(match(tmp, min(tmp)))==FALSE]
    
    # Select the T3Es belonging to the lowest CQ (syngas)
    cur.syngas.label = rbind.data.frame(cur.syngas.readcov.quant.label, cur.delta.quant.label)
    tmp = colSums(cur.syngas.label)
    syngas.tts.operon.dist[[names(ttsByOperon)[i]]] = delta.vec[is.na(match(tmp, min(tmp)))==FALSE]
    syngas.tts.operon.site[[names(ttsByOperon)[i]]] = ttsByOperon[[i]][is.na(match(tmp, min(tmp)))==FALSE]
    syngas.tts.operon.qu[[names(ttsByOperon)[i]]] = tmp[is.na(match(tmp, min(tmp)))==FALSE]
    
    # Select the T3Es belonging to the lowest CQ
    cur.label = rbind.data.frame(cur.readcov.quant.label, cur.delta.quant.label)
    tmp = colSums(cur.label)
    tts.operon.dist[[names(ttsByOperon)[i]]] = delta.vec[is.na(match(tmp, min(tmp)))==FALSE]
    tts.operon.site[[names(ttsByOperon)[i]]] = ttsByOperon[[i]][is.na(match(tmp, min(tmp)))==FALSE]
    tts.operon.qu[[names(ttsByOperon)[i]]] = tmp[is.na(match(tmp, min(tmp)))==FALSE]
    
  }
}

#
# Filter out the operons whose selected T3Es display CQs lower than the central CQ. 
#
# Q(SiteCov) is based on fructose samples
distCQ.fruct = list()
siteCQ.fruct = list()
for(i in 1 : length(fruct.tts.operon.dist)){
  
  #check if there is at least a TSS whose composite quantile fulfills the criterion
  tmp.qu = fruct.tts.operon.qu[[i]]
  lgc = array(FALSE, length(tmp.qu))
  lgc[tmp.qu <= thrs.quantile] = TRUE
  
  tmp.dist = fruct.tts.operon.dist[[i]]
  tmp.site = fruct.tts.operon.site[[i]]
  
  if( length(lgc[lgc==TRUE]) > 0 ){
    distCQ.fruct[[names(fruct.tts.operon.dist)[i]]] = tmp.dist[lgc==TRUE]
    siteCQ.fruct[[names(fruct.tts.operon.site)[i]]] = tmp.site[lgc==TRUE]
  }
  
}

# Q(SiteCov) is based on syngas samples
distCQ.syngas = list()
siteCQ.syngas = list()
for(i in 1 : length(syngas.tts.operon.dist)){
  
  #check if there is at least a TSS whose composite quantile fulfills the criterion
  tmp.qu = syngas.tts.operon.qu[[i]]
  lgc = array(FALSE, length(tmp.qu))
  lgc[tmp.qu <= thrs.quantile] = TRUE
  
  tmp.dist = syngas.tts.operon.dist[[i]]
  tmp.site = syngas.tts.operon.site[[i]]
  
  if( length(lgc[lgc==TRUE]) > 0 ){
    distCQ.syngas[[names(syngas.tts.operon.dist)[i]]] = tmp.dist[lgc==TRUE]
    siteCQ.syngas[[names(syngas.tts.operon.site)[i]]] = tmp.site[lgc==TRUE]
  }
  
}

# Q(SiteCov) is based on all samples
distCQ = list()
siteCQ = list()
for(i in 1 : length(tts.operon.dist)){
  
  #check if there is at least a TSS whose composite quantile fulfills the criterion
  tmp.qu = tts.operon.qu[[i]]
  lgc = array(FALSE, length(tmp.qu))
  lgc[tmp.qu <= thrs.quantile] = TRUE
  
  tmp.dist = tts.operon.dist[[i]]
  tmp.site = tts.operon.site[[i]]
  
  if( length(lgc[lgc==TRUE]) > 0 ){
    distCQ[[names(tts.operon.dist)[i]]] = tmp.dist[lgc==TRUE]
    siteCQ[[names(tts.operon.site)[i]]] = tmp.site[lgc==TRUE]
  }
  
}

#----------------------------------------------------------------------------------------
# Site filtering by read-coverage profile on the genomic region spanned by gene (operon).
#----------------------------------------------------------------------------------------

distCQ.fruct = distCQ.fruct[is.na(match(names(distCQ.fruct), op.unif.read.cov))==FALSE]
siteCQ.fruct = siteCQ.fruct[is.na(match(names(siteCQ.fruct), op.unif.read.cov))==FALSE]
  
distCQ.syngas = distCQ.syngas[is.na(match(names(distCQ.syngas), op.unif.read.cov))==FALSE]
siteCQ.syngas = siteCQ.syngas[is.na(match(names(siteCQ.syngas), op.unif.read.cov))==FALSE]

primary.tts = list(ops.fruct = names(siteCQ.fruct), 
                   ops.syngas = names(siteCQ.syngas),
                   primary.tts.site.fruct = siteCQ.fruct, 
                   primary.tts.dist.fruct = distCQ.fruct,
                   primary.tts.site.syngas = siteCQ.syngas, 
                   primary.tts.dist.syngas = distCQ.syngas
                   )


#-------------------------------------------------------------------------------
#
# Site clustering by physical proximity. From sites closer than 5 nt from each other, 
# the T3E closest to the stop codon of the upstream CDS is the primary site. 
#
#-------------------------------------------------------------------------------

#
# Compile the final set of sites per gene (operon) by considering representative sites per site clusters in fructose
#
final.tts.by.operon = list()
final.tts.short.dist.by.operon = c()
for(i in 1 : length(primary.tts$ops.fruct)){
  
  final.site.by.operon = c()
  
  site.by.operon = primary.tts$primary.tts.site.fruct[[i]]
  
  #
  # Consider the scenario where an operon has multiple TTSs
  #
  if(length(site.by.operon)>1){
    site.pos.vec = c()
    for(j in 1 : length(site.by.operon)){
      site.pos.vec[j] = tbseq.cap.term$term_site$TTS_pos[match(site.by.operon[j], 
                                                               tbseq.cap.term$term_site$TTS_ID)]
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
    final.tts.by.operon[[primary.tts$ops.fruct[i]]] = final.site.by.operon 
  }
  
  #
  # Consider the scenario where an operon has a single TTS
  #
  if(length(site.by.operon)==1){
    final.site.by.operon = site.by.operon
    final.tts.by.operon[[primary.tts$ops.fruct[i]]] = final.site.by.operon
  }
  
  #
  #compute the shortest distance of the final sites per operon
  #
  final.dist.pos.vec = c()
  for(j in 1 : length(final.site.by.operon)){
    
    idx = match(final.site.by.operon[j], siteCQ.fruct[[primary.tts$ops.fruct[i]]])
    final.dist.pos.vec[j] = distCQ.fruct[[primary.tts$ops.fruct[i]]][idx]
    
  }
  if(is.null(final.dist.pos.vec)==FALSE){ # parte nuova !!!!!!!
    final.tts.short.dist.by.operon[primary.tts$ops.fruct[i]] = min(final.dist.pos.vec)
  }
  if(is.null(final.dist.pos.vec)==TRUE){
    final.tts.short.dist.by.operon[primary.tts$ops.fruct[i]] = NA
  }
  
}
final.tts.by.operon.fruct = final.tts.by.operon
final.tts.short.dist.by.operon.fruct = final.tts.short.dist.by.operon


#
# Compile the final set of sites per gene (operon) by considering representative sites per site clusters in syngas
#
final.tts.by.operon = list()
final.tts.short.dist.by.operon = c()
for(i in 1 : length(primary.tts$ops.syngas)){
  
  final.site.by.operon = c()
  
  site.by.operon = primary.tts$primary.tts.site.syngas[[i]]
  #
  # Consider the scenario where an operon has multiple TTSs
  #
  if(length(site.by.operon)>1){
    site.pos.vec = c()
    for(j in 1 : length(site.by.operon)){
      site.pos.vec[j] = tbseq.cap.term$term_site$TTS_pos[match(site.by.operon[j], 
                                                               tbseq.cap.term$term_site$TTS_ID)]
    }
    
    inter.dist = c()
    for(j in 2 : length(site.pos.vec)){
      inter.dist = c( inter.dist, abs(site.pos.vec[j] - site.pos.vec[j-1]) ) 
    }
    
    #identify the indices of the sites that are distant by more than 5 nt
    retained.idx.by.operon = c()
    #
    # Consder the case where an operon displays more than two sites
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
        #consider the case where all the pairs are distant by less than 5 nt
        if( length(inter.dist[inter.dist<=5]) == length(inter.dist) ){
          retained.idx.by.operon = c(retained.idx.by.operon, length(site.pos.vec))
        }
      }
    }
    #
    # Consder the case where an operon displays two sites
    #
    if(length(inter.dist)==1){
      if(inter.dist>5){
        retained.idx.by.operon = c(retained.idx.by.operon, 1:2)
      }
      if(inter.dist<=5){
        retained.idx.by.operon = c(retained.idx.by.operon, 2)
      }
    }
    #collect sites corresponding to distinct positions
    final.site.by.operon = site.by.operon[unique(retained.idx.by.operon)]
    final.tts.by.operon[[primary.tts$ops.syngas[i]]] = final.site.by.operon 
  }
  
  #
  # Consider the scenario where an operon displays a single TTS
  #
  if(length(site.by.operon)==1){
    final.site.by.operon = site.by.operon
    final.tts.by.operon[[primary.tts$ops.syngas[i]]] = final.site.by.operon
  }
  
  #
  #compute the shortest distance of the final sites per operon
  #
  final.dist.pos.vec = c()
  for(j in 1 : length(final.site.by.operon)){
    
    idx = match(final.site.by.operon[j], siteCQ.syngas[[primary.tts$ops.syngas[i]]])
    final.dist.pos.vec[j] = distCQ.syngas[[primary.tts$ops.syngas[i]]][idx]
    
  }
  if(is.null(final.dist.pos.vec)==FALSE){ 
    final.tts.short.dist.by.operon[primary.tts$ops.syngas[i]] = min(final.dist.pos.vec)
  }
  if(is.null(final.dist.pos.vec)==TRUE){
    final.tts.short.dist.by.operon[primary.tts$ops.syngas[i]] = NA
  }

}
final.tts.by.operon.syngas = final.tts.by.operon
final.tts.short.dist.by.operon.syngas = final.tts.short.dist.by.operon



#-------------------------------------------------------------------------------
#
# Plot the results of T3Es filtering
#
#-------------------------------------------------------------------------------

# Get the genes per operon
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


jpeg("../results/TIS-T3E/T3E.jpeg", res=600, width=4500, height=3000)
par(mfrow=c(1,2))
par(mar=c(7.5,6,1,1))

#
# Display the 3' UTR length in syngas and fructose
#

# Fructose
final.tts.short.dist.by.gene.fruct = c()
for( i in 1 : length(final.tts.short.dist.by.operon.fruct)){
  tmp = gene.set.by.operon[match(names(final.tts.short.dist.by.operon.fruct)[i], names(gene.set.by.operon))][[1]]
  for(j in 1 : length(tmp)){
    final.tts.short.dist.by.gene.fruct[tmp[j]] = as.numeric(final.tts.short.dist.by.operon.fruct[i])
  }
}
final.tts.short.dist.by.gene.fruct = final.tts.short.dist.by.gene.fruct[is.na(final.tts.short.dist.by.gene.fruct)==FALSE]
counts = c()
bins = seq(0,2000,by=100)
for(i in 1 : (length(bins)-1)){
  counts[i] = length(final.tts.short.dist.by.gene.fruct[final.tts.short.dist.by.gene.fruct > bins[i] & 
                                                          final.tts.short.dist.by.gene.fruct < bins[i+1]])
}
counts.fruct = counts


# Syngas
final.tts.short.dist.by.gene.syngas = c()
for( i in 1 : length(final.tts.short.dist.by.operon.syngas)){
  tmp = gene.set.by.operon[match(names(final.tts.short.dist.by.operon.syngas)[i], names(gene.set.by.operon))][[1]]
  for(j in 1 : length(tmp)){
    final.tts.short.dist.by.gene.syngas[tmp[j]] = as.numeric(final.tts.short.dist.by.operon.syngas[i])
  }
}
final.tts.short.dist.by.gene.syngas = final.tts.short.dist.by.gene.syngas[is.na(final.tts.short.dist.by.gene.syngas)==FALSE]
counts = c()
bins = seq(0,2000,by=100)
for(i in 1 : (length(bins)-1)){
  counts[i] = length(final.tts.short.dist.by.gene.syngas[final.tts.short.dist.by.gene.syngas > bins[i] & 
                                                           final.tts.short.dist.by.gene.syngas < bins[i+1]])
}
counts.syngas = counts

mpl = t(as.matrix(cbind.data.frame(counts.fruct, counts.syngas)))

no.bins.plot = 7
b = barplot(mpl[,1:no.bins.plot], plot=FALSE)
barplot(mpl[,1:no.bins.plot], ylab="Frequency", beside=TRUE, col=c("#e08214", "#8073ac"),las=2, 
        cex.axis=1.25, cex.names=1.65)
legend("topright", legend=c("Fructose","Syngas"), fill=c("#e08214", "#8073ac"), bty="n", 
       border="white", ncol=1, cex=1.15)

binlabels = c()
for(i in 1 : (length(bins)-1)){
  binlabels[i] = paste0("[", bins[i], ",", bins[i+1], "]")
}
axis(1, at=seq(2,length(b)*3,by=3), labels = binlabels[1:no.bins.plot], las=2, cex=.85)
mtext(side=1, line=6, text="T3E distance from stop codon")


#
# Display the number of genes by the number of selected T3E sites per gene (fructose and syngas)
#
final.tts.by.gene.syngas = list()
for(i in  1 : length(final.tts.by.operon.syngas)){
  tmp = gene.set.by.operon[match(names(final.tts.by.operon.syngas)[i], names(gene.set.by.operon))][[1]]
  for(j in 1 : length(tmp)){
    final.tts.by.gene.syngas[[tmp[j]]] = final.tts.by.operon.syngas[[i]]
  }
}

final.tts.by.gene.fruct = list()
for(i in  1 : length(final.tts.by.operon.fruct)){
  tmp = gene.set.by.operon[match(names(final.tts.by.operon.fruct)[i], names(gene.set.by.operon))][[1]]
  for(j in 1 : length(tmp)){
    final.tts.by.gene.fruct[[tmp[j]]] = final.tts.by.operon.fruct[[i]]
  }
}

final.gene = union(names(final.tts.by.gene.syngas), names(final.tts.by.gene.fruct))
final.tts.by.gene = list()
for(i in 1 : length(final.gene)){
  final.tts.by.gene[[final.gene[i]]] = unique(c(final.tts.by.gene.syngas[[final.gene[i]]], 
                                                final.tts.by.gene.fruct[[final.gene[i]]]))
  if(length(final.tts.by.gene[[final.gene[i]]])==5){show(final.gene[i])}
}

gene.freq.by.tts.no.fruct = table(unlist(lapply(final.tts.by.gene.fruct, function(x) length(x))))
gene.freq.by.tts.no.syngas = table(unlist(lapply(final.tts.by.gene.syngas, function(x) length(x))))
gene.freq.by.tts.no = table(unlist(lapply(final.tts.by.gene, function(x) length(x))))
tts.no = unique(c(names(gene.freq.by.tts.no.fruct), names(gene.freq.by.tts.no.syngas), names(gene.freq.by.tts.no)))

m = matrix(0, nrow=length(tts.no), ncol=3, dimnames=list(tts.no, c("Fructose","Syngas","Either subs.")))
for(i in 1 : length(tts.no)){
  m[i,1] = as.numeric(gene.freq.by.tts.no.fruct[match(tts.no[i], names(gene.freq.by.tts.no.fruct))])
  m[i,2] = as.numeric(gene.freq.by.tts.no.syngas[match(tts.no[i], names(gene.freq.by.tts.no.syngas))])
  m[i,3] = as.numeric(gene.freq.by.tts.no[match(tts.no[i], names(gene.freq.by.tts.no))])
}

barplot(t(m[,1:2]),beside=TRUE, col=c("#e08214", "#8073ac"), ylab="Frequency", las=1, cex.axis=1.25, cex.names=1.45)
legend("topright", legend=c("Fructose","Syngas"), fill=c("#e08214", "#8073ac"), bty="n", 
       border="white", ncol=1, cex=1.15)
mtext(side=1, line=3, text="T3E no. by gene")

dev.off()


#-------------------------------------------------------------------
# Print the final T3Es
#-------------------------------------------------------------------

# Load gene (operon) annotations
foperon = read.xlsx("../data/Operon.description.gene.annot.xlsx", sheetIndex=1)
# Obtain gene sets per operon
gene.set.by.operon = list()
for(i in 1 : dim(foperon)[1]){
  LABRINI = foperon$GenesInOperon_CP110420[i]
  gene.set.by.operon[[paste0("Operon-",f$OperonID[i])]] = strsplit(LABRINI, ",")[[1]]
}


# Collect the final T3Es for each operon on either growth condition 
ops = union(names(final.tts.by.operon.syngas), names(final.tts.by.operon.fruct))

ttsByOperon = list()
ttsByGene = list()

for(i in 1 : length(ops)){
  cur.tts = c()
  
  idx = match(ops[i], names(final.tts.by.operon.fruct))
  cur.tts = c(cur.tts, unlist(final.tts.by.operon.fruct[[idx]]))
  
  idx = match(ops[i], names(final.tts.by.operon.syngas))
  cur.tts = c(cur.tts, unlist(final.tts.by.operon.syngas[[idx]]))
  
  # Assign T3Es to operon
  ttsByOperon[[ops[i]]] = unique(cur.tts)
  
  # Assign T3Es to genes
  cur.gene.by.operon = gene.set.by.operon[[ops[i]]]
  for( j in 1 : length(cur.gene.by.operon) ){
    ttsByGene[[cur.gene.by.operon[j]]] = unique(cur.tts)
  }
}

# Print T3Es by operon
out.name = names(ttsByOperon)
out.signal = as.vector(
  unlist(lapply(ttsByOperon, function(x) paste(x, collapse=",")))
)
out = cbind(out.name, out.signal)
colnames(out) = c("Operon","T3E")
write.xlsx(out, file = "../results/TIS-T3E/T3E.xlsx", row.names=FALSE, sheetName = "T3E_by_operon",append=FALSE)

# Print T3E by gene
out.name = names(ttsByGene)
out.signal = as.vector(
  unlist(lapply(ttsByGene, function(x) paste(x, collapse=",")))
)
out = cbind(out.name, out.signal)
colnames(out) = c("Gene","T3E")
write.xlsx(out, file = "../results/TIS-T3E/T3E.xlsx", row.names=FALSE, sheetName = "T3E_by_gene",append=TRUE)
