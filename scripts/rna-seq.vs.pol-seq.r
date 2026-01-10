
require(xlsx)
library(dplyr)
require(ggvenn)
require(ggpubr)
require(colorRamp2)


#-----------------------------------------------------------------------
# Load gene annotations to KO pathways according to Level 3 
# in Table S3 of Valgepea et al. (10.1128/msystems.00026-22) 
#-----------------------------------------------------------------------
fid = read.xlsx("../data/Kegg Orthology pathway annotations.xlsx", sheetName = "Kaspar curated")

#--------------------------------------------------------------
#Retrieve the NCBI annotations of Clostridium autoethanogenum
#--------------------------------------------------------------
ncbif = read.delim("../data/GCA_040166795.1/genomic.gff", header=FALSE, comment.char = "#", sep="\t")
annotf = ncbif[ncbif[,3] %in% "CDS",]
# Locus Tags
locusTags = as.vector(sapply(annotf[,9], function(x) strsplit(x,";")[[1]][2]))
locusTags = gsub("Parent=gene-", "", locusTags)
locusTags = gsub("Name=", "", locusTags)


#-------------------------------------------------------------
# Print the table of differentially transcribed (DTR) genes
#-------------------------------------------------------------

# Load edgeR results
load("../results/DTR/totrna.de.rda")

# Remove genes external to C. autoethanogenum
totrna.de = totrna.de[grepl("ERCC",rownames(totrna.de))==FALSE,]
# Remove genes that are not protein-coding 
totrna.de = totrna.de[is.na(match(rownames(totrna.de),locusTags))==FALSE,]

# Sort data according to P-value (from low to high)
totrna.ord = totrna.de[order(totrna.de$PValue),]

# Identify genes whose ranks lie in the top 5% positions 
top_n = round( 0.05 * dim(totrna.ord)[1] )
top_n_totdeg = rownames(totrna.ord)[1:top_n]
lgcTop = array(FALSE, dim(totrna.ord)[1])
lgcTop[is.na(match(rownames(totrna.ord),top_n_totdeg))==FALSE] = TRUE

# Identify genes whose q-value is lower than 0.05
qvalues = p.adjust(totrna.ord$PValue, method="BH") 
lgcQvalue = array(FALSE, dim(totrna.ord)[1])
lgcQvalue[qvalues < 0.05] = TRUE
qvalue_deg =  rownames(totrna.ord[lgcQvalue==TRUE,])

# Identify differentially transcribed (DTR) genes
dtr_genes = rownames(totrna.ord[lgcQvalue == TRUE & lgcTop == TRUE,])
setdiff(top_n_totdeg, dtr_genes)

# Identify down-/up-regulated genes
tmp = totrna.ord[1:top_n,]
tr.top.n.dn = rownames(tmp[tmp$logFC<0,])
tr.top.n.up = rownames(tmp[tmp$logFC>0,])

# Retrieve gene annotations to KO pathways
tmp = totrna.ord[1:top_n,]
gene_pathway = fid[match(rownames(totrna.ord),fid[,4]),3]
gene_pathway[is.na(gene_pathway)==TRUE] = "n/a"

# Retrieve gene names of common use
common_name = fid[match(rownames(totrna.ord),fid[,4]),7]
common_name[grepl("^RS[0-9]*", common_name)==TRUE] = "n/a"
common_name[grepl("^E[0-9].*\\.*", common_name)==TRUE] = "n/a"
common_name[is.na(common_name)==TRUE] = "n/a"

# Retrieve gene annotation according to NCBI 
protprod_LT_LAbrini_ncbi = c()
for(i in 1 : length(rownames(totrna.ord))){
  tmp.field = annotf[match(rownames(totrna.ord)[i],locusTags),9]
  tmp.field.vec = strsplit(tmp.field,";")[[1]]
  if(length(tmp.field.vec[grepl("product=",tmp.field.vec)==TRUE])==1){
    protprod_LT_LAbrini_ncbi[i] = gsub("product=","",tmp.field.vec[grepl("product=",tmp.field.vec)==TRUE])
  }
}

# Assemble and print the table of differentially transcribed genes
df = cbind.data.frame(rownames(totrna.ord), protprod_LT_LAbrini_ncbi, 
                      common_name, gene_pathway, totrna.ord$logFC, qvalues, lgcTop
                      )
colnames(df) = c("GeneID","Description","Name","KO functional category L3","log2FC","q-value","DTR")
rownames(df) = rownames(totrna.ord)
write.xlsx(df, file = "../results/DTR/dtr.xlsx", row.names = FALSE, col.names = TRUE)


#---------------------------------------------------------------------
# Volcano plot showing genes determined as differentially 
# transcribed as a subset of genes with statistically 
# different expression between syngas and fructose chemostats. 
#---------------------------------------------------------------------

x = totrna.ord$logFC
qvalues = p.adjust(totrna.ord$PValue, method="BH") 
y = -log10(qvalues)

# Plot colors
plotcolors = array("#deebf7", length(x))

# Identify genes whose q-value is lower than 0.05 
lgcQvalue = array(FALSE, dim(totrna.ord)[1])
lgcQvalue[qvalues < 0.05] = TRUE
plotcolors[lgcQvalue==TRUE] = "#92c5de"

# Identify ranked genes that lie in the top 5% positions
lgcTop = array(FALSE, dim(totrna.ord)[1])
lgcTop[is.na(match(rownames(totrna.ord),top_n_totdeg))==FALSE] = TRUE
plotcolors[lgcTop==TRUE] = "#08306b"

png("../results/DTR/volcano.plot.png", res=600, width=3000, height=3000)
plot(x,y, col=plotcolors, pch=19, xlab="", ylab="", cex.axis=1.25, las=2, axes=FALSE)
axis(1, at = seq(-10,6,by=2), labels = seq(-10,6,by=2), cex.axis=1.25)
mtext(side=1, text = expression(log[2](FC)), line=2.5, cex=1.25)
axis(2, at = seq(0, 12, by=2), labels = seq(0, 12, by=2), cex.axis=1.25, las=2)
mtext(side=2, text = expression(-log[10](q-value)), line=2.75, cex=1.25)
legend("topright", bty="n", border=NA, cex=.85,legend=c("NR","q-value < 0.05","DTR"), fill=c("#deebf7","#92c5de","#08306b"))
dev.off()



#--------------------------------------------------------------
# Load gene ranks resulting from the three metrics used to 
# determine differentially translated (DTL) genes
#--------------------------------------------------------------

f.abs = read.delim("../results/DTL/pol.rnk/abs_translation_diff_results.csv", sep=",")
f.boot = read.delim("../results/DTL/pol.rnk/bootstrap_results.csv", sep=",")
f.wa = read.delim("../results/DTL/pol.rnk/weighted_aitchison_results.csv", sep=",")

# Remove genes external to C. autoethanogenum
f.abs = f.abs[grepl("ERCC",f.abs$X)==FALSE,]
f.boot = f.boot[grepl("ERCC",f.boot$gene)==FALSE,]
f.wa = f.wa[grepl("ERCC",f.wa$gene)==FALSE,]
# Remove genes that are not protein-coding 
f.abs = f.abs[is.na(match(f.abs$X,locusTags))==FALSE,]
f.boot = f.boot[is.na(match(f.boot$gene,locusTags))==FALSE,]
f.wa = f.wa[is.na(match(f.wa$gene,locusTags))==FALSE,]

genes = unique(c(f.abs[,1], f.boot[,2], f.wa[,2]))


# Obtain normalized ranks by the TL-STATUS approach
rnk.abs = c()
round.vec = round(f.abs$abs_diff, digits=4)
unique.round.vec = unique(round.vec)
sort.unique.round.vec = unique.round.vec[rev(order(unique.round.vec))]
for(i in 1 : length(genes)){ 
  idx = match(genes[i], f.abs$X)
  rnk.abs[i] = match(round(as.numeric(f.abs$abs_diff[idx]), digits=4), sort.unique.round.vec)
}
rnk.abs.norm = 100 * rnk.abs/max(rnk.abs)
names(rnk.abs.norm) = genes

# Obtain normalized ranks by the POLY-FRACTIONS approach
rnk.boot = c()
round.vec = c()
for(i in 1 : length(f.boot$gene)){
  round.vec[i] = paste(f.boot$num_sig_phases[i], round(f.boot$max_diff[i], digits=4), sep=";")
}
unique.round.vec = unique(round.vec)
sort.unique.round.vec = unique.round.vec[rev(order(unique.round.vec))]
for(i in 1 : length(genes)){
  idx = match(genes[i], f.boot$gene)
  tmp = paste(f.boot$num_sig_phases[idx], round(f.boot$max_diff[idx], digits=4), sep=";")
  rnk.boot[i] = match(tmp, sort.unique.round.vec)
}
rnk.boot.norm = 100 * rnk.boot/max(rnk.boot)
names(rnk.boot.norm) = genes

# Obtain normalized ranks by the TL-PROFILE approach
rnk.wa = c()
round.vec = c()
for(i in 1 : length(f.wa$gene)){
  round.vec[i] = round(f.wa$weighted_distance[i], digits=4)
}
unique.round.vec = unique(round.vec)
sort.unique.round.vec = unique.round.vec[rev(order(unique.round.vec))]
for(i in 1 : length(genes)){
  idx = match(genes[i], f.wa$gene)
  rnk.wa[i] = match(round(as.numeric(f.wa$weighted_distance[idx]), digits=4), sort.unique.round.vec)
}
rnk.wa.norm = 100 * rnk.wa/max(rnk.wa)
names(rnk.wa.norm) = genes


# Compute the average rank for each gene
rnk.mat = as.matrix(cbind.data.frame(rnk.abs.norm, rnk.boot.norm, rnk.wa.norm))
rownames(rnk.mat) = genes
# Compute the arithmetic mean
rnk.mean = apply(rnk.mat,1,mean) 
# Sort genes by the computed average rank
rnk.mean.sorted = rnk.mean[order(as.numeric(rnk.mean))]
# Save ranks
pol.rnk = rnk.mean.sorted


#------------------------------------------------
# Define the top-ranking genes by each method
#------------------------------------------------

# POLY-FRACTIONS approach
boot.round.un = unique(round(rnk.boot.norm, digits=4))
boot.round.un.ord = boot.round.un[order(boot.round.un)]
idx = round( 0.05*length(boot.round.un.ord) )
rnk.th = boot.round.un.ord[idx]
top_n = length(rnk.boot.norm[rnk.boot.norm<rnk.th])
top.n.boot = names(rnk.boot.norm[rnk.boot.norm<rnk.th])

# TL-STATUS approach
abs.round.un = unique(round(rnk.abs.norm, digits=4))
abs.round.un.ord = abs.round.un[order(abs.round.un)]
idx = round( 0.05*length(abs.round.un.ord) )
rnk.th = abs.round.un.ord[idx]
top_n = length(rnk.abs.norm[rnk.abs.norm<rnk.th])
top.n.abs = names(rnk.abs.norm[rnk.abs.norm<rnk.th] )

# TL-PROFILE approach
wa.round.un = unique(round(rnk.wa.norm, digits=4))
wa.round.un.ord = wa.round.un[order(wa.round.un)]
idx = round( 0.05*length(wa.round.un.ord) )
rnk.th = wa.round.un.ord[idx]
top_n = length(rnk.wa.norm[rnk.wa.norm<rnk.th])
top.n.wa = names(rnk.wa.norm[rnk.wa.norm<rnk.th])

# CONSENSUS approach
rnk.mean.sorted.round.unique = unique(round(rnk.mean.sorted, digits=4))
idx = round( 0.05*length(rnk.mean.sorted.round.unique) )
rnk.th = rnk.mean.sorted.round.unique[idx]
top_n = length(rnk.mean.sorted[rnk.mean.sorted<rnk.th])
top.n.avg = names(rnk.mean.sorted[rnk.mean.sorted<rnk.th])


#--------------------------------------------------------------------------
# Load polysome profiles (RNA proportions by fraction and replicate) 
#--------------------------------------------------------------------------

ffru = read.delim("../data/Polysome_seq/fructose_data.txt",sep="\t")
ffru = ffru[grepl("^ERCC-", ffru[,1])==FALSE,]
geneId = ffru[,1]
ffru = ffru[,-1]
rownames(ffru) = geneId

fgas = read.delim("../data/Polysome_seq/gas_data.txt", sep="\t")
fgas = fgas[grepl("^ERCC-", fgas[,1])==FALSE,]
geneId = fgas[,1]
fgas = fgas[,-1]
rownames(fgas) = geneId

png("../results/DTL/pol-seq_reproducibility.png", res=600, width=5900, height=3500)
par(mfrow=c(2,3),mar=c(5,5,2.5,2.25))

tmp = fgas[,grepl("_1",colnames(fgas))==TRUE]
clr.tmp = matrix(NA, nrow=nrow(tmp), ncol=ncol(tmp), dimnames=list(rownames(tmp), colnames(tmp)))
for(i in 1 : dim(tmp)[1]){
  clr.tmp[i,]=as.numeric(log(tmp[i,])-0.25*sum(log(tmp[i,])))
}
clr_gas_rep_1 = clr.tmp

tmp = fgas[,grepl("_2",colnames(fgas))==TRUE]
clr.tmp = matrix(NA, nrow=nrow(tmp), ncol=ncol(tmp), dimnames=list(rownames(tmp), colnames(tmp)))
for(i in 1 : dim(tmp)[1]){
  clr.tmp[i,]=as.numeric(log(tmp[i,])-0.25*sum(log(tmp[i,])))
}
clr_gas_rep_2 = clr.tmp

tmp = fgas[,grepl("_3",colnames(fgas))==TRUE]
clr.tmp = matrix(NA, nrow=nrow(tmp), ncol=ncol(tmp), dimnames=list(rownames(tmp), colnames(tmp)))
for(i in 1 : dim(tmp)[1]){
  clr.tmp[i,]=as.numeric(log(tmp[i,])-0.25*sum(log(tmp[i,])))
}
clr_gas_rep_3 = clr.tmp

cor.test.rep1.vs.rep2 = as.numeric(cor.test(clr_gas_rep_1, clr_gas_rep_2, method="spearman")$estimate)
cor.test.rep1.vs.rep3 = as.numeric(cor.test(clr_gas_rep_1, clr_gas_rep_3, method="spearman")$estimate)
cor.test.rep2.vs.rep3 = as.numeric(cor.test(clr_gas_rep_2, clr_gas_rep_3, method="spearman")$estimate)

plot(clr_gas_rep_1, clr_gas_rep_2, xlab=expression(clr-transformed~mRNA~counts^(SYNGAS_REP_1)), 
     ylab=expression(clr-transformed~mRNA~counts^(SYNGAS_REP_2)), 
     main="", pch=19, col="#b2abd2", cex=.75, ylim=c(-4.5,3))
mtext(side=3,line=1.5, text="a", adj=0, font=2, cex=.75)
text(x=-6.75, y=2.5, cex=.85, labels=expression(paste('Spearman\'s ',rho,' = ')))
text(x=-4, y=2.5, cex=.85, labels = signif(cor.test.rep1.vs.rep2, digits=2))

plot(clr_gas_rep_1, clr_gas_rep_3, xlab=expression(clr-transformed~mRNA~counts^(SYNGAS_REP_1)), 
     ylab=expression(clr-transformed~mRNA~counts^(SYNGAS_REP_3)), 
     main="", pch=19, col="#b2abd2", cex=.75, ylim=c(-4.5,3))
mtext(side=3,line=1.5, text="b", adj=0, font=2, cex=.75)
text(x=-6.75, y=2.5, cex=.85, labels=expression(paste('Spearman\'s ',rho,' = ')))
text(x=-4, y=2.5, cex=.85, labels = signif(cor.test.rep1.vs.rep3, digits=2))

plot(clr_gas_rep_2, clr_gas_rep_3, xlab=expression(clr-transformed~mRNA~counts^(SYNGAS_REP_2)), 
     ylab=expression(clr-transformed~mRNA~counts^(SYNGAS_REP_3)), 
     main="", pch=19, col="#b2abd2", cex=.75, ylim=c(-4.5,3))
mtext(side=3,line=1.5, text="c", adj=0, font=2, cex=.75)
text(x=-2.9, y=2.5, cex=.85, labels=expression(paste('Spearman\'s ',rho,' = ')))
text(x=-1.65, y=2.5, cex=.85, labels = signif(cor.test.rep2.vs.rep3, digits=2))


tmp = ffru[,grepl("_1",colnames(ffru))==TRUE]
clr.tmp = matrix(NA, nrow=nrow(tmp), ncol=ncol(tmp), dimnames=list(rownames(tmp), colnames(tmp)))
for(i in 1 : dim(tmp)[1]){
  clr.tmp[i,]=as.numeric(log(tmp[i,])-0.25*sum(log(tmp[i,])))
}
clr_fru_rep_1 = clr.tmp

tmp = ffru[,grepl("_2",colnames(ffru))==TRUE]
clr.tmp = matrix(NA, nrow=nrow(tmp), ncol=ncol(tmp), dimnames=list(rownames(tmp), colnames(tmp)))
for(i in 1 : dim(tmp)[1]){
  clr.tmp[i,]=as.numeric(log(tmp[i,])-0.25*sum(log(tmp[i,])))
}
clr_fru_rep_2 = clr.tmp

tmp = ffru[,grepl("_3",colnames(ffru))==TRUE]
clr.tmp = matrix(NA, nrow=nrow(tmp), ncol=ncol(tmp), dimnames=list(rownames(tmp), colnames(tmp)))
for(i in 1 : dim(tmp)[1]){
  clr.tmp[i,]=as.numeric(log(tmp[i,])-0.25*sum(log(tmp[i,])))
}
clr_fru_rep_3 = clr.tmp

cor.test.rep1.vs.rep2 = as.numeric(cor.test(clr_fru_rep_1, clr_fru_rep_2, method="spearman")$estimate)
cor.test.rep1.vs.rep3 = as.numeric(cor.test(clr_fru_rep_1, clr_fru_rep_3, method="spearman")$estimate)
cor.test.rep2.vs.rep3 = as.numeric(cor.test(clr_fru_rep_2, clr_fru_rep_3, method="spearman")$estimate)

plot(clr_fru_rep_1, clr_fru_rep_2, xlab=expression(clr-transformed~mRNA~counts^(FRUCTOSE_REP_1)), 
     ylab=expression(clr-transformed~mRNA~counts^(FRUCTOSE_REP_2)), 
     main="", pch=19, col="#fdb863", cex=.75, xlim=c(-3.5,3), ylim=c(-4,3))
mtext(side=3,line=1.5, text="d", adj=0, font=2, cex=.75)
text(x=-2.65, y=2.5, cex=.85, labels=expression(paste('Spearman\'s ',rho,' = ')))
text(x=-1.25, y=2.5, cex=.85, labels = signif(cor.test.rep1.vs.rep2, digits=2))

plot(clr_fru_rep_1, clr_fru_rep_3, xlab=expression(clr-transformed~mRNA~counts^(FRUCTOSE_REP_1)), 
     ylab=expression(clr-transformed~mRNA~counts^(FRUCTOSE_REP_3)), 
     main="", pch=19, col="#fdb863", cex=.75, xlim=c(-3.5,3), ylim=c(-4,3))
mtext(side=3,line=1.5, text="e", adj=0, font=2, cex=.75)
text(x=-2.65, y=2.5, cex=.85, labels=expression(paste('Spearman\'s ',rho,' = ')))
text(x=-1.25, y=2.5, cex=.85, labels = signif(cor.test.rep1.vs.rep3, digits=2))

plot(clr_fru_rep_2, clr_fru_rep_3, xlab=expression(clr-transformed~mRNA~counts^(FRUCTOSE_REP_2)), 
     ylab=expression(clr-transformed~mRNA~counts^(FRUCTOSE_REP_3)), 
     main="", pch=19, col="#fdb863", cex=.75, xlim=c(-3.5,3), ylim=c(-4,3))
mtext(side=3,line=1.5, text="f", adj=0, font=2, cex=.75)
text(x=-2.65, y=2.5, cex=.85, labels=expression(paste('Spearman\'s ',rho,' = ')))
text(x=-1.25, y=2.5, cex=.85, labels = signif(cor.test.rep2.vs.rep3, digits=2))

dev.off()


#--------------------------------------------------------------------------
# Retrieve mRNA proportions for ranked genes
#--------------------------------------------------------------------------
rna.prop.mean.gas = list()
rna.prop.mean.fru = list()
rna.prop.sd.gas = list()
rna.prop.sd.fru = list()
for(i in 1 : length(pol.rnk)){
  
  curGene = names(pol.rnk)[i]
  
  if(is.na(match(curGene, rownames(ffru)))==FALSE & is.na(match(curGene, rownames(fgas)))==FALSE){
    row.idx.fru = match(curGene, rownames(ffru))
    row.idx.gas = match(curGene, rownames(fgas))
    
    gas.A = mean( as.numeric( fgas[row.idx.gas,is.na(match(colnames(fgas), c("A_2", "A_3")))==FALSE] ) ) 
    gas.B = mean( as.numeric( fgas[row.idx.gas,is.na(match(colnames(fgas), c("B_2", "B_3")))==FALSE] ) )  
    gas.C = mean( as.numeric( fgas[row.idx.gas,is.na(match(colnames(fgas), c("C_2", "C_3")))==FALSE] ) )
    gas.D = mean( as.numeric( fgas[row.idx.gas,is.na(match(colnames(fgas), c("D_2", "D_3")))==FALSE] ) )
    
    gas.A.sd = sd( as.numeric( fgas[row.idx.gas,is.na(match(colnames(fgas), c("A_2", "A_3")))==FALSE] ) ) 
    gas.B.sd = sd( as.numeric( fgas[row.idx.gas,is.na(match(colnames(fgas), c("B_2", "B_3")))==FALSE] ) )  
    gas.C.sd = sd( as.numeric( fgas[row.idx.gas,is.na(match(colnames(fgas), c("C_2", "C_3")))==FALSE] ) )
    gas.D.sd = sd( as.numeric( fgas[row.idx.gas,is.na(match(colnames(fgas), c("D_2", "D_3")))==FALSE] ) )
    
    fru.A = mean( as.numeric( ffru[row.idx.fru,is.na(match(colnames(ffru), c("A_1", "A_3")))==FALSE] ) )
    fru.B = mean( as.numeric( ffru[row.idx.fru,is.na(match(colnames(ffru), c("B_1", "B_3")))==FALSE] ) )
    fru.C = mean( as.numeric( ffru[row.idx.fru,is.na(match(colnames(ffru), c("C_1", "C_3")))==FALSE] ) )
    fru.D = mean( as.numeric( ffru[row.idx.fru,is.na(match(colnames(ffru), c("D_1", "D_3")))==FALSE] ) )
    
    fru.A.sd = sd( as.numeric( ffru[row.idx.fru,is.na(match(colnames(ffru), c("A_1", "A_3")))==FALSE] ) )
    fru.B.sd = sd( as.numeric( ffru[row.idx.fru,is.na(match(colnames(ffru), c("B_1", "B_3")))==FALSE] ) )
    fru.C.sd = sd( as.numeric( ffru[row.idx.fru,is.na(match(colnames(ffru), c("C_1", "C_3")))==FALSE] ) )
    fru.D.sd = sd( as.numeric( ffru[row.idx.fru,is.na(match(colnames(ffru), c("D_1", "D_3")))==FALSE] ) )
    
    rna.prop.mean.gas[[curGene]] = c(gas.A, gas.B, gas.C, gas.D)
    rna.prop.mean.fru[[curGene]] = c(fru.A, fru.B, fru.C, fru.D)
    
    rna.prop.sd.gas[[curGene]] = c(gas.A.sd, gas.B.sd, gas.C.sd, gas.D.sd)
    rna.prop.sd.fru[[curGene]] = c(fru.A.sd, fru.B.sd, fru.C.sd, fru.D.sd)
  }
}
rna.prop.mean.gas.mat = matrix(unlist(rna.prop.mean.gas), ncol = 4, byrow = TRUE)
rna.prop.mean.fru.mat = matrix(unlist(rna.prop.mean.fru), ncol = 4, byrow = TRUE)
rownames(rna.prop.mean.gas.mat) = names(rna.prop.mean.gas)
rownames(rna.prop.mean.fru.mat) = names(rna.prop.mean.fru)  

rna.prop.sd.gas.mat = matrix(unlist(rna.prop.sd.gas), ncol = 4, byrow = TRUE)
rna.prop.sd.fru.mat = matrix(unlist(rna.prop.sd.fru), ncol = 4, byrow = TRUE)
rownames(rna.prop.sd.gas.mat) = names(rna.prop.sd.gas)
rownames(rna.prop.sd.fru.mat) = names(rna.prop.sd.fru)  


#-------------------------------------------------------------------------------
#
# Print the results on differentially translated genes by the different 
# methods for manuscript supplementary information
#
#-------------------------------------------------------------------------------

genes = unique(c(f.abs[,1], f.boot[,2], f.wa[,2]))

# Logical variables stating, for each method, if a gene features a rank lying in the top 5% positions 
lgc.boot = array(FALSE, length(rnk.abs.norm))
names(lgc.boot) = genes
lgc.abs = array(FALSE, length(rnk.abs.norm))
names(lgc.abs) = genes
lgc.wa = array(FALSE, length(rnk.abs.norm))
names(lgc.wa) = genes
# Logical variable stating if a gene features a rank in the top 5% positions by avergae rank
lgc.avg = array(FALSE, length(rnk.abs.norm))
names(lgc.avg) = genes

lgc.boot[is.na(match(genes, top.n.boot))==FALSE] = TRUE
lgc.abs[is.na(match(genes, top.n.abs))==FALSE] = TRUE
lgc.wa[is.na(match(genes, top.n.wa))==FALSE] = TRUE
lgc.avg[is.na(match(genes, top.n.avg))==FALSE] = TRUE

# Retrieve detailed information from the methods used to identify DTL genes
abs.info = f.abs[match(genes, f.abs$X),2:5]
wa.info = f.wa$weighted_distance[match(genes, f.wa$gene)]
boot.info = f.boot[match(genes, f.boot$gene),3:4]

# Retrieve gene annotation according to NCBI 
protprod_LT_LAbrini_ncbi = c()
for(i in 1 : length(genes)){
  tmp.field = annotf[match(genes[i],locusTags),9]
  tmp.field.vec = strsplit(tmp.field,";")[[1]]
  if(length(tmp.field.vec[grepl("product=",tmp.field.vec)==TRUE])==1){
    protprod_LT_LAbrini_ncbi[i] = gsub("product=","",tmp.field.vec[grepl("product=",tmp.field.vec)==TRUE])
  }
}

# Retrieve gene annotations to KO pathways
gene_pathway = fid[match(genes,fid[,4]),3]

# Retrieve gene names of common use
common_name = fid[match(genes,fid[,4]),7]

common_name[is.na(common_name)==TRUE] = "n/a"
common_name[grepl("^RS[0-9]*", common_name)==TRUE] = "n/a"
common_name[grepl("^E[0-9].*\\.*", common_name)==TRUE] = "n/a"
gene_pathway[is.na(gene_pathway)==TRUE] = "n/a"

# Retrieve the mRNA proportions for each gene and polysome fraction
gas.prop = rna.prop.mean.gas.mat[match(genes, rownames(rna.prop.mean.gas.mat)),]
fru.prop = rna.prop.mean.fru.mat[match(genes, rownames(rna.prop.mean.fru.mat)),]

# Retrieve the sign of change at translational level
SignChange = array(NA, length(genes))
SignChange[abs.info$diff > 0] = "Up-regulated"
SignChange[abs.info$diff < 0] = "Down-regulated"

info.df = cbind.data.frame(genes, protprod_LT_LAbrini_ncbi, common_name, gene_pathway, 
                           gas.prop[,1:4], fru.prop[,1:4],
                           lgc.avg[match(genes, names(lgc.avg))], SignChange, 
                           as.numeric(pol.rnk[match(genes, names(pol.rnk))]),
                           as.numeric(rnk.abs.norm[match(genes, names(rnk.abs.norm))]), abs.info,
                           as.numeric(rnk.wa.norm[match(genes, names(rnk.wa.norm))]), wa.info,
                           as.numeric(rnk.boot.norm[match(genes, names(rnk.boot.norm))]), boot.info,
                           lgc.abs[match(genes, names(lgc.abs))], lgc.wa[match(genes, names(lgc.wa))], 
                           lgc.boot[match(genes, names(lgc.boot))]
                           )
rownames(info.df) = genes
colnames(info.df) = c("GeneID","Description","Name","KO functional category L3",c("Gas.A","Gas.B","Gas.C","Gas.D"),
                      c("Fructose.A","Fructose.B","Fructose.C","Fructose.D"),
                      "DTL.by.avg.rank","Sign.of.change","AVG.rank","ABS.rank","Gas.status","Fru.status","Dif","Abs.dif",
                      "WA.rank","WA", "BOOT.rank", "Num.sig.fractions","Max.dif", 
                      "DTL.by.abs.rank","DTL.by.wa.rank","DTL.by.boot.rank"
                      )
info.df.p = info.df[order(info.df$AVG.rank),]
write.xlsx(info.df.p, file = "../results/DTL/dtl.xlsx", row.names=FALSE, col.names=TRUE)


#---------------------------------------------------------------------------
# Plot Venn diagram showing overlaps between up- or down-regulated 
# differentially transcribed genes (DTR) and differentially translated genes 
# (DTL) between syngas and fructose chemostats. 
#---------------------------------------------------------------------------

# RNA-seq 
gene.tot.deg = union(tr.top.n.dn, tr.top.n.up)

# Polysome-seq
gene.pol.deg = top.n.avg

tmp = f.abs$diff[match(gene.pol.deg, f.abs$X)]
lgcUp = array(FALSE, length(gene.pol.deg))
lgcUp[tmp>0] = TRUE
lgcDown = array(FALSE, length(gene.pol.deg))
lgcDown[tmp<0] = TRUE
pol.top.n.dn = gene.pol.deg[lgcDown==TRUE]
pol.top.n.up = gene.pol.deg[lgcUp==TRUE]

# Identify genes occuring as DEGs at either level
gene.tot.deg.only = setdiff(gene.tot.deg, gene.pol.deg)
gene.pol.deg.only = setdiff(gene.pol.deg, gene.tot.deg)

# Identify concordant genes
gene.con = union( intersect(tr.top.n.dn,pol.top.n.dn), intersect(tr.top.n.up,pol.top.n.up) )
gene.tot.up.pol.up = intersect(tr.top.n.up, pol.top.n.up)
gene.tot.dn.pol.dn = intersect(tr.top.n.dn, pol.top.n.dn)

# Identify discordant genes
gene.dis = union( intersect(tr.top.n.dn,pol.top.n.up), intersect(tr.top.n.up,pol.top.n.dn) )
gene.tot.up.pol.dn = intersect(tr.top.n.up, pol.top.n.dn)
gene.tot.dn.pol.up = intersect(tr.top.n.dn, pol.top.n.up)

# Plot Venn diagram
jpeg("../results/dtr.dtl.venn.diagram.jpeg", res=500, width=3000, height=3000)
x = list(tr.up = tr.top.n.up, pol.up = pol.top.n.up, 
         tr.dn = tr.top.n.dn, pol.dn = pol.top.n.dn)
ggvenn(
  x, 
  fill_color = c("#fed9a6", "#cab2d6", "#8dd3c7", "#fb8072"),
  stroke_size = 0.35, text_size = 4.5, set_name_size=3.5
)
dev.off()


# Save DTR and DTL genes in RDA objects
save(gene.pol.deg, file="../results/DTL/dtl.rda")
save(gene.tot.deg, file="../results/DTL/dtr.rda")

