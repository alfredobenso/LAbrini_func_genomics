

require(edgeR)
library(dplyr)
library(dendextend)


# Load gene annotations to Kegg Orthology pathways according to Level 3 in ST3 of Valgepea et al. (10.1128/msystems.00026-22)
fid = read.xlsx("../data/Kegg Orthology pathway annotations.xlsx", sheetName = "Kaspar curated")

# Load Clostridium autoethanogenum gene NCBI annotations 
ncbif = read.delim("../data/GCA_040166795.1/genomic.gff", header=FALSE, comment.char = "#", sep="\t")
annotf = ncbif[ncbif[,3] %in% "CDS",]
# CDS IDs
cdsIds = as.vector(sapply(annotf[,9], function(x) strsplit(x,";")[[1]][1]))
cdsIds = gsub("ID=cds-", "", cdsIds)
# Locus tags
locusTags = as.vector(sapply(annotf[,9], function(x) strsplit(x,";")[[1]][2]))
locusTags = gsub("Parent=gene-", "", locusTags)


#-------------------------------------------------------------------------------
# Replicate consistency analysis on RNA-seq data  
#-------------------------------------------------------------------------------

# Load transcript raw read counts in LAbrini fructose and syngas chemostats
x = read.delim("../data/Gene_count_transcriptome.txt", header=TRUE, sep="\t")

group <- factor(c(1,2,1,2,1,2))
y <- DGEList(counts=x[,-1],group=group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
# Normalization
y <- normLibSizes(y)
design <- model.matrix(~group)
y <- estimateDisp(y) 


# Check the replicate consistency (gas)
cor_gas_1_vs_2 = signif(as.numeric(cor.test(y$count[,2], y$count[,4], method="spearman")$estimate), digits=2)
cor_gas_1_vs_3 = signif(as.numeric(cor.test(y$count[,2], y$count[,6], method="spearman")$estimate), digits=2)
cor_gas_2_vs_3 = signif(as.numeric(cor.test(y$count[,4], y$count[,6], method="spearman")$estimate), digits=2)

# Check the replicate consistency (fructose)
cor_fru_1_vs_2 = signif(as.numeric(cor.test(y$count[,1], y$count[,3], method="spearman")$estimate), digits=2)
cor_fru_1_vs_3 = signif(as.numeric(cor.test(y$count[,1], y$count[,5], method="spearman")$estimate), digits=2)
cor_fru_2_vs_3 = signif(as.numeric(cor.test(y$count[,3], y$count[,5], method="spearman")$estimate), digits=2)

#
# Plotting section on replicate consistency analysis
#

# relabel the columns of the count matrix for plotting purposes
mplot = y$counts
colnames(mplot) = c("FRUCTOSE_REP_1","GAS_REP_1","FRUCTOSE_REP_2","GAS_REP_2","FRUCTOSE_REP_3","GAS_REP_3")

jpeg("../results/DTR/totRNAseq_corplot_gas.jpeg", res=600,width=5500, height=2000)
par(mfrow=c(1,3), mar=c(5,5,1,1))

plot(log2(y$count[,2]), log2(y$count[,4]), col="#807dba",pch=19, cex=.85, 
     xlab=expression(log[2]~read~counts^(SYNGAS_REP_1)), ylab=expression(log[2]~read~counts^(SYNGAS_REP_2)))
text(x=5, y=18, cex=.95, labels=expression(paste("Spearman's ",rho," = ")))
text(x=9.15, y=18, cex=.95, labels=eval(cor_gas_1_vs_2))

plot(log2(y$count[,2]), log2(y$count[,6]), col="#807dba",pch=19, cex=.8, 
     xlab=expression(log[2]~read~counts^(SYNGAS_REP_1)), ylab=expression(log[2]~read~counts^(SYNGAS_REP_3)))
text(x=5, y=18, cex=.95, labels=expression(paste("Spearman's ",rho," = ")))
text(x=9.15, y=18, cex=.95, labels=eval(cor_gas_1_vs_2))

plot(log2(y$count[,4]), log2(y$count[,6]), col="#807dba",pch=19, cex=.8, 
     xlab=expression(log[2]~read~counts^(SYNGAS_REP_2)), ylab=expression(log[2]~read~countse^(SYNGAS_REP_3)))
text(x=5, y=18, cex=.95, labels=expression(paste("Spearman's ",rho," = ")))
text(x=9.15, y=18, cex=.95, labels=eval(cor_gas_1_vs_2))

dev.off()


png("../results/DTR/totRNAseq_corplot_fructose.png", res=600,width=5500, height=2000)

par(mfrow=c(1,3), mar=c(5,5,1,1))
plot(log2(y$count[,1]), log2(y$count[,3]), col="#fc8d59",pch=19, cex=.85, 
     xlab=expression(log[2]~read~counts^(SYNGAS_REP_1)), ylab=expression(log[2]~read~counts^(SYNGAS_REP_2)))
text(x=5, y=17.5, cex=.95, labels=expression(paste("Spearman's ",rho," = ")))
text(x=9.15, y=17.5, cex=.95, labels=eval(cor_fru_1_vs_2))

plot(log2(y$count[,1]), log2(y$count[,5]), col="#fc8d59",pch=19, cex=.85, 
     xlab=expression(log[2]~read~counts^(SYNGAS_REP_1)), ylab=expression(log[2]~read~counts^(SYNGAS_REP_3)))
text(x=5, y=17.5, cex=.95, labels=expression(paste("Spearman's ",rho," = ")))
text(x=9.15, y=17.5, cex=.95, labels=eval(cor_fru_1_vs_3))

plot(log2(y$count[,3]), log2(y$count[,5]), col="#fc8d59",pch=19, cex=.85,
     xlab=expression(log[2]~read~counts^(SYNGAS_REP_2)), ylab=expression(log[2]~read~counts^(SYNGAS_REP_3)))
text(x=5, y=17.5, cex=.95, labels=expression(paste("Spearman's ",rho," = ")))
text(x=9.15, y=17.5, cex=.95, labels=eval(cor_fru_2_vs_3))

dev.off()


png("../results/DTR/totRNAseq_mds.png", res=600,width=2500, height=2500)

e1 = expression(Leading~log[2]~fold-change~dim~1~("77%"))
e2 = expression(Leading~log[2]~fold-change~dim~2~("12%"))
a <- plotMDS(log2(mplot), xlim=c(-2,2), ylim=c(-2,2), pch=19, 
             col=c("#807dba","#fc8d59","#807dba","#fc8d59","#807dba","#fc8d59"),
             xlab = "",ylab="", cex=.7, var.explained=FALSE, main="MDS plot")
mtext(side=1, line=2.5, text=e1, cex=.7)
mtext(side=2, line=2.5, text=e2, cex=.7)
text(x=c(a$x[1]-.25, a$x[2]-.35, a$x[3]-.25, a$x[4]+.5, a$x[5]-.25, a$x[6]+.25), 
     y=c(a$y[1]+.15, y=a$y[2]+0.15,y=a$y[3]+0.15,y=a$y[4]-0.15,y=a$y[5]-0.15,y=a$y[6]-0.25),
     labels = c("FRUCTOSE_REP_1","SYNGAS_REP_1","FRUCTOSE_REP_2","SYNGAS_REP_2","FRUCTOSE_REP_3","SYNGAS_REP_3"), 
     col=c("#807dba","#fc8d59","#807dba","#fc8d59","#807dba","#fc8d59"),ce=.5)

dev.off()


# Fitting Hierarchical clustering Model to training dataset
distance_mat <- dist(t(log2(mplot+1)), method = 'euclidean')
png("../results/DTR/totRNAseq_cluster_dendrogram.png", res=600,width=2500, height=2750)
par(mar=c(8,5,1,1))
set.seed(240)  # Setting seed
hc <- hclust(distance_mat, method = "average")
dend <- as.dendrogram(hc)
col=c("#807dba","#807dba","#807dba","#fc8d59","#fc8d59","#fc8d59")
labels_colors(dend) <- col
labels(dend) <- c("SYNGAS_REP_2","SYNAGS_REP_1","SYNGAS_REP_3",
                  "FRUCTOSE_REP_3","FRUCTOSE_REP_1","FRUCTOSE_REP_2")
plot(dend, cex=.5, main="",sub="", xlab="", ylab="Height", ylim=c(0,100))
dev.off()


#-------------------------------------------------------------------------------
# Analysis of differentially transcribed (DTR) genes
#-------------------------------------------------------------------------------

fit <- glmQLFit(y, design)
qlf.2vs1 <- glmQLFTest(fit, coef=2)
topTags(qlf.2vs1)
summary(decideTests(qlf.2vs1, p.value=.05))
res=decideTests(qlf.2vs1, p.value=.05)
res@.Data[res@.Data==1,]
res@.Data[res@.Data==-1,]
# Get the table of ordered genes
totrna.de = qlf.2vs1$table
rownames(totrna.de) = x[keep==TRUE,1]
save(totrna.de, file="../results/DTR/totrna.de.rda")

