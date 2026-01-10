

require(xlsx)
require(ComplexHeatmap)
require("circlize")
require("GenomicRanges")

load("../data/tss.tts.rda")
load("../data/cap.term.data.rda")

#-------------------------------------------------------------------------------
# Load RNA-binding proteins (RBPs) annotations. ilter out RBPs which are 
# directly annotated to translation processes.
#-------------------------------------------------------------------------------

# Load RBP annotations
rbpf = read.delim("../results/protein-RNA/rbp.descriptions.txt")

# BPs that are not directly annotated to translation processes
rbp.retained = c("LAbrini_16815","LAbrini_08135","LAbrini_15230","LAbrini_16855","LAbrini_17015","LAbrini_15315",
                 "LAbrini_02630","LAbrini_16885","LAbrini_16975","LAbrini_11465","LAbrini_13925","LAbrini_04005",
                 "LAbrini_17055","LAbrini_08685","LAbrini_15020","LAbrini_08580","LAbrini_17670","LAbrini_00995",
                 "LAbrini_15825","LAbrini_12155","LAbrini_12395","LAbrini_10520","LAbrini_16780","LAbrini_15795",
                 "LAbrini_09970","LAbrini_18795","LAbrini_06850","LAbrini_08225","LAbrini_16090","LAbrini_07100",
                 "LAbrini_16970","LAbrini_02175","LAbrini_10855","LAbrini_16010","LAbrini_01925"
                 )

#-------------------------------------------------------------------------------
# Retrieve functional annotations of C. autoethanogenum
#-------------------------------------------------------------------------------

# Load gene annotations to Kegg Orthology pathways according to Level 3 in Table S3 of Valgepea et al. (10.1128/msystems.00026-22)
fid = read.xlsx("./LT_LAbrini_func_annot/LAbrini KEGG IDs_Kaspar curated.xlsx", sheetName = "Kaspar curated")

# Load C. autoethanogenum gene annotations according to NCBI
ncbif = read.delim("./LT_LAbrini_ncbi_dataset/ncbi_dataset/data/GCA_040166795.1/genomic.gff",
                   header=FALSE, comment.char = "#", sep="\t")
annotf = ncbif[ncbif[,3] %in% "CDS",]
# CDS IDs
cdsIds = as.vector(sapply(annotf[,9], function(x) strsplit(x,";")[[1]][1]))
cdsIds = gsub("ID=cds-", "", cdsIds)
# Locus tags
locusTags = as.vector(sapply(annotf[,9], function(x) strsplit(x,";")[[1]][2]))
locusTags = gsub("Parent=gene-", "", locusTags)

#-------------------------------------------------------------------------------
# Load results of the TL-STATUS approach used to define differentially 
# translated (DTL) genes
#-------------------------------------------------------------------------------

f.abs = read.delim("./pol.rnk/abs_translation_diff_results.csv", sep=",")

#_______________________________________________________________________________
#_______________________________________________________________________________
#
# Display RBP interactions with 5’ UTRs of differentially transcribed genes 
# (DTRs) and differentially translated genes (DTLs). Predictions are obtained 
# by catRAPID omics v2 (10.1093/nar/gkab393)
#
#_______________________________________________________________________________
#_______________________________________________________________________________

#-------------------------------------------------------------------------------
# Load predicted RBP-RNA interactions involving differentially translated genes
# and RBPs not directly involved in translation
#-------------------------------------------------------------------------------

# Load differentially translated (DTL) genes
load("../results/DTL/dtl.rda")
pol.top.n = gene.pol.deg
  
# Create a data frame containing TISs per gene (multiple TSSs per gene on multiple rows)
df = data.frame()
for(i in 1 : length(lst$tssByGene)){
  for(j in 1 : length(lst$tssByGene[[i]])){
    df = rbind.data.frame(df, cbind.data.frame(gsub("LABRINI","LAbrini",names(lst$tssByGene)[i]), lst$tssByGene[[i]][j]))
  }
}

# Retrieve TISs of differentially translated genes
pol.deg.tss = df[is.na(match(df[,1], pol.top.n))==FALSE,]
pol.deg.tss[,2] = gsub("\\_","\\.",pol.deg.tss[,2])
pol.deg.tss[,2] = gsub("\\-","\\.",pol.deg.tss[,2])
pol.deg.tss[,2] = gsub("\\+","\\.",pol.deg.tss[,2])

# Retrieve catRAPID predictions of protein-RNA interactions 
f = read.delim("catRAPID.output.genomewide.txt", header=TRUE)
f$Protein_ID = gsub( "\\.","_",f$Protein_ID)
f = f[is.na(f$rnaFrag_start)==FALSE & is.na(f$rnaFrag_end)==FALSE,]

# Retain predicted predictions in the fourth quartile 
f.keep = f[f$Ranking > as.numeric(summary(f$Ranking)[5]) & f$Z_score > as.numeric(summary(f$Z_score)[5]),]
show(summary(f.keep$Z_score))
show(summary(f.keep$Ranking))

# Select protein-RNA interactions involving DTL genes
f.keep = f.keep[is.na(match(f.keep$RNA_ID, pol.deg.tss[,2]))==FALSE,]

# Select protein-RNA interactions involving selected RBPs
f.keep = f.keep[is.na(match(f.keep$Protein_ID, rbp.retained))==FALSE,] 

# Remove interactions where the binding site does not fully fall within the UTR
utrlen = as.vector(sapply(f.keep$RNA_ID, function(x) as.numeric((strsplit(x, "\\.")[[1]][2]))))
f.keep = f.keep[f.keep$rnaFrag_start < utrlen & f.keep$rnaFrag_end < utrlen,]

# Assemble results
rbp.tss.pred = cbind.data.frame(f.keep$Protein_ID, f.keep$RNA_ID, 
                                f.keep$rnaFrag_start, f.keep$rnaFrag_end,
                                f.keep$Z_score, f.keep$RBP_Propensity, 
                                f.keep$numof.RNA_Binding_Motifs_Instances, f.keep$Ranking
                                )
rbp.tss.pred[,1] = gsub("\\.","_",rbp.tss.pred[,1])
colnames(rbp.tss.pred) = c("RBP","Target.TSS","Target.frag.start","Target.frag.end","Zscore","RBP.propensity",
                           "numof.RNA.Binding.Motifs.Instances","Ranking"
                           )
# Get RBP targets
target.tss = pol.deg.tss[match(rbp.tss.pred$Target.TSS, pol.deg.tss[,2]),1]

# Get RBP annotations according to NCBI
protprod_ncbi_rbp = c()
for(i in 1 : length(rbp.tss.pred[,1])){
  tmp.field = annotf[match(rbp.tss.pred[i,1],locusTags),9]
  tmp.field.vec = strsplit(tmp.field,";")[[1]]
  if(length(tmp.field.vec[grepl("product=",tmp.field.vec)==TRUE])==1){
    protprod_ncbi_rbp[i] = gsub("product=","",tmp.field.vec[grepl("product=",tmp.field.vec)==TRUE])
  }
}

# Get RBP target annotations according to NCBI
protprod_ncbi_target = c()
for(i in 1 : length(target.tss)){
  tmp.field = annotf[match(target.tss[i],locusTags),9]
  tmp.field.vec = strsplit(tmp.field,";")[[1]]
  if(length(tmp.field.vec[grepl("product=",tmp.field.vec)==TRUE])==1){
    protprod_ncbi_target[i] = gsub("product=","",tmp.field.vec[grepl("product=",tmp.field.vec)==TRUE])
  }
}

# Assemble RBP-mRNA target pairs
rbp.target.pred = cbind.data.frame(rbp.tss.pred[,1], 
                                   protprod_ncbi_rbp,
                                   target.tss, 
                                   protprod_ncbi_target,
                                   fid[match(target.tss, fid[,4]),3], 
                                   rbp.tss.pred[,3:8]
                                   )
colnames(rbp.target.pred) = c("RBP.locustag","RBP.fullname","Target.locustag","Target.fullname","Target.pathway",
                              "Target.frag.start","Target.frag.end",
                              "Zscore","RBP.propensity","numof.RNA.Binding.Motifs.Instances","Rank"
                              )
rbp.target.pred.names = unique(cbind.data.frame(rbp.target.pred$RBP.locustag, 
                                                rbp.target.pred$RBP.fullname, 
                                                rbp.target.pred$Target.locustag, 
                                                rbp.target.pred$Target.fullname
                                                )
                               )

# Report the unique sets of RBPs and targets
rbp.pred = unique(rbp.target.pred[,1])
target.pred = unique(rbp.target.pred[,3])

# Save the results
rbp.pol.deg.top.n = rbp.target.pred
rbp.pol.deg.top.n.names = rbp.target.pred.names

# Print protein-RBA interactions
write.xlsx(rbp.pol.deg.top.n.names, col.names=TRUE, row.names=FALSE, 
           file="../results/protein-RNA/catrapid.rbp.target.pred.xlsx", 
           sheetName = "pol.deg.top.n")
write.xlsx(rbp.pol.deg.top.n, col.names=TRUE, row.names=FALSE, 
           file="../results/protein-RNA/catrapid.rbp.target.pred.xlsx",
           append = TRUE, sheetName = "detail.pol.deg.top.n")



#-------------------------------------------------------------------------------
# Load predicted RBP-RNA interactions involving differentially transcribed genes
# and RBPs not directly involved in translation
#-------------------------------------------------------------------------------

# Load differentially transcribed (DTR) genes
load("../results/DTR/dtr.rda")
tot.deg = gene.tot.deg

# Create a data farne containing TISs per gene (multiple TSSs per gene on multiple rows)
df = data.frame()
for(i in 1 : length(lst$tssByGene)){
  for(j in 1 : length(lst$tssByGene[[i]])){
    df = rbind.data.frame(df, cbind.data.frame(gsub("LABRINI","LAbrini",names(lst$tssByGene)[i]), lst$tssByGene[[i]][j]))
  }
}

# Retrieve TISs of differentially transcribed genes
tot.deg.tss = df[is.na(match(df[,1], totdeg))==FALSE,]
tot.deg.tss[,2] = gsub("\\_","\\.",tot.deg.tss[,2])
tot.deg.tss[,2] = gsub("\\-","\\.",tot.deg.tss[,2])
tot.deg.tss[,2] = gsub("\\+","\\.",tot.deg.tss[,2])

# Retrieve catRAPID predictions of protein-RNA interactions 
f = read.delim("catRAPID.output.genomewide.txt", header=TRUE)
f$Protein_ID = gsub( "\\.","_",f$Protein_ID)
f = f[(is.na(f$rnaFrag_start)==FALSE) & (is.na(f$rnaFrag_end)==FALSE),]

# Retain predicted predictions in the fourth quartile 
f.keep = f[(f$Ranking > as.numeric(summary(f$Ranking)[5])) & 
             (f$Z_score > as.numeric(summary(f$Z_score)[5])),]
show(summary(f.keep$Z_score))
show(summary(f.keep$Ranking))

# Select protein-RNA interactions involving DTR genes
f.keep = f.keep[is.na(match(f.keep$RNA_ID, tot.deg.tss[,2]))==FALSE,]

# Select protein-RNA interactions involving selected RBPs
f.keep = f.keep[is.na(match(f.keep$Protein_ID, rbp.retained))==FALSE,] 

# Remove interactions where the binding site does not fully fall within the UTR
utrlen = as.vector(sapply(f.keep$RNA_ID, function(x) as.numeric((strsplit(x, "\\.")[[1]][2]))))
f.keep = f.keep[f.keep$rnaFrag_start < utrlen & f.keep$rnaFrag_end < utrlen,]

rbp.tss.pred = cbind.data.frame(f.keep$Protein_ID, f.keep$RNA_ID, 
                                f.keep$rnaFrag_start, f.keep$rnaFrag_end,
                                f.keep$Z_score, f.keep$RBP_Propensity, 
                                f.keep$numof.RNA_Binding_Motifs_Instances, f.keep$Ranking
                                )
rbp.tss.pred[,1] = gsub("\\.","_",rbp.tss.pred[,1])
colnames(rbp.tss.pred) = c("RBP","Target.TSS","Target.frag.start","Target.frag.end","Zscore","RBP.propensity",
                           "numof.RNA.Binding.Motifs.Instances","Ranking"
                           )
# Get RBP targets
target.tss = tot.deg.tss[match(rbp.tss.pred$Target.TSS, tot.deg.tss[,2]),1]

# Get RBP annotations according to NCBI
protprod_ncbi_rbp = c()
for(i in 1 : length(rbp.tss.pred[,1])){
  tmp.field = annotf[match(rbp.tss.pred[i,1],locusTags),9]
  tmp.field.vec = strsplit(tmp.field,";")[[1]]
  if(length(tmp.field.vec[grepl("product=",tmp.field.vec)==TRUE])==1){
    protprod_ncbi_rbp[i] = gsub("product=","",tmp.field.vec[grepl("product=",tmp.field.vec)==TRUE])
  }
}

# Get RBP target annotations according to NCBI
protprod_ncbi_target = c()
for(i in 1 : length(target.tss)){
  tmp.field = annotf[match(target.tss[i],locusTags),9]
  tmp.field.vec = strsplit(tmp.field,";")[[1]]
  if(length(tmp.field.vec[grepl("product=",tmp.field.vec)==TRUE])==1){
    protprod_ncbi_target[i] = gsub("product=","",tmp.field.vec[grepl("product=",tmp.field.vec)==TRUE])
  }
}

# Assemble RBP-mRNA target pairs
rbp.target.pred = cbind.data.frame(rbp.tss.pred[,1], 
                                   protprod_ncbi_rbp, 
                                   target.tss, 
                                   protprod_ncbi_target,
                                   fid[match(target.tss, fid[,4]),3],rbp.tss.pred[,3:8]
                                   )
colnames(rbp.target.pred) = c("RBP.locustag","RBP.fullname","Target.locustag","Target.fullname","Target.pathway",
                              "Target.frag.start","Target.frag.end",
                              "Zscore","RBP.propensity","numof.RNA.Binding.Motifs.Instances","Rank"
                              )
rbp.target.pred.names = unique(cbind.data.frame(rbp.target.pred$RBP.locustag, 
                                                rbp.target.pred$RBP.fullname, 
                                                rbp.target.pred$Target.locustag, 
                                                rbp.target.pred$Target.fullname
                                                )
                               )


# Report the unique sets of RBPs and targets
rbp.pred = unique(rbp.target.pred[,1])
target.pred = unique(rbp.target.pred[,3])

# Save the results
rbp.tot.deg.top.n = rbp.target.pred
rbp.tot.deg.top.n.names = rbp.target.pred.names

# Print the results
write.xlsx(rbp.tot.deg.top.n.names, col.names=TRUE, row.names=FALSE, 
           file="../results/protein-RNA/catrapid.rbp.target.pred.xlsx", 
           append = TRUE, sheetName = "tot.deg")
write.xlsx(rbp.tot.deg.top.n, col.names=TRUE, row.names=FALSE, 
           file="../results/protein-RNA/catrapid.rbp.target.pred.xlsx", 
           append = TRUE, sheetName = "detail.tot.deg")


#----------------------------------------------------------
# Plot the matrix of RBP by target 
#----------------------------------------------------------

rbp.target = rbind.data.frame(rbp.tot.deg.top.n, rbp.pol.deg.top.n)
targets = union(rbp.tot.deg.top.n$Target.locustag, rbp.pol.deg.top.n$Target.locustag)
rbps = union(rbp.tot.deg.top.n$RBP.locustag, rbp.pol.deg.top.n$RBP.locustag)


#
# Annotate RBP targets according to the differential analysis results by RNA-seq and Pol-seq
#

# Retrieve the LFC of targets (if any)
lfc = c()
for(i in 1 : length(targets)){
  if(is.na(match(targets[i], totdeg))==TRUE){
    lfc[i]="NONE"
  }
  if(is.na(match(targets[i], totdeg))==FALSE){
    lfc[i] = sign(totrna.de$logFC[match(targets[i], rownames(totrna.de))])
  }
}

# Retrieve the difference in ribosome occupancy between syngas and fructose according to the TL-STATUS approach
diff = c()
for(i in 1 : length(targets)){
  if(is.na(match(targets[i], pol.top.n))==TRUE){
    diff[i]="NONE"
  }
  if(is.na(match(targets[i], pol.top.n))==FALSE){
    diff[i] = sign(f.abs$diff[match(targets[i], f.abs$X)])
  }
}

# Retrieve down-/up-regulated genes by RNA-seq and Pol-seq
deg.m = matrix("NONE", nrow=2, ncol=length(lfc), dimnames=list(c("TRS","TRL"),targets))
for(i in 1 : length(lfc)){
  if(lfc[i] == 1){deg.m[1,i]="UP"}
  if(lfc[i] == -1){deg.m[1,i]="DOWN"}
  if(diff[i] == 1){deg.m[2,i]="UP"}
  if(diff[i] == -1){deg.m[2,i]="DOWN"}
}


#
# Annotate RBPs according to the differential analysis results by RNA-seq and Pol-seq
#

# Retrieve the LFC of targets (if any)
lfc = c()
for(i in 1 : length(rbps)){
  if(is.na(match(rbps[i], totdeg))==TRUE){
    lfc[i]="NONE"
  }
  if(is.na(match(rbps[i], totdeg))==FALSE){
    lfc[i] = sign(totrna.de$logFC[match(rbps[i], rownames(totrna.de))])
  }
}

# Retrieve the difference in ribosome occupancy between syngas and fructose according to the TL-STATUS approach
diff = c()
for(i in 1 : length(rbps)){
  if(is.na(match(rbps[i], pol.top.n))==TRUE){
    diff[i]="NONE"
  }
  if(is.na(match(rbps[i], pol.top.n))==FALSE){
    diff[i] = sign(f.abs$diff[match(rbps[i], f.abs$X)])
  }
}

# Retrieve down-/up-regulated genes by RNA-seq and Pol-seq
rbp.deg.m = matrix("NONE", nrow=length(lfc), ncol=2, dimnames=list(rbps, c("DTR","DTL")))
for(i in 1 : length(lfc)){
  if(lfc[i] == 1){rbp.deg.m[i,1]="UP"}
  if(lfc[i] == -1){rbp.deg.m[i,1]="DOWN"}
  if(diff[i] == 1){rbp.deg.m[i,2]="UP"}
  if(diff[i] == -1){rbp.deg.m[i,2]="DOWN"}
}


#
# Annotate RBP targets according to Kegg Orthology pathway
#
target.pathway = fid[match(targets,fid[,4]),3]
target.pathway[is.na(target.pathway)==TRUE] = "No KO ID"
target.pathway[is.na(match(target.pathway, "Function unknown"))==FALSE] = "No KO ID"
pathways = unique(target.pathway)
class.m = matrix(0, ncol=length(targets), nrow=length(pathways), dimnames=list(pathways, targets))
for(i in 1 : length(target.pathway)){
  row.idx = match(target.pathway[i], pathways)
  class.m[row.idx,i] = 1
}

# Create a matrix where cells are 0 if an RBP does not regulate a gene and a number >0 otherwise
m = matrix(0, nrow = length(rbps), ncol = length(targets), dimnames=list(rbps, targets))
for(i in 1 : dim(rbp.target)[1]){
  
  # Retrieve the maximal score for a certain RBP and target
  cur.rbp = rbp.target$RBP.locustag[i]
  cur.target = rbp.target$Target.locustag[i]
  tmp = rbp.target[is.na(match(rbp.target$RBP.locustag, cur.rbp))==FALSE & 
                     is.na(match(rbp.target$Target.locustag, cur.target))==FALSE,]  
  bestScore = max(tmp$Rank)
  
  row.idx = match(rbp.target$RBP.locustag[i], rownames(m))
  col.idx = match(rbp.target$Target.locustag[i], colnames(m))
  m[row.idx, col.idx] = bestScore
  
}


png("../results/protein-RNA/catrapid.rbp.target.tot.pol.deg.png", res=600, width=10100, height=4450) #w=8800, 3750

grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 5, widths=unit(c(84.75,2.75,16,58,10.5),"null"), 
                                           heights=unit(c(25),"null") )))

col_fun = colorRamp2(c(0,.45,.5,.55,.6,.65,.7,.75), 
                     c("#eff3ff","#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5","#084594","#08306b"))


class_col_fun = colorRamp2(c(0,1), c("#ece7f2", "#41ab5d")) 

#color map the sources of functional classes
ha_col = HeatmapAnnotation(Pathway = t(class.m), DTR = deg.m[1,], DTL = deg.m[2,],
                           col = list(DTR = c("UP" = "#e08214", "DOWN" = "#542788", "NONE" = "#ece7f2"),
                                      DTL = c("UP" = "#e08214", "DOWN" = "#542788", "NONE" = "#ece7f2"),
                                      Pathway = class_col_fun
                           ), which="column", annotation_name_gp = gpar(cex=.75), gp = gpar(cex=.65), 
                           show_legend = FALSE, simple_anno_size = unit(4, "mm")
)


ha_row = HeatmapAnnotation(DTR = rbp.deg.m[,1], DTL = rbp.deg.m[,2],
                           col = list(DTR = c("UP" = "#e08214", "DOWN" = "#542788", "NONE" = "#ece7f2"),
                                      DTL = c("UP" = "#e08214", "DOWN" = "#542788", "NONE" = "#ece7f2")
                           ), which="row", annotation_name_gp = gpar(cex=.75), gp = gpar(cex=.65),  
                           show_legend = FALSE, height = NULL, simple_anno_size = unit(2.75, "mm")
)

rowLabels = paste(rownames(m), rbpf$Name.to.use[match(rownames(m), rbpf$X)], sep=" - ")

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))

draw( Heatmap(m, col = col_fun, row_labels = gsub("LAbrini_","",rowLabels), row_names_gp = gpar(fontsize = 6.5),
              column_names_gp = gpar(fontsize = 6.5),
              column_labels = gsub("LAbrini_","", colnames(m)),
              top_annotation = ha_col, right_annotation = ha_row,
              width = ncol(m)*unit(3, "mm"), 
              height = nrow(m)*unit(3, "mm"),
              show_heatmap_legend = FALSE, 
              column_title = "RBP binding to 5'UTRs\n\nDTRs & DTLs",column_title_gp = gpar(fontsize = 7.5),
              row_title="RBPs", row_title_side = "left", row_title_gp = gpar(fontsize = 7.5), 
              row_dend_width = unit(4, "mm"), column_dend_height = unit(4, "mm")
              ), newpage=FALSE
      ) 

upViewport()

lgd = Legend(at = c(0,1), labels = c("FALSE","TRUE"), legend_gp = gpar(fontsize = 3, fill = c("#ece7f2", "#41ab5d")),
             title = "Pathway", title_position = "topcenter",
             title_gp = gpar(fontsize = 7, fontface = "bold"), grid_height = unit(3, "mm"), 
             labels_gp = gpar(fontsize = 6.5),
             direction ="vertical"
             )

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))

grid.draw(lgd)

upViewport()

lgd = Legend(at = c("UP","DOWN","NONE"), labels = c("UP","DOWN","NR"), 
             legend_gp = gpar(fontsize = 3, fill = c("#e08214", "#542788","#ece7f2")),
             title = "Reg.", title_position = "topcenter",
             title_gp = gpar(fontsize = 7, fontface = "bold"), grid_height = unit(3, "mm"), 
             labels_gp = gpar(fontsize = 6.15),
             direction ="vertical"
             )

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
grid.draw(lgd)

upViewport()


#_______________________________________________________________________________
#_______________________________________________________________________________
#
# Display RBP interactions with 3’ UTRs of differentially transcribed genes 
# (DTRs) and differentially translated genes (DTLs). Predictions are obtained 
# by catRAPID omics v2 (10.1093/nar/gkab393)
#
#_______________________________________________________________________________
#_______________________________________________________________________________


#-------------------------------------------------------------------------------
# Load predicted RBP-RNA interactions involving differentially translated genes
# and RBPs not directly involved in translation
#-------------------------------------------------------------------------------

# Load differentially translated (DTL) genes
load("../results/DTL/dtl.rda")
pol.top.n = gene.pol.deg


# Create a data frame containing T3Es per gene (multiple T3Es per gene on multiple rows)
df = data.frame()
for(i in 1 : length(lst$ttsByGene)){
  for(j in 1 : length(lst$ttsByGene[[i]])){
    df = rbind.data.frame(df, cbind.data.frame(gsub("LABRINI","LAbrini",names(lst$ttsByGene)[i]), lst$ttsByGene[[i]][j]))
  }
}

# Retrieve T3Es of differentially translated genes
pol.deg.tts = df[is.na(match(df[,1], pol.top.n))==FALSE,]
pol.deg.tts[,2] = gsub("\\_","\\.",pol.deg.tts[,2])
pol.deg.tts[,2] = gsub("\\-","\\.",pol.deg.tts[,2])
pol.deg.tts[,2] = gsub("\\+","\\.",pol.deg.tts[,2])

# Retrieve catRAPID predictions of protein-RNA interactions
f = read.delim("catRAPID.3utr.output.genomewide.txt", header=TRUE)
f$Protein_ID = gsub( "\\.","_",f$Protein_ID)
f = f[is.na(f$rnaFrag_start)==FALSE & is.na(f$rnaFrag_end)==FALSE,]

# Retain predicted predictions in the fourth quartile 
f.keep = f[f$Ranking > as.numeric(summary(f$Ranking)[5]) & 
             f$Z_score > as.numeric(summary(f$Z_score)[5]),]
show(summary(f.keep$Z_score))
show(summary(f.keep$Ranking))

# Select protein-RNA interactions involving DTL genes
f.keep = f.keep[is.na(match(f.keep$RNA_ID, pol.deg.tts[,2]))==FALSE,]

# Select protein-RNA interactions involving selected RBPs
f.keep = f.keep[is.na(match(f.keep$Protein_ID, rbp.retained))==FALSE,] 

# Remove interactions where the binding site does not fully fall within the UTR
utrlen = as.vector(sapply(f.keep$RNA_ID, function(x) as.numeric((strsplit(x, "\\.")[[1]][2]))))
f.keep = f.keep[f.keep$rnaFrag_start < utrlen & f.keep$rnaFrag_end < utrlen,]

rbp.tts.pred = cbind.data.frame(f.keep$Protein_ID, f.keep$RNA_ID, 
                                f.keep$rnaFrag_start, f.keep$rnaFrag_end,
                                f.keep$Z_score, f.keep$RBP_Propensity, 
                                f.keep$numof.RNA_Binding_Motifs_Instances, f.keep$Ranking
                                )
colnames(rbp.tts.pred) = c("RBP","Target.TTS","Target.frag.start","Target.frag.end","Zscore","RBP.propensity",
                           "numof.RNA.Binding.Motifs.Instances","Ranking"
                           )

# Get RBP targets
target.tts = pol.deg.tts[match(rbp.tts.pred$Target.TTS, pol.deg.tts[,2]),1]

# Get RBP annotations according to NCBI
protprod_ncbi_rbp = c()
for(i in 1 : length(rbp.tts.pred[,1])){
  tmp.field = annotf[match(rbp.tts.pred[i,1],locusTags),9]
  tmp.field.vec = strsplit(tmp.field,";")[[1]]
  if(length(tmp.field.vec[grepl("product=",tmp.field.vec)==TRUE])==1){
    protprod_ncbi_rbp[i] = gsub("product=","",tmp.field.vec[grepl("product=",tmp.field.vec)==TRUE])
  }
}

# Get RBP target annotations according to NCBI
protprod_ncbi_target = c()
for(i in 1 : length(target.tts)){
  tmp.field = annotf[match(target.tts[i],locusTags),9]
  tmp.field.vec = strsplit(tmp.field,";")[[1]]
  if(length(tmp.field.vec[grepl("product=",tmp.field.vec)==TRUE])==1){
    protprod_ncbi_target[i] = gsub("product=","",tmp.field.vec[grepl("product=",tmp.field.vec)==TRUE])
  }
}

rbp.target.pred = cbind.data.frame(rbp.tts.pred[,1], 
                                   protprod_ncbi_rbp, 
                                   target.tts, 
                                   protprod_ncbi_target,
                                   fid[match(target.tts, fid[,4]),3], 
                                   rbp.tts.pred[,3:8]
                                   )
colnames(rbp.target.pred) = c("RBP.locustag","RBP.fullname","Target.locustag","Target.fullname","Target.pathway",
                              "Target.frag.start","Target.frag.end",
                              "Zscore","RBP.propensity","numof.RNA.Binding.Motifs.Instances","Rank"
                              )
rbp.target.pred.names = unique(cbind.data.frame(rbp.target.pred$RBP.locustag, rbp.target.pred$RBP.fullname, 
                                                rbp.target.pred$Target.locustag, rbp.target.pred$Target.fullname
                                                )
                               )

# Report the unique sets of RBPs and targets
rbp.pred = unique(rbp.target.pred[,1])
target.pred = unique(rbp.target.pred[,3])

# Save the results
rbp.pol.deg.top.n = rbp.target.pred
rbp.pol.deg.top.n.names = rbp.target.pred.names

# Print the results
write.xlsx(rbp.pol.deg.top.n.names, col.names=TRUE, row.names=FALSE, 
           file="../results/protein-RNA/catrapid.3utr.rbp.target.pred.xlsx",
           append = TRUE, sheetName = "pol.deg.top.n")
write.xlsx(rbp.pol.deg.top.n, col.names=TRUE, row.names=FALSE, 
           file="../results/protein-RNA/catrapid.3utr.rbp.target.pred.xlsx",
           append = TRUE, sheetName = "detail.pol.deg.top.n")


#-------------------------------------------------------------------------------
# Load predicted RBP-RNA interactions involving differentially transcribed genes
# and RBPs not directly involved in translation
#-------------------------------------------------------------------------------

# Load differentially transcribed (DTR) genes
load("../results/DTR/dtr.rda")
tot.deg = gene.tot.deg

# Create a data farne containing T3Es per gene (multiple T3Es per gene on multiple rows)
df = data.frame()
for(i in 1 : length(lst$ttsByGene)){
  for(j in 1 : length(lst$ttsByGene[[i]])){
    df = rbind.data.frame(df, cbind.data.frame(gsub("LABRINI","LAbrini",names(lst$ttsByGene)[i]), lst$ttsByGene[[i]][j]))
  }
}

# Retrieve T3Es of the DTR genes
tot.deg.tts = df[is.na(match(df[,1], totdeg))==FALSE,]
tot.deg.tts[,2] = gsub("\\_","\\.",tot.deg.tts[,2])
tot.deg.tts[,2] = gsub("\\-","\\.",tot.deg.tts[,2])
tot.deg.tts[,2] = gsub("\\+","\\.",tot.deg.tts[,2])

# Retrieve catRAPID predictions of protein-RNA interactions
f = read.delim("catRAPID.3utr.output.genomewide.txt", header=TRUE)
f$Protein_ID = gsub( "\\.","_",f$Protein_ID)
f = f[is.na(f$rnaFrag_start)==FALSE & is.na(f$rnaFrag_end)==FALSE,]

# Retain predicted predictions in the fourth quartile 
f.keep = f[f$Ranking > as.numeric(summary(f$Ranking)[5]) & 
             f$Z_score > as.numeric(summary(f$Z_score)[5]),]
show(summary(f.keep$Z_score))
show(summary(f.keep$Ranking))

# Select protein-RNA interactions involving DTR genes
f.keep = f.keep[is.na(match(f.keep$RNA_ID, tot.deg.tts[,2]))==FALSE,]

# Select protein-RNA interactions involving selected RBPs
f.keep = f.keep[is.na(match(f.keep$Protein_ID, rbp.retained))==FALSE,] 

# Remove interactions where the binding site does not fully fall within the UTR
utrlen = as.vector(sapply(f.keep$RNA_ID, function(x) as.numeric((strsplit(x, "\\.")[[1]][2]))))
f.keep = f.keep[f.keep$rnaFrag_start < utrlen & f.keep$rnaFrag_end < utrlen,]

rbp.tts.pred = cbind.data.frame(f.keep$Protein_ID, 
                                f.keep$RNA_ID, 
                                f.keep$rnaFrag_start, 
                                f.keep$rnaFrag_end,
                                f.keep$Z_score, 
                                f.keep$RBP_Propensity, 
                                f.keep$numof.RNA_Binding_Motifs_Instances, 
                                f.keep$Ranking
                                )
rbp.tts.pred[,1] = gsub("\\.","_",rbp.tts.pred[,1])
colnames(rbp.tts.pred) = c("RBP","Target.TTS","Target.frag.start","Target.frag.end","Zscore","RBP.propensity",
                           "numof.RNA.Binding.Motifs.Instances","Ranking"
                           )

# Get unique RBP targets
target.tts = tot.deg.tts[match(rbp.tts.pred$Target.TTS, tot.deg.tts[,2]),1]

# Get RBP annotations to Kegg Orthology pathways
protprod_ncbi_rbp = c()
for(i in 1 : length(rbp.tts.pred[,1])){
  tmp.field = annotf[match(rbp.tts.pred[i,1],locusTags),9]
  tmp.field.vec = strsplit(tmp.field,";")[[1]]
  if(length(tmp.field.vec[grepl("product=",tmp.field.vec)==TRUE])==1){
    protprod_ncbi_rbp[i] = gsub("product=","",tmp.field.vec[grepl("product=",tmp.field.vec)==TRUE])
  }
}

# Get RBP target annotations to Kegg Orthology pathways
protprod_ncbi_target = c()
for(i in 1 : length(target.tts)){
  tmp.field = annotf[match(target.tts[i],locusTags),9]
  tmp.field.vec = strsplit(tmp.field,";")[[1]]
  if(length(tmp.field.vec[grepl("product=",tmp.field.vec)==TRUE])==1){
    protprod_ncbi_target[i] = gsub("product=","",tmp.field.vec[grepl("product=",tmp.field.vec)==TRUE])
  }
}


rbp.target.pred = cbind.data.frame(rbp.tts.pred[,1], 
                                   protprod_ncbi_rbp, 
                                   target.tts, 
                                   protprod_ncbi_target,
                                   fid[match(target.tts, fid[,4]),3],
                                   rbp.tts.pred[,3:8]
                                   )
colnames(rbp.target.pred) = c("RBP.locustag","RBP.fullname","Target.locustag","Target.fullname","Target.pathway",
                              "Target.frag.start","Target.frag.end",
                              "Zscore","RBP.propensity","numof.RNA.Binding.Motifs.Instances","Rank"
                              )
rbp.target.pred.names = unique(cbind.data.frame(rbp.target.pred$RBP.locustag, 
                                                rbp.target.pred$RBP.fullname, 
                                                rbp.target.pred$Target.locustag, 
                                                rbp.target.pred$Target.fullname
                                                )
                               )


# Report the unique sets of RBPs and targets
rbp.pred = unique(rbp.target.pred[,1])
target.pred = unique(rbp.target.pred[,3])

# Save the results
rbp.tot.deg.top.n = rbp.target.pred
rbp.tot.deg.top.n.names = rbp.target.pred.names

# Print the results
write.xlsx(rbp.tot.deg.top.n.names, col.names=TRUE, row.names=FALSE, 
           file="../results/protein-RNA/catrapid.3utr.rbp.target.pred.xlsx",
            append = TRUE, sheetName = "tot.deg")
write.xlsx(rbp.tot.deg.top.n, col.names=TRUE, row.names=FALSE, 
           file="../results/protein-RNA/catrapid.3utr.rbp.target.pred.xlsx", 
           append = TRUE, sheetName = "detail.tot.deg")


#------------------------------------------------------
# Plot the matrix of RBP by target 
#------------------------------------------------------

rbp.target = rbind.data.frame(rbp.tot.deg.top.n, rbp.pol.deg.top.n)
targets = union(rbp.tot.deg.top.n$Target.locustag, rbp.pol.deg.top.n$Target.locustag)
rbps = union(rbp.tot.deg.top.n$RBP.locustag, rbp.pol.deg.top.n$RBP.locustag)

#
# Annotate RBP targets according to the differential analysis results by RNA-seq and Pol-seq
# 
# Retrieve the LFC of targets (if any)
lfc = c()
for(i in 1 : length(targets)){
  if(is.na(match(targets[i], totdeg))==TRUE){
    lfc[i]="NONE"
  }
  if(is.na(match(targets[i], totdeg))==FALSE){
    lfc[i] = sign(totrna.de$logFC[match(targets[i], rownames(totrna.de))])
  }
}

# Retrieve the difference in ribosome occupancy between syngas and fructose according to the TL-STATUS approach
diff = c()
for(i in 1 : length(targets)){
  if(is.na(match(targets[i], pol.top.n))==TRUE){
    diff[i]="NONE"
  }
  if(is.na(match(targets[i], pol.top.n))==FALSE){
    diff[i] = sign(f.abs$diff[match(targets[i], f.abs$X)])
  }
}

deg.m = matrix("NONE", nrow=2, ncol=length(lfc), dimnames=list(c("DTR","DTL"),targets))
for(i in 1 : length(lfc)){
  if(lfc[i] == 1){deg.m[1,i]="UP"}
  if(lfc[i] == -1){deg.m[1,i]="DOWN"}
  if(diff[i] == 1){deg.m[2,i]="UP"}
  if(diff[i] == -1){deg.m[2,i]="DOWN"}
}


#
# Annotate RBPs according to the differential analysis results by RNA-seq and Pol-seq
# 
# Retrieve the LFC of RBPs (if any)
lfc = c()
for(i in 1 : length(rbps)){
  if(is.na(match(rbps[i], totdeg))==TRUE){
    lfc[i]="NONE"
  }
  if(is.na(match(rbps[i], totdeg))==FALSE){
    lfc[i] = sign(totrna.de$logFC[match(rbps[i], rownames(totrna.de))])
  }
}

# Retrieve the difference in ribosome occupancy between syngas and fructose according to the TL-STATUS approach
diff = c()
for(i in 1 : length(rbps)){
  if(is.na(match(rbps[i], pol.top.n))==TRUE){
    diff[i]="NONE"
  }
  if(is.na(match(rbps[i], pol.top.n))==FALSE){
    diff[i] = sign(f.abs$diff[match(targets[i], f.abs$X)])
  }
}

rbp.deg.m = matrix("NONE", nrow=length(lfc), ncol=2, dimnames=list(rbps, c("TRS","TRL")))
for(i in 1 : length(lfc)){
  if(lfc[i] == 1){rbp.deg.m[i,1]="UP"}
  if(lfc[i] == -1){rbp.deg.m[i,1]="DOWN"}
  if(diff[i] == 1){rbp.deg.m[i,2]="UP"}
  if(diff[i] == -1){rbp.deg.m[i,2]="DOWN"}
}

#
# Annotate RBP targets according to Kegg Orthology pathways
#
target.pathway = fid[match(targets,fid[,4]),3]
target.pathway[is.na(target.pathway)==TRUE] = "No KO ID"
target.pathway[is.na(match(target.pathway, "Function unknown"))==FALSE] = "No KO ID"
pathways = unique(target.pathway)
class.m = matrix(0, ncol=length(targets), nrow=length(pathways), dimnames=list(pathways, targets))
for(i in 1 : length(target.pathway)){
  row.idx = match(target.pathway[i], pathways)
  class.m[row.idx,i] = 1
}

# Create a matrix where cells are 0 if an RBP does not regulate a gene and a number >0 otherwise
m = matrix(0, nrow = length(rbps), ncol = length(targets), dimnames=list(rbps, targets))
for(i in 1 : dim(rbp.target)[1]){
  
  # Retrieve the maximal score for a certain RBP and target
  cur.rbp = rbp.target$RBP.locustag[i]
  cur.target = rbp.target$Target.locustag[i]
  tmp = rbp.target[is.na(match(rbp.target$RBP.locustag, cur.rbp))==FALSE & 
                     is.na(match(rbp.target$Target.locustag, cur.target))==FALSE,]  
  bestScore = max(tmp$Rank)
  
  row.idx = match(rbp.target$RBP.locustag[i], rownames(m))
  col.idx = match(rbp.target$Target.locustag[i], colnames(m))
  m[row.idx, col.idx] = bestScore
  
}


#
# Plot section
#

col_fun = colorRamp2(c(0,.45,.5,.55,.6,.65,.7,.75), 
                     c("#eff3ff","#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5","#084594","#08306b"))

#color map the sources of functional classes and if the RBP targets are regulated
ha_col = HeatmapAnnotation(Pathway = t(class.m), DTR = deg.m[1,], DTL = deg.m[2,],
                           col = list(DTR = c("UP" = "#e08214", "DOWN" = "#542788","NONE" = "#ece7f2"),
                                      DTL = c("UP" = "#e08214", "DOWN" = "#542788","NONE" = "#ece7f2"),
                                      Pathway = c("1"="#41ab5d","0"="#ece7f2")
                           ), which="column", annotation_name_gp = gpar(cex=.65), 
                           show_legend = FALSE, simple_anno_size = unit(4, "mm")
                           )

# color if the RBPs is up-/down regulared at the total RNA / polysome-bound RNA level
ha_row = HeatmapAnnotation(DTR = rbp.deg.m[,1], DTL = rbp.deg.m[,2],
                           col = list(DTR = c("UP" = "#e08214", "DOWN" = "#542788", "NONE" = "#ece7f2"),
                                      DTL = c("UP" = "#e08214", "DOWN" = "#542788", "NONE" = "#ece7f2")
                           ), which="row", annotation_name_gp = gpar(cex=.65), gp = gpar(cex=.65),
                           show_legend = FALSE, simple_anno_size = unit(2.75, "mm")
                           )

rowLabels = paste(rownames(m), rbpf$Name.to.use[match(rownames(m), rbpf$X)], sep=" - ")

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 4))

draw( Heatmap(m, col = col_fun, row_labels = gsub("LAbrini_", "", rowLabels), row_names_gp = gpar(fontsize = 6.5),
              column_names_gp = gpar(fontsize = 6.5),
              column_labels = gsub("LAbrini_", "", colnames(m)),
              top_annotation = ha_col, right_annotation = ha_row,
              width = ncol(m)*unit(3, "mm"), 
              height = nrow(m)*unit(3, "mm"),
              show_heatmap_legend = FALSE, 
              column_title = "RBP binding to 3'UTRs\n\nDTRs & DTLs",column_title_gp = gpar(fontsize = 7.5),
              row_title="RBPs", row_title_side = "left", row_title_gp = gpar(fontsize = 7.5), 
              row_dend_width = unit(4, "mm"), column_dend_height = unit(4, "mm")
              ), newpage=FALSE
      ) 

upViewport()

lgd = Legend(at = c(0,.45,.5,.55,.6,.65,.7,.75), 
             col_fun = c("#eff3ff","#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5","#084594","#08306b"), 
             title = "Score", title_position = "topcenter",
             title_gp = gpar(fontsize = 7, fontface = "bold"), legend_gp = gpar(fontsize = 4.5), 
             grid_height = unit(3, "mm"), labels_gp = gpar(fontsize = 8.5),
             direction ="vertical"
             )
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 5))
grid.draw(lgd)

upViewport()

dev.off()
