
library(pheatmap)
library(scatterplot3d) 

OriginCodes = unique(ct$OriginCode)

if (!dir.exists("validation")){dir.create("validation")}
for (code in OriginCodes){
ct_subset = subset(ct, OriginCode == code)
#find GSE
GSE = ct_subset$GSE[1]
#find GPL
GPL = ct_subset$GPL[1]
#map

#Pull out raw expression
rawValues = read.csv(paste("rawExpression/", GPL, ".RawExpressionValues.table", sep=""), sep="\t" )

#pull out gene mapping
geneMappings = read.csv(paste("geneMappings/", GPL, ".Probe2EntrezMap.table", sep=""), sep="\t" )
geneMappings = subset(geneMappings, select=c("probe", "GeneID"))
#build matrix
rawValues$probe_id = rownames(rawValues)

expressions = merge(rawValues, geneMappings, by.x="probe_id", by.y="probe")

expressions = expressions[expressions$GeneID %in% dz_sig$GeneID, colnames(expressions) %in% c("GeneID", ct_subset$GSM)]

if (nrow(expressions) < 3) {next}
expressions_by_gene = aggregate(. ~ GeneID, data= expressions, mean)
geneids = expressions_by_gene$GeneID
expressions_by_gene = expressions_by_gene[, -1]
genenames = merge(geneids, dz_sig, by.x=1, by.y="GeneID", sort=F)$Symbol
rownames(expressions_by_gene) = genenames

annotation = subset(ct_subset, select= c("GSM", "ClassCode"))
rownames(annotation) = annotation$GSM
annotation$ClassCode = as.factor(annotation$ClassCode)
annotation = subset(annotation, select=c("ClassCode"))

my.cols <- brewer.pal(9, "Blues")
pheatmap(expressions_by_gene, col = my.cols, annotation = annotation, cellwidth  = 6,  cellheight= 12, show_colnames=F, legend=F, show_rownames=T , filename=paste("validation/heatmap_", code, ".pdf", sep="")) #

if (nrow(expressions_by_gene) < ncol(expressions_by_gene)){
  pca = princomp(t(expressions_by_gene))
  sample_class = annotation[rownames(pca$scores),1]
  
  pdf(paste("validation/pca_", code, ".pdf", sep=""))
  scatterplot3d(pca$scores[,1:3], pch=20, color=sample_class) 
  dev.off()
}
}