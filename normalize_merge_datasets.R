hcc_meta <- read.csv("dz_sig_meta_analysis.txt", sep="\t")

ct<-readCodingTable("meta_input.txt")

studies = c("GPL96", "GPL570", "GPL887", "GPL1352", "GPL1708", "GPL4133", "GPL10558")
ct = subset(ct, GPL %in% studies)

OriginCodes = unique(ct$GPL)

gene_ids = read.csv("~/Desktop/gene_info_hs.csv")
expr = data.frame(gene = gene_ids$GeneID )

#merge all data into single matrix
gene_id_mapping = read.csv("~/Desktop/gene_info_hs.csv")
gene_id_mapping = gene_id_mapping[, c("GeneID","Symbol")]

for (code in OriginCodes){
  print(code)
  ct_subset = subset(ct, GPL == code)
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
  #geneMappings = merge(geneMappings, gene_id_mapping, by="GeneID")
  #geneMappings = subset(geneMappings, select=c("probe", "Symbol"))
  
  #build matrix
  rawValues$probe_id = rownames(rawValues)
  
  expressions = merge(rawValues, geneMappings, by.x="probe_id", by.y="probe")
  
  expressions = expressions[, !colnames(expressions) %in% c("probe_id")]
  
  expressions_by_gene = aggregate(. ~ GeneID, data= expressions, function(x){mean(x, na.rm=T)})
  
  write.csv(expressions_by_gene, paste("geneExpression/", code, ".csv", sep=""))
}

expression_all = gene_id_mapping
for (code in OriginCodes){
  expression  =  read.csv(paste("geneExpression/", code, ".csv", sep=""), row.names=1)
  print(dim(expression))
  # expression_all = merge(expression_all, expression, by="Symbol", all.x=T)
}

valid_sample_counts = sapply(1:nrow(expression_all), function(id){sum(is.na(expression_all[id,]))})

expression_all_subset = expression_all[valid_sample_counts < 1362, ]

write.csv(expression_all_subset, paste("geneExpression/",  "breast.csv", sep=""))

#expression_all_subset[1:3, 1:3]


