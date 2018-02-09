#merge gene expression levels
library(Rtsne)
library(sva)
library(rgl)
library(glmnet)
ct<-readCodingTable("meta_input_new.txt")
#ct <- subset(ct, valid == 1)
#ct <- subset(ct, er_pos == "er_pos" & !is.na(ClassCode) )

#studies = c("GSE18728",   "GSE50948", "GSE32646", "GSE23988", "GSE42822",  "GSE4779", "GSE22226", "GSE59515", "GSE25066") #"GSE8465",
studies = c("GPL1390") #c("GPL96", "GPL570","GPL571", "GPL887", "GPL1352",  "GPL1708", "GPL4133", "GPL5325", "GPL6884", "GPL7504", "GPL10558")
#GPL1390 has not Entrez GeneID

ct = subset(ct, !GPL %in% studies)

GSM = ct$GSM[duplicated(ct$GSM)]

ct$duplicate = duplicated(ct$GSM)

write.csv(ct, "meta_input_v2.csv")

#statistics
sort(table(ct$GPL))
sort(table(ct$GSE))
sort(table(ct$OriginCode))
sort(table(ct$pr))
sort(table(ct$her2))
sort(table(ct$grade))
sort(table(ct$stage))
sort(table(ct$pre_treatment))


OriginCodes = unique(ct$GPL)

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
  #pull out gene mapping
  
  #Pull out raw expression
  rawValues = read.csv(paste("rawExpression/", GPL, ".RawExpressionValues.table", sep=""), sep="\t" )
  
  #pull out gene mapping
  geneMappings = read.csv(paste("geneMappings/", GPL, ".Probe2EntrezMap.table", sep=""), sep="\t" )
  geneMappings = subset(geneMappings, select=c("probe", "GeneID")) #GPL570: " /// ",  #GPL5325: "|"
  #
  if (class(geneMappings$GeneID) != "integer"){
    newGeneMappings = data.frame()
    for (i in 1:nrow(geneMappings)){
      newGeneMappings = rbind(newGeneMappings, data.frame(probe = geneMappings$probe[i], GeneID = unique(na.omit(as.numeric(unlist(strsplit(unlist(geneMappings$GeneID[i]), "[^0-9]+")))))))
    }
    geneMappings = newGeneMappings
  }
  
  #print(class(geneMappings$GeneID))
  #geneMappings = merge(geneMappings, gene_id_mapping, by="GeneID")
  #geneMappings = subset(geneMappings, select=c("probe", "Symbol"))
  
  #build matrix
  rawValues$probe_id = rownames(rawValues)
  
  expressions = merge(rawValues, geneMappings, by.x="probe_id", by.y="probe")
  
  expressions = expressions[, !colnames(expressions) %in% c("probe_id")]
  
  expressions_by_gene = aggregate(. ~ GeneID, data= expressions, function(x){mean(x, na.rm=T)})
  
  write.csv(expressions_by_gene, paste("geneExpression/", code, ".csv", sep=""))
}


####
#merge into a big matrix
####
expression_all = data.frame(GeneID = gene_id_mapping$GeneID)
OriginCodes_select = OriginCodes  #[!OriginCodes %in% c("GPL96", "GPL570")] #GPL96 and 570 included a small probe set
for (code in OriginCodes_select ){
  expression  =  read.csv(paste("geneExpression/", code, ".csv", sep=""), row.names=1)
  print( dim(expression))
  print(code)
  expression_all = merge(expression_all, expression, by="GeneID", all.x=T)
}

na_sample_counts = apply(expression_all, 1, function(x){sum(is.na(x))})

expression_all_subset = expression_all[na_sample_counts < 1, ]

write.csv(expression_all_subset, paste("geneExpression/",  "breast.csv", sep=""))

#batch effect
expression_all_subset = read.csv(paste("geneExpression/",  "breast.csv", sep=""), row.names = 1)

#expression_all_subset[1:3, 1:3]
expr = t(expression_all_subset[,-1]) #first col is GeneID
expr = expr[!duplicated(expr), ] #why there are duplicated rows

#read organization
gse_meta = read.csv("datasets.csv")
ct_meta = ct # merge(ct, gse_meta, by=c("GSE", "GPL"))
ct_meta = ct_meta[!duplicated(ct_meta$GSM),]
rownames(ct_meta) = ct_meta$GSM
pheno = ct_meta[rownames(expr),]


tsne_out <- Rtsne(expr) # Run TSNE
#batch effect is profound.
plot(tsne_out$Y, col = as.numeric(as.factor(pheno$GSE))) # Plot the result

#
#remove batch effect
modcombat = model.matrix(~1, data=pheno)
batch = pheno$GSE
combat_edata = ComBat(dat=t(expr), batch=batch, mod=modcombat)

write.csv(combat_edata, paste("geneExpression/",  "breast_combat.csv", sep=""), row.names = T)


tsne_out <- Rtsne(t(combat_edata), dims = 2) # Run TSNE
#effect from studies looks removed
plot(tsne_out$Y, col = as.numeric(as.factor(pheno$prc)), pch=20, xlab="", ylab="") # Plot the result
legend()

plot(tsne_out$Y[,1], tsne_out$Y[,2], col = as.numeric(as.factor(paste( pheno$er))), pch=20) # Plot the result
plot(tsne_out$Y[,1], tsne_out$Y[,2],  col = as.numeric(as.factor(pheno$prc)), pch=20) # Plot the result

#plot3d(tsne_out$Y, col = as.numeric(as.factor(pheno$er)), pch=20) # Plot the result

#normalize data?
#select genes 
#binormial response

x = as.matrix(t(combat_edata)) #feature_columns
y = NA
y[pheno$prc == "pcr_yes"]  = 1
y[pheno$prc == "pcr_no"]  = 0

x = x[!is.na(y), ]
y = y[!is.na(y)]

#find the best alpha
alphas = seq(0.2, 0.9, 0.1)
mse = NULL
for (a in alphas){
  cvfit = cv.glmnet(x, y, alpha = a, standardize=F, parallel = TRUE, family = "binomial", type.measure = "auc")
  mse = c(mse, max(cvfit$cvm))
}
best_alpha = alphas[mse == max(mse)][1]

cvfit = cv.glmnet(x, y, alpha = best_alpha, standardize=F, parallel = TRUE, family = "binomial", type.measure = "auc")
pred_feature = as.matrix(coef(cvfit, s = "lambda.min"))
best_lambda = cvfit$lambda.min
pred_feature = data.frame(feature = rownames(pred_feature), pred_coeff = pred_feature)

feature_selected = pred_feature$feature[pred_feature[,2] != 0]

tsne_out <- Rtsne(t(combat_edata[rownames(combat_edata) %in% feature_selected, ]), dims = 2) # Run TSNE
plot(tsne_out$Y[,1], tsne_out$Y[,2], col = as.numeric(as.factor(pheno$er)), pch=20) # Plot the result

#cytolytic activity 
gzma = expression_all[expression_all$GeneID == "3001",]
prf1 = expression_all[expression_all$GeneID == "5551",]
cyto = data.frame(GSM = names(gzma), gzma = as.numeric(gzma), prf1 = as.numeric(prf1))
ct_cyto = merge(ct, cyto, by="GSM")
ct_cyto = subset(ct_cyto, !is.na(gzma) & !is.na(prf1))
ct_cyto$cyto = (ct_cyto$gzma + ct_cyto$prf1)/2
anova(lm(cyto ~ er, ct_cyto))
by(ct_cyto$cyto, ct_cyto$er, mean)
anova(lm(cyto ~ prc + er + OriginCode, ct_cyto))
by(ct_cyto$cyto, ct_cyto$prc, mean)

library(ggplot2)
ggplot(ct_cyto[ct_cyto$er == "er_pos",], aes(y=cyto,x=prc)) + geom_boxplot()
ct_cyto_subset = ct_cyto[ct_cyto$er == "er_pos",]
table(ct_cyto_subset$OriginCode, ct_cyto_subset$prc)

anova(lm(cyto ~ prc, ct_cyto_subset[ct_cyto_subset$OriginCode == "GPL1708_GSE22226",]))

ct_cyto_subset_subset = ct_cyto[ct_cyto$OriginCode == "GPL10558_GSE59515",]
by(ct_cyto_subset_subset$cyto, ct_cyto_subset_subset$prc, mean)

