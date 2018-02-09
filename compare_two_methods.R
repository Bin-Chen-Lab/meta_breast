mese <- read.delim("effectsMeta/MappedMetaEffect.Probe2EntrezMap-GeneID-OriginCodeFandROverMappedBelow.table",stringsAsFactors=F)

library(qvalue)
mese$q_values <- qvalue(mese$RandomEffectsP)$qvalues

top.genes.effect <- mese[mese$q_values < .05  & abs(mese$RandomEffects) >= 0.5,]
top.up.genes.effect <- mese[mese$q_values < .05  & (mese$RandomEffects) >= 0.5 ,]
top.down.genes.effect <- mese[mese$q_values < .05  & (mese$RandomEffects) <= -0.5,]

samr.meta <- read.delim("omnibusMeta/FisherAndLogit-SAM-OriginCode-SelectedMapped-GeneID-HumanMeta.metaPVal.fisher.table",header=T,stringsAsFactors=F)
samr.meta <- samr.meta[!is.na(samr.meta$pvalFisher),]
top.genes.samr <- samr.meta[samr.meta$pvalFisher<0.05,]
top.up.genes.samr <- top.genes.samr[top.genes.samr$a > top.genes.samr$b,]
top.down.genes.samr <- top.genes.samr[top.genes.samr$a < top.genes.samr$b,]

common.genes <- merge(top.genes.effect,top.genes.samr, by.x="ID",by.y="ID")
common.up.genes <- merge(top.up.genes.effect,top.up.genes.samr, by.x="ID",by.y="ID")
common.down.genes <- merge(top.down.genes.effect,top.down.genes.samr, by.x="ID",by.y="ID")


if (nrow(common.up.genes) > 0 & nrow(common.down.genes) >0){
  common.up.genes$up_down <- "up"
  common.down.genes$up_down <- "down"
  dz_sig <- rbind(common.up.genes,common.down.genes)
}else if (nrow(common.up.genes) > 0) {
  common.up.genes$up_down <- "up"
  dz_sig = common.up.genes
}else if (nrow(common.down.genes) > 0) {
  common.down.genes$up_down <- "down"
  dz_sig = common.down.genes
}

dz_sig <- subset(dz_sig,select=c("ID","Name.x","RandomEffects","q_values","p_values","up_down"))
names(dz_sig) <- c("GeneID","Symbol","fold.change","q_value","p_value","up_down")
write.table(dz_sig,"dz_sig_meta_analysis.txt",sep="\t",row.names=F,col.names=T,quote=F)
            
