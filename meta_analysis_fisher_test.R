results <- read.table("omnibusMeta/FisherAndLogit-SAM-OriginCode-SelectedMapped-GeneID-HumanMeta.metaPVal.table",header=T,sep="\t")
totalPositive <- sum(results$q.1)
totalPos_Neg <- sum(results$totMeas)

as <-NULL
bs <- NULL
cs <- NULL
ds <- NULL
p_values <- NULL

for (i in which(results$totMeas>0)){
  #print (i)
  a <- results$q.1[i]
  b <- results$totMeas[i] - a
  c <- totalPositive - a
  d <- totalPos_Neg - a - b - c
  
   p <- fisher.test (matrix(c(a,b,c,d),nrow=2,byrow=T))
  
  as <- c(as,a)
  bs <- c(bs,b)
  cs <- c(cs,c)
  ds <- c(ds,d)
  p_values <- c(p_values,p$p.value)
}

results <- cbind(results, as)
results <- cbind(results,bs)
results <- cbind(results,cs)
results <- cbind(results,ds)
results <- cbind(results,p_values)

write.table(results,"omnibusMeta/FisherAndLogit-SAM-OriginCode-SelectedMapped-GeneID-HumanMeta.metaPVal.fisher.table",sep="\t",row.names=F,col.names=T,quote=F)

