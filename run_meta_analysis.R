# First compile the 'MultiMeta-6.R' file - it has all the functions you need.
setwd("~/Documents/stanford/breast/meta_analysis//")
source("MultiMeta-6.R")

# set working directory for all files - change to yours
dataDir<-"~/Documents/stanford/breast/meta_analysis//"

#choose taxane related samples
ct = read.csv("BCRA_GEO_pCR_all_drug.txt", sep ="\t")
ct = ct[!is.na(ct$ClassCode), ]
ct$GSM = ct$gsm_name
ct$GSE = ct$gse_name
ct$GPL = ct$gpl_name
ct = ct[grep("drug_taxane|drug_docetaxel|drug_paclitaxel", ct$drug), ]
ct$OriginCode = sapply(1:nrow(ct), function(id){
  if (ct$er[id] == "er_pos"){
    paste(ct$OriginCode[id], "_er_pos", sep = "")
  }else if (ct$her2[id] == "her2_positive"){
    paste(ct$OriginCode[id], "_her2_pos", sep = "")
  }else{
    paste(ct$OriginCode[id], "_TNBC", sep = "")
  }
})
valid_code = data.frame(table(ct$prc, ct$OriginCode))
valid_code = valid_code[valid_code$Freq > 3, ]
valid_code = names(table(valid_code$Var2)[table(valid_code$Var2) > 1])
ct = ct[ct$OriginCode %in% valid_code, ]

# coding table (excel spreadsheet)
ct_2 = readCodingTable("meta_input_new.txt")


#ct = subset(ct, GPL %in% c("GPL571", "GPL7504")) #GPL1390 has no Entrez Gene ID
#ct <- subset(ct, er_pos == "er_pos" & !is.na(ClassCode) )
studies = c("GPL1390") #GSE34138", "GSE22226",  "GSE4779",  "GSE8465", "GSE21997") #,  "GSE59515","GSE4779", "GSE22226", "GSE59515", "GSE25066") #"GSE8465",

ct = subset(ct, !GPL %in% studies)
#table(ct$GSE, ct$ClassCode)

mappingCol <- "GeneID"
groupingCol <- "OriginCode"
mappingSpecifier <- "Probe2EntrezMap"

#ge DB connection
#source("/home/l3li/script/R/DBcon.R")
#library("RMySQL")
#mysql_drvr <-dbDriver("MySQL")
#con <- dbConnect(mysql_drvr,group="client",dbname="proj_nosology")

# Get normalized expression values - this loops through all files in your directory and reads in the expression values, as well as maps the probes to entrez ids, and writes the mapping files back to your working directory. You will need to give it 'con', the database connection since it uses the database to do the gene mappings.
#con <- dbConnect(dbDriver("MySQL"), dbname="user_l3li",username="l3li",  host="buttelab-db1", password='');
getRawDataCT(ct, dataDir, DEBUG = 3, noclobber=T)

# one of many functions that computes differential expression
computeExpressionSD(ct, dataDir, groupingCol=groupingCol, DEBUG=2)
runSAM(ct, dataDir,  groupingCol=groupingCol)
omnibusAcrossSAM(ct, dataDir, mappingSpecifier=mappingSpecifier, mappingCol=mappingCol, mname="HumanMeta-PooledTop") -> omniSAM
# use this mapped across different dataset to generate p<0.05 and <0.01 as user-defined cutoff
omnibusAcrossSAMselectivelyMapped(ct, dataDir, mappingSpecifier=mappingSpecifier, mappingCol=mappingCol, groupingCol="OriginCode", mname="HumanMeta") -> omniSAMselmapSGmap


# this does meta-analysis - it computes meta effect sizes and writes the results to a table
geneLevelMetaEffectSizeEstimate(ct, dataDir, mappingSpecifier=mappingSpecifier, mappingCol=mappingCol,  dataExtension="DifMeanLog-Welch",groupingCol="OriginCode",  resSpec="MappedMetaEffect", mname="Human", method="both", mappedName="Symbol", noclobber=TRUE, DEBUG=2)
metaEffectSizeEstimate(ct, dataDir=dataDir, mappingSpecifier=mappingSpecifier, mappingCol=mappingCol, groupingCol=groupingCol, mname="HumanMeta-PooledTop") -> mfc.OC
#metaEffectSizeEstimate(ct, dataDir, mappingSpecifier="Probe2EntrezMap", mappingCol="GeneID", groupingCol="OriginCode", DEBUG=3 )
metaEffectSizeEstimateOverMappedLevel(ct, dataDir, mappingSpecifier=mappingSpecifier, mappingCol=mappingCol,  dataExtension="DifMeanLog-Welch", groupingCol="OriginCode",  resSpec="MappedMetaEffect", mname="Human", noclobber=F, DEBUG=2)

source("meta_analysis_fisher_test.R")
source("compare_two_methods.R")
source("valid_meta_results.R")
print("done!")
q()

# this is the table written by the previous function.  You can read it in, and set thresholds for meta effect size, p-value, whatever you want
#mese <- read.table("/home/l3li/projects/NKT/meta_lupus/9set/effectsMeta/MappedMetaEffect.Probe2EntrezMap-GeneID-OriginCodeFandROverMappedBelow.table")
#mese <- read.table("effectsMeta/MappedMetaEffect.Probe2EntrezMap-GeneID-OriginCodeFandROverMappedBelow.table")

# example of thresholds, can be self defined
top.genes <- mese[mese$RandomEffectsP < .01 & mese$I2rem < .2 & abs(mese$RandomEffects) >= .5,]
#################################################
# graph for GPL
#################################################
#biomagene=c(3550)
biomagene = c(960, 1021)
#biomagene=c(3550, 6992,9349,6160,6154,64601,4793,8542,5204,83693,27132,5880,5720,5528,80344,84153,916,1991,6037,8993,6944,7388,54543,50640,9343,81928,79895,26147,348262,27005,4535,9881,252950,408050,728358,441272,653121,160428,3519)
#biomagene2=c(3429,10964,129607,3434,91543,51191,3437,11274,94240,10561,4938,54739,8638,4939,
  #           26010,5359,64135,64108,8519,3959,6241,2633,8743,91351,51056,9111,54809,983,80055,3669,4173,58472)
#biomagene2<- mese[mese$RandomEffectsP < .02 & mese$I2rem < .2 & mese$RandomEffects >= 1,]
# Make graphs
#test=top.genes$ID
eam <- createEffectsAttributesFromCT(ct, cfcol=c("OriginCode", "GSE", "GPL"))

# this makes the forest plots, and puts them into a folder by themselves

#groupFactor will be the legend, corresponding to the cfcol in eam e.g, GSE
graphEffectsModelNamedExperiment(dataDir,gplID=eam$GPL, expID=eam$OriginCode, groupFactor=eam$GPL, 
groupingCol="OriginCode", mappingSpecifier=mappingSpecifier, mappingCol=mappingCol, dataExtension="DifMeanLog-Welch", 
method="random", nameList=eam$OriginCode ,estimate=TRUE, graphDir="forestplots-gpl", outName="gpl-dengue",  graphList=biomagene, outSize=9, DEBUG=5)

#################################################
# graph for cell types
#################################################
biomagene=c(3550, 6992,9349,6160,6154,64601,4793,8542,5204,83693,27132,5880,5720,5528,80344,84153,916,1991,6037,8993,6944,7388,54543,50640,9343,81928,79895,
26147,348262,27005,4535,9881,252950,408050,728358,441272,653121,160428,3519)

#biomagene2<- mese[mese$RandomEffectsP < .02 & mese$I2rem < .2 & mese$RandomEffects >= 1,]
# Make graphs
#test=top.genes$ID
eam <- createEffectsAttributesFromCT(ct, cfcol=c("OriginCode", "GSE", "Celltype"))

# this makes the forest plots, and puts them into a folder by themselves

#groupFactor will be the legend, corresponding to the cfcol in eam e.g, GSE, nameList is the y axis
graphEffectsModelNamedExperiment(dataDir,gplID=eam$GPL, expID=eam$OriginCode, groupFactor=eam$Celltype, 
groupingCol="OriginCode", mappingSpecifier=mappingSpecifier, mappingCol=mappingCol, dataExtension="DifMeanLog-Welch", 
method="random", nameList=eam$OriginCode ,estimate=TRUE, graphDir="forestplots-celltype", outName="celltype-dengue",  graphList=biomagene, outSize=9, DEBUG=5)

#################################################
# graph for populations
#################################################
biomagene=c(3550, 6992,9349,6160,6154,64601,4793,8542,5204,83693,27132,5880,5720,5528,80344,84153,916,1991,6037,8993,6944,7388,54543,50640,9343,81928,79895,
26147,348262,27005,4535,9881,252950,408050,728358,441272,653121,160428,3519)
#biomagene2<- mese[mese$RandomEffectsP < .02 & mese$I2rem < .2 & mese$RandomEffects >= 1,]
# Make graphs
#test=top.genes$ID
eam <- createEffectsAttributesFromCT(ct, cfcol=c("OriginCode", "GSE", "Population"))

# this makes the forest plots, and puts them into a folder by themselves

#groupFactor will be the legend, corresponding to the cfcol in eam e.g, GSE, nameList is the y axis
graphEffectsModelNamedExperiment(dataDir,gplID=eam$GPL, expID=eam$OriginCode, groupFactor=eam$Population, 
groupingCol="OriginCode", mappingSpecifier=mappingSpecifier, mappingCol=mappingCol, dataExtension="DifMeanLog-Welch", 
method="random", nameList=eam$OriginCode ,estimate=TRUE, graphDir="forestplots-population", outName="pop-dengue",  graphList=biomagene, outSize=9, DEBUG=5)

#################################################
# graph for cell types and y label is yname
#################################################
biomagene2=c(3429,10964,129607,3434,91543,51191,3437,11274,94240,10561,4938,54739,8638,4939,
             26010,5359,64135,64108,8519,3959,6241,2633,8743,91351,51056,9111,54809,983,80055,3669,4173,58472)
#biomagene2<- mese[mese$RandomEffectsP < .02 & mese$I2rem < .2 & mese$RandomEffects >= 1,]
# Make graphs
#test=top.genes$ID
eam <- createEffectsAttributesFromCT(ct, cfcol=c("OriginCode", "GSE", "Celltype", "Yname"))

# I think this makes the forest plots, and puts them into a folder by themselves

#groupFactor will be the legend, corresponding to the cfcol in eam e.g, GSE, nameList is the y axis
graphEffectsModelNamedExperiment(dataDir,gplID=eam$GPL, expID=eam$OriginCode, groupFactor=eam$Celltype, 
groupingCol="OriginCode", mappingSpecifier=mappingSpecifier, mappingCol=mappingCol, dataExtension="DifMeanLog-Welch", 
method="random", nameList=eam$Yname ,estimate=TRUE, graphDir="forestplots-Filtered", outName="-sle",  graphList=biomagene2, outSize=9, DEBUG=5)

#graphEffectsModelNamedExperiment(dataDir,gplID=eam$GPL,expID=eam$OriginCode, groupFactor=eam$OriginCode,  
#groupingCol="OriginCode", mappingSpecifier=mappingSpecifier, mappingCol=mappingCol, dataExtension="DifMeanLog-Welch", 
#method="random", nameList=eam$OriginCode ,estimate=TRUE, graphDir="forestplots-Selected", outName="-SLE", graphList=genesToPlot, outSize=9, DEBUG=5)

