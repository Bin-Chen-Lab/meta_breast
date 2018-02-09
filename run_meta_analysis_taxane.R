# First compile the 'MultiMeta-6.R' file - it has all the functions you need.
setwd("~/Documents/stanford/breast/meta_analysis_cleaned///")
# set working directory for all files - change to yours
dataDir<-"~/Documents/stanford/breast/meta_analysis_cleaned//"


source("MultiMeta-6.R")

#choose taxane related samples
ct = read.csv("BCRA_GEO_pCR_all_drug.txt", sep ="\t")
ct = ct[!is.na(ct$ClassCode), ]
ct$GSM = ct$gsm_name
ct$GSE = ct$gse_name
ct$GPL = ct$gpl_name
ct = ct[grep("drug_taxane|drug_docetaxel|drug_paclitaxel", ct$drug), ]
ct$subtype = sapply(1:nrow(ct), function(id){
  if (ct$er[id] == "er_pos"){
    "ER+"
  }else if (ct$her2[id] == "her2_positive"){
    "HER2+"
  }else if (ct$her2[id] == "her2_negative" & ct$er[id] == "er_neg" & ct$pr[id] == "pr_negative"){
    "TNBC"
  }else{
    "Others"
  }
})
ct = ct[ct$subtype != "Others", ]
ct$OriginCode = paste(ct$OriginCode, ct$subtype, sep="_")
valid_code = data.frame(table(ct$prc, ct$OriginCode))
valid_code = valid_code[valid_code$Freq > 4, ]
valid_code = names(table(valid_code$Var2)[table(valid_code$Var2) > 1])
ct = ct[ct$OriginCode %in% valid_code, ]
ct$Logged = NA

#GPL1390 does not work
studies = c("GPL1390") #GSE34138", "GSE22226",  "GSE4779",  "GSE8465", "GSE21997") #,  "GSE59515","GSE4779", "GSE22226", "GSE59515", "GSE25066") #"GSE8465",
ct = subset(ct, !GPL %in% studies)

mappingCol <- "GeneID"
groupingCol <- "OriginCode"
mappingSpecifier <- "Probe2EntrezMap"

# Get normalized expression values - this loops through all files in your directory and reads in the expression values, as well as maps the probes to entrez ids, and writes the mapping files back to your working directory. You will need to give it 'con', the database connection since it uses the database to do the gene mappings.
getRawDataCT(ct, dataDir, DEBUG = 3, noclobber=T)

# one of many functions that computes differential expression
computeExpressionSD(ct, dataDir, groupingCol=groupingCol, DEBUG=2)
runSAM(ct, dataDir,  groupingCol=groupingCol)
omnibusAcrossSAM(ct, dataDir, mappingSpecifier=mappingSpecifier, mappingCol=mappingCol, mname="HumanMeta-PooledTop") #-> omniSAM
# use this mapped across different dataset to generate p<0.05 and <0.01 as user-defined cutoff
omnibusAcrossSAMselectivelyMapped(ct, dataDir, mappingSpecifier=mappingSpecifier, mappingCol=mappingCol, groupingCol="OriginCode", mname="HumanMeta") #-> omniSAMselmapSGmap

# this does meta-analysis - it computes meta effect sizes and writes the results to a table
geneLevelMetaEffectSizeEstimate(ct, dataDir, mappingSpecifier=mappingSpecifier, mappingCol=mappingCol,  dataExtension="DifMeanLog-Welch",groupingCol="OriginCode",  resSpec="MappedMetaEffect", mname="Human", method="both", mappedName="Symbol", noclobber=TRUE, DEBUG=2)
metaEffectSizeEstimate(ct, dataDir=dataDir, mappingSpecifier=mappingSpecifier, mappingCol=mappingCol, groupingCol=groupingCol, mname="HumanMeta-PooledTop") 
metaEffectSizeEstimateOverMappedLevel(ct, dataDir, mappingSpecifier=mappingSpecifier, mappingCol=mappingCol,  dataExtension="DifMeanLog-Welch", groupingCol="OriginCode",  resSpec="MappedMetaEffect", mname="Human", noclobber=F, DEBUG=2)

source("meta_analysis_fisher_test.R")
source("compare_two_methods.R")
source("valid_meta_results.R")

biomagene = dz_sig$GeneID
eam <- createEffectsAttributesFromCT(ct, cfcol=c("OriginCode", "GSE", "GPL", "subtype"))
graphEffectsModelNamedExperiment(dataDir,gplID=eam$GPL, expID=eam$OriginCode, groupFactor=eam$subtype, 
                                 groupingCol="OriginCode", mappingSpecifier=mappingSpecifier, mappingCol=mappingCol, dataExtension="DifMeanLog-Welch", 
                                 method="random", nameList=eam$OriginCode ,estimate=TRUE, graphDir="forestplots-gpl", outName="gpl-taxane",  graphList=biomagene, outSize=9, DEBUG=5)

print("done!")


