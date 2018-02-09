

library(RankProd)
library(genefilter)
#library(annotate)
library(RMySQL)  # Use dbGetQuery()  !!!
library(RColorBrewer)
library(preprocessCore) # include quantile normalization
library(rmeta)
library(gplots)
library(siggenes)
library(GEOquery)
library(ROCR)
library(meta)
#library(multicore)

## Fix the groupingCol issue, when does it make sense?


options(stringsAsFactors = FALSE)


################################
# General Utility Functions

stripNumberPrefix <-function(v, prefix) {
  # Fix this to only work if there is a prefix!  Right now, will error!
  # Definitely need to fix this up!!
  
	#return(as.numeric(as.matrix(data.frame(strsplit(v, paste("^", prefix, sep="") ))[2,] )))
        #use GSUB instead?  
  return(as.matrix(gsub(paste("^", prefix, sep="" ), "", v)))
}


stripPrefix <-function(v, prefix) {
  # Fix this to only work if there is a prefix!  Right now, will error!
  # Definitely need to fix this up!!  
	return(as.matrix(gsub(paste("^", prefix, sep="" ), "", v)))
}

stripPostfix <-function(v, postfix) {
  # Fix this to only work if there is a prefix!  Right now, will error!
  # Definitely need to fix this up!!
  #spd <- gsub(paste(postfix, "$", sep=""), "", v)
  #return(spd)
  return(as.matrix(gsub(paste(postfix, "$", sep="" ), "", v) ))
}

mergeViaRowNames <- function(a, b, aName="", bName="", all.a=F, all.b=F) {
	# This function takes in two matrices and then merges them by their rownames
	# It can take in optional arguments that pre-pend text to the column names of the parts, useful for keeping track of them
	a <- data.frame(a)
	b <- data.frame(b)
	if (aName != "") {colnames(a) <- paste(aName, colnames(a), sep="")}
	if (bName != "") {colnames(b) <- paste(bName, colnames(b), sep="")}
	a <- data.frame(mindex=rownames(a),a )
	b <- data.frame(mindex=rownames(b), b)
	m <- merge(a, b, by.x="mindex", by.y="mindex", all.x=all.a, all.y=all.b)
	rownames(m) <- m[,1]
	m <- m[,2:(dim(m)[2])]
	return(m)
}

cnr <- function(a, x=4, y=8) {
  # Show the left corner of the data
  return(a[1:min(x, dim(a)[1]), 1:min(y, dim(a)[2])])
}


extractTopRow <- function(m) {
	return(m[1,])
}

escapeText <- function(x, y=".") {
  # Escapes multiple bad characters in a row with a period (or whatever), contrast with escapeTextEach which escapes each occurrence
  gsub("[[:punct:][:blank:]]+", y, x)
}

escapeTextEach <- function(x, y=".") {
  # Escapes each occurrence of a bad text character with a period (as opposed to escaping the whole block)
  gsub("[[:punct:][:blank:]]", y, x)
}


################################
# Statistical Functions taking from Rob Tibshirani's SAMR
## Uses class labels 1 & 2, will switch to that (makes more sense and agilent can be 0 or negative)

est.s0<-function(tt,sd,s0.perc=seq(0,.99,len=40)){
	# estimate s0 (exchangeability) factor for denominator.
	# returns the actual estimate s0 (not a percentile)

	a<-cut(sd,breaks=quantile(sd,seq(0,1,len=101)),labels=F)
	a[is.na(a)]<-1
	cv.sd<-rep(0,length(s0.perc))

	for(j in 1:length(s0.perc)){
 		w<-quantile(sd,s0.perc[j])
 		w[j==1]<-0
 		tt2<-tt*sd/(sd+w)
 		sds<-rep(0,100)
		
    	for(i in 1:100){
      		sds[i]<-mad(tt2[a==i])
   		}

   		cv.sd[j]<-sqrt(var(sds))/mean(sds)
	}
	o=(1:length(s0.perc))[cv.sd==min(cv.sd)]

	s0.hat=min(sd[o])
	return(list(s0.perc=s0.perc,cv.sd=cv.sd, s0.hat= s0.hat))
}

mttest.func <- function(x,y,s0=0){
	# Modified t-test
	# x is a vector of values, y is the class labels

	n1 <- sum(y==1)
	n2 <- sum(y==2)

  	p <- nrow(x)
    m1 <- rowMeans(x[,y==1])
    m2 <- rowMeans(x[,y==2])

 	sd <- sqrt( ((n2-1) * varr(x[, y==2], meanx=m1) + (n1-1) * varr(x[, y==1], meanx=m2) )*
          	(1/n1+1/n2)/(n1+n2-2) )
	numer <-  m2 - m1
   	dif.obs <- (numer)/(sd + s0)
	return(list(tt=dif.obs,numer=numer, sd=sd))
}

varr <- function(x, meanx=NULL){
  n <- ncol(x)
  p <- nrow(x)
   Y <-matrix(1,nrow=n,ncol=1)
     if(is.null(meanx)){   meanx <- rowMeans(x)}
        ans<- rep(1, p)
                xdif <- x - meanx %*% t(Y)
                ans <- (xdif^2) %*% rep(1/(n - 1), n)
                ans <- drop(ans)
       return(ans)
}


computeCV <- function(x) {
	sdx <- sd(x, na.rm=TRUE)
	xbar <- abs(mean(x, na.rm=TRUE))
	if (is.na(sdx) || is.na(xbar) || xbar==0 ) {return(NA)} 
	return(sdx/xbar)
}



################################
# File IO

readCodingTable <- function(infile, DEBUG=0) {
  # Coding table needs to have a ClassCode (1 or 2), an Origin Code (some unique term for its identifier), GSM (particular sample ID), and GPL (the GEO platform ID)
	ct <- read.csv(infile, sep="\t", as.is=T)
	if (is.null(ct$ClassCode) | is.null(ct$OriginCode) | is.null(ct$GSM) | is.null(ct$GPL) )  {
		print("ERROR!!  Problem with coding table, needs to have columns for ClassCode (1/2 for treatment/control or a zero, negative or NA for a two color), OriginCode, GSM, and GPL.  This information is all essential")
		return(NA)		
	}
	if (is.null(ct$Logged)) {ct$Logged <- NA}
        
	return(ct)
}

getStoredTable <- function(dataDir, dataType, specifier, extension, DEBUG=0) {
	# This gets a stored data table, which is in the data directory, is of a certain type, has a specifier, and an extension
	# Usually things are stored like :  dataDir/dataType/spcifier.extension.table
	# dataTypes are usually things like 'rawData', 'effectSize', 'foldChange' or whatever
	# specifier could be something like 'GPL96', 'K' (for origin code), 'ratMuscle', etc.
	# extension could be something like 'fisherMethod', 'logitMethod'
	# This is to allow future versions to use SQL or some other serialization
	infile <- paste(c(dataDir, "/", dataType, "/", specifier, ".", extension, ".table" ), collapse="")
	if (DEBUG) {print(paste("reading:", infile))}
	st <- read.csv(infile, sep="\t", as.is=T)
}

writeStoredTable <- function(dtab, dataDir, dataType, specifier, extension, DEBUG=0, ...) {
	outpath <- paste(c(dataDir, "/", dataType), collapse="")
	if (DEBUG > 2) print(outpath)
        if (DEBUG) {dir.create(outpath, recursive=T, showWarnings=T)} else {dir.create(outpath, recursive=T, showWarnings=F)}
        if (! is.data.frame(dtab) & ! is.matrix(dtab)) {
          print("WARNINGS, call to write out something not a dataframe and not a matrix; Returning from call")
          print(head(dtab))
          warning("WARNINGS, call to write out something not a dataframe and not a matrix; Returning from call")
          return()
        }
	outfile <- paste(c(dataDir, "/", dataType, "/", specifier, ".", extension, ".table" ), collapse="")
	print(paste("writing:", outfile))
	write.table(as.matrix(dtab), outfile, sep="\t", quote=F, col.names=T, ...)
	
}

checkStoredTable <- function(dataDir, dataType, specifier, extension, DEBUG=0) {
  infile <- paste(c(dataDir, "/", dataType, "/", specifier, ".", extension, ".table" ), collapse="")
  fc <- file.exists(infile)
  if (DEBUG) {
    print(c("File Exists:", fc))
    print(c(infile))}
  return(fc)
}


getGSMfromFile <- function(gsmid, filename, probecol="probe", DEBUG=0) {
  if (DEBUG) print(filename)
  ifo <- read.csv(filename, sep="\t", as.is=T,check.names=F)
   if (DEBUG) print(gsmid)
   # assume there is a column name
   #print(ifo[1:4, 1:4])
   rf <- data.frame(gsm=gsmid, id=ifo[,probecol], val=ifo[,gsmid])
  if (DEBUG > 2) print(head(rf))
   return(rf)
   
 }

################################
# Butte Lab DB Interactions


db <- function() {
 # try(dbDisconnect(con))
  #source("/home/l3li/projects/NKT/meta_lupus/R/DBcon.R")
  library("RMySQL")
  mysql_drvr <-dbDriver("MySQL")
  con <- dbConnect(mysql_drvr,group="client",host="buttelab-db1.stanford.edu",dbname="proj_nosology")
  
}


storeHomologeneMappings <- function(ct, dataDir, con, DEBUG=0) {
	for (gpl in unique(ct$GPL)) {
        
		print(gpl)
		try({
                  egm <- getStoredTable(dataDir, "geneMappings", gpl, "Probe2EntrezMap", DEBUG=DEBUG)
		#print(head(egm)) 
                  query <- paste("select HID, GeneID from annot_gene.homologene where GeneID in (", paste(as.character(egm$GeneID), collapse=", "), ")", sep="" )
		#print(query)
                  entrez2Homologene <- dbGetQuery(con, query)
		#print(head(entrez2Homologene))
                  h2p <- merge(egm, entrez2Homologene, by.x="GeneID", by.y="GeneID", all.x=F, all.y=F)
                                        #print(c(dim(egm), dim(entrez2Homologene), dim(h2p)))
                  if (DEBUG) {print(head(h2p))}
                  dupmid <- unique(h2p$probe[duplicated(h2p$probe)])
                  h2p <- h2p[! h2p$probe %in% dupmid, c("probe", "HID", "GeneID", "Symbol")]
                  rownames(h2p) <- h2p$probe
                  writeStoredTable(h2p, dataDir, "geneMappings", gpl, "Probe2HID", DEBUG=DEBUG)
                }
                    )
	}
}

getEntrezMappings <- function(gpl, con, DEBUG=0) {
	query <- paste("select probe, GeneID, Symbol from annot_gpl.gpl_", stripNumberPrefix(gpl, "GPL"), sep="" )
	if (DEBUG) {print(query)}
	GeneMappings <- dbGetQuery(con, query)
	dupProbes <- unique(GeneMappings$probe[duplicated(GeneMappings$probe)])
	GeneMappings <- GeneMappings[! GeneMappings$probe %in% dupProbes,]
	rownames(GeneMappings) <- GeneMappings$probe
	return(GeneMappings)	
}

getGSMdata <- function(gsmList, con, logged, gsource="database", idcol="probe", DEBUG=0) {
  ## All hacky and messed up for getting material from files and so forth, need to fix this up, assume for files that there is a file that just has everything 
	expressionList <-  NULL
	#gsmList <- stripNumberPrefix(gsmList, "GSM")  # Need to do this later, only works when there is a GSM!
        #print(data.frame(gsmList, logged))
        #print(typeof(logged))
	for (i in 1:length(gsmList)) {
           if (gsource[i] == "database" | gsource[i] == "expr_geo") {
             if (DEBUG > 2) {print(gsource[i])}
             query <- paste("select gsm, id, val from expr_geo.gsm_gene_data where gsm =", stripNumberPrefix(gsmList[i], "GSM"), sep="" )
             expressionListPart <- dbGetQuery(con, query)
             if (dim(expressionListPart)[1] < 1) {
			print(gsmList[i])
			print("Error, empty db return")
                      } else {
                        expressionListPart$gsm <- paste("GSM", expressionListPart$gsm, sep="")
                      }
           } else {
             if (gsource[i] =="GEO") {
               if (DEBUG > 2) {print(paste("getGEO",gsmList[i], collapse=" "))}
               geores <- getGEO(gsmList[i])
               #print(geores)
               expressionListPart <- data.frame(gsm=gsmList[i], id=as.character(Table(geores)[,"ID_REF"]), val=as.numeric(Table(geores)[,"VALUE"]))
               if (DEBUG > 2) {print(head(expressionListPart))}
               
             } else {
               # Try to read in a file
               #print(c(i, gsmList[i], gsource[i]))
               #print(paste("Trying to load: ", gsmList[i], " from: ", gsource[i], collapse=" "))
               expressionListPart <- c()
               tryCatch( {expressionListPart <- getGSMfromFile(gsmList[i], gsource[i])}, message=paste("Errors while trying to load: ", gsource[i], " to locate the expression data for: ", gsmList[i], collapse=" "), DEBUG=DEBUG)
             }

             
           }
           if (dim(expressionListPart)[1] > 0)
                 {
                  # If missing logged information, need to guess based on range of values
                   #print(summary(expressionListPart))
                   #print(max(expressionListPart$val))
                   #print(min(expressionListPart$val))
                 if (is.na(logged[i]) ) {
                  	if ( (max(expressionListPart$val, na.rm=TRUE) > 100) | (min(expressionListPart$val, na.rm=TRUE) < -50)) {
                          if (min(expressionListPart$val, na.rm=TRUE) < 1) {
				expressionListPart$val <- expressionListPart$val + abs(min(expressionListPart$val, na.rm=TRUE)) + 1
				}
			expressionListPart$val <- log2(as.numeric(expressionListPart$val))
                        }
                } else {
		if (!logged[i]) {
			if (max(expressionListPart$val, na.rm=TRUE) < 20) {
				print(query)
				print(data.frame(i=i, GSM=gsmList[i], logged=logged[i], max=max(expressionListPart$val)))
				print("Potential Error, check to make sure the 'Logged' column is accurate.  The max log2(expression) < 20, which seems unlikely with unlogged data, particularly if Affymetrix.")
			}
			if (min(expressionListPart$val, na.rm=TRUE) < 1) {
				expressionListPart$val <- expressionListPart$val + abs(min(expressionListPart$val, na.rm=TRUE)) +90
				}
			expressionListPart$val <- log2(as.numeric(expressionListPart$val))
		} else {
			#print(data.frame(i=i, GSM=gsmList[i], logged=logged[i]))
			if (max(expressionListPart$val, na.rm=TRUE) > 50 | min(expressionListPart$val, na.rm=TRUE) < -50) {
				print(c(gsmList[i], max(expressionListPart$val, na.rm=TRUE), min(expressionListPart$val, na.rm=TRUE) ))
				print("Potential Error, check to make sure the 'Logged' column is accurate.  It looks like there is a log2(expression) magnitude over 50 which seems infeasible.")
			}
		}}
                
		if (is.null(expressionList)) {expressionList <- expressionListPart}
		else {expressionList <- rbind(expressionList, expressionListPart)}
		}
	}
        if (DEBUG > 1) { print(head(expressionList))}
        if (dim(expressionList)[1] > 10) {
          emat <- reshape(expressionList, v.names="val", idvar="id", timevar="gsm", direction="wide")
          expressionMatrix <- normalize.quantiles(as.matrix(emat[,2:dim(emat)[2]]))
          rownames(expressionMatrix) <- emat$id
                                        #print(head(emat))	
          cn <- stripPrefix(colnames(emat[,2:dim(emat)[2]]), "val.")

          colnames(expressionMatrix) <- cn
          return(expressionMatrix)
        } else {return(NULL)}
}


getHomologeneSpeciesMap <- function(taxid, con, DEBUG=0) {
		query <- paste("select HID, GeneID, Symbol from annot_gene.homologene where tax_id=", taxid, sep="" )
		entrez2Homologene <- dbGetQuery(con, query)
                return(entrez2Homologene)
}





################################
#  Stage 1, getting the data

#  View simple version
makeExperimentSummaryTable <- function(ct, cols, DEBUG=0) {
	# Takes in a coding table, and returns the columns from the table where the OriginCode is unique
  #  Add in some other things that do things like calculate summary statistics by the groupingColumn factor
#  TODO!
  
  mr <- by(ct[,cols], ct$OriginCode, extractTopRow)
  return(do.call("rbind", mr))
}


getRawDataCT <- function(ct, dataDir, getEntrez=TRUE, DEBUG=0, noclobber = TRUE) {
	# This gets the expression data and the entrex mappings
	# codingTableFile provides all ht einformation about what is going on
	# dataDir is the working directory where files will be stored
	# con is the DB connection
	dups <- ct$GSM[duplicated(ct$GSM)]
	if (length(dups) > 0) {
		print("Duplicated GSM entries")
		print(dups)
	}
	for (gpl in unique(ct$GPL)) {
		if (DEBUG > 0) {print(gpl)}
		ctS <- ct[ct$GPL == gpl & !is.na(ct$OriginCode), c("GSM", "ClassCode", "Logged", "OriginCode", "Source") ]
		print(ctS[1:2, 1:4])
		if (nrow(ctS) >= 4) {
                  if (!noclobber | !checkStoredTable(dataDir, "rawExpression", gpl, "RawExpressionValues", DEBUG=DEBUG)) {
                    if(DEBUG) {print(c("Getting Data for gpl:", gpl, noclobber,checkStoredTable(dataDir, "rawExpression", gpl, "RawExpressionValues") ))}
                    #print(ctS$GSM)
                    eMat <- getGSMdata(ctS$GSM, con, ctS[, "Logged"], ctS[, "Source"], DEBUG=DEBUG)
                    if(!is.null(eMat)){writeStoredTable(eMat, dataDir, "rawExpression", gpl, "RawExpressionValues", DEBUG=DEBUG)}
                  }
                  
                  if (getEntrez & (!noclobber | !checkStoredTable(dataDir, "geneMappings", gpl, "Probe2EntrezMap"))) { 
                    egm <- getEntrezMappingsFromGEO(gpl) #getEntrezMappings(gpl, con) 
			              writeStoredTable(egm, dataDir, "geneMappings", gpl, "Probe2EntrezMap", DEBUG=DEBUG)
		                }
		} else {print("Error for this GPL, less than 4 rows of data for it.  Skipping it, and the information needs te be removed from the Coding Table.  Otherwise, there will be errors trying to analyze this data.")}
	}
	#return(ct)

}




################################
#  Simple significance tests

simpleTT <- function(diffline, cl, alternative, DEBUG=0) {
  #  If there is any error, the results are suppressed and a .5 is returned for a one sided test and a .99999999 for a two sided test
  #  Usually the errors are for a constant value or not enough values, etc.  This is the same as a non-informative result.  At this point, better to have a value that says the test did not show anything than anything else.  Something a fudge factor, but maybe close enough
	if (sum(!is.na(diffline[cl==1])) > 1  & sum(!is.na(diffline[cl==2])) > 1 ) {
           tryCatch({         
             t <- t.test(diffline[cl==2], diffline[cl==1], alternative=alternative, na.action='na.omit')
             return(t$p.value)
           }, error = function(e) {
             #print(e)
             # Suppress printing errors, it happens too much
             ifelse(alternative=="two.sided", return(0.9999999999), return(0.5))
           })
         }
	else {ifelse(alternative=="two.sided", return(0.9999999999), return(0.5))}
}


simpleSignificanceTest <- function(ct, dataDir, method, groupingCol="OriginCode", DEBUG=0, noclobber=T) {	
	# This function runs a significance test for each experiment, it groups all across GPL for rankproduct or does each origin separately for t-test and modified t-test

	# method={"TT"}   for method of rank products, t-test, will include other tests later, maybe modified t-test, etc.
	# extension is if I want to use normalized data
	# writes out the up/down p-values for each one sided hypothesis test
	
	### For each gpl, open up the raw data file, for each origin, do the test down each experiment
	for (gpl in unique(ct$GPL)) {
		print(gpl)
		ctS <- ct[ct$GPL == gpl & !is.na(ct[,groupingCol]), c("GSM", "ClassCode", groupingCol) ]
		eMat <- getStoredTable(dataDir, "rawExpression", gpl, "RawExpressionValues")
		for (oc in unique(ctS[,groupingCol])) {
		  if (!noclobber | !checkStoredTable(dataDir, "significanceValues", paste(gpl, oc, groupingCol, method,  sep="-"), "PValues")) {
                  tryCatch({         
			ctC <- ctS[ctS[,groupingCol] == oc, ]
			gsmList <- ctC$GSM
			sm <- as.matrix(eMat[,gsmList])
			cl <- ctC$ClassCode
                        if (sum(cl != 1 & cl != 2) >= 1) {
                          print("Potential problem with class coding, expecting 1 & 2")
                        }
			#print(cl)
			if (method=="TT") {sf <- simpleTT}
			up <- apply(sm, 1, sf, cl=cl, alternative="greater" )
			dp <- apply(sm, 1, sf, cl=cl, alternative="less" )
			tp <- apply(sm, 1, sf, cl=cl, alternative="two.sided" )
			if(length(up) != length(dp)) 
        {print("Error, the length of the up p-values does not equal the down p-values.")} 
      else {
                          
        if(length(up) == sum(is.na(up))) {
          print(c("Error, the values are all NA for up", gpl, oc))
        } else {
          
          pv <- data.frame(pvalUp=up, pvalDown=dp, pvalTwoSided=tp)
          writeStoredTable(pv, dataDir, "significanceValues", paste(gpl, oc, groupingCol, method, sep="-"), "PValues")
        }
          }
                    }, error = function(e) {
                        print(c(gpl, oc))
                        print(str(e))
                         }
                           )
		}
	}
	}
	
}

computeRankProdSig <- function(ct, dataDir, groupingCol="OriginCode", outName="", DEBUG=0, noclobber=TRUE) {
	for (gpl in unique(ct$GPL)) {
		print(gpl)
		ctS <- ct[ct$GPL == gpl & !is.na(ct[,groupingCol]), c("GSM", "ClassCode", groupingCol) ]
                #print(ctS)
		if (max(ctS$ClassCode) == 2 & min(ctS$ClassCode) == 1) {
			if (!noclobber | !checkStoredTable(dataDir, "significanceValues", gpl, paste("RankProd-PValues", outName, sep=""))) {
				eMat <- getStoredTable(dataDir, "rawExpression", gpl, "RawExpressionValues")
                                #print(head(eMat))
                                print(eMat[1:3, 1:4])
				ofac <- factor(ctS[,groupingCol])
				olabels <- as.numeric(ofac)
				classLabels <- ctS$ClassCode - 1   # this is the 1,2 case, turn to 0 & 1 for RP
                                if (sum(classLabels != 0 & classLabels != 1) >= 1) {
                                  print("Potential problem with class coding, expecting 1 & 2")
                                }
                                
				#RP.out <- RP(eMat[,ctS$GSM], classLabels,olabels, num.perm = 100, logged = T, na.rm = TRUE, plot = FALSE, rand = 123)
				
				RP.out <- RPadvance(eMat[,ctS$GSM], classLabels,olabels, num.perm=100, logged=TRUE, rand=123 )
				rptab <- cbind(RP.out$pval, RP.out$AveFC, RP.out$RPrank, RP.out$pfp )
				colnames(rptab) <- c("pvalDown", "pvalUp", "AveFC", "RPRDown", "RPRUp", "PFPDown", "PFPUp")
				writeStoredTable(rptab, dataDir, "significanceValues", paste(gpl, "RankProd", sep="-"), "PValues")
				
				}
			}
		# Other cases are for two color arrays, but these I guess should be coded with -1
	}
}




################################
#  Stage 3, Omnibus statistical procedures

calcOmnibus <- function(pvList, method, minp=1e-10, DEBUG=0) {
	# takes in a p-value list and computes the mea p-values
		# method={'fisher'|'logit'}
	# minp is a little fudge factor to help when there is a zero p-value, keeps from taking log(0)
        pvList <- pvList[!is.na(pvList)]
	pvList <- pvList + minp
	pvList[pvList >= 1.0] <- .99999999  # keep it from getting to be one either

	if (method =="logit") {
		L <- sum(log(pvList/(1-pvList)))
		k <- length(pvList)
		Lstar <- abs(L)*sqrt(.3*(5*k+4)/(k*(5*k+2)))
		mp  <- 2*pt(Lstar, 5*k + 4, lower.tail=F)  # get the p-value from the T-distribution for 5k+4 degrees of freedom
	}
	if (method =="fisher") {
			ss<- sum(-2*log(pvList))
			mp <- pchisq(ss, df=length(pvList)*2, lower.tail=FALSE)
	}
	if (method=="mean") {
		
		mp <- mean(pvList)
	}
	return(mp)
}


omnibusAcrossMappedGrouped <- function(specList, platList, mappingSpecifier, dataDir, method, extension = "PValues", mappingCol,  sig.method = "TT", nameCol="Symbol", resSpec="PValues", groupingCol="OriginCode", mname="metaPvalue",  DEBUG=0) {
	# specList is a list of specifiers to do the analysis over, these are used to determine which files to load and then merge over the table loaded for the mappingSpecifier, which contains the probeID and then a mapped identifier which is used
	# For an experimentwise version: specList <- paste(ct$GPL[!duplicated(ct$OriginCode)], ct$OriginCode[!duplicated(ct$OriginCode)], "TT", sep="-")
	# platList <- ct$GPL[!duplicated(ct$OriginCode)]
  # if using across platform RankProduct, then specList=unique(ct$GPL) & platList=unique(ct$GPL)
  # Example for single species only RankProduct, across platform results:  omnibusAcrossMappedGrouped(unique(ct$GPL), unique(ct$GPL), "Probe2EntrezMap", dataDir, "fisher",  "RankProd-PValues", "GeneID", "RP2GPL2Entrez") 
  	# Takes the mapping column info and then maps
	# method={'fisher'|'logit'}
  if (sig.method == "RankProd"){
    specList <- paste(unique(ct$GPL),  "RankProd", sep="-")
  }else{
    specList <- paste(ct$GPL[!duplicated(ct$OriginCode)], ct$OriginCode[!duplicated(ct$OriginCode)], groupingCol, "TT", sep="-")
  }
  platList <- ct$GPL[!duplicated(ct$OriginCode)]
  method <- 'fisher'
  
	mt <- data.frame()
	
	for (i in 1:length(specList)) {
		oi <- specList[i]
		plat <- platList[i]
          tryCatch({
            print(c(oi, plat))
            mf <- getStoredTable(dataDir, "geneMappings", plat, mappingSpecifier, DEBUG=DEBUG)
            mf <- mf[!duplicated(mf$probe), ]
            rownames(mf) <- mf$probe
            #print(head(mf))
            nt <- getStoredTable(dataDir, "significanceValues", oi, extension, DEBUG=DEBUG) 
            #print(head(nt))
            
            twoSide = !is.null(nt$pvalTwoSided)
            if (twoSide) {
                nt <- nt[,c("pvalDown", "pvalUp", "pvalTwoSided")]
            } else {
                nt <- nt[,c("pvalDown", "pvalUp")]
            }
            nt <- mergeViaRowNames(mf, nt, all.a=F, all.b=F)
            #print(head(nt))
            mt <- rbind(mt, nt)
            #print(head(mt))
          },
                   error = function(e) print(str(e)) )
                   
	}
	#print(dim(mt))
  
	midList <- unique(mt[,mappingCol])
	metaRes <- data.frame()
        print("Calculating omnibus Statistics")
	#print(head(midList))
	for (mid in midList) {		
		#print(mid)
		symbol <- mt[mt[,mappingCol] == mid, nameCol][1]
		
		pListD <- mt[mt[,mappingCol] == mid, "pvalDown"]
		mpD <- calcOmnibus(pListD, method)
		pListU <- mt[mt[,mappingCol] == mid, "pvalUp"]
		mpU <- calcOmnibus(pListU, method)
		if (twoSide) {
		  pListTS <- mt[mt[,mappingCol] == mid, "pvalTwoSided"]
		  mpTS <- calcOmnibus(pListTS, method)
		}
              
		#mr <- list(ID=mid, pvalDown=mpD, pvalUp=mpU)
   if (twoSide) {mr <-  list(ID=mid, Symbol=symbol, pvalDown=mpD, pvalUp=mpU, pvalTwoSided=mpTS)} else {mr <-  list(ID=mid, Symbol=symbol, pvalDown=mpD, pvalUp=mpU)}
		metaRes <- rbind(metaRes, mr)
  }
      
  metaRes$qvalDown = p.adjust(metaRes$pvalDown, method="fdr")
  metaRes$qvalUp = p.adjust(metaRes$pvalUp, method="fdr")
  
	writeStoredTable(metaRes, dataDir, "omnibusMeta", resSpec, paste(sig.method, method, mappingCol, groupingCol, mname, sep="-"))
	return(metaRes)
}



################################
#  Stage 4  Effect Sizes


cipv <- function(x,y) {
  # p value from confidence interval, given a value and it's sd
  a <- ci(x, y)
  return(a$p)
}




computeRowExpressionSD <- function(x, y, DEBUG=0) {
	# x is the matrix of values
	# y is the class labels (1 or 2)
  #print(sum(!is.na(x[y==1])))
  #print(sum(!is.na(x[y==2])))
  	n1 <- sum(!is.na(x[y==1]))
	n2 <- sum(!is.na(x[y==2]))
	if (n1 >= 2  & n2 >= 2 ) {
	
	    m1 <- mean(x[y==1], na.rm=TRUE)
	    m2 <- mean(x[y==2], na.rm=TRUE)
                sd1 <- sd(x[y==1], na.rm=TRUE)
                sd2 <- sd(x[y==2], na.rm=TRUE)
               
	 	sd <- sqrt(var(x[y==1], na.rm=TRUE)/n1 + var(x[y==2], na.rm=TRUE)/n2    )
                revval <- c(m1, m2, sd, m2-m1, n1, n2, sd1, sd2)
		
	}
	else {revval <- c(NA, NA, NA, NA, n1, n2, NA, NA)}
  #print(x)
  #print(y)
  #print(revval)
  #blah()
  return(revval)
	
}


computePairedRowExpressionSD <- function(x, DEBUG=0) {
  # Assumes the data is an even number length vector, takes second half away from first half and computes the variance in the result
  lv <- length(x)
  if ((lv %%2 ) != 0 ) {
    print(c("ERROR, paired expression vector not even length", lv))
          print(x)
          return(NA)
          }
  a1 <- x[1:(lv/2)]
  a2 <- x[(1+(lv/2)):lv]
  de <- a2-a1
  n <- sum(!is.na(de))*2
 
  if (n >=2 ) {
    m <- mean(de, na.rm=TRUE)
    sdv <- sd(de, na.rm=TRUE)
    revval <- c(m, sdv, n)
  } else {
    revval <- c(NA, NA, n)
  }
  if (DEBUG > 6) {
      print("computePairedRowExpressionSD: x, lv, revval")
      print(x)
      print(lv)
      print(revval)
  }
  return(revval)
}


computeExpressionSD <- function(ct, dataDir, groupingCol="OriginCode", pairID=NA, DEBUG=0, noclobber=TRUE) {
	# assuming logged expression data, looks at difference of means and computes SD
  # If pairID is not NA, it will try to look for a column that has pair information in it, and then use that to calculate a paired effect size
  if (is.na(pairID)) {
    print("Doing unpaired effect size calculations in computeEpressionSD")
  
	for (gpl in unique(ct$GPL)) {
		print(gpl)
		ctS <- ct[ct$GPL == gpl & !is.na(ct[,groupingCol]), c("GSM", "ClassCode", groupingCol) ]
		eMat <- getStoredTable(dataDir, "rawExpression", gpl, "RawExpressionValues", DEBUG=DEBUG)
                if (DEBUG) {print(cnr(eMat))}
		for (oc in unique(ctS[,groupingCol])) {
                  if (!noclobber | !checkStoredTable(dataDir, "expressionVariation", paste(gpl, oc, groupingCol, sep="-"), "DifMeanLog-Welch", DEBUG=DEBUG)) {
                    tryCatch({
                      if (DEBUG) {print(oc)}
                      ctC <- ctS[ctS[,groupingCol] == oc, ]
                      if (DEBUG > 2) {print(cnr(ctC))}
                      gsmList <- ctC$GSM
                      if (sum(gsmList %in% colnames(eMat)) != length(gsmList)) {
                        print(c("Error message, the gsmList does not match the colnames in the expression data table"))
                        print(rbind(gsmList, gsmList %in% colnames(eMat)))
                        print(colnames(eMat))
                        
                      }
                      sm <- as.matrix(eMat[,gsmList])
                      cl <- ctC$ClassCode
                                if (DEBUG) {print(head(ctC))}
                                if (sum(cl != 1 & cl != 2) >= 1) {
                                  print("Potential problem with class coding, expecting 1 & 2")
                                }
                                #print(head(sm))
				ev <- t(apply(sm, 1, computeRowExpressionSD, y=cl, DEBUG=DEBUG))
				#ev <- data.frame(ev)
				colnames(ev) <- c("mean1", "mean2", "sd", "DiffMeans", "n1", "n2", "sd1", "sd2")
                                if (DEBUG)  print(ev[1:3,])
				writeStoredTable(ev, dataDir, "expressionVariation", paste(gpl, oc, groupingCol, sep="-"), "DifMeanLog-Welch")
                              }, 
                                       error = function(e) print(str(e)) )
			}
		}
	}
  
  } else {
    # pairID is not NA
    print("Trying PAIRED effect size calculations in computeEpressionSD")
    for (gpl in unique(ct$GPL)) {
		  print(gpl)
	  	ctS <- ct[ct$GPL == gpl & !is.na(ct[,groupingCol]), c("GSM", "ClassCode", groupingCol, pairID) ]
      eMat <- getStoredTable(dataDir, "rawExpression", gpl, "RawExpressionValues", DEBUG=DEBUG)
                if (DEBUG) {print(cnr(eMat))}
		for (oc in unique(ctS[,groupingCol])) {
                  if (!noclobber | !checkStoredTable(dataDir, "expressionVariation", paste(gpl, oc, groupingCol, sep="-"), "PairedDifMeanLog", DEBUG=DEBUG)) {
                      tryCatch({
                        if (DEBUG) {print(oc)}
                        ctC <- ctS[ctS[,groupingCol] == oc, ]
                        ctC$xpair <- factor(ctC[,pairID])
                        
                        ctC <- ctC[ctC$xpair %in% ctC$xpair[duplicated(ctC$xpair)],] # remove anything not paired up
                        ctC$xpair <- as.numeric(factor(ctC$xpair))
                        ctC <- ctC[order(ctC$ClassCode, ctC$xpair),]
                        gsmList <- ctC$GSM  # Line the GSM list up in order, first half is ClassCode = 1, Second half is ClassCode = 2, then ordered by their pair number, can add one to the other
            
                      if (DEBUG > 2) {print(cnr(ctC))
                        print(rbind(colnames(eMat), colnames(eMat) %in% gsmList ))
                        print(rbind(gsmList, gsmList %in% colnames(eMat) ))                        
                        }
                      sm <- as.matrix(eMat[,gsmList])
                			ev <- t(apply(sm, 1, computePairedRowExpressionSD, DEBUG=DEBUG))
				colnames(ev) <- c("DiffMeans", "sd", "n")
                                if (DEBUG)  print(ev[1:3,])
				writeStoredTable(ev, dataDir, "expressionVariation", paste(gpl, oc, groupingCol, sep="-"), "PairedMeanLog")
                              }, error = function(e) print(str(e)) )
                  }      
                  }
    }
  }
}




metaEffectSizeEstimate <- function(ct, dataDir, mappingSpecifier, mappingCol, groupingCol="OriginCode",  ...) {
  ### metaEffectSizeEstimate(ct, dataDir, mappingSpecifier="Probe2EntrezMap", mappingCol="GeneID", groupingCol="GSE", DEBUG=3 )
  uniqOC <- !duplicated(ct[,groupingCol])
  computeEffectsModel(dataDir,gplID=ct$GPL[uniqOC], expID=ct[,groupingCol][uniqOC], groupingCol=groupingCol, mappingSpecifier=mappingSpecifier, mappingCol=mappingCol, ...)


  
}


  
computeEffectsModel <- function(dataDir, gplID, expID, mappingSpecifier, mappingCol, dataExtension="DifMeanLog-Welch", groupingCol="", resSpec="metaEffect", mname="", method="both", mappedName="Symbol", DEBUG=0)  {
	# dataDir
	# gplID -- list of GPL identifiers, one item for each list in the groupID list, eg (c("GPL81", "GPL81", "GPL97"))
	# expID -- a grouping colum identifier, one item for each group ID, is the ID used in the file for essentially each origin which is examples, eg c("A", "B", "K")
	# dataExtension -- what is the file extension/class of things that are being merge/meta-analyzed, eg "DifMeanLog-Welch"
	# resSpec -- the specifier for the results file
  # method is essentially ignored, as it does them all anyway
  # micro/macroaveraged results too
  outName <- mname
  
    if (DEBUG > 0) print(match.call())

  
	masterTab <- data.frame(stringsAsFactors = FALSE)
	for (gpl in unique(gplID)) {
		gb <- gplID==gpl  # gpl boolean, which ones are relevant
		mf <- getStoredTable(dataDir, "geneMappings", gpl, mappingSpecifier, DEBUG=DEBUG)
		mf <- mf[!duplicated(mf$probe), ]
		rownames(mf) <- mf$probe
		#  Get the mapping we want
		for (oc in expID[gb]) {
			#gf <- groupFactor[gb & expID== oc][1]
			nt <- getStoredTable(dataDir, "expressionVariation", paste(gpl, oc, groupingCol, sep="-"), dataExtension, DEBUG=DEBUG) 
			mr <- mergeViaRowNames(mf, nt, all.a=F, all.b=F)
			#mr <- cbind(mr, Exp=paste(gpl, oc, sep="-"), GroupFactor=gf)
			mr <- cbind(mr, Exp=paste(gpl, oc, sep="-"), OC=oc)
			if (DEBUG > 2) {print(head(mr))}
			masterTab <- rbind(masterTab, mr)
		}
	}
	if (DEBUG) print(head(masterTab))
	finalResult <- data.frame(stringsAsFactors = FALSE)
	for (mid in unique( masterTab[,mappingCol]) ) {
          ## Figure out how to best do the micro and macro averaging and those separate metanalyses here
          ## TODO!
          
          subTab <- masterTab[ !is.na(masterTab[,mappingCol]) & masterTab[,mappingCol]==mid ,]
          if (sum(!is.na(subTab$sd2) &  !is.na(subTab$sd1)) > 0) {
            subTab <- subTab[!is.na(subTab$sd2) &  !is.na(subTab$sd1), ]
                                        #print(subTab)
            tn1 <- tapply(subTab$n1, subTab[,mappingCol], sum)
            tn2 <- tapply(subTab$n2, subTab[,mappingCol], sum)
            #kbar1 <- mean(tn1)
            #kbar2 <- mean(tn2)
            M1 <-  tapply(subTab$n1, subTab[,mappingCol], max)
            M2 <-  tapply(subTab$n2, subTab[,mappingCol], max)
            expN <- length(unique(subTab$OC))
            
            mcr <- metacont(subTab$n2, subTab$mean2, subTab$sd2, subTab$n1, subTab$mean1, subTab$sd1)
            mcrS <- summary(mcr)
            frm <-  data.frame(ID=mid, Name=subTab[1,mappedName], FixedEffects=mcr$TE.fixed, FixedEffectsSE=mcr$seTE.fixed, FixedEffectsP=mcrS$fixed$p, FixedEffectsZ=mcrS$fixed$z, RandomEffects=mcr$TE.random, RandomEffectsSE=mcr$seTE.random, RandomEffectsP=mcrS$random$p, RandomEffectsZ=mcrS$random$z, k=mcr$k, nsum=sum(subTab$n1, subTab$n2), Q=mcr$Q, tau=mcr$tau, I2=mcrS$I2[[1]], I2.lower=mcrS$I2[[2]], I2.upper=mcrS$I2[[3]], H=mcrS$H[[1]], H.lower=mcrS$H[[2]], H.upper=mcrS$H[[3]], k1=sum(tn1), k2=sum(tn2), M1=M1, M2=M2, expN=expN   )
                  #print(frm)
            finalResult <- rbind(finalResult, frm)
                }
        }
	
	writeStoredTable(finalResult, dataDir, "effectsMeta", resSpec, paste(paste(mappingSpecifier, method, groupingCol, sep="-"), outName, sep=""), DEBUG=DEBUG)
	return(finalResult)
}

metaEffectEstimateOnGroup <- function(subTab, mappingCol, mappedName, DEBUG=0) {
  #  Can use rma.uni instead of metaconf if want more flexibility, need to check the speed and options for rma.uni
  
    if (DEBUG > 3) {print(subTab)}
    
    
  if (sum(!is.na(subTab$sd2) &  !is.na(subTab$sd1)) > 0) {
    id <- subTab[,mappingCol][1]
    subTab <- subTab[!is.na(subTab$sd2) &  !is.na(subTab$sd1), ]
    tn1 <- tapply(subTab$n1, subTab[,mappingCol], sum)
    tn2 <- tapply(subTab$n2, subTab[,mappingCol], sum)
    kbar1 <- mean(tn1)
    kbar2 <- mean(tn2)
    M1 <-  tapply(subTab$n1, subTab[,mappingCol], max)
    M2 <-  tapply(subTab$n2, subTab[,mappingCol], max)
    mcr <- metacont(subTab$n2, subTab$mean2, subTab$sd2, subTab$n1, subTab$mean1, subTab$sd1)
    mcrS <- summary(mcr)
    if (DEBUG > 3) {print(mcrS)}
    frm <-  data.frame(ID=id, Name=subTab[1,mappedName], FixedEffects=mcr$TE.fixed, FixedEffectsSE=mcr$seTE.fixed, FixedEffectsP=mcrS$fixed$p, FixedEffectsZ=mcrS$fixed$z, RandomEffects=mcr$TE.random, RandomEffectsSE=mcr$seTE.random, RandomEffectsP=mcrS$random$p, RandomEffectsZ=mcrS$random$z, kMapped=mcr$k, nsum=sum(subTab$n1, subTab$n2), Q=mcr$Q, tau=mcr$tau, I2=mcrS$I2[[1]], I2.lower=mcrS$I2[[2]], I2.upper=mcrS$I2[[3]], rep1=sum(tn1), rep2=sum(tn2), M1=M1, M2=M2   )
    if (DEBUG > 3) print(frm)
    return(frm)
  } else{
     if (DEBUG > 2) print(c(subTab$sd1, subTab$sd2))
    return()
  }
}
  
geneLevelMetaEffectSizeEstimate <- function(ct, dataDir, mappingSpecifier, mappingCol,  dataExtension="DifMeanLog-Welch",groupingCol="OriginCode",  resSpec="MappedMetaEffect", mname="", method="both", mappedName="Symbol", parallelVers=FALSE, noclobber=TRUE, DEBUG=0){
#dataDir, gplID, expID, mappingSpecifier, mappingCol, dataExtension="DifMeanLog-Welch", groupingCol="", resSpec="metaEffect", mname="", method="both", mappedName="Symbol", DEBUG=0
# ct, dataDir, mappingSpecifier, mappingCol, groupingCol="OriginCode",
  
    if (DEBUG > 0) print(match.call())

  uniqOC <- !duplicated(ct[,groupingCol])
  gplID=ct$GPL[uniqOC]
  expID=ct[,groupingCol][uniqOC]
  outName <- mname
  if (DEBUG > 2) {
    print(expID)
    print(gplID)
  }
  for (gpl in unique(gplID)) {
    print(gpl)
    gb <- gplID==gpl  # gpl boolean, which ones are relevant
    mf <- getStoredTable(dataDir, "geneMappings", gpl, mappingSpecifier, DEBUG=DEBUG)
    mf <- mf[!duplicated(mf$probe), ]
    rownames(mf) <- mf$probe
    if (DEBUG > 1) print(expID[gb])
    for (oc in expID[gb]) {
      if(DEBUG) print(oc)
      if (!noclobber | !checkStoredTable(dataDir, "effectsMeta", resSpec, paste(paste(mappingSpecifier, oc, groupingCol, sep="-"), outName, sep=""))) {
        nt <- getStoredTable(dataDir, "expressionVariation", paste(gpl, oc, groupingCol, sep="-"), dataExtension, DEBUG=DEBUG)
        if (DEBUG) {print(head(nt))}
        mr <- mergeViaRowNames(mf, nt, all.a=F, all.b=F)
        mr <- cbind(mr, Exp=paste(gpl, oc, sep="-"), OC=oc)
        if (DEBUG) {print(head(mr))}
        mr <- mr[!is.na(mr$sd1) & !is.na(mr$sd2) & !is.na(mr$mean1) & !is.na(mr$mean2) & mr$sd1 > 0 & mr$sd2 > 0,]
        metaRes <- by(mr,mr[,mappingCol], metaEffectEstimateOnGroup, mappingCol=mappingCol, mappedName=mappedName, DEBUG=DEBUG)
        metaResTM <- metaRes[as.numeric(lapply(metaRes, length)) != 0]
        if(DEBUG) {print(head(metaRes))}
        metaResM <- data.frame(t(sapply(metaResTM, rbind)))
        if(DEBUG) {print(head(metaResM))}
        writeStoredTable( metaResM, dataDir, "effectsMeta", resSpec, paste(paste(mappingSpecifier, oc, groupingCol, sep="-"), outName, sep=""), DEBUG=DEBUG)
      }
    }
  }
}


## DO NOT USE:
    geneLevelMESEoverMultiLevels <- function(ct, dataDir, mappingSpecifier, mappingCol,  dataExtension="DifMeanLog-Welch", groupingCol="OriginCode", higherGroupingCol="StudyGroup", resSpec="MappedMetaEffect", mname="", method="both", mappedName="Symbol", pairID=NA, parallelVers=FALSE, noclobber=TRUE, DEBUG=0){
## DO NOT USE
  # higherGroupingCol takes in a high level group (such as representing multiple arrays with different origin codes really mapping to a single study)
  # The meta effects estimate than operates over this
  # Replaced by the function "metaEffectSizeEstimateOverMappedLevel"
  #  If pairID is not set to NA and is given a name, it assumes that values are paired and does a meta-analysis over paired values
    
    if (DEBUG > 0) print(match.call())

  hgOCL <- unique(ct[,higherGroupingCol])
  outName = mname
  if (DEBUG) print(hgOCL)
  for (hgoc in hgOCL) {
      if (!noclobber | !checkStoredTable(dataDir, "effectsMeta", resSpec, paste(paste(mappingSpecifier, hgoc, higherGroupingCol, sep="-"), outName, sep=""), DEBUG=DEBUG)) {
        emT <- data.frame(stringsAsFactors = FALSE)
        if (DEBUG > 1) {print(hgoc)}
        ctSM <- ct[ct[,higherGroupingCol] == hgoc,]
        sgplL <- unique(ctSM$GPL)
        for (gpl in sgplL) {
          ctVM <- ctSM[ctSM$GPL == gpl,]
          ocSL <- unique(ctVM[,groupingCol])
          for (oc in ocSL) {
                mf <- getStoredTable(dataDir, "geneMappings", gpl, mappingSpecifier, DEBUG=DEBUG)
                mf <- mf[!duplicated(mf$probe), ]
                rownames(mf) <- mf$probe
                nt <- getStoredTable(dataDir, "expressionVariation", paste(gpl, oc, groupingCol, sep="-"), dataExtension, DEBUG=DEBUG)
                if (DEBUG > 2) { print(cnr(mf)); print(cnr(nt))}
                mr <- mergeViaRowNames(mf, nt, all.a=F, all.b=F)
                mr <- cbind(mr, Exp=paste(gpl, oc, sep="-"), OC=oc)
                emT <- rbind(emT, mr)
                if (DEBUG > 2) {
                  print(c(gpl, oc, dim(emT)))
                  print(cnr(emT, 2, 3))}
              }}
        if (DEBUG > 1) {
          print(dim(emT))
          print(cnr(emT, 2, 3))
        }
        metaRes <- by(emT,emT[,mappingCol], metaEffectEstimateOnGroup, mappingCol=mappingCol, mappedName=mappedName, DEBUG=DEBUG)
        metaResTM <- metaRes[as.numeric(lapply(metaRes, length)) != 0]
        metaResM <- data.frame(t(sapply(metaResTM, rbind)))
        #if (length(dim(metaResM)) < 2) {return(metaResM)}
        
        if(DEBUG) {print(cnr(metaResM))}
        writeStoredTable(metaResM, dataDir, "effectsMeta", resSpec, paste(paste(mappingSpecifier, hgoc, higherGroupingCol, sep="-"), outName, sep=""), DEBUG=DEBUG)
      }
    }
}


metaEffectSizeEstimateOverMappedLevel <- function(ct, dataDir, mappingSpecifier, mappingCol,  dataExtension="DifMeanLog-Welch", groupingCol="OriginCode", resSpec="MappedMetaEffect", mname="FandROverMappedBelow", outName="FandROverMappedBelow", method="both", mappedName="Symbol", outDir="effectsMeta", parallelVers=FALSE, noclobber=TRUE, DEBUG=0){
  # resSpec is part of the file name for the files which represent a study below this that is at the level of the genes (or mappingSpecifier) to a single value for expression
  if (DEBUG > 0) print(match.call())
  
    if (noclobber & checkStoredTable(dataDir, outDir, resSpec, paste(paste(mappingSpecifier, mappingCol,  groupingCol, sep="-"), outName, sep=""), DEBUG=(DEBUG) ) ) {
      return(getStoredTable(dataDir, outDir, resSpec, paste(paste(mappingSpecifier, mappingCol,  groupingCol, sep="-"), outName, sep=""), DEBUG=(DEBUG)) )
    }
  masterTab <- data.frame(stringsAsFactors = FALSE)
  masterRes <- data.frame(stringsAsFactors = FALSE)
  uniqOC <- !duplicated(ct[,groupingCol])
  gplID=ct$GPL[uniqOC]
  expID=ct[,groupingCol][uniqOC]
  if (DEBUG > 2) {
    print(expID)
    print(gplID)
  }
  for (gpl in unique(gplID)) {
    if (DEBUG) print(gpl)
    gb <- gplID==gpl  # gpl boolean, which ones are relevant
    if (DEBUG > 1) print(expID[gb])
    for (oc in expID[gb]) {
      if(DEBUG) print(oc)
    
        mt <- getStoredTable(dataDir, "effectsMeta", resSpec, paste(paste(mappingSpecifier, oc, groupingCol, sep="-"), mname, sep=""), DEBUG=DEBUG)
      mt$GPL <- gpl
      mt$OriginCode <- oc
        if (DEBUG > 2) print(c(gpl, oc))
        if (DEBUG > 2) print(head(mt))
        masterTab <- rbind(masterTab, mt)
    }
  }
  for (id in unique(masterTab$ID)) {
    st <- masterTab[masterTab$ID == id,]
    rem <- metagen(st$RandomEffects, st$RandomEffectsSE, comb.random=T, comb.fixed=F)
    fem <- metagen(st$FixedEffects, st$FixedEffectsSE, comb.random=F, comb.fixed=T)
    if (DEBUG > 3) {print(st); print(rem)}
    srem <- summary(rem)
    sfem <- summary(fem)
    rres <- data.frame(ID=id, Name=st[1,"Name"], RandomEffects=rem$TE.random, RandomEffectsSE=rem$seTE.random, RandomEffectsP=srem$random$p, RandomEffectsLogP=-log10(srem$random$p), I2rem=srem$I2[[1]],
                       FixedEffects=fem$TE.fixed, FixedEffectsSE=fem$seTE.fixed, FixedEffectsP=sfem$fixed$p, FixedEffectsLogP=-log10(sfem$fixed$p), I2fem=sfem$I2[[1]], rep1=sum(st$rep1), rep2=sum(st$rep2), M1=sum(st$M1), M2=sum(st$M2), kMapped=sum(st$kMapped), nsum=sum(st$nsum), mesStudies=dim(st)[1])
    masterRes <- rbind(masterRes, rres)
  }
  if (DEBUG) print(head(masterRes))
  writeStoredTable(masterRes, dataDir, outDir, resSpec, paste(paste(mappingSpecifier, mappingCol, groupingCol, sep="-"), outName, sep=""), DEBUG=DEBUG)
    #  Some issues in naming changes here, modify later maybe - try to fix this up using previous conventions
  return(masterRes)
}


leaveOneOutMetaES.OverMappedLevel <- function(ct, dataDir, mappingSpecifier, mappingCol,  dataExtension="DifMeanLog-Welch", groupingCol="OriginCode",  resSpec="MappedMetaEffect",mname="Human", method="both", mappedName="Symbol", outDir="effectsMetaLOO", inDir="effectsMeta", parallelVers=FALSE, p.thresh=.0001, isq.thresh=.0001, nstud.thresh=(length(unique(ct$OriginCode))),  noclobber=TRUE, DEBUG=1){
    # This does leave one out meta-analysis for effect size, leaving out one study/exp (as indicated by groupingCol) at a time
  #print("leaveOneOutMetaES.OverMappedLevel()")
    
  if (DEBUG > 0) print(match.call())
  
    figCells <- ceiling(sqrt(unique(ct[,c(groupingCol)])))
    par(mfrow=c(figCells, figCells))
    
  for (sg in unique(ct[,c(groupingCol)]) ) {
    print(c("sg:", sg))
    ctSS <- ct[ct[,c(groupingCol)] != sg,]    
    outName.loo <- paste("LOO", sg, sep="-")
    print(outName.loo)
    mes <- metaEffectSizeEstimateOverMappedLevel(ctSS, dataDir=dataDir, mappingSpecifier=mappingSpecifier, mappingCol=mappingCol,  dataExtension=dataExtension, groupingCol="StudyCode", resSpec=resSpec, mname=mname, outName=outName.loo, , outDir=outDir, noclobber=noclobber, DEBUG=DEBUG)
    colnames(mes) <- paste(colnames(mes), "meloo", sep=".")
    
    mes <- mes[mes$RandomEffectsP.meloo <= p.thresh & mes$I2rem.meloo <= isq.thresh & mes$mesStudies.meloo >= nstud.thresh, ]  #
   mt <- getStoredTable(dataDir, inDir, resSpec, paste(paste(mappingSpecifier, sg, groupingCol, sep="-"), mname, sep=""), DEBUG=DEBUG)
  colnames(mt) <- paste(colnames(mt), "singular", sep=".")    
  mz <- mergeViaRowNames(mes, mt)

    
  #kendall.es <- cor.test(mes$RandomEffects.meloo, mt$RandomEffects.singular, method="kendall")
   # print(kendall.es)
  #  kendall.es.tau <- kendall.es$statistic
    
  }
}
                       
twoLevelMappedMetaAnalysisPairwiseMappable <- function(ct, dataDir, mappingSpecifier, mappingCol, groupingCol, higherGroupingCol="StudyGroup", pairID=NULL, outName="twoLevelMappedMAPM", noclobber=TRUE, DEBUG=0) {
  #  Added Sept 7, 2011
  # Does a multilevel meta-analysis based on effect size, allows effect size to be based on paired samples
    # First meta-analysis goes across all platforms (GPL) and all groups/arms (groupingCol) of a study and gets a mapped (ie gene) level expression arcross all the platforms (GPL) in that study
  #  Second meta-analysis goes across all the studies and gets an overall meta-estimate for those studies
  outDir = "effectsMeta"      
  if (noclobber & checkStoredTable(dataDir, outDir, outName, paste(mappingSpecifier, mappingCol,  groupingCol, higherGroupingCol, sep="-"), DEBUG=(DEBUG) ) ) {
  return(getStoredTable(dataDir, outDir, outName, paste(mappingSpecifier, mappingCol,  groupingCol, higherGroupingCol, sep="-"), DEBUG=(DEBUG)) )
    } else {
      allSubstudyMetaRes <- data.frame()
      for (hloc in unique(ct[,higherGroupingCol])){ 
        if (DEBUG) {print(c("Top of Studies Arm, HLOC:", hloc))}
          if (noclobber & checkStoredTable(dataDir, outDir, "WholeStudyMappedMA", paste(hloc, higherGroupingCol, mappingCol, mappingSpecifier,  sep="-"), DEBUG=(DEBUG) ) ) {
              metaResM <- getStoredTable(dataDir, outDir, "WholeStudyMappedMA", paste(hloc, higherGroupingCol, mappingCol, mappingSpecifier,  sep="-"), DEBUG=(DEBUG)) 
          } else {
              for (oc in unique(ct[ct[,higherGroupingCol]==hloc,groupingCol])) {
              ocED = data.frame()
              for (gpl in unique(ct[ct[,groupingCol]==oc & ct[,higherGroupingCol]==hloc, "GPL"])) {
                if (DEBUG) {print (c("OC:", oc, "GPL:", gpl))}    
                mf <- getStoredTable(dataDir, "geneMappings", gpl, mappingSpecifier, DEBUG=DEBUG)
                mf <- mf[!duplicated(mf$probe), ]
                rownames(mf) <- mf$probe
                if (is.null(pairID)) {
                  nt <- getStoredTable(dataDir, "expressionVariation", paste(gpl, oc, groupingCol, sep="-"), "DifMeanLog-Welch" , DEBUG=DEBUG)} else {
                  nt <- getStoredTable(dataDir, "expressionVariation", paste(gpl, oc, groupingCol, sep="-"), "PairedMeanLog" , DEBUG=DEBUG)} 
                mr <- mergeViaRowNames(mf, nt, all.a=F, all.b=F)
                mr <- cbind(mr)
                ocED <- rbind(ocED, mr)
                if (DEBUG > 2) { print(c(gpl, oc));print(c("dim", dim(ocED)));print(head(ocED))}  
              }
            }
              metaResSS <- by(ocED,ocED[,mappingCol], singleLevelMetaAnalysisPairwiseMapped, mappingCol=mappingCol, mappedName="Symbol", pairID=pairID, DEBUG=DEBUG)
            metaResTM <- metaResSS[as.numeric(lapply(metaResSS, length)) != 0] # emove emply lines
            metaResM <- data.frame(t(sapply(metaResTM, rbind)))   # Put it into a data.frame, from a list of mixed values
            if (DEBUG >2) {print("metaResM:"); print(head(metaResM))}
            writeStoredTable(metaResM, dataDir, outDir, "WholeStudyMappedMA", paste(hloc, higherGroupingCol, mappingCol, mappingSpecifier,  sep="-"), DEBUG=(DEBUG))          
      }
      metaResM <- cbind(metaResM, hloc)
      allSubstudyMetaRes <- rbind(allSubstudyMetaRes, metaResM)
      }                            
      # Do the higher level meta-analysis over the subgroups
      if (DEBUG) {print("tail of the allSubstudyMetaRes"); print(tail(allSubstudyMetaRes))}
      OAmetaRes <- by(allSubstudyMetaRes,allSubstudyMetaRes[,"ID"], secondLevelMetaAnalysisMapped, mappingCol=mappingCol, mappedName="Name", DEBUG=DEBUG)
      browser(); browser(); browser(); browser(); browser(); browser()  # Calls to debugger
      OAmetaResTM <- OAmetaRes[as.numeric(lapply(OAmetaRes, length)) != 0] # remove emply lines
      OAmetaResM <- data.frame(t(sapply(OAmetaResTM, rbind)))                      
      if (DEBUG) {print("head of OAmetaResM");print(head(OAmetaResM)) }
      #write the overall meta table
      writeStoredTable(OAmetaResM, dataDir, outDir, outName, paste(mappingSpecifier, mappingCol,  groupingCol, higherGroupingCol, sep="-"), DEBUG=(DEBUG) )  
      return(OAmetaResM)
  }
}

secondLevelMetaAnalysisMapped <- function(subTab, mappingCol, mappedName="Symbol", pairID=NULL, DEBUG=0, ...) {
  # Does an inverse variance weighted meta-estimate using previous meta-analysis of studies (does not do fixed over random as that doesn't make much sense)
  # Fixed and Random giving the same results, need to check why, is it metagen instead of metacont?

if (DEBUG > 4) {
      print("secondLevelMetaAnalysisMapped"); print(subTab);  print(c(mappingCol, mappedName))}
       if (sum(!is.na(subTab$RandomEffectsSE)) > 0) {

      mFE <- metagen(subTab$FixedEffects, subTab$FixedEffectsSE, comb.random=T, comb.fixed=T)
      mFESummary <- summary(mFE)
      mRE <- metagen(subTab$RandomEffects, subTab$RandomEffectsSE, comb.random=T, comb.fixed=F)
      mRESummary <- summary(mRE) 
      meanQ <- mean(subTab$Q, na.rm=T)
      meanI2 <- mean(subTab$I2, na.rm=T)
      frm <- data.frame(ID=subTab$ID[1], Name=subTab[1,mappedName],
            FixedOverFixed=mFE$TE.fixed, FixedOverFixedSE=mFE$seTE.fixed,  FixedOverFixedP=mFESummary$fixed$p, RandomOverFixed=mFE$TE.random, RandomOverFixedSE=mFE$seTE.random, RandomOverFixedP=mFESummary$random$p, RandomOverRandom=mRE$TE.random, RandomOverRandomSE=mRE$seTE.random,RandomOverRandomP=mRESummary$random$p, numberStudies=mFE$k, sumKMapped = sum(subTab$kMapped), nsum=(sum(subTab$nsum)), QoverF=mFE$Q, QoverR=mRE$Q, I2overF=mFESummary$I2[[1]], I2overR=mRESummary$I2[[1]], rep1=sum(subTab$rep1), rep2=sum(subTab$rep2), M1=sum(subTab$M1), M2=sum(subTab$M2))
      if (DEBUG > 4) print(frm)
    return(frm)
    }else{
       if (DEBUG) {print("Error, missing values in secondLevelMetaAnalysisMapped:"); print(head(subTab))}
          return()
  }
}      
        
singleLevelMetaAnalysisPairwiseMapped <- function(subTab, mappingCol, mappedName="Symbol", pairID=NULL, DEBUG=0 ){
  # Meant to be called by multiLevelMappedMetaAnalysisPairwiseMappable
  # Does a meta-analysis across the values
  #  If results are from a paired analysis, there is not a separate sd, just single deviations
      
    if (DEBUG > 3) {
      print("singleLevelMetaAnalysisPairwiseMapped")
      print(subTab)
      print(c(mappingCol, mappedName))}
  if (is.null(pairID)) {   
  if (sum(!is.na(subTab$sd2) &  !is.na(subTab$sd1)) > 0) {
    id <- subTab[,mappingCol][1]
    subTab <- subTab[!is.na(subTab$sd2) &  !is.na(subTab$sd1), ]
    tn1 <- tapply(subTab$n1, subTab[,mappingCol], sum)
    tn2 <- tapply(subTab$n2, subTab[,mappingCol], sum)
    kbar1 <- mean(tn1)
    kbar2 <- mean(tn2)
    M1 <-  tapply(subTab$n1, subTab[,mappingCol], max)
    M2 <-  tapply(subTab$n2, subTab[,mappingCol], max)
    mcr <- metacont(subTab$n2, subTab$mean2, subTab$sd2, subTab$n1, subTab$mean1, subTab$sd1)
    mcrS <- summary(mcr)
    if (DEBUG > 3) {print(mcrS)}
    frm <-  data.frame(ID=id, Name=subTab[1,mappedName], FixedEffects=mcr$TE.fixed, FixedEffectsSE=mcr$seTE.fixed, FixedEffectsP=mcrS$fixed$p, FixedEffectsZ=mcrS$fixed$z, RandomEffects=mcr$TE.random, RandomEffectsSE=mcr$seTE.random, RandomEffectsP=mcrS$random$p, RandomEffectsZ=mcrS$random$z, kMapped=mcr$k, nsum=sum(subTab$n1, subTab$n2), Q=mcr$Q, I2=mcrS$I2[[1]], I2.lower=mcrS$I2[[2]], I2.upper=mcrS$I2[[3]], rep1=sum(tn1), rep2=sum(tn2), M1=M1, M2=M2   )
    if (DEBUG > 3) print(frm)
    return(frm)
  } else{
     if (DEBUG > 2) print(c(subTab$sd1, subTab$sd2))
    return()
  }
  } else{
    
     if (sum(!is.na(subTab$sd)) > 0) {
      mcr <- metagen(subTab$DiffMeans, subTab$sd, comb.random=T, comb.fixed=T)
      mcrS <- summary(mcr)
      frm <- data.frame(ID=subTab[1,mappingCol], Name=subTab[1,mappedName],FixedEffects=mcr$TE.fixed, FixedEffectsSE=mcr$seTE.fixed, FixedEffectsP=mcrS$fixed$p, FixedEffectsZ=mcrS$fixed$z, RandomEffects=mcr$TE.random, RandomEffectsSE=mcr$seTE.random, RandomEffectsP=mcrS$random$p, RandomEffectsZ=mcrS$random$z, kMapped=mcr$k, nsum=(sum(subTab$n)*2), Q=mcr$Q, tau=mcr$tau,I2=mcrS$I2[[1]], I2.lower=mcrS$I2[[2]], I2.upper=mcrS$I2[[3]], M1=max(subTab$n), M2=max(subTab$n) )
      if (DEBUG > 3) print(frm)
        return(frm)
      } else{
          if (DEBUG) {print("Error, missing values:"); print(head(subTab))}
          return() 
        }
  }
}  
  
  

########################################
  
computePairwiseDistances <- function(dataDir, gplID, expID, groupFactor, dataNames, mappingSpecifier, method="spearman", graphPrefix="Dist", mappingCol="ID", metaMethod="fixed", effectsCol= "MetaEffect", effectsSD="MetaSD", graphDir="distancePlots", DEBUG=0)  {
		# dataDir
		# gplID -- list of GPL identifiers, one item for each list in the groupID list, eg (c("GPL81", "GPL81", "GPL97"))
		# expID -- one item for each group ID, is the ID used in the file for essentially each origin which is examples, eg c("A", "B", "K")
		# groupFactor -- goes along with the gplID & expID, is some text that groups the experiment/origin, could be species or tissue or both, eg c("Mouse", "Mouse", "Human")
		# dataNames are the short names of the datasets (typically longer than the expID's, but can be the same)
	# mappingSpecifier is what was used, in many cases "Probe2HID"
	# mappingCol is the name of the column in the mapped data that has the data
	
	mt <- data.frame(stringsAsFactors = FALSE)
	for (gpl in unique(gplID)) {
		gb <- gplID==gpl  # gpl boolean, which ones are relevant
		for (oc in expID[gb]) {
			print(c(gpl, oc))
			gf <- groupFactor[gb & expID== oc][1]
		
			nt <- getStoredTable(dataDir, "mappedEffectsMeta", paste(gpl, oc, sep="-"), paste(mappingSpecifier, metaMethod, "metaEffect", sep="-"))
			#print(head(nt))
			#print(c(mappingCol, effectsCol, effectsSD))
			nt <- nt[,c(mappingCol, effectsCol, effectsSD)]
			#print(head(nt))
			ont <- cbind(expID=rep(oc, dim(nt)[1]), nt)
			mt <- rbind(mt, ont )
		}
	}
	#print(head(mt))
	
	n <- length(expID)
	sMat <- matrix(0, n, n)
	rownames(sMat) <- dataNames
	colnames(sMat) <- dataNames

	
	for (i in 1:(n-1)) {	
		exp1 <- expID[i]
		es1 <- mt[expID==exp1,c(mappingCol, effectsCol, effectsSD)]
		for (j in (i+1):n) {		
			exp2 <- expID[j]
			#print(c(i, j, exp1, exp2))
			es2 <- mt[expID==exp2,c(mappingCol, effectsCol, effectsSD)]
			mexp <- merge(es1, es2, by.x=mappingCol, by.y=mappingCol, all.x=F, all.y=F)
			#print(head(mexp))
			size <- dim(mexp)[1]
			if (method == "spearman") {d <- 1- cor(mexp[,2], mexp[,4], method="spearman")}
			if (method == "intersect" | method=="tanimoto")  {
				un <- length(unique(c(es1[,mappingCol],es2[,mappingCol] )))
				isect <- dim(mexp)[1]
				d <- 1 - ((un - isect)/un)
			}
			if (method == "weightedspearman") {
				w <- 1/sqrt( (es1[,effectsSD])^2 + (es2[,effectsSD])^2 )
				d <- 1-  (cov.wt(cbind((rank(es1[,mappingCol])), rank(es2[,mappingCol])), w, cor=TRUE)$cor[1,2]) 
				}
			if (method == "weightedpearson") {
				w <- 1/sqrt( (es1[,effectsSD])^2 + (es2[,effectsSD])^2 )
				d <-  1- (cov.wt(cbind(es1[,mappingCol], es2[,mappingCol]), w, cor=TRUE)$cor[1,2]) 
				}
			sMat[i,j] <- d
			sMat[j,i] <- d
		}
	}
	return(sMat)

}



################################
# Graphing Functions


graphExpressionDensity <- function(emat, type="l", na.rm=TRUE, leg=NA, ...) {
	    # Makes smoothed histograms of expression Density
	    x.density <- apply(emat, 2, density, na.rm = na.rm)
	    all.x <- do.call("cbind", lapply(x.density, function(x) x$x))
	    all.y <- do.call("cbind", lapply(x.density, function(x) x$y))
		colNum <- dim(emat)[2]
		if (colNum < 12) {
			col <- brewer.pal(11, "Spectral")
		}  else {
			bcpal <- brewer.pal(11, "Spectral")
			col <- colorRampPalette(bcpal)(colNum)
		}
	    matplot(all.x, all.y, col = col, lty=1, type=type, xlab="Expression Level", ylab="Density", ...)
	    #invisible(list(all.x = all.x, all.y = all.y))
		if (is.na(leg)) {leg <- colnames(emat)}
		legend("topright", leg, col=col, cex=0.5, inset=0.1, lty=1)
}



createEffectsAttributesFromCT <- function(ct, cfcol=c("Species") ) {
  # The most common way you would want to do a meta effects model
	eam <- data.frame()
	for (oc in unique(ct$OriginCode)) {
		sl <- ct[ct$OriginCode==oc, c("GPL", "OriginCode", cfcol)][1,]
		eam <- rbind(eam, sl)
	}
	return(eam)
}



graphEffectsModelNamedExperiment <- function(dataDir, gplID, expID, groupFactor, mappingSpecifier, mappingCol, dataExtension, groupingCol=NULL, nameList=NULL,  outName="ForestPlot",  method="fixed", graphDir="forestPlots", graphList=NULL , colorLegend=TRUE, mappedName="Symbol", estimate=FALSE, outSize=8, DEBUG=0)  {
  if (method != "fixed") {method <- "random"}  #  Just set it
  print(method)
  gcolorfl <- levels(factor(groupFactor))
  if (length(gcolorfl) == 1) {
    colors <- c("black")
  } else if (length(gcolorfl) < 3) {colors <- c("red", "blue")} else if (length(gcolorfl) <= 9) {
    colors <- brewer.pal(length(gcolorfl), "Set1")
  }  else {
    bcpal <- brewer.pal(11, "Spectral")
    colors <- colorRampPalette(bcpal)(length(gcolorfl))
  }
  masterTab <- data.frame(stringsAsFactors = FALSE)
  for (gpl in unique(gplID)) {
    gb <- gplID==gpl  # gpl boolean, which ones are relevant
    mf <- getStoredTable(dataDir, "geneMappings", gpl, mappingSpecifier, DEBUG)
    mf <- mf[!duplicated(mf$probe), ]
    rownames(mf) <- mf$probe
                                        #  Get the mapping we want
    for (oc in expID[gb]) {
      gf <- groupFactor[gb & expID== oc][1]
      if (is.null(nameList)) {  name = oc} else{name <- nameList[gb & expID == oc][1] }
      if (is.null(groupingCol)) { nt <- getStoredTable(dataDir, "expressionVariation", paste(gpl, oc, sep="-"), dataExtension, DEBUG=DEBUG+2)  } else{nt <- getStoredTable(dataDir, "expressionVariation", paste(gpl, oc, groupingCol, sep="-"), dataExtension, DEBUG=DEBUG)  }

      mr <- mergeViaRowNames(mf, nt, all.a=F, all.b=F)
      mr <- cbind(mr, Exp=paste(gpl, oc, sep="-"), GroupFactor=gf, ExpName=name)
                                        #print(head(mr))
      masterTab <- rbind(masterTab, mr[])
    }
  }
  #print(gcolorfl)
  #print(colors)
  for (mid in graphList) {
    indMeta <- data.frame(stringsAsFactors = FALSE)  # individual meta for a single exp, for a single mapping ID
    subTab <- masterTab[masterTab[,mappingCol]==mid ,]
    subTab <- subTab[!is.na(subTab$sd1) & !is.na(subTab$sd2), ] # Remove the missing deviation rows
    if (dim(subTab)[1] > 0) {
      Name <- subTab[1,mappedName]
      if (DEBUG > 1) {print(c(mid, "subTab", dim(subTab) ))
                      print(head(subTab))
                    }
      for (oc in unique(subTab$Exp)) {
        sst <- subTab[subTab$Exp == oc, ]
        if (sum(!is.na(sst$sd2) &  !is.na(sst$sd1)) > 0) {
          sst <- sst[!is.na(sst$sd2) &  !is.na(sst$sd1), ]
                                        #print(sst)
          N = sum(sum(sst$n1), sum(sst$n2))
          mcolor <- colors[gcolorfl == sst$GroupFactor[1]]
          if (dim(sst)[1] == 1) {
            imr <- list(EffectSize=sst$DiffMeans[1], EffectSD=sst$sd[1], N=N, Exp=oc, GroupFactor=sst$GroupFactor[1], ExpName=sst$ExpName[1], mcolor=mcolor)
          } else {
            mcr <- metacont(sst$n2, sst$mean2, sst$sd2, sst$n1, sst$mean1, sst$sd1)
        
            if (method == "fixed") {
              effectsize = mcr$TE.fixed
              standarderror = mcr$seTE.fixed
            } else {
              effectsize = mcr$TE.random
              standarderror = mcr$seTE.random
            }
            imr <- list(EffectSize=effectsize, EffectSD=standarderror, N=N, Exp=oc, GroupFactor=sst$GroupFactor[1], ExpName=sst$ExpName[1], mcolor=mcolor)
          }
          indMeta <- rbind(indMeta, imr)
                                        #print(indMeta)
        }
      }
      if (DEBUG) {print(head(indMeta))}
      if (nrow(indMeta) > 0) {
				#print(dim(indMeta))
        indMeta <- indMeta[!is.na(indMeta$EffectSD) & !is.na(indMeta$EffectSize),]
        indMeta <- indMeta[order(indMeta$EffectSize, decreasing=T),]
      
        outpath <- paste(dataDir, graphDir, sep="/")
        dir.create(outpath, recursive=T, showWarnings=F)
        outname <- paste(outpath, "/", Name,  "-", as.character(mid),"-", outName, "-",  method,".pdf", sep="" )
        if (DEBUG) {print(outname)}
        if (!estimate) {
      
          pdf(outname, height=outSize, width=outSize)
          metaplot(indMeta$EffectSize, indMeta$EffectSD, indMeta$N, labels=indMeta$ExpName, colors=meta.colors(box=indMeta$mcolor), xlab="Log Fold Change", ylab="Experiment", main=subTab[,mappedName][1] )
          if (!is.null(colorLegend) & colorLegend) {
            legend("bottomright", legend=gcolorfl, col=colors, pch=19)
          }
          dev.off()
        } else {
          EST <- subTab
          N = sum(sum(EST$n1), sum(EST$n2))
          if (DEBUG) {
            print(head(EST))
            print(c("N", N))
          }
          mcr <- metacont(EST$n2, EST$mean2, EST$sd2, EST$n1, EST$mean1, EST$sd1)                                
          if (method == "fixed") {
            Seffectsize = mcr$TE.fixed
            Sstandarderror = mcr$seTE.fixed
          } else {
            Seffectsize = mcr$TE.random
            Sstandarderror = mcr$seTE.random
          }
          if (DEBUG) {print(mcr)}
          pdf(outname, height=outSize, width=outSize)
          metaplot(indMeta$EffectSize, indMeta$EffectSD, indMeta$N, labels=indMeta$ExpName, colors=meta.colors(box=indMeta$mcolor), xlab="Log Fold Change", ylab="Experiment", main=subTab[,mappedName][1], summn=Seffectsize, sumse=Sstandarderror, sumnn=N, summlabel="Meta-Estimate")
          if (!is.null(colorLegend) & colorLegend) {
            legend("bottomright", legend=gcolorfl, col=colors, pch=19, cex=.7)
          }
          dev.off()
                                
        }
      }	
    }	
  }
}




################################
#   Pairs plot graphics


     panel.hist <- function(x, col="grey",...)
     {
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(usr[1:2], 0, 1.5) )
         h <- hist(x, plot = FALSE)
         breaks <- h$breaks; nB <- length(breaks)
         y <- h$counts; y <- y/max(y)
         rect(breaks[-nB], 0, breaks[-1], y, col=col, ...)
     }
   panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
     {
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(0, 1, 0, 1))
         r <- cor(x[!is.na(x) & !is.na(y)], y[!is.na(x) & !is.na(y)])
         ar <- abs(r)
         txt <- format(c(r, 0.123456789), digits=digits)[1]
         txt <- paste(prefix, txt, sep="")
         if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
         text(0.5, 0.5, txt, cex = cex * ar)
     }

panel.cor.pvals  <- function(x, y, digits=2, prefix="", cex.cor, ...)
     {
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(0, 1, 0, 1))
         c <- cor.test(x[!is.na(x) & !is.na(y)], y[!is.na(x) & !is.na(y)], method="pearson")
         r <- c$estimate
         pval <- c$p.value
         rtxt <- format(c(r, 0.123456789), digits=digits)[1]
         ptxt <- format(c(pval, 0.123456789), digits=digits)[1]
         txt <- paste(prefix, "r=", rtxt, ", p=", ptxt,sep="")
         if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
         text(0.5, 0.5, txt, cex = cex)
     }


panel.dot <- function(x, y, col="black", bg=NA, pch=".", cex=1, ...) {
      points(x, y, pch = pch, col = col, bg = bg, ...)
}




################################
#  ROC


addIndividualUpDownTestsToMatrix <- function(om, labList, specList, platList, mappingSpecifier, dataDir, method, extension, mappingCol, resSpec, nameCol="Symbol", DEBUG=0) {
	#specList <- paste(ct$GPL[!duplicated(ct$OriginCode)], ct$OriginCode[!duplicated(ct$OriginCode)], "TT", sep="-")
	#platList <- ct$GPL[!duplicated(ct$OriginCode)]
	
	j <- dim(om)[2]
	for (i in 1:length(specList)) {
		j <- j +1
		oi <- specList[i]
		plat <- platList[i]
		print(c(i, oi, plat))
          tryCatch({          
            mf <- getStoredTable(dataDir, "geneMappings", plat, mappingSpecifier, DEBUG=DEBUG)
            mf <- mf[!duplicated(mf$probe), ]
            rownames(mf) <- mf$probe
            nt <- getStoredTable(dataDir, "significanceValues", oi, extension, DEBUG=DEBUG) 
            nt <- mergeViaRowNames(mf, nt, all.a=F, all.b=F)
			nt$minp <- apply(nt[, c("pvalUp",  "pvalDown")], 1, min)
			#print(tail(nt))
			#print(head(om))
			
			om <- merge(om, nt[,c(mappingCol, "minp")], by.x=mappingCol, by.y=mappingCol, all.x=T)
			colnames(om)[j] <- labList[i]
			#print(head(om))
		}, 
                   error = function(e) print(str(e)) )
	}
	return(om)
}

plotROCBasedOnLabels <- function(om, plotL, truthC, labels=NULL, legendB=TRUE,cols=NA, DEBUG=0, ...) {
	# Plots ROC curves for the columns listed in plotL, vs the truth in the truth column, truthC
	#  This lets one plot multiple predictions against the same truth values
  # om is a matrix/dataframe
  # plotL is a list of column headings for columns containing predictions
  # truthC is the column heading of the labels of truth 
  if(sum(plotL %in% colnames(om)) != length(plotL) | length(plotL) < 1) {
    print("Error in call to plot ROCR, showing columns in plotting data")
    print(data.frame(plotL=plotL, inOM=c(plotL %in% colnames(om))))
    return()
  }
  
	nc <- length(plotL)
  if (is.na(cols)) {
	if (nc < 3) {cols <- c("darkred", "darkblue")} else if (nc <= 8) {
			cols <- brewer.pal(nc, "Dark2")
	}  else {
			bcpal <- brewer.pal(8, "Dark2")
			cols <- colorRampPalette(bcpal)(nc)
	}
      }
        
	tM <- c()
	for (i in 1:nc) {
		tM <- cbind(tM, om[,c(truthC)])
	}
        if (is.null(labels)) { labels <- plotL}
  if (DEBUG) {print(plotL)}
  pred.GA <- prediction(om[,plotL], tM)
	auc.GA <- as.numeric(performance(pred.GA, 'auc')@y.values)
  if (DEBUG) {print(auc.GA)}
	perf.GA <- performance(pred.GA, 'tpr', 'fpr')
	plot(perf.GA, col=as.list(cols), ...)
  if (legendB) {
	legend("bottomright", paste(labels, sprintf("%.3f", auc.GA), sep=": "), col=cols, pch=15)}
  
        rl <- data.frame(predictor=plotL, auc=auc.GA)
  if (DEBUG) print(rl)
	return(rl)
}





################################
#   Heatmaps

heatmapIndicatedGenes <- function(ct, dataDir, groupingCol, idList, mappingSpecifier, mappingCol, inName, exGroupFactor=NULL, geneGroupFactor=NULL, DEBUG=0,  mapObj, mappedName="Symbol", dispCol="RandomEffects", resSpec="MappedMetaEffect", ... ) {
  #  The ct is a coding table (or subset of a coding table), indicating which datasets to look at
  # groupingCol is the column name for the grouping, usually OriginCode or maybe GSE of a single experiment
  # idList is a list of the ID's of what to show on the heatmap, usually Entrez id's or Homologene ID's, depending
  # mappingSpecifier is
  # dispCol is the column which contains the data which is actually displayed, this is usually the effect size point estimates
  # mappingCol is which is the column to map from
  # expGroupFactor is something to decide how to group things, if null doesn't include this on the experiment side, this is for putting a colored bar above the heatmap
  # geneGroupFactor is something that allows putting a colored bar on the side
  # inName is a note about what is included "MixedSpecies" is one example
  # mapObj is something that has a column which corresponds to the mappingCol (e.g. HID) and another with
  # mappedName -- such as "HumanSymbol" or "Symbol"

  # Need to figure out to map the id's to the symbols, what is the best way?

  uniqOC <- !duplicated(ct[,groupingCol])
  gplID=ct$GPL[uniqOC]
  expID=ct[,groupingCol][uniqOC]
  
  mergedExpressTab <- data.frame(stringsAsFactors = FALSE)
  for (gpl in unique(gplID)) {
    gb <- gplID==gpl  # gpl boolean, which ones are relevant
    for (oc in expID[gb]) {
      evinex <- paste(paste(mappingSpecifier, oc, groupingCol, sep="-"), inName, sep="") 
      if (DEBUG > 1) print(evinex)
       
      ev <- getStoredTable(dataDir, "effectsMeta", resSpec, evinex, DEBUG=DEBUG)
      if (DEBUG) print(head(ev, 2))
      ev <- ev[ev$ID %in% idList,]
      if (dim(ev)[1] > 1) {
        ev$OC <- oc
        mergedExpressTab <- rbind(mergedExpressTab, ev[,c("ID", "OC", dispCol)])
      }
      }
  }
  if (DEBUG) {print("mergedExpressTab") ; print(head(mergedExpressTab))}
  mergedExpressMat <- reshape(mergedExpressTab, v.names=dispCol, idvar="ID", timevar="OC", direction="wide")
  mappedMEM <- merge(mapObj[, c(mappingCol, mappedName)], mergedExpressMat, by.x=mappingCol, by.y="ID", all.x=F, all.y=F)
  rownames(mappedMEM) <- paste(mappedMEM[,c(mappedName)],mappedMEM[,c(mappingCol)], sep=".")
  if(DEBUG) print(head(mergedExpressMat))
  mappedMEM <- mappedMEM[, !colnames(mappedMEM) %in% c(mappedName, mappingCol)]
  colnames(mappedMEM) <- stripNumberPrefix(colnames(mappedMEM), "RandomEffects.")
  try(heatmap(as.matrix(mappedMEM), scale="none", cexCol=2, Colv=NA, Rowv=NA))
  return(mergedExpressMat)
}




################################
#   Disease Associations


getDiseaseAssociatedEntrezOptra <- function(con, phenotypeName, DEBUG=0) {
  # Uses 'like' so the phenotype name can include a wildcard, for prostate cancer something like "Prostat%"
	dbSNPquery <- paste("select distinct(dbSNP)  from var_disease_snp2.Association_Data where Broad_phenotype like '", phenotypeName, "'", sep="")
        if (DEBUG) {print(dbSNPquery)}
        dbsIDs <- dbGetQuery(con, dbSNPquery)
        dbsNums <- as.numeric(data.frame(strsplit(dbsIDs[,1], "rs"))[2,])
        gidquery <- paste("select distinct(GeneID) from var_dbSNP.hs_GeneID where dbSNP in (", paste(dbsNums[!is.na(dbsNums)], collapse=", ") ,")", sep = "")
        #print(gidquery)
        associatedGeneIDs <- dbGetQuery(con, gidquery)$GeneID
        return(associatedGeneIDs)
        
      }







################################
#  Merge expression data within chip, get id expression values

averageToGeneLevelData <- function(nadf, groupID, avef=median) {
  # This function takes in a dataframe of expression values with one column that is a mapping target, usually a gene or a homologene group, and then it simply averages each probe for that group ID (gene) and returns that simplified matrix
  # This assumes that the expression levels are normalized in some way, or are ranks or something like that
  groupList <- list(groupID=factor(nadf[,c(groupID)]))
  names(groupList) <- groupID
    a <- aggregate(data.frame(nadf[,! colnames(nadf) == groupID]), groupList, FUN=avef,na.rm=T)    
    rownames(a) <- a[,groupID]
    a <- a[,!colnames(a) == groupID]
    return(a)
}

relrank <- function(x) {
    rx <- rank(x, na.last="keep")
    rxrel <- rx/length(rx)
    return(rxrel)
}

rankAndScaleAverageToGeneLevelData <- function(anedf, groupID, avef=median) {
  # This function takes in a data.frame of expression values (anedf) with grouping annotations (such as gene level or homologene group) in one column named groupID, and then averages the relative expression within a study, and then scales and centers the genes, assumes here at the rows are the genes and the columns are the experiments
  ma <- data.frame(apply(anedf[,! colnames(anedf) == groupID], 2, relrank))
  ma$xnames <- anedf[,c(groupID)]
  avld <- averageToGeneLevelData(ma, groupID="xnames", avef=avef)
  t(scale(t(avld))) 
}

################################
#  SAM & RankProd

runSAM <- function(ct, dataDir, groupingCol="OriginCode", noclobber=T, DEBUG=0) {
	for (gpl in unique(ct$GPL)) {
		print(gpl)
		ctS <- ct[ct$GPL == gpl & !is.na(ct[,groupingCol]), c("GSM", "ClassCode", groupingCol) ]
		eMat <- getStoredTable(dataDir, "rawExpression", gpl, "RawExpressionValues")
		for (oc in unique(ctS[,groupingCol])) {
			if (!noclobber | !checkStoredTable(dataDir, "significanceValues",  paste(gpl, oc, "SAM", groupingCol, sep="-"), "QValues")) {
                           print(paste(gpl, oc, "SAM", sep="-"))
                            tryCatch({          
				print(c(gpl, oc))
				ctC <- ctS[ctS[,groupingCol] == oc, ]
				gsmList <- ctC$GSM
                                sm <- as.matrix(eMat[,gsmList])
                                if (DEBUG) {
                                  print(dim(sm))
                                  print(head(sm))}
				cl <- ctC$ClassCode-1
                                if (DEBUG) {
                                  print(length(cl))
                                  print(cl)}
				sam.out <- sam(sm, cl, gene.names=rownames(sm))
				sl <- (summary(sam.out, 1e-100))@mat.sig
                                if (DEBUG) {print(head(sl))}
                                if (length(sl$rawp) != sum(is.na(sl$rawp)))
                                  {
                                    writeStoredTable(sl, dataDir, "significanceValues", paste(gpl, oc, "SAM", groupingCol, sep="-"), "QValues")
                                  } else {  print("Error, The SAM results are all NA") }
                                },    error = function(e) print(str(e)) )
}}}}

runRankProd <- function(ct, dataDir, groupingCol="OriginCode", outName="", DEBUG=0, noclobber=TRUE) {
  #compare to computeRankProdSig, here only one side is used
  
  for (gpl in unique(ct$GPL)) {
    print(gpl)
    ctS <- ct[ct$GPL == gpl & !is.na(ct[,groupingCol]), c("GSM", "ClassCode", groupingCol) ]
    #print(ctS)
    if (max(ctS$ClassCode) == 2 & min(ctS$ClassCode) == 1) {
      if (!noclobber | !checkStoredTable(dataDir, "significanceValues", gpl, paste("RankProd-PValues", outName, sep=""))) {
        eMat <- getStoredTable(dataDir, "rawExpression", gpl, "RawExpressionValues")
        #print(head(eMat))
        print(eMat[1:3, 1:4])
        ofac <- factor(ctS[,groupingCol])
        olabels <- as.numeric(ofac)
        classLabels <- ctS$ClassCode - 1   # this is the 1,2 case, turn to 0 & 1 for RP
        if (sum(classLabels != 0 & classLabels != 1) >= 1) {
          print("Potential problem with class coding, expecting 1 & 2")
        }
        
        RP.out <- RP(eMat[,ctS$GSM], classLabels, num.perm = 10, logged = T, na.rm = TRUE, plot = FALSE, rand = 123)
        
        siggenes <- topGene(RP.out,cutoff=1, method="pfp",logged=FALSE)
        
        RP.out <- RPadvance(eMat[,ctS$GSM], classLabels,olabels, num.perm=100, logged=TRUE, rand=123 )
        rptab <- cbind(RP.out$pval, RP.out$AveFC, RP.out$RPrank, RP.out$pfp )
        colnames(rptab) <- c("pvalDown", "pvalUp", "AveFC", "RPRDown", "RPRUp", "PFPDown", "PFPUp")
        writeStoredTable(rptab, dataDir, "significanceValues", gpl, paste("RankProd-PValues", outName, sep=""))
      }
    }
    # Other cases are for two color arrays, but these I guess should be coded with -1
  }
}

runSAMpaired <-  function(ct, dataDir, groupingCol="OriginCode", pairID="pair_id" ,DEBUG=0) {
  if (DEBUG) {print(c("Running paired SAM"))}
  ctR <- ct[!is.na(ct[,c(pairID)]),]
  for (gpl in unique(ctR$GPL)) {
    if(DEBUG) {print(gpl)}
    ctS <- ctR[ctR$GPL == gpl & !is.na(ctR[,groupingCol]), c("GSM", "ClassCode", groupingCol, pairID) ]
    for (oc in unique(ctS[,groupingCol])) {
      if (!checkStoredTable(dataDir, "significanceValues", paste(gpl, oc, "SAMpaired", groupingCol, sep="-"), "QValues")) {
        eMat <- getStoredTable(dataDir, "rawExpression", gpl, "RawExpressionValues")
        print(paste(gpl, oc, "SAM", sep="-"))
        tryCatch({          
          print(c(gpl, oc))
          ctC <- ctS[ctS[,groupingCol] == oc, ]
          
          ctC$xpair <- factor(ctC[,pairID])
          ctC <- ctC[ctC$xpair %in% ctC$xpair[duplicated(ctC$xpair)],] # remove anything not paired up
          ctC$xpair <- as.numeric(factor(ctC$xpair))
          ctC[ctC$ClassCode == 1, c("xpair")] <- - ctC[ctC$ClassCode == 1, c("xpair")]
          gsmList <- ctC$GSM
          sm <- as.matrix(eMat[,gsmList])
          if (DEBUG) {
            print(dim(sm))
            print(head(sm))}

          cl <-  ctC[, c("xpair")]
          if (DEBUG) {
            print(length(cl))
            print(cl)}
          sam.out <- sam(sm, cl, gene.names=rownames(sm))
          sl <- (summary(sam.out, 1e-100))@mat.sig
          if (DEBUG) {print(head(sl))}
          if (length(sl$rawp) != sum(is.na(sl$rawp)))
            {
              writeStoredTable(sl, dataDir, "significanceValues", paste(gpl, oc, "SAMpaired", groupingCol, sep="-"), "QValues")
            } else {  print("Error, The SAM results are all NA") }
        },    error = function(e) print(str(e)) )
}}}}



turnPOneSided <- function(pval, fc) {
  # takes in a list of p-vaules (pvals) and a list of log fold changes (-inf, +inf) and then turns them into the p-values which might have been obtained had it been a one sided test, so that 0.95 is the same as 0.05, just for genes fold changed in the 'opposite direction'
  # if f > 0, newp = p*.5 ; if f == 0, newp=p; if f < 0, newp = 1 - .5*p
  if (fc < 0) { return(1-pval/2)}
  if (fc >= 0) { return(pval/2)}
}


omnibusAcrossSAM <- function(ct, dataDir, outName="", mappingSpecifier="Probe2EntrezMap", mappingCol="GeneID", groupingCol="OriginCode", mname="metaV",  DEBUG=0) {
	mt <- data.frame()
	for (gpl in unique(ct$GPL)) {
		print(gpl)
		ctS <- ct[ct$GPL == gpl & !is.na(ct[,groupingCol]), c("GSM", "ClassCode", groupingCol) ]
		mf <- getStoredTable(dataDir, "geneMappings", gpl, mappingSpecifier, DEBUG=DEBUG)
		mf <- mf[!duplicated(mf$probe), ]
		rownames(mf) <- mf$probe
		for (oc in unique(ctS[,groupingCol])) {
      tryCatch({          
			samM <- getStoredTable(dataDir, "significanceValues", paste(gpl, oc, "SAM", groupingCol, sep="-"), "QValues", DEBUG=DEBUG)
			nt <- mergeViaRowNames(mf, samM, all.a=F, all.b=F)
      nt[,groupingCol] = oc
			mt <- rbind(mt, nt)   },
      error = function(e) print(str(e)) )
		}
	}
  
  if (DEBUG) print(head(mt))
  
	midList <- unique(mt[,mappingCol])
  
	metaRes <- data.frame()
	print(head(midList))
	for (mid in midList) {	
      mtS <- mt[mt[,mappingCol] == mid,]
                
      if(DEBUG > 3) {
        print(mid)
        print(dim(mtS))
        print(head(mtS))
      }
                
		pList <- mtS[, "rawp"]
    qList <- mtS[, "q.value"]
		symbol <- mtS[, "Symbol"][1]
		mpF <- calcOmnibus(pList, method="fisher")
		mpL <- calcOmnibus(pList, method="logit")
      
      # Thresholded values, at the experiment level:
      t.25 <- length(unique(mtS[,groupingCol][qList < .25]))
      t.1 <- length(unique(mtS[,groupingCol][qList < .1]))
      t.001 <- length(unique(mtS[,groupingCol][qList < .001]))
      t.denom <- length(unique(mtS[,groupingCol]))
      #originalPList <- mt[mt[,mappingCol] == mid, "rawp"]
      fcList <-  mt[mt[,mappingCol] == mid, "d.value"]
                                        #  Very important to log the fcList!
      OSpList <- mapply(turnPOneSided, pList, fcList)
      OSmpL <- calcOmnibus(OSpList, method="logit")
      mr <- list(ID=mid, Name=symbol, pvalFisher=mpF, pvalLogit=mpL, pvalOSLogit=OSmpL, t.25=t.25, t.1=t.1, t.001=t.001, t.denom=t.denom)
      if (DEBUG >3) print(data.frame(mr))
		metaRes <- rbind(metaRes, mr)
	}
	writeStoredTable(metaRes, dataDir, "omnibusMeta", paste("FisherAndLogit", "SAM", groupingCol, mname, sep="-"), paste("metaPVal", outName, sep=""))
	return(metaRes)
}



################################
#  SAM, doing MHC at the gene level and taking the min


## Note, problem, what does this mean when the p-value is one sided?  What is the result?, need to get the MHC correction right!  It's okay, taking the smaller value and expanding is okay from a SAM perspective
mappedSAMselector <-function(msm, mappingCol, DEBUG=0) {
  # takes in the column of mapped signficance values from sam and then selects a single sig value by do mt <- data.frame() # where all the values end up
  mt <- data.frame()
  if (DEBUG > 5) {print(head(msm))}

  for (ms in unique(msm[,mappingCol])) {
    msss <- msm[msm[,mappingCol]==ms,]
    if (DEBUG > 5) {print(ms);print(msss)}
    if (dim(msss)[1] == 1) {
      mt <- rbind(mt,msss )
    } else {
      msr <- msss[msss$rawp == min(msss$rawp),][1,]
      msr$rawp <- min(msr$rawp * dim(msss)[1], .9999)
      mt <- rbind(mt, msr)
      if (DEBUG > 5) {print(msr)}
    }

  }
  if (DEBUG) {print(head(mt))}
  return(mt)
}


  

omnibusAcrossSAMselectivelyMapped <- function(ct, dataDir, outName="", mappingSpecifier="Probe2EntrezMap", mappingCol="GeneID", groupingCol=NA, higherGroupingCol=NA, lowerGrouping="OriginCode", mname="metaV", tname="SAM", noclobber=TRUE, DEBUG=0) {
  # groupingCol is just to be backwards compatible, higherGroupingCol is the prefered name
  
  if (is.na(higherGroupingCol)) {higherGroupingCol=groupingCol}
    
  
  if (!noclobber | !checkStoredTable(dataDir, "omnibusMeta", paste("FisherAndLogit", tname, higherGroupingCol, "SelectedMapped", mappingCol, mname, sep="-"), paste("metaPVal", outName, sep=""), DEBUG=DEBUG)) {
    mt <- data.frame()
    for (hgoc in unique(ct[,higherGroupingCol])) {
      if (checkStoredTable(dataDir, "significanceValuesMetaMapped", paste(tname, higherGroupingCol, hgoc, mappingCol, sep="-"), "SelectedMappedRawPValue", DEBUG=DEBUG)) {
        ntSelect <- getStoredTable(dataDir, "significanceValuesMetaMapped", paste(tname, higherGroupingCol, hgoc, mappingCol, sep="-"), "SelectedMappedRawPValue", DEBUG=DEBUG)
      } else {
        ctSM <- ct[ct[,c(higherGroupingCol)] ==hgoc,]
        ot <- data.frame()
        for (gpl in unique(ctSM$GPL)) {
          mf <- getStoredTable(dataDir, "geneMappings", gpl, mappingSpecifier, DEBUG=DEBUG)
          mf <- mf[!duplicated(mf$probe), ]
          rownames(mf) <- mf$probe
          for (oc in unique(ctSM[ctSM$GPL == gpl,c(lowerGrouping)])) {
            samM <- getStoredTable(dataDir, "significanceValues", paste(gpl, oc, tname, lowerGrouping, sep="-"), "QValues", DEBUG=DEBUG)
            nt <- mergeViaRowNames(mf, samM, all.a=F, all.b=F)
            nt[,higherGroupingCol] = oc
            ot <- rbind(ot, nt)
          }
        }
        ntSelect <- mappedSAMselector(ot, mappingCol,  DEBUG=DEBUG)
        writeStoredTable(ntSelect, dataDir, "significanceValuesMetaMapped", paste(tname, higherGroupingCol, hgoc, mappingCol, sep="-"), "SelectedMappedRawPValue", DEBUG=DEBUG)
      }        
      mt <- rbind(mt, ntSelect) 
    }
    
  if (DEBUG) {print(head(mt))}
    
  midList <- unique(mt[,mappingCol])
  metaRes <- data.frame()
  for (mid in midList) {		
    mtS <- mt[mt[,mappingCol] == mid,]
    if(DEBUG > 3) {
      print(mid)
      print(dim(mtS))
      print(head(mtS))
    }
    pList <- mtS[, "rawp"]
    qList <- mtS[, "q.value"]
    symbol <- mtS[, "Symbol"][1]
    mpF <- calcOmnibus(pList, method="fisher")
    mpL <- calcOmnibus(pList, method="logit")
    fcList <-  mt[mt[,mappingCol] == mid, "d.value"]
    OSpList <- mapply(turnPOneSided, pList, fcList)
                                        #if (DEBUG) {print("oneside p values"); print(head(cbind(originalPList, fcList, pList)))}
    OSmpL <- calcOmnibus(OSpList, method="logit")
    mr <- list(ID=mid, Name=symbol, pvalFisher=mpF, pvalLogit=mpL, pvalOSLogit=OSmpL, q.1=sum(qList < .1),q.01=sum(qList < .01),  q.001=sum(qList < .001), p.001=sum(pList<.001), p.00001=sum(pList<.00001), p.000001=sum(pList<.000001),totMeas=length(pList))
    if (DEBUG >3) print(data.frame(mr))
    metaRes <- rbind(metaRes, mr)
  }
    
  if (DEBUG) {print(head(metaRes))}
  writeStoredTable(metaRes, dataDir, "omnibusMeta", paste("FisherAndLogit", tname, higherGroupingCol, "SelectedMapped", mappingCol, mname, sep="-"), paste("metaPVal", outName, sep=""))
  }
  
  metaRes <- getStoredTable(dataDir, "omnibusMeta", paste("FisherAndLogit", tname, higherGroupingCol, "SelectedMapped", mappingCol, mname, sep="-"), paste("metaPVal", outName, sep=""))
  return(metaRes)
}


omnibusAcrossTTselectivelyMapped <- function(ct, dataDir, outName="", mappingSpecifier="Probe2EntrezMap", mappingCol="GeneID", groupingCol=NA, higherGroupingCol=NA, lowerGrouping="OriginCode", mname="metaV", tname="TT", noclobber=TRUE, DEBUG=0) {
  # groupingCol is just to be backwards compatible, higherGroupingCol is the prefered name
  
  if (is.na(higherGroupingCol)) {higherGroupingCol=groupingCol}
  
  
  if (!noclobber | !checkStoredTable(dataDir, "omnibusMeta", paste("FisherAndLogit", tname, higherGroupingCol, "SelectedMapped", mappingCol, mname, sep="-"), paste("metaPVal", outName, sep=""), DEBUG=DEBUG)) {
    mt <- data.frame()
    for (hgoc in unique(ct[,higherGroupingCol])) {
      if (checkStoredTable(dataDir, "significanceValuesMetaMapped", paste(tname, higherGroupingCol, hgoc, mappingCol, sep="-"), "SelectedMappedRawPValue", DEBUG=DEBUG)) {
        ntSelect <- getStoredTable(dataDir, "significanceValuesMetaMapped", paste(tname, higherGroupingCol, hgoc, mappingCol, sep="-"), "SelectedMappedRawPValue", DEBUG=DEBUG)
      } else {
        ctSM <- ct[ct[,c(higherGroupingCol)] ==hgoc,]
        ot <- data.frame()
        for (gpl in unique(ctSM$GPL)) {
          mf <- getStoredTable(dataDir, "geneMappings", gpl, mappingSpecifier, DEBUG=DEBUG)
          mf <- mf[!duplicated(mf$probe), ]
          rownames(mf) <- mf$probe
          for (oc in unique(ctSM[ctSM$GPL == gpl,c(lowerGrouping)])) {
            samM <- getStoredTable(dataDir, "significanceValues", paste(gpl, oc,  lowerGrouping, tname, sep="-"), "PValues", DEBUG=DEBUG)
            nt <- mergeViaRowNames(mf, samM, all.a=F, all.b=F)
            nt[,higherGroupingCol] = oc
            ot <- rbind(ot, nt)
          }
        }
        ntSelect <- mappedSAMselector(ot, mappingCol,  DEBUG=DEBUG)
        writeStoredTable(ntSelect, dataDir, "significanceValuesMetaMapped", paste(tname, higherGroupingCol, hgoc, mappingCol, sep="-"), "SelectedMappedRawPValue", DEBUG=DEBUG)
      }        
      mt <- rbind(mt, ntSelect) 
    }
    
    if (DEBUG) {print(head(mt))}
    
    midList <- unique(mt[,mappingCol])
    metaRes <- data.frame()
    for (mid in midList) {  	
      mtS <- mt[mt[,mappingCol] == mid,]
      if(DEBUG > 3) {
        print(mid)
        print(dim(mtS))
        print(head(mtS))
      }
      pList <- mtS[, "rawp"]
      qList <- mtS[, "q.value"]
      symbol <- mtS[, "Symbol"][1]
      mpF <- calcOmnibus(pList, method="fisher")
      mpL <- calcOmnibus(pList, method="logit")
      fcList <-  mt[mt[,mappingCol] == mid, "d.value"]
      OSpList <- mapply(turnPOneSided, pList, fcList)
      #if (DEBUG) {print("oneside p values"); print(head(cbind(originalPList, fcList, pList)))}
      OSmpL <- calcOmnibus(OSpList, method="logit")
      mr <- list(ID=mid, Name=symbol, pvalFisher=mpF, pvalLogit=mpL, pvalOSLogit=OSmpL, q.1=sum(qList < .1),q.01=sum(qList < .01),  q.001=sum(qList < .001), p.001=sum(pList<.001), p.00001=sum(pList<.00001), p.000001=sum(pList<.000001),totMeas=length(pList))
      if (DEBUG >3) print(data.frame(mr))
      metaRes <- rbind(metaRes, mr)
    }
    
    if (DEBUG) {print(head(metaRes))}
    writeStoredTable(metaRes, dataDir, "omnibusMeta", paste("FisherAndLogit", tname, higherGroupingCol, "SelectedMapped", mappingCol, mname, sep="-"), paste("metaPVal", outName, sep=""))
  }
  
  metaRes <- getStoredTable(dataDir, "omnibusMeta", paste("FisherAndLogit", tname, higherGroupingCol, "SelectedMapped", mappingCol, mname, sep="-"), paste("metaPVal", outName, sep=""))
  return(metaRes)
}



################################
#  Add Annotations

addAnnotatedColumnBooleansByNameColumn <- function(sm, dbConn, annList, DEBUG=0) {
	for (geneSet in annList) {
		gsQuery <- paste("select Symbol from user_alexmo.broad_molecular_signatures where GeneSet = '", geneSet, "'", sep="", collapse="")
	
		SymbolList <- dbGetQuery(dbConn, gsQuery)$Symbol
		Boole <- rep(1, length(SymbolList))
		a <- data.frame(SymbolList, Boole)
		colnames(a) <- c("Symbol", geneSet)
		sm <- merge(sm, a, by.x="Name", by.y="Symbol", all.x=T)
		sm[is.na(sm[,geneSet]),geneSet] <- 0
		#print(cumsum(sm[,geneSet]))
	}
	return(sm)
	
	
}


addAnnotatedColumnsForGOcodes <- function(sm, dbConn, goCodeList, colname="ID", DEBUG=0) {
  #  GO codes are like "GO:0008150" but a list of text forms
	for (geneSet in goCodeList) {
		gsQuery <- paste("select GeneID, GO_term from annot_gene.gene2go where GO_ID = '", geneSet, "'", sep="", collapse="")
                print(gsQuery) #  This helps keep track of when the queries are taking place, can be kind of a pain, but better than nothing
                dbres <- dbGetQuery(dbConn, gsQuery)
		IDList <- dbres$GeneID
                goname <-  gsub("\\(", ".",  gsub("\\)", ".",  gsub(",", ".",  gsub("\\+", ".",  gsub(":", ".",  gsub("\\-", ".", gsub(" ", ".", gsub("_", ".", dbres$GO_term[1]))))))))
                sm[,c(goname)] <- sm[,c(colname)] %in% IDList
	}
	return(sm)
	
	
}


getAnnotColumnForName<- function(sm, dbConn, geneSet, mapName="Name", DEBUG=0) {
	gsQuery <- paste("select Symbol from alexmo.broad_molecular_signatures where GeneSet = '", geneSet, "'", sep="", collapse="")
	SymbolList <- dbGetQuery(dbConn, gsQuery)$Symbol
	Boole <- rep(1, length(SymbolList))
	a <- data.frame(SymbolList, Boole)
	colnames(a) <- c("Symbol", geneSet)
	sm$Name <- toupper(sm[,mapName])
	sm <- merge(sm, a, by.x=mapName, by.y="Symbol", all.x=T)
	sm[is.na(sm[,geneSet]),geneSet] <- 0
	return(sm)
}


process_raw_data <- function(filePath="", compressed = F , dest_path = "rawExpression",  method="expresso", normalize_method="quantiles",bgcorrect_method="rma",pmcorrect_method="pmonly",summary_method="avgdiff"){
  
  library(affy) #
  
  mydata <- ReadAffy(celfile.path = filePath, compress = compressed) #pattern: one treatment vs one or more vechicles
  
  if (method == "rma"){
    eset <- rma(mydata)
  }else if (method == "gcrma"){
    library(gcrma)
    eset <- gcrma(mydata)
  }else if (method == "mas5"){
    eset <- mas5(mydata)
  }else if (method == "dchip"){
    eset <- expresso(mydata, normalize.method= "invariantset",
                     bgcorrect.method= "none", pmcorrect.method= "subtractmm",
                     summary.method= "liwong") 
  }else {
    eset <- expresso(mydata, normalize.method= normalize_method,
                     bgcorrect.method= bgcorrect_method, pmcorrect.method= pmcorrect_method,
                     summary.method= summary_method) #method "wong" would not coverge    
  }
  
  exprs_mat <- exprs(eset)
  
  #write.table(exprs_mat, dest_path, sep="\t", quote=F, row.names=T, col.names=T)
  return(exprs_mat)
}

getGene2Unigene <- function(con){
  gene2unigene <- dbReadTable(con, "annot_gene.gene2unigene")
  return(gene2unigene)
}

getEntrezMappingsFromGEO <- function(gpl_id){
  library(Biobase)
  library(GEOquery)
  gpl <- getGEO(gpl_id, destdir=".")
  if (gpl_id == "GPL4133"){
    geneMappings <- Table(gpl)[,c("ID","GENE","GENE_SYMBOL","GENE_SYMBOL")]   # the name of GPL4133 is messed up     
  }else if (gpl_id == "GPL10558"){
    geneMappings <- Table(gpl)[,c("ID","Entrez_Gene_ID","Symbol","Definition")]        
  }else if (gpl_id == "GPL10687"){
    geneMappings <- Table(gpl)[,c("ID","EntrezGeneID","GeneSymbol","GeneSymbol")]   
 # }else if (gpl_id == "GPL571" | gpl_id == "GPL3921"){ #this two gpls' annotation is not right from the website, need to re-processed, so we use from our server
    #geneMappings <- Table(gpl)[,c("ID","ENTREZ_GENE_ID","Gene Symbol","Gene Title")]
  #  geneMappings <- getEntrezMappings(gpl_id, con)
 #   geneMappings$Name <- geneMappings$Symbol # name does not return
  }else if (gpl_id == "GPL10687"){
    geneMappings <- Table(gpl)[,c("ID","EntrezGeneID","GeneSymbol","GeneSymbol")]
  }else if (gpl_id == "GPL96" | gpl_id == "GPL1352"){
    geneMappings <- Table(gpl)[,c("ID","ENTREZ_GENE_ID","Gene Symbol",	"Gene Symbol")]
   }else if (gpl_id == "GPL8177"){
    geneMappings <- Table(gpl)[,c("ID","UGCluster","Symbol","Name")]
    gene2unigene <- getGene2Unigene(con)
    geneMappings <- merge(geneMappings, gene2unigene, by.x = "UGCluster", by.y="Unigene" )
    geneMappings <- subset(geneMappings, select = c("ID", "GeneID", "Symbol", "Name"))
  }else {
    #try to guess column names of entrez gene id and symbol name
    putative_fields <- grep("ENTREZ|GENE|Symbol", names(Table(gpl)), ignore.case=T)
    house_keeping_gene_ids = c("23521", "7316", "567")
    house_keeping_gene_symbols = c("RPL13A", "UBC", "B2M")
    gene_field_name = ""
    symbol_field_name = ""
    for (field in putative_fields){
      if (sum(Table(gpl)[,field] %in% house_keeping_gene_ids) > 0 ){
        gene_field_name = names(Table(gpl))[field]
      }
    }
 
    for (field in putative_fields){
      if (sum(Table(gpl)[,field] %in% house_keeping_gene_symbols) > 0 ){
        symbol_field_name = names(Table(gpl))[field]
      }
    }     
    
    if (symbol_field_name != "" & gene_field_name != ""){
       geneMappings <- Table(gpl)[,c("ID", gene_field_name, symbol_field_name,symbol_field_name)]
    }
  }

  names(geneMappings) <- c("probe","GeneID","Symbol","NAME")
  geneMappings <- subset(geneMappings, !is.na(GeneID) & GeneID != "")
  dupProbes <- unique(geneMappings$probe[duplicated(geneMappings$probe)])
  geneMappings <- geneMappings[! geneMappings$probe %in% dupProbes,] #remove duplicated probes
  rownames(geneMappings) <- geneMappings$probe
  
  #some mapping files, GENEID composed by muliple gene ids, should split.. this may cause problems for the symbol and name
  if (class(geneMappings$GeneID) != "integer"){
    newGeneMappings = data.frame()
    for (i in 1:nrow(geneMappings)){
      newGeneMappings = rbind(newGeneMappings, data.frame(probe = geneMappings$probe[i], 
                                                          GeneID = unique(na.omit(as.numeric(unlist(strsplit(unlist(geneMappings$GeneID[i]), "[^0-9]+"))))),
                                                          Symbol = geneMappings$Symbol[i], 
                                                          NAME = geneMappings$NAME[i]))
    }
    geneMappings = newGeneMappings
  }
  
  return(geneMappings)
}

getEntrezMappings_fromAE <- function(array.id, con, DEBUG=0) {
  query <- paste("select a.probe, a.GeneID, b.Symbol from expr_array_express.Expt_probe_Id1 as a join annot_gene.gene_info_hs as b using (GeneID)  where ExptID=\"", array.id,"\"", sep="" )
  if (DEBUG) {print(query)}
  GeneMappings <- dbGetQuery(con, query)
  dupProbes <- unique(GeneMappings$probe[duplicated(GeneMappings$probe)])
  GeneMappings <- GeneMappings[! GeneMappings$probe %in% dupProbes,]
  rownames(GeneMappings) <- GeneMappings$probe
  return(GeneMappings)  
}


post.process.eset <- function(emat){
  comparison_frame = as.data.frame(emat)
  values_logged <- is_logged(comparison_frame)
  if (!values_logged) {
    comparison_frame <- log2(comparison_frame+0.001) # add 1 to avoid log2(0)
    values_logged <- TRUE
  }
  
  #normalization
  comparison_frame1 <- normalize.quantiles(as.matrix(comparison_frame ))
  row.names(comparison_frame1) <- row.names(comparison_frame)
  colnames(comparison_frame1) <- colnames(comparison_frame)
  comparison_frame <- comparison_frame1
  
  return(comparison_frame)
}

is_logged <- function(sample_frame) {
  mean_sample_var <- mean(apply(sample_frame,2,var,na.rm=T),na.rm=T)
  ifelse(mean_sample_var > 10,FALSE,TRUE)
}








