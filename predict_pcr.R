#predict response
#using clincal feature + gene expression + cell type

expr = read.csv( paste("geneExpression/",  "breast_combat.csv", sep=""), row.names = 1)
expr = t(expr)

#clincal feature
meta = read.delim("meta_input.txt", sep="\t", stringsAsFactors = F)

#cell type
cell = read.csv("sample_cell_type.csv", stringsAsFactors = F)
a = apply(cell[, c("drug_cisplatin", "drug_capecitabine", "drug_anthracycline", "drug_cyclophosphamide")], 1, function(x){paste(unique(x), collapse=",")})

drugs = c("drug_cisplatin", "drug_capecitabine", "drug_cyclophosphamide" ) #colnames(cell)[grep("drug_", colnames(cell))]
clinic = c("er","her2","pr","grade","stage")
cell_types = colnames(cell)[which(colnames(cell) == "aDC"):ncol(cell)]

train_data = merge(data.frame(expr, GSM=rownames(expr)), cell[, c("GSM", clinic, cell_types, "prc")], by="GSM")
pcr = train_data$prc
train_data = train_data[, -1] #remove GSM

train_data = as.data.frame(unclass(train_data), stringsAsFactors=T)

x_train <- model.matrix( ~ . -(prc), data = train_data,  contrasts.arg = lapply(train_data[, colnames(train_data) %in% c(clinic)], contrasts, contrasts=FALSE) )

x = as.matrix(x_train) #feature_columns
y = NA
y[pcr == "pcr_yes"]  = 1
y[pcr== "pcr_no"]  = 0

x = x[!is.na(y), ]
y = y[!is.na(y)]


cvfit = cv.glmnet(x, y, alpha = 0.2, standardize=F, parallel = TRUE, family = "binomial", type.measure = "auc")
max(cvfit$cvm)
pred_feature = as.matrix(coef(cvfit, s = "lambda.min"))
best_lambda = cvfit$lambda.min
pred_feature = data.frame(feature = rownames(pred_feature), pred_coeff = pred_feature)
pred_feature = subset(pred_feature, abs(X1) > 0)
pred_feature = pred_feature[order((pred_feature$X1)), ]
