Predictor <-
function (formula, offset=NULL, LMatrix = NULL, ZLMatrix = NULL) {
  if (substr((as.character(formula[2])),1,5)=="cbind") { 
       positives <- strsplit(strsplit(as.character(formula[2]),"\\(")[[1]][2],",")[[1]][1] ## names of variables
       negatives <- strsplit(strsplit(as.character(formula[2]),",")[[1]][2],"\\)")[[1]][1] ## names of variables
       ## cbind syntax for binomial model => is converted to alternative syntax; cbind would fail in HLframes as marked there (chnage HLframessome day ?)
       oriFormula <- formula
       formula <- as.formula(paste(positives,"~",as.character(formula[3])))
       BinDenForm <- paste(positives,"+",negatives)       
  } else {BinDenForm <-NULL;oriFormula <- NULL }
  res <- list(formula=formula,offset=offset, LMatrix=LMatrix, ZLMatrix=ZLMatrix,BinDenForm=BinDenForm,oriFormula=oriFormula)
  return(res)
}
