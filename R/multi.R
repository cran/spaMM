multi <-
function(binResponse=c("npos","nneg"),binfamily=binomial(),input="types",...) {
  return(c(list(family="multi",binResponse=binResponse,binfamily=binfamily,input=input),list(...)))
}
