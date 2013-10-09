spMMexpandSlash <-
function (bb) {
    if (!is.list(bb)) 
        return(spMMexpandSlash(list(bb)))
    unlist(lapply(bb, function(x) {
        if (length(x) > 2 && is.list(trms <- slashTerms(x[[3]]))) {
          mess <- pastefrom("Check code for interactions.",prefix="(!) From ")
          stop(mess)
## 'original' code was          
#            return(lapply(unlist(makeInteraction(trms)), function(trm) substitute(foo | 
#                bar, list(foo = x[[2]], bar = trm))))
## but this requires lem4 and was never tested
        }        
        x
    }))
}
