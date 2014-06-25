`^.blockDiag` <-
function(scal,bdobject) {
  bd <- lapply(1:length(bdobject),function(v){scal^bdobject[[v]]})
  attributes(bd) <- attributes(bdobject)
  bd
}
