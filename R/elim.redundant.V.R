elim.redundant.V <-
function(vertices) { ## removes redundant vertices
  if (nrow(vertices)<=ncol(vertices)) { ## convhulln crashes!
    minimalvertices<-vertices ## stupid case, should never occur
  } else if (ncol(vertices)==1) {
    minimalvertices<-array(c(min(vertices),max(vertices)),dim=c(2,1))
  } else {
    minimalvertices<-vertices[unique(as.numeric(convhulln(vertices, "Pp"))),] ## removes redundant vertices
  }
  colnames(minimalvertices)<-colnames(vertices)
  minimalvertices
}
