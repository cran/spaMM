generateFileName <-
function(base="tmp",ext="") { ## for a file
   pattern <- paste(base,"*",ext,sep="")
   allmatches <- dir(pattern=pattern)
   allremainders <- substring(allmatches,nchar(base)+1)
   allremainders <- unlist(strsplit(allremainders,ext)) ## removes the extension from the remainder 
   allremainders <- as.numeric(allremainders[which( ! is.na(as.numeric(allremainders )))  ]) ## as.numeric("...")
   if (length(allremainders) == 0) {
            num <- 0
   } else num <- max(allremainders)+1
   validFileName <-paste ( base , num , ext,sep="") 
   return(validFileName)
}
