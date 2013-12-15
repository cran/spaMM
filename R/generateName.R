generateName <-
function(base="tmp") { ## for an object within the .GlobalEnv ##... not used in distributable code
   pattern <- paste(base,"*",sep="")
   allmatches <- ls(pattern=pattern,pos=".GlobalEnv")
   allremainders <- substring(allmatches,nchar(base)+1) 
   allremainders <- as.numeric(allremainders[which( ! is.na(as.numeric(allremainders )))  ]) ## as.numeric("...")
   if (length(allremainders) == 0) {
            num <- 0
   } else num <- max(allremainders)+1
   validname <-paste ( base , num , sep="") 
   return(validname)
}
