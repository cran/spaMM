pastefrom <-
function(..., callstack=sys.calls(),prefix="From ",full.stack=.SpaMM$MESSAGES.FULL.STACK) {  ## callstack is a list of *calls*
  if (full.stack) {
    cs <- clean_cs(callstack) ## cs is a string made of the names of function called in the call stack
  } else {
    lcs <- length(callstack)
    cs <- paste(callstack[[lcs-1]])[1] ## [[lcs-1]] because [[lcs]] is call to pastefrom.
  }
  paste(prefix,cs,": ", ...,sep="")
}
