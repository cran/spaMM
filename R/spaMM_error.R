.hack_options_error <- local({
  .R_errorfn_or_expr <- NULL ## keep original error-handling function, "stop" or "recover" or rstudio stuff or NULL
  .message <- NULL ## the message to show when an error occurs
  function(message=NULL) {
    if ( ! is.null(message)) {
      if (is.null(.message)) { ## Do not call options("error"=errorWrap) twice !
        # see ?options for error:
        # either [not set or] a function or an expression governing the handling of non-catastrophic errors such as those generated by 
        # stop as well as by signals and internally detected errors. If the option is a function, 
        # a call to that function, with no arguments, is generated as the expression.
        orig_call_or_expr <- getOption("error") ## original error-handling function call [or NULL]: save it to be able to restore it,
        if ( inherits(orig_call_or_expr,"call")) { # The messy thing is that if one assigns error= as function, getOption() return a call, whic his improper to assign again error=...
          orig_errorfn_or_expr <- orig_call_or_expr[[1]]
        } else { # "NULL", or quoted expression (tested by inherits(orig_call_or_expr,"{")) in which case .R_errorfn[[1L]] does not work.
          orig_errorfn_or_expr <- orig_call_or_expr
        }
        .R_errorfn_or_expr <<- orig_errorfn_or_expr 
        # Local copy 'orig_errorfn_or_expr' is for R CMD check not complaining about the 'errorWrap' definition.
        # I might call it '.R_errorfn_or_expr' and assign .R_errorfn_or_expr <<- .R_errorfn_or_expr and use .R_errorfn_or_expr' in 'errorWrap' def
        # But this would be less clear as it would use the same name for variables in two envirs, local and parent env.  
        errorWrap <- function() {
          message(paste0("An error occurred, possibly due to the following earlier issue:\n", .message))
          options(error=orig_errorfn_or_expr) 
          # cleans everything; otherwise non-null .message will interfere with calls to .hack_options_error() even after an error has been signalled 
          .R_errorfn_or_expr <<- NULL
          .message <<- NULL
          eval(orig_call_or_expr) ## stop or recover or quoted expression
        }
        options("error"=errorWrap) 
      }
      .message <<- message ## possibly updating the message.
    } else if ( ! is.null(.message)) { ## cleans everything
      options("error"=.R_errorfn_or_expr) ## [[]] gets the original error-handling function (without the trailing '()') from the call and restores it
      .R_errorfn_or_expr <<- NULL
      .message <<- NULL
    }
  }
})