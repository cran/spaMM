noReForm <-
function(re.form) {
  (!is.null(re.form) && !is(re.form,"formula") && is.na(re.form)) ||
    (is(re.form,"formula") && length(re.form)==2 && identical(re.form[[2]],0))
}
