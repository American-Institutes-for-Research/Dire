skip_on_cran()
context("Parsing")
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 
## Initial Parse an existing dct file
dct <- system.file("trialData", "NAEPprimer.dct", package="Dire")
parsed <- quiet(parseNAEPdct(dct)) 
## write that dct file back out
to_dct <- system.file("trialData", "test.dct", package="Dire")
# read that dct now back
parsedBack <- quiet(parseNAEPdct(to_dct))
# see if the two are the same
expect_equal(nrow(setdiff(parsed$dichotParamTab, parsedBack$dichotParamTab)), 0)
expect_equal(nrow(setdiff(parsed$polyParamTab, parsedBack$polyParamTab)), 0)
expect_equal(nrow(setdiff(parsed$testDat, parsedBack$testDat)), 0)
expect_equal(nrow(setdiff(parsed$surveyParameters, parsedBack$surveyParameters)), NULL)
