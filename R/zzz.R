# @author Paul Bailey
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("Dire v", utils::packageDescription("Dire")$Version, "\n"))
}

globalVariables(c("dopari", "testID", "subtestID", "nodes"))
