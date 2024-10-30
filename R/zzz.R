GAStartupMessage <- function()
{
  # Startup message obtained as
  # > figlet GA
  msg <- c(paste0(
    "  ____    _
 / ___|  / \\     Genetic
| |  _  / _ \\    Algorithms
| |_| |/ ___ \\
 \\____/_/   \\_\\  version ", utils::packageVersion("TGGATES")),
"\nType 'citation(\"TGGATES\")' for citing this R package in publications.")
  return(msg)
}

.onAttach <- function(libname, pkgname)
{
  # unlock .ga.default variable allowing its modification
#  unlockBinding(".get_TGGATES.default", asNamespace("TGGATES"))
  # startup message
  msg <- GAStartupMessage()
  if(!interactive())
    msg[1] <- paste("Package 'TGGATES' version", utils::packageVersion("TGGATES"))
  packageStartupMessage(msg)
  invisible()
}

.onLoad <- function(libname, pkgname) {
  data("comp_data", package=pkgname, envir=parent.env(environment()))
  invisible()
}
