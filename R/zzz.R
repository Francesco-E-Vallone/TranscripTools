.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "Welcome to TranscripTools!\n",
    "\n",
    "Documentation: https://francesco-e-vallone.github.io/TranscripTools/"
  )
}

utils::globalVariables(c("Gene", ".group", "Expression"))
