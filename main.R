library(tercen)
library(dplyr, warn.conflicts = FALSE)

library(drda)
library(nplr)
source( 'fit_functions.R' )


# =======================================================
#
#               Operator entry point
#
# =======================================================
ctx = tercenCtx()

result <- curvefit( ctx )
ctx$save( result )


