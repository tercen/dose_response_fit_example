library(tercen)
library(tercenApi)


#http://127.0.0.1:5402/admin/w/22ae949dc1a3dd3daf96768225009600/ds/d72f5099-6ecf-4944-8f33-99d0ef0e8909
options("tercen.workflowId" = "22ae949dc1a3dd3daf96768225009600")
options("tercen.stepId"     = "d72f5099-6ecf-4944-8f33-99d0ef0e8909")

source( 'fit_functions.R' )

# =======================================================
#
#               Operator entry point
#
# =======================================================
ctx = tercenCtx()

result <- curvefit( ctx )
ctx$save( result )

