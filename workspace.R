library(tercen)
library(dplyr, warn.conflicts = FALSE)

library(drda)
library(nplr)

# http://127.0.0.1:5402/admin/w/22ae949dc1a3dd3daf96768225009600/ds/d72f5099-6ecf-4944-8f33-99d0ef0e8909
options("tercen.workflowId" = "22ae949dc1a3dd3daf96768225009600")
options("tercen.stepId"     = "d72f5099-6ecf-4944-8f33-99d0ef0e8909")

.drda_fit <- function(df, npar=5){
  x <- df$.x
  y <- df$.y
  
  maxY <- max(y)
  y <- y / maxY
  
  stdX <- x[df$`Sample type` == 'Standard']
  stdY <- y[df$`Sample type` == 'Standard']
  
  qcX <- x[df$`Sample type` == 'QC']
  qcY <- y[df$`Sample type` == 'QC']
  
  
  
  mdl <- drda(y ~ x, data = data.frame( x=stdX, y=stdY ), 
              mean_function = 'logistic5',
              is_log=TRUE)
  
  
  w <- 1 / ( (stdY)**2  )
  mdlWeight <- drda(y ~ x, data = data.frame( x=stdX, y=stdY ), 
                    mean_function = 'logistic5',
                    weights = w, 
                    is_log=TRUE)
  
  coeff <- unlist( mdl$coefficients, use.names = FALSE)
  
  qcYp <- ((coeff['alpha'] + (coeff['beta'] -coeff['alpha']  ) / 
              (1 + coeff['nu'] * exp(-coeff['eta'] * (qcX - coeff['phi'])))^(1 / coeff['nu'])  ) ) * maxY
  
  coeff <- unlist( mdlWeight$coefficients, use.names = FALSE)
  qcYwp <- ((coeff['alpha'] + (coeff['beta'] -coeff['alpha']  ) / 
               (1 + coeff['nu'] * exp(-coeff['eta'] * (qcX - coeff['phi']) ))^(1 / coeff['nu'])  ) ) * maxY
  
  rowIdx <- rep( unique(df$.ri)[1], length(qcYp) )
  colIdx <- rep( unique(df$.ci)[1], length(qcYp) )
  
  outDf <- data.frame(
    .ri=rowIdx,
    .ci=colIdx,
    responseU=qcYp,
    responseW=qcYwp,
    logConcentration=qcX,
    diff=(1-(qcYwp/qcYp))*100,
    npar=npar
  ) 
  
  return(outDf)
}






.nplr_fit <- function(df, npar=5){
  x <- df$.x
  y <- df$.y
  
  maxY <- max(y)
  y <- y / maxY
  
  stdX <- x[df$`Sample type` == 'Standard']
  stdY <- y[df$`Sample type` == 'Standard']
  
  qcX <- x[df$`Sample type` == 'QC']
  qcY <- y[df$`Sample type` == 'QC']
  
  # 4PL
  # formula y = B + (T-B) / (1 + exp(b*(xmid-x)))
  
  # 5PL
  # formula y = B + (T-B) / (1 + exp(b*(xmid-x)))^S
  # where B and T are the bottom and top asymptotes, and b, xmid and s are the Hill slope, the x-coordinate
  # at the inflexion point and an asymetric coefficient, respectively.
  
  mdlU <- nplr(stdX, stdY, npars=5, useLog=FALSE, silent = TRUE)
  
  coeff <- getPar(mdlU  )$params
  
  if(npar == 5){
    qcYp <- (coeff[['bottom']] + 
               (coeff[['top']] - coeff[['bottom']])/
               ((1 + exp(coeff[['scal']]*(coeff[['xmid']]-qcX) ) )^coeff[['s']])) * maxY  
  }else if(npar == 4){
    qcYp <- (coeff[['bottom']] + 
               (coeff[['top']] - coeff[['bottom']])/
               ((1 + exp(coeff[['scal']]*(coeff[['xmid']]-qcX) ) ))) * maxY
  }
  
  
  
  #getFitValues(mdlU)
  # getInflexion(mdlU)
  mdlW <- nplr(stdX, stdY, npars=5, useLog=FALSE, silent = TRUE,
               method='gw', LPweight=2)
  
  if(npar == 5){
    qcYwp <- (coeff[['bottom']] + 
                (coeff[['top']] - coeff[['bottom']])/
                ((1 + exp(coeff[['scal']]*(coeff[['xmid']]-qcX) ) )^coeff[['s']])) * maxY  
  }else if(npar == 4){
    qcYwp <- (coeff[['bottom']] + 
                (coeff[['top']] - coeff[['bottom']])/
                ((1 + exp(coeff[['scal']]*(coeff[['xmid']]-qcX) ) ))) * maxY
  }
  
  rowIdx <- rep( unique(df$.ri)[1], length(qcYp) )
  colIdx <- rep( unique(df$.ci)[1], length(qcYp) )
  
  outDf <- data.frame(
    .ri=rowIdx,
    .ci=colIdx,
    responseU=qcYp,
    responseW=qcYwp,
    logConcentration=qcX,
    diff=(1-(qcYwp/qcYp))*100,
    npar=npar
  ) 
  
  return(outDf)
}

do.curvefit <- function(df, lib){
  
  if( lib == 'drda'){
    outDf <- .drda_fit(df)
  }else if( lib == 'nplr' ){
    outDf <- .nplr_fit(df)
  }
  
  
  return(outDf)
}


ctx = tercenCtx()

rowNames <- ctx$rnames
colorNames <- ctx$colors

if( !("Assay ID" %in% rowNames) ){
  error("Row 'Assay ID' is mandatory."  )
}

if( !("Sample type" %in% colorNames) ){
  error("Color 'Sample type' is mandatory."  )
}


lib  <- 'nplr'
npar <- 5

ctx %>%
  select( .y, .x, .ci, .ri, 'Sample type'  ) %>%
  group_by(.ri) %>%
  do( do.curvefit(., lib, npar) ) %>%
  ctx$addNamespace() %>%
  ctx$save()