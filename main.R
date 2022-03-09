library(tercen)
library(dplyr, warn.conflicts = FALSE)

library(drda)
library(nplr)


.drda_fit <- function(df){
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
    npar=5
  ) 
  
  return(outDf)
}

.nplr_fit <- function(df){
  x <- df$.x
  y <- df$.y
  
  maxY <- max(y)
  y <- y / maxY
  
  stdX <- x[df$`Sample type` == 'Standard']
  stdY <- y[df$`Sample type` == 'Standard']
  
  qcX <- x[df$`Sample type` == 'QC']
  qcY <- y[df$`Sample type` == 'QC']
  
  # formula y = B + (T-B) / (1 + exp(b*(xmid-x)))^S
  # where B and T are the bottom and top asymptotes, and b, xmid and s are the Hill slope, the x-coordinate
  # at the inflexion point and an asymetric coefficient, respectively.
  
  mdlU <- nplr(stdX, stdY, npars=5, useLog=FALSE, silent = TRUE)
  
  coeff <- getPar(mdlU  )$params
  
  qcYp <- (coeff[['bottom']] + 
             (coeff[['top']] - coeff[['bottom']])/
             ((1 + exp(coeff[['scal']]*(coeff[['xmid']]-qcX) ) )^coeff[['s']])) * maxY
  
  
  #getFitValues(mdlU)
  # getInflexion(mdlU)
  mdlW <- nplr(stdX, stdY, npars=5, useLog=FALSE, silent = TRUE,
               method='gw', LPweight=2)
  
  qcYwp <- (coeff[['bottom']] + 
              (coeff[['top']] - coeff[['bottom']])/
              ((1 + exp(coeff[['scal']]*(coeff[['xmid']]-qcX) ) )^coeff[['s']])) * maxY
  
  rowIdx <- rep( unique(df$.ri)[1], length(qcYp) )
  colIdx <- rep( unique(df$.ci)[1], length(qcYp) )
  
  outDf <- data.frame(
    .ri=rowIdx,
    .ci=colIdx,
    responseU=qcYp,
    responseW=qcYwp,
    logConcentration=qcX,
    diff=(1-(qcYwp/qcYp))*100,
    npar=5
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


lib <- 'nplr'

ctx %>%
  select( .y, .x, .ci, .ri, 'Sample type'  ) %>%
  group_by(.ri) %>%
  do( do.curvefit(., lib) ) %>%
  ctx$addNamespace() %>%
  ctx$save()