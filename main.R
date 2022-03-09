library(tercen)
library(dplyr, warn.conflicts = FALSE)

library(drda)
library(nplr)


.drda_fit <- function(df, npar=5){
  x <- df$.x
  y <- df$.y
  
  maxY <- max(y)
  y <- y / maxY
  
  stdX <- x[df$`Sample type` == 'Standard']
  stdY <- y[df$`Sample type` == 'Standard']
  
  qcX <- x[df$`Sample type` == 'QC']
  qcY <- y[df$`Sample type` == 'QC']
  
  
  if(npar==5){
    meanFunc <- 'logistic5'
  }else if(npar == 4){
    meanFunc <- 'logistic4'
  }
  
  
  mdl <- drda(y ~ x, data = data.frame( x=stdX, y=stdY ), 
              mean_function = meanFunc,
              is_log=TRUE)
  
  
  w <- 1 / ( (stdY)**2  )
  mdlWeight <- drda(y ~ x, data = data.frame( x=stdX, y=stdY ), 
                    mean_function = meanFunc,
                    weights = w, 
                    is_log=TRUE)
  
  coeff <- unlist( mdl$coefficients, use.names = FALSE)
  
  if(npar==5){
    qcYp <- ((coeff['alpha'] + (coeff['beta'] -coeff['alpha']  ) / 
                (1 + coeff['nu'] * exp(-coeff['eta'] * (qcX - coeff['phi'])))^(1 / coeff['nu'])  ) ) * maxY
  }else if(npar == 4){
    qcYp <- ((coeff['alpha'] + (coeff['beta'] -coeff['alpha']  ) / 
                (1  * exp(-coeff['eta'] * (qcX - coeff['phi']))) ) ) * maxY
  }
  
  
  coeff <- unlist( mdlWeight$coefficients, use.names = FALSE)
  if(npar==5){
    qcYwp <- ((coeff['alpha'] + (coeff['beta'] -coeff['alpha']  ) / 
                (1 + coeff['nu'] * exp(-coeff['eta'] * (qcX - coeff['phi'])))^(1 / coeff['nu'])  ) ) * maxY
  }else if(npar == 4){
    qcYwp <- ((coeff['alpha'] + (coeff['beta'] -coeff['alpha']  ) / 
                (1  * exp(-coeff['eta'] * (qcX - coeff['phi']))) ) ) * maxY
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






.nplr_fit <- function(df, npar=5){
  x <- df$.x
  y <- df$.y
  
  maxY <- max(y)
  y <- y / maxY
  
  stdX <- x[df$`Sample type` == 'Standard']
  stdY <- y[df$`Sample type` == 'Standard']
  
  qcX <- x[df$`Sample type` == 'QC']
  qcY <- y[df$`Sample type` == 'QC']
  

  mdlU <- nplr(stdX, stdY, npars=npar, useLog=FALSE, silent = TRUE)
  
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
  
  
  
  # getFitValues(mdlU)
  # getInflexion(mdlU)
  mdlW <- nplr(stdX, stdY, npars=npar, useLog=FALSE, silent = TRUE,
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


do.curvefit <- function(df, lib, npar=5){
  
  if( lib == 'drda'){
    outDf <- .drda_fit(df, npar)
  }else if( lib == 'nplr' ){
    outDf <- .nplr_fit(df, npar)
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