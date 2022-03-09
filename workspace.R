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
  
  stdX <- x[df$`Sample type` == 'Standard']
  stdY <- y[df$`Sample type` == 'Standard']
  
  maxY <- max(stdY)
  minY <- min(stdY)
  
  stdY <- stdY / maxY
  
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
  x <- df$log
  y <- df$.y
  
  stdX <- x[df$`Sample type` == 'Standard']
  stdY <- y[df$`Sample type` == 'Standard']
  
  maxY <- max(stdY)
  minY <- min(stdY)
  
  stdY <- (stdY - minY)/ (maxY-minY)
  
  qcX <- x[df$`Sample type` == 'QC']
  qcY <- y[df$`Sample type` == 'QC']
  
  
  mdlU <- nplr(stdX, stdY, npars=npar, useLog=FALSE, silent = TRUE)
  
  coeff <- getPar(mdlU  )$params
  
  if(npar == 5){
    qcYp <- (coeff[['bottom']] + 
               (coeff[['top']] - coeff[['bottom']])/
               ((1 + 10^(coeff[['scal']]*(coeff[['xmid']]-qcX) ) )^coeff[['s']])) * maxY  
  }else if(npar == 4){
    qcYp <- (coeff[['bottom']] + 
               (coeff[['top']] - coeff[['bottom']])/
               ((1 + 10^(coeff[['scal']]*(coeff[['xmid']]-qcX) ) ))) * maxY
  }
  
  
  
  # getFitValues(mdlU)
  # getInflexion(mdlU)
  
  mdlW <- nplr(stdX, stdY, npars=npar, useLog=FALSE, silent = TRUE,
               method='gw', LPweight=2)
  coeff <- getPar(mdlW  )$params
  if(npar == 5){
    qcYwp <- (coeff[['bottom']] + 
                (coeff[['top']] - coeff[['bottom']])/
                ((1 + 10^(coeff[['scal']]*(coeff[['xmid']]-qcX) ) )^coeff[['s']])) * maxY  
  }else if(npar == 4){
    qcYwp <- (coeff[['bottom']] + 
                (coeff[['top']] - coeff[['bottom']])/
                ((1 + 10^(coeff[['scal']]*(coeff[['xmid']]-qcX) ) ))) * maxY
  }
  

  outDf <- data.frame(
    .ri=df$.ri[df$`Sample type` == 'QC'],
    .ci=df$.ci[df$`Sample type` == 'QC'],
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
    outDf <- .drda_fit(df, npar=4)
    outDf <- rbind( outDf, .drda_fit(df, npar=5) )
  }else if( lib == 'nplr' ){
    outDf <- .nplr_fit(df, npar=4)
    outDf <- rbind( outDf, .nplr_fit(df, npar=5) )
  }
  
  outDf <- outDf %>%
    mutate(across(npar, as.integer)) %>%
    as_tibble()
  
  return(outDf)
}

# =======================================================
#
#               Operator entry point
#
# =======================================================
ctx = tercenCtx()

rowNames <- ctx$rnames

cnames.with.ns <- ctx$cnames


cnames.without.ns <- unlist(lapply( cnames.with.ns, function(x) {
  if( !startsWith(x, '.') ){
    x<-strsplit(x, '\\.')[[1]][2]
  }else{
    x
  }
}))


colorNames <- ctx$colors

if( !("Assay ID" %in% rowNames) ){
  stop("Row 'Assay ID' is mandatory."  )
}


if( !("log" %in% cnames.without.ns) ){
  stop("Column 'log' is mandatory."  )
}

if( !("Sample type" %in% colorNames) ){
  stop("Color 'Sample type' is mandatory."  )
}

# Read in operator parameters
lib  <- 'nplr'

operatorProps <- ctx$query$operatorSettings$operatorRef$propertyValues

for( prop in operatorProps ){
  if (prop$name == "Fitting Library"){
    lib <- prop$value
  }
}

col <- ctx$cselect() 

names(col) <- cnames.without.ns

col <- mutate(col, .ci=seq(0,nrow(col)-1))

ctx %>%  
  select( .y, .ci, .ri, 'Sample type'  ) %>%  
  left_join(col, by=".ci") %>%
  group_by(.ri) %>%
  do( do.curvefit(., lib ) ) %>%
  ctx$addNamespace() %>%
  ctx$save()