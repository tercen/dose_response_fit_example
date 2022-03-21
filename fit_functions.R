library(drda)
library(nplr)



.drda_fit <- function(df, npar=5){
  x <- df$log
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
  
  #coeff <- unlist( mdl$coefficients, use.names = FALSE)
  #browser()
  

  qcXp  <- unlist(lapply(qcY, function(yhat) approx(x = mdl$fitted.values, y = stdX, xout = yhat/maxY)$y ))
  qcXwp <- unlist(lapply(qcY, function(yhat) approx(x = mdl$fitted.values, y = stdX, xout = yhat/maxY)$y ))

  
  rowIdx <- rep( unique(df$.ri)[1], length(qcXp) )
  colIdx <- rep( unique(df$.ci)[1], length(qcXp) )
  
  outDf <- data.frame(
    .ri=rowIdx,
    .ci=colIdx,
    concentrationU=qcXp,
    concetrationW=qcXwp,
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
  
  
  browser()
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
    #mutate(across(npar, as.integer)) %>%
    as_tibble()
  
  return(outDf)
}

