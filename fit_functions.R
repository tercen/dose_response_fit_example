library(drda)
library(nplr)



.drda_fit <- function(df, npar=5, summary_index=1){
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
  
  
  
  

  qcXp  <- unlist(lapply(qcY, function(yhat) approx(x = mdl$fitted.values, y = stdX, xout = yhat/maxY)$y ))
  qcXwp <- unlist(lapply(qcY, function(yhat) approx(x = mdlWeight$fitted.values, y = stdX, xout = yhat/maxY)$y ))
 
  x_prediction <- 10^(seq(-1,1,0.05) )
  
  x_prediction <- append(x_prediction, qcXp  )
  x_prediction <- append(x_prediction, qcXwp  )
  x_prediction <- sort(x_prediction)
  
  x_prediction<- x_prediction[x_prediction>=min(stdX) & x_prediction<=max(stdX) ]
  
  y_prediction <- predict(mdl, x_prediction)*maxY
  
  
  conc_u <- y_prediction * NA
  conc_w <- y_prediction * NA
  resp_u <- y_prediction * NA
  resp_w <- y_prediction * NA
  
  
  for( i in seq(1, length(x_prediction))){
    if( x_prediction[i] %in% qcXp ){
      conc_u[i] <- x_prediction[i]
      resp_u[i] <- y_prediction[i]
    }
    
    if( x_prediction[i] %in% qcXwp ){
      conc_w[i] <- x_prediction[i]
      resp_w[i] <- y_prediction[i]
    }
  }
  
  
  rowIdx <- rep( unique(df$.ri)[1], length(conc_u) )
  colIdx <- rep( unique(df$.ci)[1], length(conc_u) )
  
  outDf <- data.frame(
    .ri=rowIdx,
    .ci=colIdx,
    concentrationU=conc_u,
    concentrationW=conc_w,
    responseU=resp_u,
    responseW=resp_w,
    x_predicted=x_prediction,
    y_predicted=y_prediction,
    diff=(1-(conc_w/conc_u))*100,
    npar=npar
  ) 
  
  coeffs <- unlist( mdl$coefficients, use.names = TRUE)
  coeffsW <- unlist( mdlWeight$coefficients, use.names = TRUE)
  
  sumDf <- data.frame(
    summary_index= summary_index,
    aucU = nauc( mdl ),
    aucW = nauc( mdlWeight ),
    gofU = 1-(sigma(mdl)/sd(stdY)),
    gofW = 1-(sigma(mdlWeight)/sd(stdY)),
    fitpar=npar,
    paramsU=paste(names(coeffs),coeffs,sep="=",collapse=", " ),
    paramsW=paste(names(coeffsW),coeffsW,sep="=",collapse=", " )
  )
  
  
  
  return( list(outDf, sumDf) )
}




.nplr_fit <- function(df, npar=5, summary_index=0){
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
  
  mdlW <- nplr(stdX, stdY, npars=npar, useLog=FALSE, silent = TRUE,
               method='gw', LPweight=2)
  
  
  qcXp  <- getEstimates(mdlU,qcY/maxY)$x
  qcXwp <- getEstimates(mdlW,qcY/maxY)$x
  

  x_prediction <- 10^(seq(-1,1,0.05) )
  
  x_prediction <- append(x_prediction, qcXp  )
  x_prediction <- append(x_prediction, qcXwp  )
  x_prediction <- sort(x_prediction)
  
  x_prediction<- x_prediction[x_prediction>=min(stdX) & x_prediction<=max(stdX) ]
  
  
  coeff <- getPar(mdlU  )$params
  
  if(npar == 5){
    y_prediction <- (coeff[['bottom']] + 
               (coeff[['top']] - coeff[['bottom']])/
               ((1 + 10^(coeff[['scal']]*(coeff[['xmid']]-x_prediction) ) )^coeff[['s']])) * maxY  
  }else if(npar == 4){
    y_prediction <- (coeff[['bottom']] + 
               (coeff[['top']] - coeff[['bottom']])/
               ((1 + 10^(coeff[['scal']]*(coeff[['xmid']]-x_prediction) ) ))) * maxY
  }
  
  
  
  conc_u <- y_prediction * NA
  conc_w <- y_prediction * NA
  resp_u <- y_prediction * NA
  resp_w <- y_prediction * NA
  
  
  for( i in seq(1, length(x_prediction))){
    if( x_prediction[i] %in% qcXp ){
      conc_u[i] <- x_prediction[i]
      resp_u[i] <- y_prediction[i]
    }
    
    if( x_prediction[i] %in% qcXwp ){
      conc_w[i] <- x_prediction[i]
      resp_w[i] <- y_prediction[i]
    }
  }
  
  
  rowIdx <- rep( unique(df$.ri)[1], length(conc_u) )
  colIdx <- rep( unique(df$.ci)[1], length(conc_u) )

  outDf <- data.frame(
    .ri=rowIdx,
    .ci=colIdx,
    concentrationU=conc_u,
    concentrationW=conc_w,
    responseU=resp_u,
    responseW=resp_w,
    x_predicted=x_prediction,
    y_predicted=y_prediction,
    diff=(1-(conc_w/conc_u))*100,
    npar=npar
  ) 
  
  
  
  sumDf <- data.frame(
    summary_index= summary_index,
    aucU = getAUC( mdlU ),
    aucW = getAUC( mdlW ),
    gofU = getGoodness(mdlU),
    gofW = getGoodness(mdlW),
    fitpar=npar,
    paramsU=paste(names(getPar(mdlU)),getPar(mdlU),sep="=",collapse=", " ),
    paramsW=paste(names(getPar(mdlW)),getPar(mdlW),sep="=",collapse=", " )
  )
  
  
  
  return( list(outDf, sumDf) )
}

do.curvefit <- function(df, lib){
  
  if( lib == 'drda'){
    outDf <- .drda_fit(df, npar=4, summary_index=1)
    outDf2 <- .drda_fit(df, npar=5, summary_index=2)

    data <- rbind( outDf[[1]], outDf2[[1]] )
    summary <- rbind( outDf[[2]], outDf2[[2]] )
    
    outDf <- cbind(data,summary)
    
  }else if( lib == 'nplr' ){
    outDf <- .nplr_fit(df, npar=4, summary_index=1)
    outDf2 <- .nplr_fit(df, npar=5, summary_index=2)

    data <- rbind( outDf[[1]], outDf2[[1]] )
    summary <- rbind( outDf[[2]], outDf2[[2]] )
    
    outDf <- cbind(data,summary)
  }

  return(outDf)
}

