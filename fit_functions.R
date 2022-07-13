library(tercen)
library(tercenApi)
library(dplyr, warn.conflicts = FALSE)

library(drda)
library(nplr)

curvefit <- function( ctx, lib='nplr' ){
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
  operatorProps <- ctx$query$operatorSettings$operatorRef$propertyValues
  
  for( prop in operatorProps ){
    if (prop$name == "Fitting Library"){
      lib <- prop$value
    }
  }
  
  
  logCol <- ctx$cselect() 
  names(logCol) <- cnames.without.ns
  logCol <- mutate(logCol, .ci=seq(0,nrow(logCol)-1))
  
  res <- ctx %>%  
    select( .y, .ci, .ri, 'Sample type'  ) %>%  
    left_join(logCol, by=".ci") %>%
    group_by(.ri) %>%
    do( do.curvefit(., lib ) ) %>%
    ctx$addNamespace() 
  
  
  colNames <- names(res)
  
  # NOTE
  # do() yields an error when returning a list of data frames.
  # To circumvent that, the resulting tables are merged with cbind(), which duplicates the rows in the sumamry table
  # This issue is addressed by issuing an unique id as the summary's first column (summary_index)
  # The code below separates the table and removes duplicates from the results
  
  # NOTE 2
  # The amount of columns is hard-coded. 
  # IF something changes in the resulting tables, these numbers need to be changed here as well
  summaryDf <- res %>% select(append('.ri', colNames[11:length(colNames)] ) ) 
  summaryDf <- summaryDf[!duplicated(summaryDf[,colNames[11]]),] %>% 
    select(-colNames[11])
  
  dataDf <- res %>% select( colNames[1:10]  ) 
  
  result = OperatorResult$new()
  result$tables = list(tercen::dataframe.as.table(dataDf), tercen::dataframe.as.table(summaryDf))
  return(result)
}

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
  
  
  conc_u<-approx( qcY, qcXp, xout = y_prediction )
  conc_w<-approx( qcY, qcXwp, xout = y_prediction )
  
  resp_u <- conc_u$x
  conc_u <- conc_u$y
  
  resp_w <- conc_w$x
  conc_w <- conc_w$y
  
  
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
  
  
  coeff   <- getPar(mdlU  )$params
  coeff_w <- getPar(mdlW  )$params
  
  if(npar == 5){
    y_prediction <- (coeff[['bottom']] + 
                       (coeff[['top']] - coeff[['bottom']])/
                       ((1 + 10^(coeff[['scal']]*(coeff[['xmid']]-x_prediction) ) )^coeff[['s']])) * maxY  
    y_prediction_w <- (coeff_w[['bottom']] + 
                         (coeff_w[['top']] - coeff_w[['bottom']])/
                         ((1 + 10^(coeff_w[['scal']]*(coeff_w[['xmid']]-x_prediction) ) )^coeff_w[['s']])) * maxY  
    
  }else if(npar == 4){
    y_prediction <- (coeff[['bottom']] + 
                       (coeff[['top']] - coeff[['bottom']])/
                       ((1 + 10^(coeff[['scal']]*(coeff[['xmid']]-x_prediction) ) ))) * maxY
    
    y_prediction_w <- (coeff_w[['bottom']] + 
                         (coeff_w[['top']] - coeff_w[['bottom']])/
                         ((1 + 10^(coeff_w[['scal']]*(coeff_w[['xmid']]-x_prediction) ) ))) * maxY
    
  }
  
  conc_u<-approx( qcY, qcXp, xout = y_prediction )
  conc_w<-approx( qcY, qcXwp, xout = y_prediction )
  
  resp_u <- conc_u$x
  conc_u <- conc_u$y
  
  resp_w <- conc_w$x
  conc_w <- conc_w$y
  
  
  rowIdx <- rep( unique(df$.ri)[1], length(conc_u) )
  colIdx <- rep( unique(df$.ci)[1], length(conc_u) )

  outDf <- data.frame(
    concentrationU=conc_u,
    concentrationW=conc_w,
    responseU=resp_u,
    responseW=resp_w,
    x_predicted=x_prediction,
    y_predicted=y_prediction,
    diff=(1-(conc_w/conc_u))*100,
    npar=npar,
    .ri=rowIdx,
    .ci=colIdx
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

