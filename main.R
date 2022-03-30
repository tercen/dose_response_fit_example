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

ctx$save( result )

