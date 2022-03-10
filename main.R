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


ctx %>%  
  select( .y, .ci, .ri, 'Sample type'  ) %>%  
  left_join(logCol, by=".ci") %>%
  group_by(.ri) %>%
  do( do.curvefit(., lib ) ) %>%
  ctx$addNamespace() %>%
  ctx$save()