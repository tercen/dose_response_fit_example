library(tercen)
library(tercenApi)
library(dplyr, warn.conflicts = FALSE)
library(tim)
library(jsonlite)

library(stringi)
# Function definitions
source('fit_functions.R')



# http://127.0.0.1:5402/admin/w/22ae949dc1a3dd3daf96768225009600/ds/d72f5099-6ecf-4944-8f33-99d0ef0e8909/inTable/main
wkfId <- "22ae949dc1a3dd3daf96768225009600"
options("tercen.workflowId"=wkfId)
stepIdList <- c("d72f5099-6ecf-4944-8f33-99d0ef0e8909")



propDictList <- list()
# Steps with properties
# propDictList <- list("Default"=list(stepId="a18c19cb-52e7-427c-a1fa-ad7af76dce95"),
                     # "Mult2"=list(stepId="a18c19cb-52e7-427c-a1fa-ad7af76dce95",
                                  # multiplicator=2),
                     # "AllProp"=list(stepId="a18c19cb-52e7-427c-a1fa-ad7af76dce95",
                                    # multiplicator=4,
                                    # power1="1",
                                    # power2="2")
# )


for( i in seq(1, length(stepIdList))){
  options("tercen.stepId"=stepIdList[i])  
  ctx = tercenCtx()
  
  # Get step name
  step_name <- get_step_name(ctx, wkfId, stepIdList[i])
  
  print(step_name)
  
  
  has_property_list <- lapply(  propDictList, function(x){
    x$stepId == stepIdList[i]
  })
  
  if( any(unname(unlist(has_property_list)))){
    
    prop_list_idx <- which(unname(unlist(has_property_list)))
    for(j in prop_list_idx){
      
      props <- propDictList[j]
      test_suff <- unlist(unname(names(props)))[1]
      props <- props[[test_suff]]
      props$stepId <- NULL
      
      plist <- props
      pnames <- unlist(unname(names(props)))
      
      
      step_name_ex <- paste0(step_name, "_", test_suff)
      result <- curvefit( ctx, !!!setNames(plist,pnames) )
      
      
      tim::build_test_data( result, ctx, paste0(step_name_ex, "_r2"),
                            version = '1.2.0',
                            r2 = 0.9,
                            skipCols=c(),
                            props=setNames(plist,pnames))
    }
  }else{
    result <- curvefit( ctx )
    
    tim::build_test_data( result, ctx, paste0(step_name, "_r2"),
                          version = '1.2.0',
                          r2 = 0.9, gen_schema=TRUE,
                          skipCols=c())
  }
  
  
}


