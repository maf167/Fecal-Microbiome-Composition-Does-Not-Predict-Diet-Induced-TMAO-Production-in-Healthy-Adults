#######################################
#
#  fitModels
#   -- supporting fxn: makeSets
#
#######################################

getpnm <- function(aa){
  bb <- gsub("\\.", " ", aa)
  (gsub("p", "+", bb))
}



fitModels <- function(dat, pred, y){
  # x.base  - observation matrix
  # pred    - named list of predictor sets
  # y       - matrix of outcomes (cols)
  
  ModelList <- structure(list(), class="VEList")
  
  ## Fit Models
  for(i in 1:length(pred)){
    sets <- makeSets(dat, pred[[i]], y)
      
    mod <- try(
      if(length(pred[[i]]) > 1){
          varExplained(sets, 
                     pred.nm = getpnm(names(pred)[i]), 
                     outcome.nm=colnames(dat)[y])
              }
      else{
          sets[[1]] <- sets[[1]][,1]
          varExplainedLM(sets, pred.nm = getpnm(names(pred)[i]), 
                         outcome.nm=colnames(dat)[y])
      }
           )
    if(class(mod) == "try-error") mod <- NA
    
    ModelList[[length(ModelList)+1]] <- mod
    
  }

  (ModelList)
    
}  
  
fitModels2 <- function(x.base, x.lst, y, PP.base, PP.lst, sets){
  # x.base  - observation matrix
  # x.lst   - named list of obs matrices
  # y       - matrix of outcomees (cols)
  # PP.base - preprocessing function for base matrix
  # PP.lst  - named list of preprocessing functions (names must match x.lst)
  
  ### Check input args
  
  class.lst <- sapply(list(x.base, x.lst, y, PP.base, PP.lst,sets), 
                      function(x) {class(x)[1]})
  if(class.lst != c("matrix", "list", "matrix", "function", "list", "list")) {
    stop("Bad input")
  }
  ## Fit Base Model

  for( i in 1:ncol(y)){
    dat <- makeSets(y[,i,drop=FALSE], filepaths, PP.base, x.base)
    OUT[[length(OUT)+1]] <-
      varExplained(dat, pred.nm = "Base", outcome.nm=colnames(y)[i])
  }

  ## Fit x.lst models alone
  for(i in 1:ncol(y)){
    for(k in 1:length(x.lst)){
      dat <- makeSets(y[,i,drop=FALSE], filepaths, PP.base=NA,
                      x.base=NA, PP.lst=PP.lst, x.lst=list(x.lst[[k]]))

      if(ncol(dat[[1]]) > 1){
        OUT[[length(OUT)+1]] <-
          varExplained(dat, pred.nm = names(x.lst)[k], outcome.nm=colnames(y)[i])
      }
      else{
        OUT[[length(OUT)+1]] <-
          varExplainedLM(dat, pred.nm = names(x.lst)[k], outcome.nm=colnames(y)[i])
      }

    }
  }

  ## Fit Base + Each x.lst

  for(i in 1:ncol(y)){
    for(k in 1:length(x.lst)){
      dat <- makeSets(sets[[i]], y[,i], filepaths, PP.base,
                      x.base, PP.lst=PP.lst, x.lst=list(x.lst[[k]]))
      OUT[[length(OUT)+1]] <-
        varExplained(dat, pred.nm = paste("Base +",names(x.lst)[k]),
                     outcome.nm=colnames(y)[i])
    }
  }

  ## Fit Full Model

  for( i in 1:ncol(y)){
    dat <- makeSets(sets[[i]], y[,i], filepaths, PP.base, x.base,
                    PP.lst=PP.lst, x.lst=x.lst)
    OUT[[length(OUT)+1]] <-
      varExplained(dat, pred.nm = "Full", outcome.nm=colnames(y)[i])
  }

  (OUT)

}


