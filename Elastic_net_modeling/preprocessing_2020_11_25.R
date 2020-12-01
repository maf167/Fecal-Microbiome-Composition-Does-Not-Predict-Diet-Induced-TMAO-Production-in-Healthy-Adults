#######################################
#
#  Preprocessing functions
# 
#######################################


mergeByRows <- function(a){
  # a is a list of matrix-like objects
  o <- a[[1]]
  for(i in 2:length(a)){
    o <- merge(o,a[[i]],by=0,all=TRUE)
    rownames(o) <- o[,"Row.names"]
    o <- o[-1]
  }
  
  (o[order(rownames(o)),])
}


ppCLIN <- function(CLINi){
  num_fields <- c(
    "Age.At.Diet.Start",
    "Diet.Start.Date",
    "TG.Ave.",
    "TC.Ave.",
    "LDL.Ave.",
    "HDLC.Ave.",
    "ApoA1.Ave.",
    "ApoB.Ave.",
    "BMI",
    "Per_BF",
    "SBP.Avg.",
    "DBP.Avg.",
    "Endopat_AI",
    "Endopat_AI_75Per",
    "Hip..cm.",
    "Waist..cm.",
    "WaistIC..cm.",
    "Wt...Lbs..Avg.",
    "Height..cm.",
    "Diet..Compliance"
  )
  suppressWarnings(
    cat_fields <- sapply(1:ncol(CLINi), function(x){
      all(unique(CLINi[,x])[order(unique(CLINi[,x]))] == c(0,1))
    }))
  
  CLIN_cat <- CLINi[,cat_fields]
  
  CLIN_num <- sapply(num_fields, function(x){
    eval(parse(text=sprintf("CLINi[,'%s'] <- scale(CLINi[,'%s'])",x,x)))
  })
  
  (cbind(CLIN_num, CLIN_cat))
  
}

pp <- function(ct, test=FALSE){
  ## ct is a gene count table normalized by gene length and lib size or
  ## a formated clinical table (dates and categorical variables re-formatted)
  
  # Filter fields if absence > 30% unless for testing set
  ## testing set predictors should match training set exactly
  cat_fields <- sapply(1:ncol(ct), function(x){
    all(unique(ct[,x])[order(unique(ct[,x]))] == c(0,1))
  })
  ct.cont <- ct[,which(!cat_fields),drop=FALSE]
  ct.cat  <- ct[,which(cat_fields),drop=FALSE]
  
  
  
  if(!test){
    ab <- apply(ct.cont, 2, function(x) (sum(x==0,na.rm=TRUE) / length(x)))

    ct.cont.fil <- ct.cont[,ab <= 0.3, drop=FALSE]
  }
  else{
    ct.cont.fil <- ct.cont
  }
  
  ct.final <- scale(ct.cont.fil)
  rownames(ct.final) <- rownames(ct.cont.fil)
  colnames(ct.final) <- colnames(ct.cont.fil)
  
  ( cbind(ct.final, ct.cat))
  
}


my.complete.cases <- function(l){
  # l is a named list of data frames
  # returns same list with only rows containing no missing data
  
  # make a list with the rownames of each table
  l <- l[!is.na(l)]
  allSamps <- lapply(l, rownames)
  
  # find the rownames shared among all tables
  sharedSamps <- intersect(allSamps[[1]], allSamps[[2]])
  if(length(l) > 2){
    for(i in 3:length(allSamps)){ sharedSamps <- intersect(sharedSamps, allSamps[[i]]) }
  }
  
  sharedSamps <- sharedSamps[order(sharedSamps)]
  
  # make tables containing only shared samps
  ll <- lapply(l, function(x) x[sharedSamps,])
  
  # 1. concatenate all tables
  # 2. find subset of shared samps with no missing data
  completeSamps <- sharedSamps[complete.cases(do.call("cbind", ll))]
  
  # return list of tables with matching rownames and complete data
  ( lapply(l, function(x) x[completeSamps,]) )
  
}


#######################################
#
#  startFrac
# 
#######################################

stratFrac <- function(s, y){
  
  ## s is a named vector of values between 0 and 1 that sum to 1
  ## y is a vector of the observed values being fractionated
  
  if(sum(s)!= 1){"s must sum to 1"}
  else{
    OUT <- list()
    yind <- 1:length(y)
    
    # Round down and distribute values according to s
    for( i in 1:length(s)){
      if(length(yind) > 1) n <- sample(yind, size=floor(s[i]*length(y))) else n <- yind[1]
      OUT[[i]] <- n
      yind <- yind[!(yind %in% n)]
    }
    
    # Distribute remaining values equally and randomly
    i = 1
    OUTind <- 1:length(OUT)
    while(length(yind) > 0){
      if(length(yind) > 1) n <- sample(yind, size=1) else n <- yind[1]
      Pn <- sample(OUTind, size=1)
      OUT[[OUTind[Pn]]] <- c(OUT[[OUTind[Pn]]], n)
      yind <- yind[!(yind %in% n)]
      if(i < length(OUT)) i = i + 1 else i <- 1
    }
    
    for( i in OUT){ message(length(i)/length(y))}
    names(OUT) <- names(s)
    ( OUT )
  }
}

#######################################
#
#  makeSets
# 
#######################################

makeSets <- function(dat, x, y){
  
  samp.xy <- complete.cases(dat[,c(x,y)])
  
  dat.x <- dat[samp.xy,x,drop=FALSE]
  dat.y <- dat[samp.xy,y]
  
  set <- stratFrac(c(train=0.7,test=0.3), dat.y)
  
  x.train <- pp(dat.x[set$train,,drop=FALSE])
  x.test  <- pp(dat.x[set$test,,drop=FALSE], test=TRUE)
  x.test  <- x.test[,colnames(x.train),drop=FALSE]
  
  y.train      <- dat.y[set$train]
  y.test       <- dat.y[set$test]
  
  (list(x.train,x.test,y.train,y.test))
}



