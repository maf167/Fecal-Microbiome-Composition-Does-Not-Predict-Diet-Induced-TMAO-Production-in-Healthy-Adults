#
# VarExplained and container objects
# 
# Author: Marc Ferrell

#######################################
#
#  varExplained - Train EN model and test on held out data
#    - dat:        list of x.train (matrix), x.test(matrix), 
#                  y.train (vector), y.test (vector)
#    - pred.nm:    names of predictor group for later search
#    - outcome.nm: names of outcome for later search
#
#  CVsummary    - Plot training results as a function of alpha and lambda 
#
#######################################

varExplained <- function(dat, pred.nm, outcome.nm){ 
  require(caret)
  
  x.train <- dat[[1]]
  x.test  <- dat[[2]]
  y.train <- dat[[3]]
  y.test  <- dat[[4]]
  ## x is a matrix of predictors (rows are samples), y is a vector of outcome values
  
  # Search grid for alpha and lambda
  alpha.grid <- seq(0,1,0.1)
  lambda.grid <- seq(0,12, length=120)
  search.grid <- expand.grid(.alpha = alpha.grid, .lambda = lambda.grid)
  
  ## CV method is 10-fold
  cv.method <- trainControl(
    method = "cv",
    number = 10
  )
  
  # Cross validate values for lambda and alpha
  my.train <- train(x.train, y.train, 
                    method = "glmnet",
                    metric = "RMSE",
                    tuneGrid = search.grid,
                    trControl = cv.method,
                    maxit = 5000
  )
  
  # Use the trained model to predict the test data
  # my.pred <- predict(my.train$finalModel, newx = x.test, 
  # s=my.train$bestTune$lambda)
  
  betas <- coef(my.train$finalModel, s=my.train$bestTune$lambda)
  
  my.pred <- apply(x.test,1,function(x) betas[1] + sum(betas[2:length(betas)] * x,na.rm=TRUE) )
  
  # Return the model coefs, rsq, and RMSE
  suppressWarnings(
    t <- summary(lm(y.test ~ my.pred))
  )
  
  pred.p <- ifelse(length(unique(my.pred)) < 2, 1, t$coefficients[2,4])
  plt = matrix(c(my.pred,y.test ),ncol=2)
  colnames(plt) <- c("Predicted","Actual")
  ( structure(
    list(coefs = betas,
         rsq = c(rsq = t$adj.r.squared, p = pred.p),
         RMSE = sqrt(mean((my.pred - y.test)^2, na.rm=TRUE)), model=my.train$finalModel,s=my.train$bestTune$lambda,
         myPlot = plt,
         paramResults = my.train$results,
         Ntest  = length(y.test),
         Ntrain = length(y.train)),
    class      ="varExplained",
    predictors = pred.nm,
    outcome    = outcome.nm)
  )
  # }
}

CVsummary <- function(x){
  if(class(x) != "varExplained") stop("Object must be of class varExplained")
  require(ggplot2)
  require(reshape2)
  require(gridExtra)
  dat <- x$paramResults
  dat <- dat[,-grep("SD",colnames(dat))]
  
  datm <- melt(dat, id.vars=c("alpha","lambda"))
  
  p1 <- ggplot(datm[datm$variable == "RMSE",], aes(x=lambda, y=alpha)) +
    geom_raster(aes(fill = value)) + theme_minimal() + xlab("") +
    ggtitle(paste(attr(x, "predictors"), " Predicting ", attr(x, "outcome"))) + 
    scale_fill_gradientn(colors=c("blue","white","red"), name="RMSE") +
    theme(plot.title=element_text(hjust = 0.5), legend.title.align=0.5)

  p2 <- ggplot(datm[datm$variable == "Rsquared",], aes(x=lambda, y=alpha)) +
    geom_raster(aes(fill = value)) + theme_minimal() + xlab("") +
    scale_fill_gradientn(colors=c("blue","white","red"), name="Rsquared") +
    theme(legend.title.align=0.5)
  
  p3 <- ggplot(datm[datm$variable == "MAE",], aes(x=lambda, y=alpha)) +
    geom_raster(aes(fill = value)) + theme_minimal() + 
    scale_fill_gradientn(colors=c("blue","white","red"), name="MAE") +
    theme(legend.title.align=0.5)
  
  grid.arrange(p1,p2,p3, nrow=3)
    
}

varExplainedLM <- function(dat, pred.nm, outcome.nm){
  x.train <- dat[[1]]
  x.test  <- dat[[2]]
  y.train <- dat[[3]]
  y.test  <- dat[[4]]
  
  model <- lm(y.train ~ x.train)
  
  betas <- coef(model)
  
  my.pred <- apply(x.test,1,function(x) betas[1] + sum(betas[2:length(betas)] * x,na.rm=TRUE) )
  
  # Return the model coefs, rsq, and RMSE
  suppressWarnings(
    t <- summary(lm(y.test ~ my.pred))
  )
  ( structure(
    list(coefs = betas,
         rsq = c(rsq = ifelse(t$adj.r.squared > 0, t$adj.r.squared, 0), p = t$coefficients[2,4]),
         RMSE = sqrt(mean((my.pred - y.test)^2, na.rm=TRUE)),
         model=model,s=NA,
         myPlot = matrix(c(my.pred,y.test ),ncol=2),
         paramResults = NA),
    class      ="varExplained",
    predictors = pred.nm,
    outcome    = outcome.nm)
  )
  
  
}




#######################################
#
#  VEList - a list of VarExplained objects with the attribute class="VEList"
#    - VEQuery: method to query VE objects by attributes for error stats
#        - vel:     VEList
#        - pred:    regex to search vel pred attributes
#        - outcome: regex to search vel outcome attributes
#        - returns indices 
#    - VESum:   method to generate an error summary table from a VEList
#        - vel:     VEList
#        - returns long-format data frame
#    - geom_VESum: ggplot settings to plot VESum table
#
#######################################

VEQuery <- function(vel, pred=".", outcome="."){
  
  out = c()
  
  for( i in 1:length(vel)){
    p.match = grepl(pred, attr(vel[[i]], "predictors"))
    o.match = grepl(outcome, attr(vel[[i]], "outcome"))
    
    if(all(c(p.match,o.match))) out <- c(out, i)
    
  }
  
  (out)
  
}

VESum <- function(vel){
  out <- data.frame(Outcome = character(),Pred = character(),
                    Err=numeric(), Err_Type = factor(levels=c("RMSE","Rsq")), stringsAsFactors = FALSE)
  
  for( i in 1:length(vel)){
    veli <- vel[[i]]
    for( j in 1:length(veli)){
      if(!is.na(veli[[j]])){
        my.o <- attr(veli[[j]], "outcome")
        my.p <- attr(veli[[j]], "predictors")
        
        
        out[nrow(out)+1,1:ncol(out)] <- c(my.o,my.p,veli[[j]]$RMSE,"RMSE")
        out[nrow(out)+1,1:ncol(out)] <- c(my.o,my.p,veli[[j]]$rsq["rsq"],"Rsq")
        
      }
    }
   
  }
  out[,"Err"] <- as.numeric(out[,"Err"])
  # out[is.na(out$Err), "Err"] <- 0
  (out)
}

geom_VESum <- function(tit){
  list(
    geom_violin(scale="width"),
    geom_boxplot(width=0.3),
    theme_bw(base_size = 20),
    ggtitle(tit),
    ylab(bquote("Adjusted"~R^2~"(Actual vs. Predicted)")),
    xlab("Model"),
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          axis.text.x = element_text(angle = 50,hjust = 1,size=rel(1.3)),
          axis.text.y = element_text(size=rel(1.3)))
  )
  
  
  
  
}

#######################################
#
# Example Plot
#
# ggplot(SumTab[fullTab.S$Lab == "Plasma.TMAO...µM." & 
#                    fullTab.S$Err_Type =="Rsq",], 
#        aes(x=Pred, y=Err, group=Pred)) + 
#   geom_VESum("")
#
#######################################

