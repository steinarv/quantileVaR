CAViaR_optim <- function(y, Model=1, prob=0.05){
  
  if(prob>=0 & prob<=1)  cavObj <- new( CAViaR, y, prob)
  else	stop("Probability level must be between 0 and 1")
  ###############
  cavObj$VAR[1] <- quantile(y, prob)
  ##############
  if(Model==1){
    tryVals <- matrix(runif(10000*3,0,1),ncol=3)
    RQscores <- cavObj$testInitials(tryVals,1)
    bestVals <- tryVals[order(RQscores), ][1:10, ]
  }else if(Model==2){
    tryVals <- matrix(runif(100000*4,0,1),ncol=4)
    RQscores <- cavObj$testInitials(tryVals,2)
    bestVals <- tryVals[order(RQscores), ][1:15, ]
  }else if(Model==3){
    tryVals <- matrix(runif(10000*3,0,1),ncol=3)
    RQscores <- cavObj$testInitials(tryVals,3)
    bestVals <- tryVals[order(RQscores), ][1:10, ]
  }else if(Model==4){
    tryVals <- matrix(runif(10000*1,0,1),ncol=1)
    RQscores <- cavObj$testInitials(tryVals,4)
    bestVals <- matrix(tryVals[order(RQscores), ][1:5], ncol=1)
  }else if(Model==5){
    tryVals <- matrix(runif(10000*4,0,1),ncol=4)
    RQscores <- cavObj$testInitials(tryVals,5)
    bestVals <- tryVals[order(RQscores), ][1:10, ]
  }else 
    stop("Unknown input: Model")
  #############
  RQoptim <- cbind(NA, bestVals,NA)
  if(Model==4){
  met <- "Brent"; low <- -10; up <- 10;
  }else{
  met <- "Nelder-Mead"; low <- -Inf; up <- Inf;
  }
  
  for(i in 1:nrow(bestVals)){
    opt <- optim(bestVals[i,], cavObj$runModel, Model=Model, method=met, lower=low, upper=up)
    RQoptim[i, 1] <- opt$value
    RQoptim[i, 2:(ncol(tryVals)+1)] <- opt$par
    
    for(j in 1:5){
      opt <- optim(RQoptim[i, 2:(ncol(tryVals)+1)], cavObj$runModel, Model=Model, method="BFGS")
      opt <- optim(opt$par, cavObj$runModel, Model=Model, 
                     method=met, lower=low, upper=up)
      if(abs(RQoptim[i, 1] - opt$value)>10000000000){ #Convergence test
          RQoptim[i, 1] <- opt$value
          RQoptim[i, (ncol(tryVals)+1)] <- opt$par      
      }else{
        RQoptim[i,ncol(RQoptim)] <- j
        break
      }
    }
  }
  
    bestPar <- RQoptim[order(RQoptim[,1]), ][1, 2:(ncol(tryVals)+1)]
	cavObj$runModel(bestPar, Model) #Object needs to be updated with the "best" parameters
  
	lOut <- list(	bestVals = RQoptim,
					bestPar = bestPar,
					VAR = cavObj$VAR,
					bestRQ = cavObj$RQ)
				
				
}