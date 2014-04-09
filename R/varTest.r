varTest <- function(prob, y, VaR){
	
	nP <- length(prob)
	if(!is.null(dim(VaR))){ if(ncol(VaR)!=nP) stop("Number of VaR columns do not mach number of probabilities.") }
					else { if(nP!=1) stop("Number of VaR columns do not mach number of probabilities.") }	
	
	mOut <- matrix(NA, ncol=5, nrow=nP)
	mOut[ , 1] <- prob
	colnames(mOut) <- c("Prob", "Hit(%)", "P-value(DQ)", "P-value(uncond)", "P-value(cond)")
  
	for(i in 1:length(prob)){
		VaR.i <- VaR[, i]
		if(length(y)!=length(VaR.i)) stop("y and VaR do not have the same length! ")
    
		Coverage <- Christoffersen1998( prob[i], y, VaR.i)
    
		mOut[i, 2] <- mean(y<VaR.i)
		mOut[i, 3] <- DQ(Y=y, VAR=VaR.i, prob=prob[i])$p_val
		mOut[i, 4] <- Coverage$p.uncond
		mOut[i, 5] <- Coverage$p.cond
	}
  
	mOut
}