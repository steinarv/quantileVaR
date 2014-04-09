DQ <- function(Y, VAR, prob, intercept=T, nVAR=0, nHIT=c(1,2,3,4), nY=F){
  N <- length(Y)
  n <- N-max(nVAR, nHIT, nY)
  HIT <- (Y<VAR)-prob
  
  X <- matrix(ifelse(intercept,1,0),ncol=1,nrow=n)
  

  if(any(is.numeric(nVAR))){
  xVAR=matrix(NA, n, length(nVAR))
  for(i in 1:length(nVAR)){ 
    xVAR[,i]=VAR[(N-nVAR[i]-n+1):(N-nVAR[i])]
  }
  X <- cbind(X, xVAR)}

  if(any(is.numeric(nHIT))){
  xHIT=matrix(NA, n, length(nHIT))
  for(i in 1:length(nHIT)){
    xHIT[,i]=HIT[(N-nHIT[i]-n+1):(N-nHIT[i])]
  }
  X <- cbind(X, xHIT)}
  
  if(any(is.numeric(nY))){
  xY=matrix(NA, n, length(nY))
  for(i in 1:length(nY)){
    xY[,i]=Y[(N-nY[i]-n+1):(N-nY[i])]
  }
  X <- cbind(X, xY)}
  
  DQ <- fastDQ(HIT[(N-n+1):N], X, prob)

  DQ[["p_val"]] <- pchisq(DQ$DQ,length(DQ$coefficients),lower.tail=F)
  
  DQ
}


