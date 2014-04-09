Christoffersen1998 <- function(prob, y, var){
  n <- length(y)
  I <- (y<var)
  n1 <- sum(I) 
  n0 <- n-n1
  n00 <- 0; n01 <- 0; n10 <- 0; n11 <- 0;
  
  for(i in 2:length(y)){
    if(I[i-1]==0 & I[i]==0)n00 <- n00+1
    else if(I[i-1]==0 & I[i]==1)n01 <- n01+1
    else if(I[i-1]==1 & I[i]==0)n10 <- n10+1
    else n11 <- n11+1
  }
  
  
  pi01 <- ifelse((n00+n01)>0, n01/(n00+n01), 0) 
  pi11 <- ifelse((n10+n11)>0, n11/(n10+n11), 0) 
  pi2 <- (n01+n11)/(n00+n10+n01+n11)
  
  t.uncond <- -2*( n0*log(1-prob) + n1*log(prob) - ( n0*log(1-(n1/n)) +n1*log(n1/n) ) )
  
  t.indep <- -2*log(
      ( ( (1-pi2)^((n00+n10)/2) / (1-pi01)^n00 ) * ( (1-pi2)^((n00+n10)/2) / pi01^n01 ) ) * 
      ( ( pi2^((n01+n11)/2) / (1-pi11)^n10 ) * ( pi2^((n01+n11)/2) / pi11^n11 )   )
  ) #Factorized for avoiding to low values for OS
  
  if(is.na(t.indep)) #If to low values for OS
  t.indep <- -2*( (n00+n10)*log(1-pi2)+(n01+n11)*log(pi2) - ( n00*log(1-pi01) + n01*log(pi01) +
                                                                n10*log(1-pi11) + n11*log(pi11) ) )
  #Does not handle n00=0 or n11=0 and is therfore not the default formula for t.indep
                           
  t.cond <- t.uncond+t.indep
  
  list(t.uncond = t.uncond, t.indep=t.indep, t.cond=t.cond, 
       p.uncond = pchisq(t.uncond, 1, lower.tail=FALSE),
       p.indep = pchisq(t.indep, 1, lower.tail=FALSE),
       p.cond = pchisq(t.cond, 2, lower.tail=FALSE) )
  
}