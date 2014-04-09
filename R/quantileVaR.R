loadModule("cppMod", TRUE)

fastDQ <- function(y, x, prob){
	.Call( "fastDQ", y, x, prob, PACKAGE = "quantileVaR" )
}