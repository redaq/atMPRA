getPowerOne <-
function(resultAll)	
{		#### collect results
		message("\nComputing the power/type I error rate for the dataset... \n")
		power = apply(as.matrix(resultAll[,-c(1:2)]), 2, function(x) sum(x<0.05 & x>=0, na.rm=T)/sum(!is.na(x)&x>=0, na.rm=T))
		print(power)
		return((power))

}
