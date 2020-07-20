getPowerOne <- function(resultAll, method="fdr", significance=0.05)	
{		#### collect results
		message("\nComputing the power/type I error rate for the dataset using ", method, " at ", significance, ": \n")
		resultsAlladj = apply(as.matrix(resultAll[,-c(1:2)]), 2, function(x) p.adjust(x, method=method))
		power = apply(as.matrix(resultsAlladj), 2, function(x) sum(x<significance & x>=0, na.rm=T)/sum(!is.na(x)&x>=0, na.rm=T))
		print(power)
		return((power))

}


