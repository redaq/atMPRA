run_ttest <-
function(oligo, nsim, nrep)
{		ratios = apply(oligo[,-c(1:2)], 1, function(x) (x[(nrep+1):(2*nrep)]+1)/(mean(as.numeric(x[1:(nrep)]))+1))
		ratios2= data.frame(t(ratios))
		ratios2$allele= oligo$allele
		ratios2$simN = oligo$simN
		combos = unique(oligo$simN)
		result = sapply(combos, function(x) ttest_each(ratios2, x, nrep))
		return(data.frame(simN=combos, ttest_paired=as.numeric(result[1,]), ttest = as.numeric(result[2,])))

}
