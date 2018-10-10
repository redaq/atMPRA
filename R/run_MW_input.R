run_MW_input <-
function(x, datt, ntag, nsim)
{	  tryCatch(
	  {
		temp <- datt[datt$simN==x, ]
		if( nrow(temp) < 70)
		{
			testInput = wilcox_test(formula=temp$input~as.factor(temp$allele), paired=FALSE, alternative="two.sided", distribution="exact")
		}else
		{	testInput = wilcox_test(formula=temp$input~as.factor(temp$allele), paired=FALSE, alternative="two.sided", distribution="asymptotic")
		}
		return(pvalue(testInput))
	  },error=function(e)
	  {	return(NA)
	  })
}
