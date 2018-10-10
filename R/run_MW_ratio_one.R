run_MW_ratio_one <-
function(x, datt)
{ 
	  tryCatch(
	  {
		temp <- datt[datt$simN==x & datt$keep, ]	
		if(nrow(temp) < 70)
		{
			   wt = wilcox_test(formula=temp$ratio~as.factor(temp$allele), paired=FALSE, alternative="two.sided", distribution="exact")
		}else
		{	   wt = wilcox_test(formula=temp$ratio~as.factor(temp$allele), paired=FALSE, alternative="two.sided", distribution="asymptotic")
		}
		return(pvalue(wt))
	  },error = function(e)
	  { return(NA)
	  })
}
