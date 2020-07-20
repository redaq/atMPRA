run_ttest <-
function(oligo, nsim, nrep)
{		notTwoAllele = names(table(oligo$simN)[table(oligo$simN)<2])
        oligo2 = oligo[!oligo$simN %in% notTwoAllele,]
		if(nrow(oligo2)==0)
		{	return(NA)
		}
		ratios = log((oligo2[,(3+nrep):(2+2*(nrep))]+1)/(oligo2[,(3:(2+nrep))]+1))
		ratios2= data.frame((ratios))
		ratios2$allele= oligo2$allele
		ratios2$simN = oligo2$simN
		
		combos = unique(oligo2$simN)
		result = sapply(combos, function(x) ttest_each(ratios2, x, nrep))
		return(data.frame(simN=combos, ttest_paired=as.numeric(result[1,]), ttest = as.numeric(result[2,])))

}
