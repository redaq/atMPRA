run_fisherExact <-
function(oligo, ntag, nsim, nrep)
{	  
		oligo$input_sum = apply(as.matrix(oligo[,(3:(2+nrep))]), 1, mean, na.rm=T)
		notTwoAllele = names(table(oligo$simN)[table(oligo$simN)<2])
        oligo2 = oligo[!oligo$simN %in% notTwoAllele,]
		if(nrow(oligo2)==0)
		{	return(NA)
		}
		oligo2$input_sum = apply(as.matrix(oligo2[,(3:(2+nrep))]), 1, sum, na.rm=T)
		oligo2$output_sum =  apply(as.matrix(oligo2[,((3+nrep):(2+2*nrep))]), 1, sum, na.rm=T)
		index = seq(1, nrow(oligo2), 2)
		oligo3= merge(oligo2[grepl("Ref", oligo2$allele), c("allele", "simN", "input_sum", "output_sum")], oligo2[grepl("Mut", oligo2$allele),c("allele", "simN", "input_sum", "output_sum")], by.x="simN", by.y="simN", sort=F)
		combos = (unique(oligo3$simN))
		results = sapply(combos, function(x) fisherExactEach(oligo3, x))
  		return(data.frame(simN=combos, fisherPvalue = results))
}
