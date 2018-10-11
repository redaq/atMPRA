fisherExactEach= function(oligo3, x)
{	
		oligoLine = oligo3[oligo3$simN==x,]
		if(nrow(oligoLine)>0)
		{   
			temp = as.integer(oligoLine[6:7])
			names(temp) = colnames(oligoLine[6:7])
			tableS = rbind(as.integer(oligoLine[3:4]), temp)
			if(tableS[1,1]>500000 & sum(tableS)>1000000)
			{	warning("Counts are too large, Fisher's Exact test will take a long time. \n")
			}
			result=fisher.test(tableS)
			return(result$p.value)
		}else
		{	return(NA)
		}
}


