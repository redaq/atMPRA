run_matching_one <-
function(x, datt, width=5)
{ 
		  tryCatch(
		  {
		  temp <- datt[datt$simN==x, ]
		  alleles = unique(temp$allele)
		  temp$treat = ifelse(grepl(alleles[1], temp$allele), 1, 0)
		  m.out = matchit(treat~input, data=temp, method="subclass", sub.by="all", subclass=ceiling(nrow(temp)/width))
		  m.data=match.data(m.out)
		  temp2 =  m.data
		  infoTable <- table(temp2$subclass, temp2$allele)
		  ignoreBin <- rownames(infoTable)[which(infoTable[,1]==0|infoTable[,2]==0)]
		  temp22 <- temp2[!temp2$subclass %in% ignoreBin,]
		  infoTable2 <- (table(temp22$subclass, temp22$allele))
		  infoTable3 = as.matrix(infoTable2)
		  mintemp = sapply(rownames(infoTable3), function(x) min(temp22$input[temp22$subclass==x]))
		  maxtemp = sapply(rownames(infoTable3), function(x) max(temp22$input[temp22$subclass==x]))
		  tempData = data.frame(bin=rownames(infoTable3), A=infoTable3[,1], T=infoTable3[,2],min=mintemp, max=maxtemp)
		  if(length(unique(temp22$allele))==2)
		  {  if(nrow(temp22)<2)
				{
				return (list(pvalue=-1))
				}else
				{   temp22$subclass <- as.factor(temp22$subclass)
						temp22$allele2 <- as.factor(temp22$allele)

						if(nrow(temp22)>=70)
						{
						wt <- wilcox_test(formula=output~allele2|subclass, data=temp22, distribution="asymptotic")
						}else
						{
						wt <- wilcox_test(formula=output~allele2|subclass, data=temp22, distribution="exact")
						}

						### number of tags for each allele
						#return (list(pvalue = pvalue(wt)[1], nbins = length(unique(temp22$subclass)), dataLost = nrow(temp)-nrow(temp22) , propDataLost = (nrow(temp)-nrow(temp22))/nrow(temp), binInfo = tempData ))
						return( list(pvalue = pvalue(wt)[1]))
		  }
		  }else
		  {return  (list(pvalue=NA))
		  }
		  },error = function(e)
		  { return (list(pvalue=NA))
		  })

}
