calEnrichment <-
function(resultAll, threshold=0.05)
{		message("Compute enrichment in allele-imbalanced SNPs among significant results... \n")
		enrichment=c()
		for(i in 3:ncol(resultAll))
		{	dataFrame = resultAll[resultAll$resInput >=0 & resultAll[,i] >= 0 & !is.na(resultAll$resInput) & !is.na(resultAll[,i]),]
			ddtemp = dataFrame[dataFrame[,i]<threshold,]
			MWP = data.frame(method = colnames(resultAll)[i], q=sum(ddtemp$resInput<threshold), m=nrow(ddtemp), enrichP=phyper(q=sum(ddtemp$resInput<threshold)-1, m=nrow(ddtemp), n=nrow(dataFrame)-nrow(ddtemp), k=sum(dataFrame$resInput<threshold), lower.tail=F))
			enrichment =rbind(enrichment, MWP)
		}
		return(enrichment)
}
