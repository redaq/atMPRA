run_edgeR <-
function(oligo, ntag, nsim, nrep)
{	    
		notTwoAllele = names(table(oligo$simN)[table(oligo$simN)<2])
        oligo2 = oligo[!oligo$simN %in% notTwoAllele,]
		if(nrow(oligo2)<=2)
		{	return(NA)
		}
		RNAdata = oligo2[,c(1:2, (3+nrep):(2+2*nrep))]
		DNAdata = oligo2[,c(1:(2+nrep))]
		
		resRNAMatrix = merge(RNAdata[grepl("Ref", RNAdata$allele),], RNAdata[grepl("Mut", RNAdata$allele),], by="simN", sort=F)
		resDNAMatrix = merge(DNAdata[grepl("Ref", DNAdata$allele),], DNAdata[grepl("Mut", DNAdata$allele),], by="simN", sort=F)
		resRNAMatrix2=resRNAMatrix[,c(3:(2+nrep), (4+nrep):(3+2*nrep))]
		resDNAMatrix2=resDNAMatrix[,c(3:(2+nrep), (4+nrep):(3+2*nrep))]
		colnames(resRNAMatrix2) = colnames(resDNAMatrix2) = c(paste(paste0("output_rep",1:nrep), "A", sep="_"), paste(paste0("output_rep",1:nrep), "T", sep="_"))
		rownames(resRNAMatrix2) = rownames(resDNAMatrix2) = resRNAMatrix$simN
		
		groups = c(rep("A", nrep), rep("T", nrep))
		cds = DGEList(ceiling(resRNAMatrix2), group=groups)
		#cds = DGEList((ratioMatrix2), group=groups)
		design = model.matrix(~groups)
		cds$offset = as.matrix(log(resDNAMatrix2+1))
		cds = estimateDisp(cds, design)
		fit = glmFit(cds, design)
		result = glmTreat(fit)
		resultTable = result$table
		resultTable$simN = resRNAMatrix$simN   
		return(resultTable)
				
		
}
