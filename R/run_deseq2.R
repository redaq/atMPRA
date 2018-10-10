run_deseq2 <-
function(oligo, ntag, nsim, nrep)
{		oligo$input_sum = apply(as.matrix(oligo[,(3:(2+nrep))]), 1, mean, na.rm=T)
		notTwoAllele = names(table(oligo$simN)[table(oligo$simN)<2])
        oligo2 = oligo[!oligo$simN %in% notTwoAllele,]
	    if(nrow(oligo2)<=2)
		{	return(NA)
		}
		meanInput = apply(oligo2[,c(3:(2+nrep))], 1, mean, na.rm=T)
		RNAdata = oligo2[,c(1:2, (3+nrep):(2+2*nrep))]
		DNAdata = oligo2[,c(1:(2+nrep))]
		for(col in 3:ncol(DNAdata))
		{	DNAdata[,col] = meanInput
		}
		resRNAMatrix = merge(RNAdata[grepl("Ref", RNAdata$allele),], RNAdata[grepl("Mut", RNAdata$allele),], by="simN", sort=F)
		resDNAMatrix = merge(DNAdata[grepl("Ref", DNAdata$allele),], DNAdata[grepl("Mut", DNAdata$allele),], by="simN", sort=F)
		resRNAMatrix2=resRNAMatrix[,c(3:(2+nrep), (4+nrep):(3+2*nrep))]
		resDNAMatrix2=resDNAMatrix[,c(3:(2+nrep), (4+nrep):(3+2*nrep))]
		colnames(resRNAMatrix2) = colnames(resDNAMatrix2) = c(paste(paste0("output_rep",1:nrep), "A", sep="_"), paste(paste0("output_rep",1:nrep), "T", sep="_"))
		rownames(resRNAMatrix2) = rownames(resDNAMatrix2) = resRNAMatrix$simN
		
		ratioMatrix2 = (as.matrix(resRNAMatrix2)+1)/(as.matrix(resDNAMatrix2)+1)
		coldata = data.frame(group=c(rep("A", nrep), rep("T", nrep)))
		rownames(coldata) = colnames(ratioMatrix2)
		dds <- DESeqDataSetFromMatrix(countData = ceiling(resRNAMatrix2),
                              colData = coldata,
                              design = ~ group)
        normalFactor = as.matrix(ceiling(resDNAMatrix2))
        normalizationFactors(dds) =normalFactor/exp(rowMeans(log(normalFactor)))
        
        result = tryCatch(
        {       DESeq(dds)
        },error=function(e)
        {       warning("Error in DESeq2 fitting, using edgeR tagwise estimates here...\n")
		### DESEQ gene estimates are really weird, so hacking in edgeR tagwise estimates
		 		groups = c(rep("Ref", nrep), rep("Mut", nrep))
                cds = DGEList(ceiling(resRNAMatrix2), group=groups)
                design = model.matrix(~groups)
                cds$offset = as.matrix(log(resDNAMatrix2))
                cds = estimateDisp(cds, design)
        	#	dds <- estimateDispersionsGeneEst(dds)
                dispersions(dds) <- cds$tagwise.dispersion
                nbinomWaldTest(dds)
        })
        restable=results(result)
		restable$simN = resRNAMatrix$simN
		return(restable)
}
