estimateMPRA <- function(datt, nrepIn, rnaCol, nrepOut, nsim, ntag, plotFigure=FALSE, plotName="simData")
{	  if(any(colnames(datt)[1:2] != c("allele", "simN")))
	  {	stop("Please format the count data (datt) such that the first two columns are named as 'allele' and 'simN'.\n")
	  }

	### input
	if(nrepIn>1)
	{
	coldata = data.frame(group=rep("Input", nrepIn))
	rownames(coldata) = colnames(datt[,c(3:(2+nrepIn))])
	dds <- DESeqDataSetFromMatrix(countData = datt[,c(3:(2+nrepIn))],
						colData = coldata,
						design = ~ 1)
	dds <- estimateSizeFactors(dds)
	norFactor_input =sizeFactors(dds)
	input_norm <- counts(dds, normalized=TRUE)			
	temp1 = estimateDispersions(dds, fitType="parametric")
	dispFunc_input = dispersionFunction(temp1)
	}else
	{   warning("There is only one input replicate, cannot compute dispersion function.\n")
		dispFunc_input = NA
		norFactor_input = 1
	}
	  
 	### output
 	if(nrepOut>1)
 	{
	coldata1 = data.frame(group=c(rep("Output", nrepOut)))
	rownames(coldata1) = colnames(datt[,rnaCol:(rnaCol+nrepOut-1)])
	dds <- DESeqDataSetFromMatrix(countData = datt[,rnaCol:(rnaCol+nrepOut-1)],
						colData = coldata1,
						design = ~ 1)
	dds <- estimateSizeFactors(dds)
	norFactor_output =sizeFactors(dds)                   
	output_norm <- counts(dds, normalized=TRUE)
	temp2 = estimateDispersions(dds, fitType="parametric")
	dispFunc_output = dispersionFunction(temp2)
	}else
	{	warning("There is only one output replicate, cannot compute dispersion function.\n")
		dispFunc_output = NA
		norFactor_output = 1

	}
	## estimate inputProp using goodTuringProportions
	propEst = apply(as.matrix(datt[,3:(2+nrepIn)]), 2, goodTuringProportions)
	if(nrepIn>1)
	{	propMean = apply(as.matrix(propEst), 1, mean, na.rm=T)
	}else
	{	propMean = propEst
	}
	cdfest = ecdf(sort(propMean))
	ptiles = cdfest(sort(propMean))
	inputProp = approxfun(x=ptiles, y=sort(propMean), rule=2)
	
	### estimate the slope
	### estimate of b from normalization using total count
	dattNorm1 = t(t(datt[,-c(1:2)]+1)/colSums(datt[,-c(1:2)], na.rm=T)) 
	dattNorm = data.frame(datt[,1:2], dattNorm1)
	oligo =aggregate(dattNorm[,3:(rnaCol+nrepOut-1)], by=list(dattNorm$allele, dattNorm$simN), FUN=sum, na.rm=T)
	meanDNA  = apply(as.matrix(oligo[,3:(2+nrepIn)]), 1, mean, na.rm=T)
	ratio = (oligo[, rnaCol:(rnaCol+nrepOut-1)])/(meanDNA)
	meanRatio = apply(as.matrix(ratio), 1, mean, na.rm=T)
	cdfest2 = ecdf(sort(meanRatio))
	ptiles2 = cdfest2(sort(meanRatio))
	transEff = approxfun(x=ptiles2, y=sort(meanRatio), rule=2)

	if(plotFigure)
	{
		png(paste0("plot_dispersion_", plotName, ".png"), width=2000, height=2000, res=300)
		message("Plotting for the data...\n")
		par(mfrow=c(2,2))
		plotDispEsts(temp1)
		plotDispEsts(temp2)
		plot(density(propMean), main="input proportions")
		plot(density(meanRatio), main="transfection efficiency")
		dev.off()
	}
	return(list(dispFunc_input=dispFunc_input, dispFunc_output=dispFunc_output, inputProp=inputProp, transEff=transEff, sizeFactor_input = norFactor_input, sizeFactor_output = norFactor_output))
	
}

