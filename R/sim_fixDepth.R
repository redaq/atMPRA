sim_fixDepth <- function(inputProp, ntag, nsim, nrepIn,  nrepOut, slope, inputDispFunc=NULL, outputDispFunc=NULL, sampleDepth=NULL, inputDispParam=NULL, outputDispParam=NULL,  meanDepth=NULL) 
{  
	  if(is.null(slope)|(length(slope)!=2*ntag*nsim&length(slope)!=1))
	  {	stop("The 'slope' parameter needs to be of length 1 or 2*ntag*nsim to specify the transfection efficiency for all tags.\n")
	  }
	  if(is.null(inputProp)|(length(inputProp)!=2*ntag*nsim))
	  {	stop("The 'inputProp' parameter needs to be of length 2*ntag*nsim.\n")
	  }
	  if(!is.null(inputDispFunc) & !is.function(inputDispFunc))
	  {	stop("'inputDispFunc' needs to be a function obtained using DESeq2 for dispersion function estimation.\n")
	  }
	  if(!is.null(outputDispFunc) &!is.function(outputDispFunc))
	  {	stop("'outputDispFunc' needs to be a function obtained using DESeq2 for dispersion function estimation.\n")
	  }
	  if(is.null(inputDispFunc) & !is.null(inputDispParam) & length(inputDispParam)!=3)
	  {	stop("'inputDispParam' should be a vector of three parameters.\n")
	  }
	  if(is.null(outputDispFunc) & !is.null(outputDispParam) & length(outputDispParam)!=3)
	  {	stop("'inputDispParam' should be a vector of three parameters.\n")
	  }

	  if(is.null(inputDispFunc) & is.null(inputDispParam))
	  {	message("Both 'inputDispFunc' and 'inputDispParam' were not specified.\nUsing default values...\n")
	  	inputDispParam = c(0, 4.37, 0.25)
	  }
	  if(is.null(outputDispFunc) & is.null(outputDispParam))
	  {	message("Both 'outputDispFunc' and 'outputDispParam' were not specified.\nUsing default values...\n")
	  	outputDispParam = c(0.54, 12, 0.25)
	  }

	  if(!is.null(meanDepth) & is.null(sampleDepth))
	  {	message("Sample depth is not specified, using meanDepth.\n")
	  	if(length(meanDepth)==1)
	  	{	
	  		sampleDepth = rep(meanDepth*ntag*nsim*2, nrepIn+nrepOut)
	  	}else
	  	{	sampleDepth = meanDepth*ntag*nsim*2
	  	}
	  }else if(!is.null(sampleDepth))
	  {	if(length(sampleDepth)==1)
	  	{	sampleDepth = rep(sampleDepth, nrepIn+nrepOut)
	  	}
	  }else
	  {		message("The sample depth is not specified by the user.\n Using default sample depth...\n")
	  		sampleDepth = rep(15000000, nrepIn+nrepOut)
	  }
	  allele <-c(rep("Ref", ntag*nsim), rep("Mut", ntag*nsim))
	  simN <- c(rep(1:nsim, each=ntag), rep(1:nsim, each=ntag))
	  datt <- data.frame(allele, simN)
	  
	  datt1 = datt
	  for(rep in 1:nrepIn)
	  {  muinput = sampleDepth[rep] * inputProp/sum(inputProp)
	  	 datt1[,paste0("input_rep", rep)] = muinput
	  } 
	  ## normalize 
	  coldata = data.frame(group=rep("Input", nrepIn))
	  rownames(coldata) = colnames(datt1[,c(3:(2+nrepIn))])
	  suppressWarnings(dds <- DESeqDataSetFromMatrix(countData = ceiling(datt1[,c(3:(2+nrepIn))]),
						colData = coldata,
						design = ~ 1))
	  dds <- estimateSizeFactors(dds)
	  input_norm1 <- counts(dds, normalized=TRUE)			
 	  mean_input = apply(input_norm1, 1, mean, na.rm=T)
 	  if(is.function(inputDispFunc))
 	  {
	  if(attr(inputDispFunc, "fitType" )=="parametric")
		{		sigma2_DNA_a0 = attr(inputDispFunc, "coefficients" )["asymptDisp"]
				sigma2_DNA_a1 = attr(inputDispFunc, "coefficients" )["extraPois"]
				sigma2_DNA_0 = sigma2_DNA_a0 + sigma2_DNA_a1/mean_input
				sigma2_DNA = exp(rnorm(length(muinput), mean=log(sigma2_DNA_0), sd=sqrt(attr(inputDispFunc, "dispPriorVar" ))))
		}else
		{	    mean_input_0 = mean_input
				mean_input_0[mean_input_0==0] = 0.5
				sigma2_DNA = inputDispFunc(mean_input_0)
			  
		}
	 }else ## using inputDispParam
	 { 	sigma2_DNA_a0 = inputDispParam[1]
	 	sigma2_DNA_a1 = inputDispParam[2]
	 	sigma2_DNA_0 = sigma2_DNA_a0 + sigma2_DNA_a1/mean_input
		sigma2_DNA = exp(rnorm(length(muinput), mean=log(sigma2_DNA_0), sd=sqrt(inputDispParam[3])))
	 }

	  for(rep in 1:nrepIn)
	  {  muinput = sampleDepth[rep]*(inputProp/sum(inputProp)) 
	  	 inputNew =rnbinom(length(muinput), mu=muinput , size=1/sigma2_DNA)
		 datt[,paste0("input_rep", rep)] = inputNew
	  }
  
	  ### generate output
	  datt2 = datt
	  for(rep in 1:nrepOut)
	  {		mean_RNA = slope*(inputProp/sum(inputProp))* sampleDepth[nrepIn+rep]
	  		datt2[,paste0("output_rep", rep)] = mean_RNA
	  }
	  ## normalize 
	  coldata = data.frame(group=rep("Output", nrepOut))
	  rownames(coldata) = colnames(datt2[,c((2+nrepIn+1):(2+nrepIn+nrepOut))])
	  suppressWarnings(ddds <- DESeqDataSetFromMatrix(countData = ceiling(datt2[,c((2+nrepIn+1):(2+nrepIn+nrepOut))]),
						colData = coldata,
						design = ~ 1))
	  dds <- estimateSizeFactors(dds)
	  output_norm1 <- counts(dds, normalized=TRUE)			

	  mean_output = apply(output_norm1, 1, mean, na.rm=T)
	  if(is.function(outputDispFunc))
	  {
		  if(attr(outputDispFunc, "fitType" )=="parametric")
		  {		sigma2_RNA_a0 = attr(outputDispFunc, "coefficients" )["asymptDisp"]
				sigma2_RNA_a1 = attr(outputDispFunc, "coefficients" )["extraPois"]
				sigma2_RNA_0 = sigma2_RNA_a0 + sigma2_RNA_a1/mean_output
				sigma2_RNA = exp(rnorm(length(muinput), mean=log(sigma2_RNA_0), sd=sqrt(attr(outputDispFunc, "dispPriorVar" ))))
		  }else
		  {	    mean_output_0 = mean_output
				mean_output_0[mean_output_0==0] = 0.5
		  		sigma2_RNA = outputDispFunc(mean_output_0)
		  }
	  }else
	  {			sigma2_RNA_a0 = outputDispParam[1]
				sigma2_RNA_a1 = outputDispParam[2]
				sigma2_RNA_0 = sigma2_RNA_a0 + sigma2_RNA_a1/mean_output
				sigma2_RNA = exp(rnorm(length(muinput), mean=log(sigma2_RNA_0), sd=sqrt(outputDispParam[3])))

	  }

	  for(rep in 1:nrepOut)
	  {		mean_RNA = slope*(inputProp/sum(inputProp))*sampleDepth[nrepIn+rep]
	  		output = rnbinom(length(muinput), mu = mean_RNA , size=1/sigma2_RNA)
			datt[,paste0("output_rep", rep)] = output
	  }
	  datt[is.na(datt)]=0
	  	  
	  return(datt)
}

