sim_fixTotalD <- function(datMean, totalDepth=total, sigma2_DNA_a0, sigma2_DNA_a1, sigma2_RNA_a0, sigma2_RNA_a1, ntag=ntagsV, nsim=nelements, nrep=5, measurementERROR=TRUE, slope) 
{  
	  if(length(slope)!= (2*ntag*nsim)  & length(slope)!=1)
	  {	stop("Length of slope parameter needs to be 2*ntag*nsim, with the first half correpsonding to the slopes for the reference allele, and second half corresponding to the mutant allele.\n")
	  }
	  allele <-c(rep("Ref", ntag*nsim), rep("Mut", ntag*nsim))
	  simN <- c(rep(1:nsim, each=ntag), rep(1:nsim, each=ntag))
	  datt <- data.frame(allele, simN)
  
	  ## varirance for the mean
	  meanDNA = sample(datMean, size=nsim*ntag*2, replace=T)
	  muinput = totalDepth*(meanDNA/sum(meanDNA))
	  inputRNAv = c()
	  for(rep in 1:nrep)
	  { if (measurementERROR)
		{   sigma2_DNA = sigma2_DNA_a0 + sigma2_DNA_a1/muinput
			inputNew = rnbinom(length(muinput), mu=muinput , size=1/sigma2_DNA)
		}else
		{	inputNew = muinput     
		} 
		datt[,paste0("input_rep", rep)] = inputNew
	  }
  
	  ### generate output
	  for(rep in 1:nrep)
	  { mean_RNA = slope*muinput
		sigma2_RNA = sigma2_RNA_a0+sigma2_RNA_a1/mean_RNA
		output = rnbinom(length(muinput), mu = mean_RNA , size=1/sigma2_RNA)
		datt[,paste0("output_rep", rep)] = output
	  }
	  return(datt)
}

