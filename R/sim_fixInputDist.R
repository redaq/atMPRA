sim_fixInputDist <- function(mean_A, mean_B, std_A, std_B,  sigma2_DNA_a0, sigma2_DNA_a1, sigma2_RNA_a0, sigma2_RNA_a1, ntag=33, nsim=100, nrep=5, measurementERROR=TRUE, slope) 
### Slope can be 1, or a vector of length ntag*nsim*2 , where every ntag is one allele 
{
	  allele <-c(rep("Ref", ntag*nsim), rep("Mut", ntag*nsim))
	  simN <- c(rep(1:nsim, each=ntag), rep(1:nsim, each=ntag))
	  datt <- data.frame(allele, simN)
  
	  ## varirance for the mean
	  paramA = getParam(meanNeed=mean_A, sdNeed=std_A)
	  input1 = rnbinom(ntag*nsim, prob=paramA$p, size=paramA$n)
	  paramB = getParam(meanNeed=mean_B, sdNeed=std_B)
	  input2 = rnbinom(ntag*nsim, prob=paramB$p, size=paramB$n)

	  muinput <- c(input1, input2)
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

