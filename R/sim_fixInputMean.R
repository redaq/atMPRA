sim_fixInputMean <- function(mean_A, mean_B,  ntag, nsim, nrepIn, nrepOut, slope, inputDist=NULL, std_A=mean_A, std_B=mean_B, inputDispFunc=NULL, outputDispFunc=NULL, inputDispParam=NULL, outputDispParam=NULL)
{
          if(length(slope)!= (2*ntag*nsim)  & length(slope)!=1)
          {     stop("Length of slope parameter needs to be 2*ntag*nsim, with the first half correpsonding to the slopes for the reference allele, and second half corresponding to the mutant allele.\n")
		  }
		  ### make input distribution
		  ## varirance for the mean
		  
		  if(!is.null(inputDist))
		  {
		  	if(length(inputDist)<ntag)
		  	{	stop("'inputDist' has less than ntag values, which will be used to sample the tag counts without replacement.\n")
		  	}
		  	input1 = rep(NA, nsim*ntag)
		  	for(i in 1:nsim)
		  	{ tempProp = sample(inputDist, ntag, replace=FALSE)
		  	  totaltemp = mean_A/(mean(tempProp/sum(tempProp)))
		  	  input1[((i-1)*ntag+1):(i*ntag)] = totaltemp * tempProp/sum(tempProp)
		  	}
		  	input2 = rep(NA, nsim*ntag)
		  	for(i in 1:nsim)
		  	{ tempProp = sample(inputDist, ntag, replace=FALSE)
		  	  totaltemp = mean_B/(mean(tempProp/sum(tempProp)))
		  	  input2[((i-1)*ntag+1):(i*ntag)] = totaltemp * tempProp/sum(tempProp)
		  	}

		  }else
		  {	 message("'inputDist' is not provided, using std_A and std_B to specify the input distribution for each SNP/simulation.\n")
		  	 paramA = getParam(meanNeed=mean_A, sdNeed=std_A)
          	 input1 = rnbinom(ntag*nsim, prob=paramA$p, size=paramA$n)
          	 paramB = getParam(meanNeed=mean_B, sdNeed=std_B)
          	 input2 = rnbinom(ntag*nsim, prob=paramB$p, size=paramB$n)

		  }
		 
		  muinput <- c(input1, input2)
		  inputProp = muinput/sum(muinput)
		  sampleDepth = sum(muinput)
		  datt = sim_fixDepth(inputProp, ntag=ntag, nsim=nsim, nrepIn=nrepIn,  nrepOut=nrepOut, slope=slope, inputDispFunc=inputDispFunc, outputDispFunc=outputDispFunc, inputDispParam=inputDispParam, outputDispParam=outputDispParam, sampleDepth=sampleDepth) 
		  return(datt)
}
		  
		  

