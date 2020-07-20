getPower <- function(nsim = 100, ntag = 10, nrepIn = 2, nrepOut=2, slope = 1, scenario=c("fixInputDist", "fixTotalDepth", "fixMeanDepth"), method=c("MW", "Matching", "Adaptive", "Fisher", "QuASAR", "T-test", "mpralm", "edgeR", "DESeq2"), fixInput  = c(20, 100), fixMeanD = 70 , fixTotalD= 20000000, std_A=mean_A, std_B=mean_B, inputDist=NULL, inputDispFunc=NULL, outputDispFunc=NULL, inputDispParam=NULL, outputDispParam=NULL, cutoff=-1, cutoffo=-1, p.adjust.method="fdr", significance=0.05)
{		if(all(is.na(slope))|(length(slope)!=2*ntag*nsim&length(slope)!=1))
	  	{	stop("The 'slope' parameter needs to be of length 1 or 2*ntag*nsim long to specify the transfection efficiency for all tags.\n")
	  	}

		if(is.null(inputDist))
		{	warning("The input distribution 'inputDist' is not specified, it should be of length nsim*ntag*2, using default distribution.\n")
			inputDist=GSE70531_params[[3]](runif(nsim*ntag*2))
	
		}
		if(scenario=="fixInputDist")
		{	    message("\n\nSimulating data using fixed input distributions...\n")
				datt = sim_fixInputMean(mean_A=fixInput[1], mean_B=fixInput[2],  ntag=ntag, nsim=nsim, nrepIn=nrepIn, nrepOut=nrepOut, slope=slope, inputDist=inputDist, std_A=std_A, std_B=std_B, inputDispFunc=inputDispFunc, outputDispFunc=outputDispFunc, inputDispParam=inputDispParam, outputDispParam=outputDispParam)
		}else if(scenario=="fixTotalDepth")
		{		message("\n\nSimulating data using given input distributions...\n")
				
				datt=sim_fixDepth(inputProp=inputDist, ntag=ntag, nsim=nsim, nrepIn=nrepIn, nrepOut=nrepOut, slope=slope, inputDispFunc=inputDispFunc, outputDispFunc=outputDispFunc, sampleDepth=fixTotalD, inputDispParam=inputDispParam, outputDispParam=outputDispParam) 
		}else if(scenario=="fixMeanDepth")
		{		message("\n\nSimulating data using given input distributions...\n")
				
				datt=sim_fixDepth(inputProp=inputDist, ntag=ntag, nsim=nsim, nrepIn=nrepIn, nrepOut=nrepOut, slope=slope, inputDispFunc=inputDispFunc, outputDispFunc=outputDispFunc, meanDepth=fixMeanD, inputDispParam=inputDispParam, outputDispParam=outputDispParam) 
		}else
		{		stop("\n\nUnknown input for parameter 'scenario'.\n")
		}
		
		results = analyzeMPRA(datt, nrepIn, rnaCol=2+nrepIn+1, nrepOut, nsim, ntag, method, cutoff, cutoff)
		return(list(Power=getPowerOne(results, p.adjust.method, significance), simData=datt, results=results))

}


