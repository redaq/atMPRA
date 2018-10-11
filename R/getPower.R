getPower <- function(nsim = 100, ntag = 10, nrep = 1, slope = NA, scenario=c("fixInputDist", "fixTotalDepth", "fixMeanDepth"), method=c("MW", "Matching", "Adaptive", "Fisher", "QuASAR", "T-test", "mpralm", "edgeR", "DESeq2"), sigma2_DNA_a0, sigma2_DNA_a1, sigma2_RNA_a0, sigma2_RNA_a1, fixInput  = NA, fixMeanD = 70 , fixTotalD= 20000000, datMean=NA, cutoff=-1, cutoffo=-1)
{		if(all(is.na(slope))|length(slope)!=2*ntag*nsim)
	  	{	stop("The 'slope' parameter needs to be 2*ntag*nsim long to specify the transfection efficiency for all tags.\n")
	  	}

		if(scenario=="fixInputDist")
		{	    message("\n\nSimulating data using fixed input distributions...\n")
				datt = sim_fixInputDist(mean_A=fixInput[1], mean_B=fixInput[2], std_A=fixInput[3], std_B=fixInput[4], sigma2_DNA_a0=sigma2_DNA_a0, sigma2_DNA_a1=sigma2_DNA_a1, sigma2_RNA_a0=sigma2_RNA_a0, sigma2_RNA_a1=sigma2_RNA_a1,  ntag=ntag, nsim=nsim, nrep=nrep, slope=slope)
		}else if(scenario=="fixTotalDepth")
		{		message("\n\nSimulating data using given input distributions...\n")
				if(all(is.na(datMean)))
				{	stop("The input distribution 'datMean' needs to be specified.\n")
				}		
				datt = sim_fixTotalD(datMean=datMean, totalDepth=fixTotalD, sigma2_DNA_a0=sigma2_DNA_a0, sigma2_DNA_a1=sigma2_DNA_a1, sigma2_RNA_a0=sigma2_RNA_a0, sigma2_RNA_a1=sigma2_RNA_a1, ntag=ntag, nsim=nsim, nrep=nrep, slope=slope)
		}else if(scenario=="fixMeanDepth")
		{		message("\n\nSimulating data using given input distributions...\n")
				if(all(is.na(datMean)))
				{	stop("The input distribution 'datMean' needs to be specified.\n")	
				}
				datt = sim_fixMeanD(datMean=datMean, meanDepth=fixMeanD, sigma2_DNA_a0=sigma2_DNA_a0, sigma2_DNA_a1=sigma2_DNA_a1, sigma2_RNA_a0=sigma2_RNA_a0, sigma2_RNA_a1=sigma2_RNA_a1, ntag=ntag, nsim=nsim, nrep=nrep, slope=slope)
		}else
		{		stop("\n\nUnknown input for parameter 'scenario'.\n")
		}
		
		results = analyzeMPRA(datt, nrep, rnaCol=2+nrep+1, nrep, nsim, ntag, method, cutoff, cutoff)
		return(list(Power=getPowerOne(results), simData=datt, results=results))

}

