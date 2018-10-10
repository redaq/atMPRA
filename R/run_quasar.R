run_quasar <-
function(datt2, ntag, nsim, nrep)
{		
		datt2 = datt2[order(datt2$allele),]
		sumall = aggregate(datt2[,(3:(2+2*nrep))], by=list(datt2$allele, datt2$simN), FUN= sum, na.rm=T)
		### remove oligos that do not have two alleles 
		notTwoAllele = names(table(sumall$Group.2)[table(sumall$Group.2)<2])
		sumall2 = sumall[!sumall$Group.2 %in% notTwoAllele,]
		refIndex = which(grepl("Ref", sumall2$Group.1))
		
		DNAprop = as.matrix(sumall2[refIndex,3:(2+nrep)]/(sumall2[refIndex,3:(2+nrep)]+sumall2[-refIndex,3:(2+nrep)]))
		output_ref = as.matrix(sumall2[refIndex, (2+nrep+1):(2+2*nrep)])
		output_alt = as.matrix(sumall2[-refIndex, (2+nrep+1):(2+2*nrep)])
		
		if(sum(sumall2[refIndex,]$Group.2 == sumall2[-refIndex,]$Group.2)!= length(refIndex))
		{ stop("not matching simN")
		}
		
        results <- list()
        
        tryCatch(
        {
        inputSum = apply(as.matrix(sumall2[, (3:(2+nrep))]), 1, sum, na.rm=T)
        rl = inputSum[refIndex]/(inputSum[-refIndex]+inputSum[refIndex])
        
        ###
        if(nrep > 1)
        {
        for(i in 1:nrep)
        {	results[[i]] = fitQuasarMpra(ref=as.numeric(output_ref[,i]), alt=as.numeric(output_alt[,i]), as.numeric(DNAprop[,i]))
        	results[[i]]$weight = 1/((results[[i]]$betas_se)^2) 
        	results[[i]]$betaWeighted = results[[i]]$betas.beta.binom*results[[i]]$weight
        	
        }
        ### fixed effect meta
        weightsMatrix = sapply(1:nrep, function(x) results[[x]]$weight)
        weightsSum = apply(weightsMatrix, 1, sum) 
        SE = sqrt(1/weightsSum)
        betaWeightedMat = sapply(1:nrep, function(x) results[[x]]$betaWeighted)
        betaSum = apply(betaWeightedMat, 1, sum)/weightsSum
        Zs = (betaSum - log(rl/(1-rl)))/SE
        pvalue = 2*pnorm(-abs(Zs))
        	resultAll = data.frame(simN = sumall2[refIndex,]$Group.2, Z=Zs, DNAprop = rl, pvalue=pvalue)
        	return(resultAll)
        }else
        {	results = fitQuasarMpra(ref=as.numeric(output_ref[,1]), alt=as.numeric(output_alt[,1]), as.numeric(DNAprop[,1]))
       		resultAll = data.frame(simN = sumall2[refIndex,]$Group.2, Z=results$betas_z, DNAprop = rl, pvalue=results$pval3)
        	return(resultAll)
        }	

        },error=function(e)
        {	warning("Error in running QuASAR-MPRA.n")
        	resultAll = data.frame(simN =sumall2[refIndex,]$Group.2, pvalue=NA)
        	return(resultAll)
        })
       
        
}
