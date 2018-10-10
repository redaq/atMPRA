analyzeMPRA <- function(datt, nrepIn, rnaCol, nrepOut, nsim, ntag, method=c("MW", "Matching", "Adaptive", "Fisher", "QuASAR", "T-test", "mpralm", "edgeR", "DESeq2"), cutoff=-1, cutoffo=-1)
{		message("Start analyzing MPRA data: \n")
		
		### normalization
		message("==== Normalizing samples using the maximum library depth...\n")
		nrep = nrepOut
		normFactor = max(colSums(datt[,3:(2+2*nrep)], na.rm=T))
		counts = data.frame(t( t(datt[,-c(1:2)])/colSums(datt[,-c(1:2)], na.rm=T))*normFactor)
		meanInput = apply(as.matrix(counts[, 1:nrepIn]), 1, mean, na.rm=T)
		newDNA =t(do.call(rbind, replicate(nrepOut, meanInput, simplify=FALSE)) )
		colnames(newDNA) = paste0("input_rep", 1:nrepOut)
		
		datt4 = cbind(datt[,1:2], newDNA, counts[, (rnaCol-2):(rnaCol-2-1+nrepOut)])
		colnames(datt4) = colnames(datt)

		#### filtering
		message("==== Filtering tags with DNA counts less than or equal to ", cutoff, ", and with RNA counts less than or equal to ", cutoffo, ".\n")
		if(cutoff==-1 & cutoffo==-1)
		{	keep=rep(TRUE, nrow(datt))  
			datt2 = datt4
			naindex = which(is.na(datt2), arr.ind=T)
			datt2[naindex] = 0
		}else
		{
			keep = apply((datt[,3:(2+2*nrep)]), 1, function(x) all(!is.na(as.numeric(x[1:nrep]))&as.numeric(x[1:nrep])>cutoff, na.rm=T)&all(!is.na(as.numeric(x[(nrep+1):(2*nrep)]))&as.numeric(x[(nrep+1):(2*nrep)])>cutoffo, na.rm=T))
			datt2 = datt4[keep,]
		}
		### This is the data frame after filtering
		datt3 = datt4
		datt3[!keep, 3:(2+2*nrep)] = NA
		
		### Data frame for MW test		
		meanInput = apply(as.matrix(datt2[,3:(2+nrep)]), 1, mean, na.rm=T)
		meanOutput = apply(as.matrix(datt2[,(3+nrep):(2+nrep*2)]), 1, mean, na.rm=T)
		ratio = (meanOutput+1)/(meanInput+1)
		dattNew = cbind(datt2[, c("allele", "simN")], meanInput, meanOutput, ratio)		
		colnames(dattNew)[3:4] = c("input" , "output")
		dattNew$keep =  !is.na(dattNew$ratio)

		#### oligo counts
		oligo = aggregate(datt2[,c(3:(2+2*nrep))], by=list(datt2$allele, datt2$simN), FUN=sum, na.rm=T)
		colnames(oligo)[1:2] = c("allele", "simN")
		
		#### apply tests
		combos = unique(dattNew$simN)
		input_MW = sapply(combos, run_MW_input, dattNew, ntag, nsim)
		resultAll = data.frame(simN = combos, resInput = input_MW)

		if("MW" %in% method | "Adaptive" %in% method)
		{	message("==== Applying the Mann-Whitney test... \n")
			result_MW <- run_MW_ratio(dattNew, ntag, nsim)
			resultAll$resMW = result_MW$res_MW
		}
		if("Matching" %in% method | "Adaptive" %in% method)
		{	message("==== Applying the matching test... \n")
			result_matching <- run_matching(dattNew, ntag, nsim=length(combos), input_MW)
			resultAll$resMatching = result_matching$res_matching
		}
		if("Adaptive" %in% method)
		{	message("==== Applying the adaptive test... \n")
			result_adaptive <- run_adaptive(combos, 0.05, ntag, nsim=length(combos), result_MW$res_MW, result_matching$res_matching, input_MW)
			resultAll$resAdaptive = result_adaptive$res_adaptive
		}
		### other methods work for both cases
		if("QuASAR" %in% method)
		{	message("==== Applying the QuASAR-MPRA method... \n")
			result_quasar = run_quasar(datt2, ntag, nsim, nrep )
			#### collect results
			if(!all(is.na(result_quasar)))
    		{	resultAll = merge(resultAll, result_quasar[,c("simN", "pvalue")], by="simN", sort=F, all.x=T)
    			colnames(resultAll)[ncol(resultAll)] = "QuASAR"
    		}
		}
		if("Fisher" %in% method)
		{	message("==== Applying the Fisher's Exact test... \n")
			result_fisher = run_fisherExact(oligo, ntag, nsim, nrep)
			resultAll = merge(result_all, result_fisher, by="simN", sort=F, all.x=T)
		}

		
		if(nrep > 1)
		{	if("T-test" %in% method)
			{	message("==== Applying the T-test... \n")
				result_t = run_ttest(oligo, nsim, nrep)
				resultAll = merge(resultAll, result_t, by="simN", sort=F, all.x=T)
			}
			if("mpralm" %in% method)
			{	message("==== Applying the mpralm methods... \n")
				result_mpralm = run_mpralm(datt3, ntag, nsim, nrep)
				if(all(!is.na(result_mpralm[[1]])))
				{	result_mpralm[[1]]$simN = rownames(result_mpralm[[1]])
					result_mpralm[[2]]$simN = rownames(result_mpralm[[2]])
					result_mpralm2 = merge(result_mpralm[[1]], result_mpralm[[2]], by="simN", sort=F, all.x=T)
					result_mpralm3 = result_mpralm2[, c("simN", "P.Value.x", "P.Value.y")]
					colnames(result_mpralm3) = c("simN", "mpralm_mean", "mpralm_sum")
					resultAll = merge(resultAll, result_mpralm3, by="simN", sort=F, all.x=T)
				}
			}
			if("edgeR" %in% method)
			{	message("==== Applying the edgeR method... \n")
				result_edgeR = run_edgeR(oligo, ntag, nsim, nrep)
				resultAll = merge(resultAll, result_edgeR[,c("simN", "PValue")], by="simN", sort=F, all.x=T)
				colnames(resultAll)[ncol(resultAll)] = "edgeR"

			}
			if("DESeq2" %in% method)
			{	message("==== Applying the DESeq2 method... \n")
				result_deseq2 = run_deseq2(oligo, ntag, nsim, nrep)
				resultAll = merge(resultAll, as.data.frame(result_deseq2[,c("simN", "pvalue")]), by="simN", sort=F, all.x=T)
				colnames(resultAll)[ncol(resultAll)] = "DESeq2"
			}
		}
		
    	return(resultAll)
}

