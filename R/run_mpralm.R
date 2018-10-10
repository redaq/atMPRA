run_mpralm <-
function(datt3, ntag, nsim, nrep)
{
		datt3 = datt3[order(datt3$allele),]
		indexd = which(grepl("Ref", datt3$allele))
		ref = datt3[indexd,]
		refDNAmean = apply(ref[, c(3:(2+nrep))], 1, mean, na.rm=T)
		alt = datt3[-indexd,]
		altDNAmean = apply(alt[, c(3:(2+nrep))], 1, mean, na.rm=T)
	
		newdna = t(do.call(rbind, replicate(nrep, refDNAmean, simplify=FALSE)) )
		newdna2 = t(do.call(rbind, replicate(nrep, altDNAmean, simplify=FALSE)) )
		newdna3 = cbind(newdna, newdna2)

		newData = cbind(ref, alt)
		if(sum(newData[,2]==newData[,(3+nrep*2+1)])!=nrow(newData))
		{	stop("not matched simN")
		}
		eid = as.character(newData$simN)
		rna = cbind(ref[,(2+nrep+1):(2+nrep+nrep)], alt[,(2+nrep+1):(2+nrep+nrep)])
		dna2 = newdna3
		rownames(dna2) = rownames(rna) = paste(newData$simN, rep(1:ntag, nsim), sep="_")
		colnames(dna2) = colnames(rna) = c(paste0("rep_A_", 1:nrep),paste0("rep_T_", 1:nrep))
		mpraob = MPRASet(DNA=dna2, RNA=rna, eid=eid, eseq=NULL, barcode=NULL)
	
		design=data.frame(intcpt=1, grepl("A", colnames(mpraob)))
		toptab=toptab2=NA
		tryCatch(
		{
		mprafit = mpralm(object=mpraob, design=design, aggregate="mean", normalize=F, model_type="indep_groups", plot=F)
		toptab= topTable(mprafit, coef=2, number=Inf)
		},error=function(e)
		{	toptab=NA
		})

		tryCatch(
		{
		mprafit2 = mpralm(object=mpraob, design=design, aggregate="sum", normalize=F, model_type="indep_groups",  plot=F)
		toptab2= topTable(mprafit2, coef=2, number=Inf)
		},error=function(e)
		{	toptab2 = NA
		})
		return(list(averageEst_corr = toptab, aggregateEst_corr= toptab2))

}
