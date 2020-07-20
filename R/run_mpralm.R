run_mpralm <-
function(datt3, ntag, nsim, nrep, paired=FALSE)
{
		datt3 = datt3[order(datt3$allele, datt3$simN),]
		indexd = which(grepl("Ref", datt3$allele))
		ref = datt3[indexd,]
		alt = datt3[-indexd,]
	
		## the assumption is both alleles have the same number of barcodes
		newData = cbind(ref, alt)
		if(sum(newData[,2]==newData[,(3+nrep*2+1)])!=nrow(newData))
		{	stop("not matched simN")
		}
		eid = as.character(newData$simN)
		rna = cbind(ref[,(2+nrep+1):(2+nrep+nrep)], alt[,(2+nrep+1):(2+nrep+nrep)])
		dna2 = newData[,c(3:(2+nrep), (2+nrep*2+3):(2+nrep*2+2+nrep))]
		rownames(dna2) = rownames(rna) = paste(newData$simN, rep(1:ntag, nsim), sep="_")
		colnames(dna2) = colnames(rna) = c(paste0("rep_A_", 1:nrep),paste0("rep_T_", 1:nrep))
		mpraob = MPRASet(DNA=dna2, RNA=rna, eid=eid, eseq=NULL, barcode=NULL)
	
		design=data.frame(intcpt=1, grepl("A", colnames(mpraob)))
		toptab=toptab2=NA
		if(paired)
		{	model_type ="corr_groups"
			block_vector =  rep(1:nrep, 2)
		}else
		{	model_type="indep_groups"
			block_vector=NULL
		}
		tryCatch(
		{
		mprafit = mpralm(object=mpraob, design=design, aggregate="mean", normalize=F, model_type=model_type, block=block_vector, plot=F)
		toptab= topTable(mprafit, coef=2, number=Inf)
		},error=function(e)
		{	toptab=NA
		})

		tryCatch(
		{
		mprafit2 = mpralm(object=mpraob, design=design, aggregate="sum", normalize=F, model_type=model_type, block=block_vector, plot=F)
		toptab2= topTable(mprafit2, coef=2, number=Inf)
		},error=function(e)
		{	toptab2 = NA
		})
		return(list(averageEst_corr = toptab, aggregateEst_corr= toptab2))

}
