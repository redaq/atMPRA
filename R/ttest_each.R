ttest_each = function(ratios2, sim, nrep)
{   	#cat(sim)
		temp = ratios2[ratios2$simN==sim,1:nrep]
		if(nrow(temp)>1 & all(!is.na(temp)) & ((sum(mean(as.numeric(temp[1,]))==temp[1,])<nrep) & (sum(mean(as.numeric(temp[2,]))==temp[2,])<nrep)))
		{
			res1=t.test(x=as.numeric(temp[1,]), y=as.numeric(temp[2,]),  paired=TRUE, var.equal=FALSE)
			res2=t.test(x=as.numeric(temp[1,]), y=as.numeric(temp[2,]),  var.equal=FALSE)
			return(list(res1=res1$p.value, res2=res2$p.value))
		}else
		{	return(list(res1=-1, res2=-1))
		}	
}

