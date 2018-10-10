run_matching <-
function(datt, ntag, nsim, input_MW, width=5)
{   	combos = unique(datt$simN)
		matchingTest = sapply(combos, run_matching_one, datt, width, simplify=F)
		testR = data.frame(input=input_MW, matching=sapply(1:length(matchingTest), function(x) matchingTest[[x]]$pvalue))
		## check enrichement
		testTemp  = testR[testR$matching >=0 & testR$input>=0 &!is.na(testR$matching) &!is.na(testR$input),]
		ddtemp = testTemp[testTemp$matching<0.05,]
		enrichP = dhyper(x=sum(ddtemp$input<0.05), m = nrow(ddtemp), n=nrow(testTemp)-nrow(ddtemp), k = sum(testTemp$input<0.05)) + phyper(q=sum(ddtemp$input<0.05), m = nrow(ddtemp), n=nrow(testTemp)-nrow(ddtemp), k = sum(testTemp$input<0.05), lower.tail=F)
		while(enrichP < 0.05 & width>1)
		{	width = max(1, width-1)
			#cat(width, ":")
			matchingTest = sapply(combos, run_matching_one, datt, width, simplify=F)
			testR = data.frame(input=input_MW, matching=sapply(1:length(matchingTest), function(x) matchingTest[[x]]$pvalue))
			## check enrichement
			testTemp  = testR[testR$matching >=0 & testR$input>=0&!is.na(testR$matching) &!is.na(testR$input),]
			ddtemp = testTemp[testTemp$matching<0.05,]
			enrichP = dhyper(x=sum(ddtemp$input<0.05), m = nrow(ddtemp), n=nrow(testTemp)-nrow(ddtemp), k = sum(testTemp$input<0.05)) + phyper(q=sum(ddtemp$input<0.05), m = nrow(ddtemp), n=nrow(testTemp)-nrow(ddtemp), k = sum(testTemp$input<0.05), lower.tail=F)
			#cat(enrichP, "\n")
	
		}
		message("Final width: ", width, ".\n")
		return(data.frame(simN= combos, res_matching=sapply(1:length(matchingTest), function(x) matchingTest[[x]]$pvalue)))
	
	
	
}
