run_MW_ratio <-
function(datt, ntag, nsim)
{	  combos = unique(datt$simN)
	  res_MW = sapply(combos, run_MW_ratio_one, datt)
	  return(data.frame(simN=combos, res_MW = res_MW))
}
