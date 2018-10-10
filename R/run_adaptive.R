run_adaptive <-
function(combos, inputCutoff, ntag, nsim, res_MW, res_matching, inputTest)
{ 		result = data.frame(combos, res_adaptive=NA)
   		for(x in combos)
   		{
  			if(inputTest[x] >= inputCutoff)
  			{		result$res_adaptive[x] = res_MW[x]
  			}else
  			{		result$res_adaptive[x] = res_matching[x]
  			}
   		}
   		return (result)
}
