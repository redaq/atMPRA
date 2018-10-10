getParam <-
function(meanNeed, sdNeed)
{ 		p = meanNeed/((sdNeed)^2)
	  	n = meanNeed*p/(1-p)
	  	return(list(p=p, n=n))
}
