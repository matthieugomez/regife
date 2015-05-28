/***************************************************************************************************

***************************************************************************************************/
program define regife, sortpreserve
	version 12
	if regexm("`0'", "(.*),(.*)"){
		local anything = regexs(1)
		local option = regexs(2)
	}
	local 0 , `option'
	syntax [, reps(string) Factors(string) CLuster(string) *]
	if "`reps'" == "" {
		di as text "Use the option reps() to compute correct Standard Errors"
		innerregife `anything', factors(`factors') `options'
	}
	else{
		if "`cluster'" == ""{
			bootstrap, reps(`reps') : ///
			innerregife `anything', factors(`factors')  `options'
		}
		else{
			tempvar clusterid
			tsset, clear
			bootstrap, reps(`reps') cluster(`cluster') idcluster(`clusterid'): ///
			innerregife `anything', factors(`factors') `options'
		}
	}
end


