/***************************************************************************************************

***************************************************************************************************/
program define regife, sortpreserve
	version 12
	syntax anything, Factors(string) [CLuster(string) reps(string) *]

	if "`reps'" == "" {
		di as text "Use the option reps() to compute errors by bootstrap"
		innerregife `anything', factors(`factors') `options'
	}
	else{
		if "`cluster'" == ""{
			bootstrap, reps(`reps') : ///
			innerregife `anything', factors(`factors') `options'
		}
		else{
			tempvar clusterid
			tsset, clear
			bootstrap, reps(`reps') cluster(`cluster') idcluster(`clusterid'): ///
			innerregife `anything', factors(`factors') `options'
		}
	}
end



