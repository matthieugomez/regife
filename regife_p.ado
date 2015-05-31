program define regife_p, eclass sortpreserve
	syntax newvarname [if] [in] , [RES RESF XB XBF F]

	* Does not currently support previous estimation using fwl or partial
	if "`e(absorb)'" != ""  {
		di in r "predict not supported after regife with absorb option"
		error 499
	}

	local option `res' `resf' `xb' `xbf' `f'

	marksample touse,  novarlist
	if "`option'"=="xb"{
		predict `varlist' if `touse' == 1
	}
	else if "`option'" == "resf"{
		qui _predict `varlist' if `touse' == 1
		replace `varlist' = `=e(depvar)' - `varlist'
	}
	else{
		replace `touse' = 0 if e(sample) == 0
		if ("`=e(f1)'"=="" | "`=e(f2)'"=="") {
			di as error "In order to predict, the factors need to be saved with the ife option (#`g' was not)" _n "For instance, instead of {it:ife(firm year)}, set ife(fe1=firm fe2=year)"
			exit 112
		}
		tempvar resf
		gen `resf' = 0 if `touse'
		forv r = 1/`=e(d)'{
			qui replace `resf' = `resf' + `=e(f1)'`r' * `=e(f2)'`r'
		}
		if "`option'" == "f"{
			qui gen `varlist' = `resf'
		}
		else if "`option'" == "res"{
			qui _predict `varlist' if `touse' == 1
			qui replace `varlist' =`=e(depvar)' - `resf' - `varlist'
		}
		else if "`option'" == "xbf"{
			qui _predict `varlist' if `touse' == 1
			qui replace `varlist' = `varlist' + `resf'
		}
	}
end





