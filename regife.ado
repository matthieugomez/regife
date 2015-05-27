/***************************************************************************************************

***************************************************************************************************/
program define regife, eclass sortpreserve

	version 12
	syntax varlist(min=1 numeric fv ts) [if] [in] [aweight fweight pweight iweight], /// 
	Factors(string)  Dimension(int) ///
	[ccep ccemg  Absorb(string) noCONS TOLerance(real 1e-6) MAXIterations(int 10000) VERBose *]


	/* tempname */
	tempvar res res2 y2 g1 g2
	tempname b V

	/* syntax */




	/* syntax factors */
	local factors = trim("`factors'")
	if regexm("`factors'", "(.+[^=]) ([^=].+)"){
		local id1 = regexs(1)
		local id2 = regexs(2)
	}

	if regexm("`id1'", "(.*)=(.*)"){
		local id1gen `= regexs(1)'
		forval i = 1/`dimension'{
			confirm new variable `id1gen'_`i'
		}
		local id1 = regexs(2)
	}
	confirm var `id1'

	if regexm("`id2'", "([^ ]*)=([^ ]*)"){
		local id2gen `= regexs(1)'
		forval i = 1/`dimension'{
			confirm new variable `id2gen'_`i'
		}
		local id2 = regexs(2)
	}
	confirm var `id2'


	if ("`weight'"!=""){
		local wt [`weight'`exp']
		local wtv = subinstr("`exp'","=","", .)
		local wtv2 "*`wtv'"
		display as text "Weight are used for regressions, and are summed accross variable `id2' for the PCA"
		local sumwt [aw`exp']
	}



	/* touse */
	marksample touse
	markout `touse' `id1' `id2' `wtv', strok


	/*syntax varlist  */
	fvrevar `varlist' if `touse'
	local varlist = r(varlist)
	foreach v in `varlist'{
		local tvar : char `v'[tsrevar]
		local fvar : char `v'[fvrevar]
		if "`tvar'`fvar'" ~=""{
			local namelist `namelist'  `tvar'`fvar'
		}
		else{
			local namelist `namelist' `v'
		}
	}
	tokenize `varlist'
	local y `1'
	macro shift 
	local x `*'
	tokenize `namelist'
	local yname `1'
	macro shift 
	local xname `*'



	/* absorb */
	tempname df_a
	if "`absorb'" ~= ""{
		cap which hdfe.ado
		if _rc {
			di as error "hdfe.ado required when using multiple absorb variables:"
			di as error "{stata cap ado uninstall hdfe}"
			di as error "{stata net install hdfe, from(https://raw.githubusercontent.com/sergiocorreia/reghdfe/master/package)}"
		}
		tempvar sample
		tempname prefix
		hdfe `y' `x'  `wt' if `touse', a(`absorb') gen(`prefix') sample(`sample')
		scalar `df_a' = r(df_a)
		qui replace `touse' = 0 if `sample' == 0
		tempvar `prefix'`y'
		qui gen  ``prefix'`y'' = `prefix'`y' 
		drop `prefix'`y'
		local py ``prefix'`y''
		foreach v in `x'{
			tempvar `prefix'`v'
			qui gen  ``prefix'`v'' = `prefix'`v'
			drop `prefix'`v'
			local px `px' ``prefix'`v''
		}
		local cons nocons
	}
	else{
		local px `x'
		local py `y'
		scalar `df_a' = 0
	}
	/* count number of observations (after hdfe since it reads the absorb syntax) */
	sort `touse'
	qui count if `touse' 
	local touse_first = _N - r(N) + 1
	local touse_last = _N
	local obs = `touse_last'-`touse_first' + 1


	/* create group for i and t */
	sort `touse' `id1'
	qui by `touse' `id1': gen long `g1' = _n == 1 if `touse'
	qui replace `g1' = sum(`g1') if `touse'
	local N = `g1'[_N]

	sort `touse' `id2'
	qui by `touse' `id2': gen long `g2' = _n == 1 if `touse'
	qui replace `g2' = sum(`g2') if `touse'
	local T = `g2'[_N]


	cap assert `T' >= `dimension'
	if _rc{
		di as error "The dimension `dimension' should be lower than the number of distinct values for `id2'"
		exit 0
	}

	/* nocons */
	if "`cons'" == ""{
		tempvar vcons
		qui gen `vcons' = 1
		local pcx `px' `vcons'
		local xname `xname' _cons
	}
	else{
		local pcx `px'
	}

	/* algorithim */
	* first reg  
	qui _regress `py' `pcx' `wt' if `touse', nocons
	matrix `b' = e(b)
	qui predict `res' if `touse'



	* iterate 
	mata: iteration("`py'","`res'", "`pcx'", "`wtv'", "`g1'", "`g2'", `N', `T', `dimension', `tolerance', `maxiterations', "`b'", `touse_first', `touse_last', "`id1gen'", "`id2gen'", "`verbose'")
	local iter = r(N)
	tempname error
	scalar `error' = r(error)


	* last reg
	qui gen `res2' = `py' - `res'
	qui reg `res2' `px' `wt' if `touse', `cons' `options'

	/* return results */
	tempname df_m
	scalar `df_m' = e(df_m) + (`N'+`T')* `dimension' + `df_a'
	tempname df_r
	scalar `df_r' = `obs' - `df_m'

	mat `b' = e(b)
	mat colnames `b' =`xname'
	mat `V' = e(V) * e(df_r)/ `df_r'
	mat colnames `V' =`xname'
	mat rownames `V' =`xname'


	tempname r2w 
	scalar `r2w' = e(r2)
	qui sum `y' `sumwt'
	local rtotal = r(sd)^2
	qui sum `res2' `sumwt'
	local rpartial = r(sd)^2
	tempname r2
	scalar `r2' = (`r2w' * `rpartial' + `rtotal'- `rpartial') / `rtotal'

	di
	if `iter' == `maxiterations'{
		display as text "The algorithm did not converge : error is" in ye %4.3gc `error' in text " (higher than tolerance" in ye %4.3gc `tolerance' in text")"
		display as text "Use the maxiterations options to increase the amount of iterations"
	}

	ereturn post `b' `V', depname(`yname') obs(`obs') esample(`touse') 
	ereturn scalar df_r = `df_r'
	ereturn scalar df_m = `df_m'
	ereturn scalar r2_within = `r2w'
	ereturn scalar r2 = `r2'
	ereturn scalar error = `error'


	ereturn local absorb `absorb'
	ereturn local id1 `id1'
	ereturn local id2 `id2'
	ereturn local f1 `id1'
	ereturn local f2 `id2gen'
	ereturn local d `dimension'
	ereturn local predict "regife_p"

	ereturn display
	display as text "{lalign 26:Number of obs = }" in ye %10.0fc `obs'
	display as text "{lalign 26:R-squared  = }" in ye %10.3fc `r2'
	display as text "{lalign 26:Within R-sq  = }" in ye %10.3fc `r2w'
	display as text "{lalign 26:Number of iterations = }" in ye %10.0fc `iter'

end

/***************************************************************************************************
helper functions
***************************************************************************************************/
set matastrict on
mata:


	void iteration(string scalar y, string scalar res, string scalar x, string scalar w, string scalar id, string scalar time, real scalar N, real scalar T, real scalar d, real scalar tolerance, maxiterations, string scalar bname, real scalar first, real scalar last, string scalar idgen, string scalar timegen, string scalar verbose){
		real matrix Y 
		real matrix X
		real matrix tY
		real matrix M
		real matrix Ws
		real scalar iindex
		real scalar tindex

		real scalar iter
		real scalar obs
		real scalar col

		real scalar error
		real matrix U
		real matrix V
		real colvector s
		real matrix R
		real colvector b1
		real colvector b2

		string scalar name

		iindex = st_varindex(id)
		tindex = st_varindex(time)
		Y = st_data((first::last), y)
		b1 = st_matrix(bname)'
		X = st_data((first::last), x)

		Ws = J(N, T, 1)
		if (strlen(w) > 0) {
			W = st_data((first::last), w)
			M = invsym(cross(X, W, X)) * X'* diag(W)
			for (obs = first; obs <= last ; obs++) {    
				Ws[_st_data(obs, iindex), _st_data(obs, tindex)] = W[obs - first + 1, 1]
			}
			Ws= rowsum(Ws)
			Ws = sqrt(Ws:/mean(Ws))
		}
		else{
			M = invsym(cross(X, X)) * X'
		}

		R = J(N, T, 0)
		V = J(T, T, .)
		iter = 0
		error = 1
		while (((maxiterations == 0) | (iter < maxiterations)) & (error >= tolerance)){
			if (strlen(verbose) > 0){
				if ((mod(iter, 100)==0) & (iter > 0)){
					if (iter == 100){
						stata(`"display as text "each .=100 iterations""')
					}
					stata(`"display "." _c"')
				}
			}
			iter = iter + 1
			tY = Y :-  X * b1
			for (obs = first; obs <= last ; obs++) {    
				R[_st_data(obs, iindex), _st_data(obs, tindex)] = tY[obs - first + 1, 1]
			}
			if (strlen(w) > 0) {
				R = Ws :* R
			}
			_svd(R, s, V)
			U = R[.,(1::d)] * diag(s[1::d]) 
			if (strlen(w) > 0) {
				 U = U :/ Ws
			}
			R = U *  V[(1::d),.]
			for (obs = first; obs <= last ; obs++) {    
				tY[obs - first + 1,1] = R[_st_data(obs, iindex),_st_data(obs, tindex)] 
			}	
			b2 = M * (Y :- tY)
			error = max(abs(b2 :- b1))
			b1 = b2
		}
		st_store(first::last, res, tY)
		if (strlen(idgen) > 0){
			for (col = 1; col <= d; col++){
				name =  idgen +  "_" + strofreal(col)
				st_addvar("float", name)
				for (obs = first; obs <= last ; obs++) { 
					st_store(obs, 	name, U[_st_data(obs, iindex), col])
				} 
			}
		}
		if (strlen(timegen) > 0){
			for (col = 1; col <= d; col++){
				name =  timegen + "_" + strofreal(col)
				st_addvar("float", name)
				for (obs = first; obs <= last ; obs++) { 
					st_store(obs, name, V[col, _st_data(obs, tindex)])
				} 
			}
		}
		st_numscalar("r(N)", iter)
		st_numscalar("r(error)", error)
	}
end

/***************************************************************************************************
wls in mata
set matastrict on
cap mata: mata drop wols()
mata:
	void wols(string scalar y, string scalar x, string scalar w, numeric scalar first, numeric scalar last, string scalar newres){
		real matrix Y 
		real matrix X
		real matrix W
		real matrix M
		real scalar b
		real matrix res
		Y = st_data((first::last), y)
		X = st_data((first::last), x)
		if (strlen(w) > 0) {
			W = st_data((first::last), w)
			/*W = sqrt(W/mean(W))  */
			M = invsym(cross(X, W, X)) * X'* diag(W)
		}
		else{
			M = invsym(cross(X, X)) * X'
		}
		b = M * Y 
		b
		res = Y :-  X * b
		st_addvar("float", newres)
		st_store(first::last, newres, res)

	}
end
***************************************************************************************************/



