program define regife, eclass sortpreserve

	version 12
	syntax varlist(min=1 numeric fv ts) [if] [in] [aweight fweight pweight iweight], /// 
	Factors(string)   ///
	[ccep ccemg Dimension(string) Absorb(string) noCONS TOLerance(real 1e-6) MAXIterations(int 10000) *]


	/* tempname */
	tempvar res res2 y2 g1 g2
	tempname b V

	/* syntax */
	if ("`weight'"!=""){
		local wt [`weight'`exp']
		local wtv = subinstr("`exp'","=","", .)
		local wtv2 "*`wtv'"
		display as text "Weight are used only for regressions, not for the PCA"
		local sumwt [aw`exp']
	}



	/* syntax factors */
	local factors = trim("`factors'")
	if regexm("`factors'", "(.+[^=]*) ([^=]*.+)"){
		local id = regexs(1)
		local time = regexs(2)
	}

	if regexm("`id'", "(.*)=(.*)"){
		local idgen `= regexs(1)'
		forval i = 1/`dimension'{
			confirm new variable `idgen'_`i'
		}
		local id = regexs(2)
	}
	confirm var `id'

	if regexm("`time'", "([^ ]*)=([^ ]*)"){
		local timegen `= regexs(1)'
		forval i = 1/`dimension'{
			confirm new variable `timegen'_`i'
		}
		local time = regexs(2)
	}
	confirm var `time'




	/* touse */
	marksample touse
	markout `touse' `id' `time' `wtv', strok


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
	/* count number of observations. since I don't look at syntax of absorb, I need to count after hdfe redefines touse */
	sort `touse'
	qui count if `touse' 
	local touse_first = _N - r(N) + 1
	local touse_last = _N
	local obs = `touse_last'-`touse_first' + 1


	/* create group for i and t */
	sort `touse' `id'
	qui by `touse' `id': gen long `g1' = _n == 1 if `touse'
	qui replace `g1' = sum(`g1') if `touse'
	local id `g1'
	local N = `id'[_N]

	sort `touse' `time'
	qui by `touse' `time': gen long `g2' = _n == 1 if `touse'
	qui replace `g2' = sum(`g2') if `touse'
	local time `g2'
	local T = `time'[_N]


	cap assert `T' >= `dimension'
	if _rc{
		di as error "The dimension `dimension' should be lower than the number of distinct values for `time'"
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
	mata: iteration("`py'","`res'", "`pcx'", "`id'", "`time'", `N', `T', `dimension', `tolerance', `maxiterations', "`b'", `touse_first', `touse_last', "`idgen'", "`timegen'")
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
	ereturn local id1 `id'
	ereturn local id2 `time'
	ereturn local f1 `idgen'
	ereturn local f2 `timegen'
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
	void meanvar(real matrix b, real matrix w, string scalar sb1, string scalar sV1){
		b1 = J(1, cols(b), .)
		V1 =  J(1, cols(b), .)
		for(i=1;i<= cols(b);++i){
			b1[1, i] = mean(b[.,i], w)
			V1[1, i] = variance(b[.,i], w)
		}
		st_matrix(sb1, editmissing(b1,0))
		st_matrix(sV1, editmissing(diag(V1),0))
	}
end


mata:
	void iteration(string scalar y, string scalar res, string scalar x, string scalar id, string scalar time, real scalar N, real scalar T, real scalar d, real scalar tolerance, maxiterations, string scalar bname, real scalar first, real scalar last, string scalar idgen, string scalar timegen){
		real matrix Y 
		real matrix X
		real matrix tY
		real matrix M
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
		tY = st_data((first::last), res)
		b1 = st_matrix(bname)'
		X = st_data((first::last), x)
		M = invsym(cross(X, X)) * X'

		R = J(N, T, 0)
		V = J(T, T, .)
		iter = 0
		error = 1
		while (((maxiterations == 0) | (iter < maxiterations)) & (error >= tolerance)){
			iter = iter + 1
			tY = Y :-  X * b1
			for (obs = first; obs <= last ; obs++) {    
				R[_st_data(obs, iindex), _st_data(obs, tindex)] = tY[obs - first + 1, 1]
			}
			_svd(R, s, V)
			U = R[.,(1::d)] * diag(s[1::d]) 
			R = U *  V[(1::d),.]
			for (obs = first; obs <= last ; obs++) {    
				tY[obs - first + 1,1] = R[_st_data(obs, iindex),_st_data(obs, tindex)] 
			}	
			b2 = M * (Y :- tY)
			error = max(abs(b2 :- b1))
			b1 = b2
		}
		st_store(first::last, res, tY)
		st_numscalar("r(N)", iter)
		st_numscalar("r(error)", error)
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
	}
end




