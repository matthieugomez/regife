cap program drop regife
program define regife, eclass sortpreserve
	version 13
	syntax anything [if] [in] [aweight/], Id(string) Time(string) Dimension(integer)  [Absorb(varlist) noCONS convergence(real 0.000001) MAXiteration(int 500) gen(string)]

	/* syntax  */
	fvrevar `anything'
	local anything = r(varlist)
	tokenize `anything'
	local y `1'
	macro shift 
	local x `*'

	if "`gen'" ~= ""{
		confirm new variable `gen'
	}

	if regexm("`id'", "(.*)=(.*)"){
		local idgen `= regexs(1)'
		cap confirm new variable `idgen'
		local id = regexs(2)
	}
	confirm var `id'

	if regexm("`time'", "(.*)=(.*)"){
		local timegen `= regexs(1)'
		cap confirm new variable `timegen'
		local time = regexs(2)
	}
	confirm var `time'

	marksample touse
	markout `touse' `id' `time' `absorb'



	tempname df_a
	if "`absorb'" ~= ""{
		local oldx `x'
		local oldy `y'
		local x ""
		local y ""
		tempvar sample
		tempname prefix
		hdfe `oldy' `oldx' if `touse', a(`absorb') gen(`prefix') sample(`sample')
		scalar `df_a' = r(df_a)
		qui replace `touse' = 0 if `sample' == 0
		tempvar `prefix'`oldy'
		gen  ``prefix'`oldy'' = `prefix'`oldy' 
		drop `prefix'`oldy'
		local y ``prefix'`oldy''
		foreach v in `oldx'{
			tempvar `prefix'`v'
			gen  ``prefix'`v'' = `prefix'`v'
			drop `prefix'`v'
			local x `x' ``prefix'`v''
		}
		local nocons cons
	}
	else{
		scalar `df_a' = 0
	}

	/* count after potential redefinition by absorb */
	qui count if `touse'
	local touse_first = _N - r(N) + 1
	local touse_last = _N
	local obs = `touse_last'-`touse_first'


	/* tempname */
	tempvar res res2 y2 g1 g2
	tempname b V




	/* create group for i and t */
	sort `touse' `id'
	qui by `touse' `id': gen `g1' = _n == 1 if `touse'
	qui replace `g1' = sum(`g1') if `touse'
	local id `g1'
	local N = `id'[_N]

	sort `touse' `time'
	qui by `touse' `time': gen `g2' = _n == 1 if `touse'
	qui replace `g2' = sum(`g2') if `touse'
	local time `g2'
	local T = `time'[_N]


	if "`nocons'" == ""{
		tempvar cons
		gen `cons' = 1
		local xc `x' `cons'
		local oldx `oldx' _cons
	}
	else{
		local xc `x'
	}

	/* reg  */
	qui _regress `y' `xc', nocons
	matrix `b' = e(b)
	qui predict `res', res

	mata: iteration("`y'","`res'", "`xc'", "`id'", "`time'", `N', `T', `dimension', `convergence', `maxiteration', "`b'", `touse_first', `touse_last', "`idgen'", "`timegen'")
	local iter = r(N)
	tempname error
	scalar `error' = r(error)

	gen `res2' = `y' - `res'
	if "`gen'" ~= ""{
		rename `res' `gen'
	}
	qui reg `res2' `xc'

	mat `b' = e(b)
	tempname df_m
	scalar `df_m' = e(df_m) + (`N'+`T')* `dimension' + `df_a'
	tempname df_r
	scalar `df_r' = `obs' - `df_m'
	mat `V' = e(V) * e(df_r)/ `df_r'

	mat colnames `b' =`oldx'
	mat colnames `V' =`oldx'
	mat rownames `V' =`oldx'


	tempname r2w
	scalar `r2w' = e(r2)
	ereturn post `b' `V', depname(`yname') obs(`obs') esample(`touse') 
	ereturn display

	qui sum `y'
	local rtotal = r(sd)^2
	qui sum `res2'
	local rpartial = r(sd)^2
	tempname r2
	scalar `r2' = (`r2w' * `rpartial' + `rtotal'- `rpartial') / `rtotal'

	ereturn scalar df_r = `df_r'
	ereturn scalar df_m = `df_m'
	ereturn scalar r2_within = `r2w'
	ereturn scalar r2 = `r2'
	ereturn scalar error = `error'

	display as text "{lalign 26:Number of obs = }" in ye %10.0fc `obs'
	display as text "{lalign 26:R-squared  = }" in ye %10.3fc `r2'
	display as text "{lalign 26:Within R-sq  = }" in ye %10.3fc `r2w'
	display as text "{lalign 26:Number of iterations = }" in ye %10.0fc `iter'
	display as text "{lalign 26:Convergence error = }" in ye %10.3gc `error'
	
end

/***************************************************************************************************
mata helper
***************************************************************************************************/
set matastrict on
mata:
	void iteration(string scalar y, string scalar res, string scalar x, string scalar id, string scalar time, real scalar N, real scalar T, real scalar d, real scalar convergence, maxiteration, string scalar bname, real scalar first, real scalar last, string scalar idgen, string scalar timegen){
		real matrix Y 
		real matrix X
		real matrix tY
		real matrix M
		real scalar iindex
		real scalar tindex
		real scalar obs
		real scalar i
		real scalar iter
		real scalar error
		real matrix V
		real colvector s
		real matrix R
		iindex = st_varindex(id)
		tindex = st_varindex(time)
		Y = st_data((first::last), y)
		X = st_data((first::last), x)
		tY = st_data((first::last), res)
		R = J(N, T, 0)
		M = invsym(cross(X, X)) * X'
		b1 = st_matrix(bname)'
		iter = 0
		while ((iter < maxiteration) & (error >= convergence)){
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
				name =  idgen + strofreal(col)
				st_addvar("float", name)
				for (obs = first; obs <= last ; obs++) { 
					st_store(obs, idgen +   strofreal(col), U[_st_data(obs, iindex), col])
				} 
			}
		}
		if (strlen(timegen) > 0){
			for (col = 1; col <= d; col++){
				name =  timegen + strofreal(col)
				st_addvar("float", name)
				for (obs = first; obs <= last ; obs++) { 
					st_store(obs, name, V[col, _st_data(obs, tindex)])
				} 
			}
		}
	}
end



