cap program drop regife
program define regife, eclass sortpreserve

	version 12
	syntax anything [if] [in], Factors(string) Dimension(integer)  [Absorb(varlist) noCONS convergence(real 0.000001) MAXiteration(int 10000) GENerate(string)]


	tokenize `anything'
	local yname `1'
	macro shift 
	local xname `*'

	fvrevar `anything', tsonly 
	local anything = r(varlist)
	tokenize `anything'
	local y `1'
	macro shift 
	local x `*'


	if "`gen'" ~= ""{
		confirm new variable `gen'
	}

	local factors = trim("`factors'")
	if regexm("`factors'", "(.+[^=]*) ([^=]*.+)"){
		local id = regexs(1)
		local time = regexs(2)
	}

	if regexm("`id'", "(.*)=(.*)"){
		local idgen `= regexs(1)'
		cap confirm new variable `idgen'
		local id = regexs(2)
	}
	confirm var `id'

	if regexm("`time'", "([^ ]*)=([^ ]*)"){
		local timegen `= regexs(1)'
		cap confirm new variable `timegen'
		local time = regexs(2)
	}
	confirm var `time'

	marksample touse
	markout `touse' `id' `time' `y' `x'



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
		qui gen  ``prefix'`oldy'' = `prefix'`oldy' 
		drop `prefix'`oldy'
		local y ``prefix'`oldy''
		foreach v in `oldx'{
			tempvar `prefix'`v'
			qui gen  ``prefix'`v'' = `prefix'`v'
			drop `prefix'`v'
			local x `x' ``prefix'`v''
		}
		local nocons cons
	}
	else{
		local oldx `x'
		local oldy `y'
		scalar `df_a' = 0
	}



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


	cap assert `T' >= `dimension'
	if _rc{
		di as error "The dimension `dimension' should be lower than the number of distinct values for `time'"
		exit 0
	}

	if "`nocons'" == ""{
		tempvar cons
		qui gen `cons' = 1
		local xc `x' `cons'
		local oldx `oldx' _cons
		local xname `xname' _cons
	}
	else{
		local xc `x'
	}

	/* reg  */
	qui _regress `y' `xc' if `touse', nocons
	matrix `b' = e(b)
	qui predict `res' if `touse', res

	/* count after potential redefinition by absorb */
	sort `touse'
	qui count if `touse' 
	local touse_first = _N - r(N) + 1
	local touse_last = _N
	local obs = `touse_last'-`touse_first' + 1


	mata: iteration("`y'","`res'", "`xc'", "`id'", "`time'", `N', `T', `dimension', `convergence', `maxiteration', "`b'", `touse_first', `touse_last', "`idgen'", "`timegen'")
	local iter = r(N)
	tempname error
	scalar `error' = r(error)


	qui gen `res2' = `y' - `res'
	if "`generate'" ~= ""{
		rename `res' `generate'
	}
	qui reg `res2' `xc', nocons



	/* results */
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
	qui sum `y'
	local rtotal = r(sd)^2
	qui sum `res2'
	local rpartial = r(sd)^2
	tempname r2
	scalar `r2' = (`r2w' * `rpartial' + `rtotal'- `rpartial') / `rtotal'


	ereturn post `b' `V', depname(`yname') obs(`obs') esample(`touse') 
	ereturn scalar df_r = `df_r'
	ereturn scalar df_m = `df_m'
	ereturn scalar r2_within = `r2w'
	ereturn scalar r2 = `r2'
	ereturn scalar error = `error'

	ereturn display
	display as text "{lalign 26:Number of obs = }" in ye %10.0fc `obs'
	display as text "{lalign 26:R-squared  = }" in ye %10.3fc `r2'
	display as text "{lalign 26:Within R-sq  = }" in ye %10.3fc `r2w'
	display as text "{lalign 26:Number of iterations = }" in ye %10.0fc `iter'

	if `iter' == `maxiteration'{
		display as error "{lalign 26:Convergence error = }" in ye %10.3gc `error'
	}

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

