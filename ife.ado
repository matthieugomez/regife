program define ife, eclass sortpreserve
	version 13
	syntax varname [if] [in], Factors(string) Dimension(integer) GENerate(string) [TOLerance(real 1e-6) MAXIterations(int 10000)]

	local y `varlist'
	if "`generate'" ~= ""{
		confirm new variable `generate'
	}

	local factors = trim("`factors'")
	if regexm("`factors'", "(.+[^=]*) ([^=]*.+)"){
		local id = regexs(1)
		local time = regexs(2)
	}

	if regexm("`id'", "(.*)=(.*)"){
		local idgen `= regexs(1)'
		confirm new variable `idgen'
		local id = regexs(2)
	}
	confirm var `id'

	if regexm("`time'", "([^ ]*)=([^ ]*)"){
		local timegen `= regexs(1)'
		confirm new variable `timegen'
		local time = regexs(2)
	}
	confirm var `time'

	marksample touse
	markout `touse' `id' `time' `y', strok




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

	/* count after potential redefinition by absorb */
	sort `touse'
	qui count if `touse' 
	local touse_first = _N - r(N) + 1
	local touse_last = _N
	local obs = `touse_last'-`touse_first' + 1

	tempvar res
	gen `res' = `y'

	mata: iterationf("`res'", "`id'", "`time'", `N', `T', `dimension', `tolerance', `maxiterations', `touse_first', `touse_last', "`idgen'", "`timegen'")
	local iter = r(N)
	tempname error
	scalar `error' = r(error)
	display as text "{lalign 26:Number of iterations = }" in ye %10.0fc `iter'
	if `iter' == `maxiterations'{
		display as error "{lalign 26:Convergence error = }" in ye %10.3gc `error'
	}
	rename `res' `generate'
end

/***************************************************************************************************
mata helper
***************************************************************************************************/
set matastrict on
mata:
	void iterationf(string scalar y, string scalar id, string scalar time, real scalar N, real scalar T, real scalar d, real scalar convergence, maxiterations, real scalar first, real scalar last, string scalar idgen, string scalar timegen){
		real matrix Y 
		
		real scalar iindex
		real scalar tindex
		real scalar index

		real scalar iter
		real scalar obs
		real scalar col

		
		real scalar error
		real matrix U
		real matrix V
		real colvector s
		real matrix R1
		real matrix R2
		real matrix na

		string scalar name

		index = st_varindex(y)
		iindex = st_varindex(id)
		tindex = st_varindex(time)


		Y = J(N, T, .)
		R1 = J(N, T, 0)

		for (obs = first; obs <= last ; obs++) {    
			Y[_st_data(obs, iindex), _st_data(obs, tindex)] = _st_data(obs, index)
		}
		na = Y :==.
		Y = editmissing(Y, 0)
		error = 1
		iter = 0
		while (((maxiterations == 0) | (iter < maxiterations)) & (error >= convergence)){
			iter = iter + 1
			R2 = Y :+ na :* R1
			_svd(R2, s, V)
			U = R2[.,(1::d)] * diag(s[1::d]) 
			R2 = U *  V[(1::d),.]
			error = sqrt(sum((R2:-R1):^2)/ (cols(R1)*rows(R1)))
			R1 = R2
		}
		for (obs = first; obs <= last ; obs++) {    
			st_store(obs, index, R2[_st_data(obs, iindex), _st_data(obs, tindex)])
		}	


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

