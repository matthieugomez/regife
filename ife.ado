program define ife, eclass sortpreserve
	version 12.0
	syntax varname [if] [in] [aweight fweight pweight iweight], Factors(string)   [Absorb(string) RESiduals(string) TOLerance(real 1e-6) MAXIterations(int 10000) VERBose]

	/***************************************************************************************************
	check syntax
	***************************************************************************************************/

	/* syntax absorb */
	while (regexm("`absorb'", "[ ][ ]+")) {
		local absorb : subinstr local absorb "  " " ", all
	}
	local absorb : subinstr local absorb " =" "=", all
	local absorb : subinstr local absorb "= " "=", all
	

	tokenize `absorb'
	while "`1'" ~= ""{	 
		if regexm("`1'", "(.+)=(.+)"){
			cap confirm new variable `=regexs(1)'
			if _rc{
				if _rc == 198{
					di as error "variable `=regexs(1)' is not a valid name"
				}
				else if _rc == 110{
					di as error "variable `=regexs(1)' already defined"
				}
				exit 198
			}
			local savefe "yes"
			local absorbvars `absorbvars' `=regexs(2)'
		}
		else{
			local absorbvars `absorbvars' `1'
		}
		macro shift
	}

	/* syntax factors */
	if regexm("`factors'", "(.*),(.*)"){
		local factors  = regexs(1)
		local dimension = regexs(2)
	}
	else{
		di as error "dimensions should be specified within the option factors"
		exit
	}
	while (regexm("`factors'", "[ ][ ]+")) {
		local factors : subinstr local factors "  " " ", all
	}
	local factors : subinstr local factors " =" "=", all
	local factors : subinstr local factors "= " "=", all

	cap assert `: word count `factors'' == 2
	if _rc{
		di as error "There must be exactly two variables in the option factors"
		exit 
	}


	forv i = 1/2{	 
		local f : word `i' of `factors'
		if regexm("`f'", "(.+)=(.+)"){
			local id`i'gen `=regexs(1)'
			local id`i' = regexs(2)
			forval d = 1/`dimension'{
				cap confirm new variable `id`i'gen'`d'
				if _rc{
					if _rc == 198{
						di as error "variable `id`i'gen'`d' is not a valid name"
					}
					else if _rc == 110{
						di as error "variable `id`i'gen'`d' already defined"
					}
					exit 198
				}
			}
		}
		else{
			local id`i' `f'
		}
		confirm var `id`i''
	}


	if "`residuals'" ~= "" {
		confirm new variable `residuals'
	}

	if "`residuals'`id1gen'`id2gen'" == ""{
		di as error "Nothing to return. Either save the loading / factors (using the factors option),  or the residuals (using the residuals option)"
		exit 0
	}

	if ("`weight'"!=""){
		local wtype `weight'
		local wvar  `=subinstr("`exp'","=","", .)'
		local sumwt [aw=`wvar']
	}


	/***************************************************************************************************
	prepare the data
	***************************************************************************************************/

	/* touse */
	marksample touse
	markout `touse' `id1' `id2' `wvar', strok

	count if `touse'
	if r(N) == 0{
		di as error "There are no observations in the sample"
		exit 0
	}
	/* tempname */
	tempvar res res2 y2 g1 g2 y
	tempname b V


	/* residual */
	if "`absorb'" ~= ""{
		cap which reghdfe.ado
		if _rc {
			di as error "reghdfe.ado required when using multiple absorb variables: {stata ssc install reghdfe}"
			exit 111
		}
		qui reghdfe `varlist' `wt' if `touse', a(`absorbvars') residuals(`res')
		if "`savefe'" == "yes"{
			qui reghdfe `varlist' `wt' if `touse', a(`absorb') 
		}
		qui replace `touse' = e(sample)
	}
	else{
		display as text "Variable `varlist' is not demeaned. Use absorb(`id1') (resp. absorb(`id2')) if you want to demean with respect to `id1' (resp. `id2')"
		gen `res' = `varlist'
	}

	/* create group for i and t */
	if "`id1'" ~= "`: char _dta[_IDpanel]'" | "`id2'" == "`: char _dta[_TStvar]'"{
		sort `touse' `id1' `id2'
		cap bys `touse' `id1' `id2' : assert _N == 1 if `touse'
		if _rc{
			di as error "repeated observations for `id2' within `id1'"
			exit 451
		}
	}

	sort `touse' `id1'
	qui by `touse' `id1': gen double `g1' = _n == 1 if `touse'
	qui replace `g1' = sum(`g1') if `touse'
	local N = `g1'[_N]

	sort `touse' `id2'
	qui by `touse' `id2': gen double `g2' = _n == 1 if `touse'
	qui replace `g2' = sum(`g2') if `touse'
	local T = `g2'[_N]


	qui count if `touse' 
	local touse_first = _N - r(N) + 1
	local touse_last = _N
	local obs = `touse_last'-`touse_first' + 1



	cap assert `T' >= `dimension'
	if _rc{
		di as error "The factor structure dimension should be lower than the number of distinct values of the time variable"
		exit 0
	}


	/* mata program to find factors */
	mata: iterationf("`res'", "`g1'", "`g2'", "`wvar'", `N', `T', `dimension', `tolerance', `maxiterations', `touse_first', `touse_last', "`id1gen'", "`id2gen'", "`residuals'", "`verbose'")

	/* display results */
	local iter = r(N)
	tempname error
	scalar `error' = r(error)
	display as text "{lalign 26:Number of iterations = }" in ye %10.0fc `iter'
	if `iter' == `maxiterations'{
		display as text "The algorithm did not converge"
		display as text "{lalign 26:Convergence error = }" in ye %4.3gc `error' 
	}
end


/***************************************************************************************************
mata helper
***************************************************************************************************/
set matastrict on
mata:

	void iterationf(string scalar y, string scalar id, string scalar time, string scalar w, real scalar N, real scalar T, real scalar d, real scalar tolerance, real scalar maxiterations, real scalar first, real scalar last, string scalar id1gen, string scalar id2gen, string scalar residuals, string scalar verbose){
		

		
		real scalar iindex, tindex, windex, index, idx, iter, obs, col, error
		real matrix Y, U, V, Ws, Wm, R1, R2, na
		real colvector s
		string scalar name
		index = st_varindex(y)
		iindex = st_varindex(id)
		tindex = st_varindex(time)
		Y = J(N, T, .)
		R1 = J(N, T, 0)
		factorsfull = J(T, T, .)
		variance = J(T, T, .)
		variance2 = J(T, T, .)
		S = J(T, 1, .)

		if (strlen(w) > 0) {
			Ws = J(N, T, .)
			windex = st_varindex(w)
			for (obs = first; obs <= last ; obs++) {    
				Y[_st_data(obs, iindex), _st_data(obs, tindex)] = _st_data(obs, index)
				Ws[_st_data(obs, iindex), _st_data(obs, tindex)] =  _st_data(obs, windex)
			}
			Wm = rowsum(Ws :!=.)
			Wm= rowsum(editmissing(Ws, 0)):/ Wm
			Wm = sqrt(Wm)		
		}
		else{
			for (obs = first; obs <= last ; obs++) {  
				Y[_st_data(obs, iindex), _st_data(obs, tindex)] = _st_data(obs, index)  
			}
		}


		na = Y :==.
		Y = editmissing(Y, 0)
		error = 1
		iter = 0


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
			R = Y :+ na :* R1
			/* do PCA of residual */
			if (strlen(w) > 0) {
				variance = cross(R, Ws, R)
			}
			else{
				variance = cross(R, R)
			}
			_symeigensystem(variance, factorsfull, S)
			variance2 = cross(factorsfull[., 1::d]', factorsfull[., 1::d]')
			R2 = R * variance2
			error = max(abs(R2:-R1))
			R1 = R2
		}
		factors = factorsfull[., 1::d] :* sqrt(T)
		loadings = (R * factors) :/ T


		st_numscalar("r(N)", iter)
		st_numscalar("r(error)", error)

		if (strlen(residuals) > 0){
			residuals
			idx = st_addvar("double", residuals)
			for (obs = first; obs <= last ; obs++) {    
				st_store(obs, idx, Y[_st_data(obs, iindex), _st_data(obs, tindex)] - R2[_st_data(obs, iindex), _st_data(obs, tindex)])
			}	
		}


		if (strlen(id1gen) > 0){
			for (col = 1; col <= d; col++){
				name =  id1gen  + strofreal(col)
				idx = st_addvar("double", name)
				for (obs = first; obs <= last ; obs++) { 
					st_store(obs, idx, loadings[_st_data(obs, iindex), col])
				} 
			}
		}

		if (strlen(id2gen) > 0){
			for (col = 1; col <= d; col++){
				name =  id2gen + strofreal(col)
				idx = st_addvar("double", name)
				for (obs = first; obs <= last ; obs++) { 
					st_store(obs, idx, factors[_st_data(obs, tindex), col])
				} 
			}
		}


	}
end

