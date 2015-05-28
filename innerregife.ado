
/***************************************************************************************************

***************************************************************************************************/
program define innerregife, eclass 
	version 12
	syntax , Dimension(int) [ /// 
	id1(string) id2(string) id1gen(string) id2gen(string) /// 
	y(string) x(string) xname(string) yname(string) /// 
	touse(string)  wvar(string) wtype(string)  ///
	Absorb(string) noCONS TOLerance(real 1e-6) MAXIterations(int 10000) VERBose]


	/* tempname */
	tempvar res res2 y2 g1 g2
	tempname b V


	if ("`wtype'"!=""){
		local wt [`wtype'=`wvar']
		local wvar = "`wvar'"
		local sumwt [aw=`wvar']
	}


	/* absorb */
	tempname df_a
	if "`absorb'" ~= ""{
		cap which hdfe.ado
		if _rc {
			di as error "hdfe.ado required when using multiple absorb variables: {stata ssc install hdfe}"
			exit 111
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
		di as error "The factor structure dimension should be lower than the number of distinct values of the time variable"
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
	mata: iteration("`py'","`res'", "`pcx'", "`wvar'", "`g1'", "`g2'", `N', `T', `dimension', `tolerance', `maxiterations', "`b'", `touse_first', `touse_last', "`id1gen'", "`id2gen'", "`verbose'")
	local iter = r(N)
	tempname error
	scalar `error' = r(error)


	* last reg
	qui gen `res2' = `py' - `res'
	qui reg `res2' `px' `wt' if `touse', `cons' 

	/* return results */
	tempname df_m
	scalar `df_m' = e(df_m) + (`N'+`T')* `dimension' + `df_a'
	tempname df_r
	scalar `df_r' = `obs' - `df_m'

	mat `b' = e(b)
	mat `V' = e(V)
	mat colnames `b' =`xname'
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
		display as text "The algorithm did not converge : convergence error is" in ye %4.3gc `error' in text " (tolerance" in ye %4.3gc `tolerance' in text")"
		display as text "Use the maxiterations options to increase the amount of iterations"
	}

	ereturn post `b' `V', depname(`yname') obs(`obs') esample(`touse') 
	ereturn scalar df_r = `df_r'
	ereturn scalar df_m = `df_m'
	ereturn scalar r2_within = `r2w'
	ereturn scalar r2 = `r2'
	ereturn scalar convergence_error = `error'


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


	void iteration(string scalar y, string scalar res, string scalar x, string scalar w, string scalar id, string scalar time, real scalar N, real scalar T, real scalar d, real scalar tolerance, real scalar maxiterations, string scalar bname, real scalar first, real scalar last, string scalar id1gen, string scalar id2gen, string scalar verbose){
		real matrix Y 
		real matrix X
		real matrix tY
		real matrix M
		real matrix Ws
		real scalar iindex
		real scalar tindex
		real scalar windex


		real scalar iter
		real scalar obs
		real scalar col
		real scalar idx

		real scalar error
		real matrix U
		real matrix V
		real matrix W
		real matrix Ws
		real matrix Wm
		real colvector s
		real matrix R
		real colvector b1
		real colvector b2

		string scalar name

		iindex = st_varindex(id)
		tindex = st_varindex(time)
		windex = st_varindex(w)

		Y = st_data((first::last), y)
		b1 = st_matrix(bname)'
		X = st_data((first::last), x)

		if (strlen(w) > 0) {
			W = st_data((first::last), w)
			M = invsym(cross(X, W, X)) * X'* diag(W)
			/* define vector weight for each N (as sum of individual weight) */
			Ws = J(N, T, m)
			for (obs = first; obs <= last ; obs++) {    
				Ws[_st_data(obs, iindex), _st_data(obs, tindex)] = W[obs - first + 1, 1]
			}
			Wm = rowsum(Ws :==.)
			Wm= rowsum(editmissing(Ws, 0)) :/ Wm
			Wm =sqrt(Wm)
		}
		else{
			M = invsym(cross(X, X)) * X'
		}

		R = J(N, T, 0)
		V = J(T, T, .)
		iter = 0
		error = 1
		while (((maxiterations == 0) | (iter < maxiterations)) & (error >= tolerance)){
			iter = iter + 1
			if (strlen(verbose) > 0){
				if ((mod(iter, 100)==0) & (iter > 0)){
					if (iter == 100){
						stata(`"display as text "each .=100 iterations""')
					}
					stata(`"display "." _c"')
				}
			}
			/* construct residual matrix */
			tY = Y :-  X * b1
			for (obs = first; obs <= last ; obs++) {    
				R[_st_data(obs, iindex), _st_data(obs, tindex)] = tY[obs - first + 1, 1]
			}
			if (strlen(w) > 0) {
				R = R :* Wm
			}
			/* do PCA of residual */
			_svd(R, s, V)
			U = R[.,(1::d)] * diag(s[1::d]) 

			if (strlen(w) > 0) {
				U = U :/ Wm
			}
			
			R = U *  V[(1::d),.]
			for (obs = first; obs <= last ; obs++) {    
				tY[obs - first + 1,1] = R[_st_data(obs, iindex),_st_data(obs, tindex)] 
			}	
			/* estimate coefficient of (Y- PCA(RES)) on b */
			b2 = M * (Y :- tY)
			error = max(abs(b2 :- b1))
			b1 = b2
		}

		/* store factors if asked. note normalization in Bai different from svd in stata that has VV" = I rather than VV'/T = I */
		st_store(first::last, res, tY)
		if (strlen(id1gen) > 0){
			for (col = 1; col <= d; col++){
				name =  id1gen +  "_" + strofreal(col)
				idx = st_addvar("float", name)
				U = U :/ sqrt(T)
				for (obs = first; obs <= last ; obs++) { 
					st_store(obs, 	idx , U[_st_data(obs, iindex), col])
				} 
			}
		}
		if (strlen(id2gen) > 0){
			for (col = 1; col <= d; col++){
				name =  id2gen + "_" + strofreal(col)
				V = V:* sqrt(T)
				idx = st_addvar("float", name)
				for (obs = first; obs <= last ; obs++) { 
					st_store(obs, idx , V[col, _st_data(obs, tindex)])
				} 
			}
		}
		st_numscalar("r(N)", iter)
		st_numscalar("r(error)", error)
	}
end
/***************************************************************************************************

***************************************************************************************************/