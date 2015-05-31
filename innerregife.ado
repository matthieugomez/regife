
/***************************************************************************************************

***************************************************************************************************/
program define innerregife, eclass 
	version 12
	syntax , Dimension(int) [ /// 
	id1(string) id2(string) id1gen(string) id2gen(string) /// 
	y(string) x(string) xname(string) yname(string) /// 
	touse(string)  wvar(string) wtype(string)  ///
	fast ///
	Absorb(string) TOLerance(real 1e-6) MAXIterations(int 10000) VERBose partial * ///
	]


	/* tempname */
	tempvar res res2 y2 g1 g2
	tempname b V


	if ("`wtype'"!=""){
		local wt [`wtype'=`wvar']
		local wvar = "`wvar'"
		local sumwt [aw=`wvar']
	}

	if "`fast'" == ""{
		forval d = 1/`dimension'{
			tempvar tfactor`d'
			local tfactorlist `tfactorlist' `tfactor`d''
		}
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
		local touse `sample'
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
	}
	else{
		local px `x'
		local py `y'
		scalar `df_a' = 0
	}
	/* count number of observations (after hdfe since it reads the absorb syntax) */

	if "`fast'" == ""{
		if "`id1'" ~= "`: char _dta[_IDpanel]'" | "`id2'" == "`: char _dta[_TStvar]'"{
			sort `touse' `id1' `id2'
			cap bys `touse' `id1' `id2' : assert _N == 1 if `touse'
			if _rc{
				di as error "repeated observations for `id2' within `id1'"
				exit 451
			}
		}
	}


	sort `touse' `id1'
	qui count if `touse' 
	local touse_first = _N - r(N) + 1
	local touse_last = _N
	local obs = `touse_last'-`touse_first' + 1


	/* create group for i and t */
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

	if "`absorb'" == ""{
		tempvar cons
		qui gen `cons' = 1
		local xname `xname' _cons
	}


	/* algorithim */
	* first reg  
	qui _regress `py' `px' `cons' `wt' in `touse_first'/`touse_last', nocons
	*qui ccemg `py' `pcx' `wt' if `touse', f(`id1' `id2')
	matrix `b' = e(b)
	qui predict `res' in `touse_first'/`touse_last'



	* iterate 
	mata: info = iteration("`py'","`res'", "`px' `cons'", "`wvar'", "`g1'", "`g2'", `N', `T', `dimension', `tolerance', `maxiterations', "`b'", `touse_first', `touse_last', "`id1gen'", "`tfactorlist'", "`verbose'")
	local iter = r(N)
	tempname error
	scalar `error' = r(error)



	if `iter' == `maxiterations'{
		display as text "The algorithm did not converge : convergence error is" in ye %4.3gc `error' in text " (tolerance" in ye %4.3gc `tolerance' in text")"
		display as text "Use the maxiterations options to increase the amount of iterations"
	}


	ereturn local id1 `id1'
	ereturn local id2 `id2'
	ereturn local f1 `id1gen'
	ereturn local f2 `id2gen'
	ereturn local d `dimension'

	if "`fast'" == ""{

		/*  Errors are computed using the fact B_{Asy} has the same error distribution than the real b. beta_{ASI} can be estimated by reghdfe or by partialing dependent and regressor first. reghdfe looks simpler to me and give same result*/ 

		if "`partial'" == ""{
			qui reghdfe `y' `x' `cons' `wt' in `touse_first'/`touse_last', a(`absorb' `id1'#c.(`tfactorlist')) `options'
			tempname df_a
			scalar `df_a' = e(df_a) + `T' * `dimension'
			tempname df_r
			scalar `df_r' = e(df_r) - `T' * `dimension'

			tempname b V
			matrix `b' = e(b)
			matrix `V' = e(V) * e(df_r) / `df_r'
			mat colnames `b' =`xname'
			mat colnames `V' =`xname'
			mat rownames `V' =`xname'
			tempvar esample
			gen `esample' = `touse'
			ereturn post `b' `V', depname(`yname') obs(`obs') esample(`esample') 

			ereturn scalar convergence_error = `error'
			tempname df_m
			scalar `df_m' = e(df_m)
			ereturn scalar df_m = `df_m'
			tempname N
			scalar `N' = e(N)
			ereturn scalar N = `N'
			tempname rss
			scalar `rss' = e(rss)
			ereturn scalar rss = `rss'
			tempname F
			scalar `F' = e(F)
			ereturn scalar F = `F'
			tempname rank
			scalar `rank' = e(rank)
			ereturn scalar rank = `rank'
			tempname N_hdfe
			scalar `N_hdfe' = e(N_hdfe)
			ereturn scalar N_hdfe = `N_hdfe'
			tempname N_hdfe_extended
			scalar `N_hdfe_extended' = e(N_hdfe_extended)
			ereturn scalar N_hdfe_extended = `N_hdfe_extended'
			tempname mobility
			scalar `mobility' = e(mobility)
			ereturn scalar mobility = `mobility' 
			tempname M_due_to_nested
			scalar `M_due_to_nested' = e(M_due_to_nested)
			ereturn scalar M_due_to_nested = `M_due_to_nested'
			tempname df_a
			scalar `df_a' = e(df_a)
			ereturn scalar df_a = `df_a'
			tempname M1
			scalar `M1' = e(M1)
			ereturn scalar M1 = `M1'
			tempname K1
			scalar `K1' = e(K1)
			ereturn scalar K1 = `K1'
			tempname M2
			scalar `M2' = e(M2)
			ereturn scalar M2 = `M2'
			tempname K2
			scalar `K2' = e(K2)
			ereturn scalar K2 = `K2'
			tempname tss
			scalar `tss' = e(tss)
			ereturn scalar tss = `tss'
			tempname mss
			scalar `mss' = e(mss)
			ereturn scalar mss = `mss'
			tempname rmse
			scalar `rmse' = e(rmse)
			ereturn scalar rmse = `rmse'
			tempname F_absorb
			scalar `F_absorb' = e(F_absorb)
			ereturn scalar F_absorb = `F_absorb'
			tempname r2_a_within
			scalar `r2_a_within' = e(r2_a_within)
			ereturn scalar r2_a_within = `r2_a_within'
			tempname r2_a
			scalar `r2_a' = e(r2_a)
			ereturn scalar r2_a = `r2_a'
			tempname r2_within
			scalar `r2_within' = e(r2_within)
			ereturn scalar r2_within = `r2_within'
			tempname r2
			scalar `r2' = e(r2)
			ereturn scalar r2 = `r2'
			tempname ll_0
			scalar `ll_0' = e(ll_0)
			ereturn scalar ll_0 = `ll_0'
			tempname ll
			scalar `ll' = e(ll)
			ereturn scalar ll = `ll'
			ereturn scalar df_m = `df_m'
			ereturn scalar df_r = `df_r'
		ereturn local predict "regife_p"



			ereturn display
			display as text "{lalign 26:Number of obs = }" in ye %10.0fc `obs'
			display as text "{lalign 26:R-squared  = }" in ye %10.3fc `r2'
			display as text "{lalign 26:Within R-sq  = }" in ye %10.3fc `r2_within'
			display as text "{lalign 26:Number of iterations = }" in ye %10.0fc `iter'
		}
		else{
			* alternatively, I can residualize everything with respect to factors. Interestingly, I get same option
			local newlist
			foreach var in `py' `px' `cons'{
				tempvar new`var'
				mata: transform(info, "`var'", "`new`var''", `N', `T', `touse_first', `touse_last')
				local newlist  `newlist' `new`var''
			}
			qui reg `newlist' `wt' in `touse_first'/`touse_last', nocons `options'
			scalar `df_a' =`df_a' + (`N'+`T')*`dimension'
			tempname df_r
			scalar `df_r' = e(df_r) - `df_a'
			tempname b V
			matrix `b' = e(b)
			matrix `V' = e(V)* e(df_r) / `df_r'
			mat colnames `b' =`xname'
			mat colnames `V' =`xname'
			mat rownames `V' =`xname'
			tempvar esample
			gen `esample' = `touse'
			ereturn post `b' `V', depname(`yname') obs(`obs') esample(`esample') 
			tempname df_m
			scalar `df_m' = e(df_m)
			tempname df_r
			scalar `df_r' = e(N) - `df_m'
			tempname r2
			scalar `r2' = e(r2)
			ereturn scalar r2 = `r2'
			ereturn scalar convergence_error = `error'
			ereturn display
			display as text "{lalign 26:Number of obs = }" in ye %10.0fc `obs'
			display as text "{lalign 26:R-squared  = }" in ye %10.3fc `r2'
			display as text "{lalign 26:Number of iterations = }" in ye %10.0fc `iter'
		}
	}
	else{
		/* does not give righ errors, even after degree of freedom because I don't residualize Y and X correctly */
		tempvar y2
		qui gen `y2' = `py' - `res'
		qui _regress `y2' `px' `cons' `wt'  in `touse_first'/`touse_last', nocons
		matrix `b' = e(b)
		matrix `V' = e(V)
		mat colnames `b' =`xname'
		mat colnames `V' =`xname'
		mat rownames `V' =`xname'
		tempvar esample
		gen `esample' = `touse'
		ereturn post `b' `V', depname(`yname') obs(`obs') esample(`esample') 
		ereturn display
	}
	if "`id2gen'"~=""{
		forval d = 1/`dimension'{
			gen `id2gen'`d' = `tfactor'`d'
		}
	}

end

/***************************************************************************************************
helper functions
***************************************************************************************************/
set matastrict on
mata:


	pointer iteration(string scalar y, string scalar res, string scalar x, string scalar w, string scalar id, string scalar time, real scalar N, real scalar T, real scalar d, real scalar tolerance, real scalar maxiterations, string scalar bname, real scalar first, real scalar last, string scalar id1gen, string scalar tfactorlist, string scalar verbose){

		real matrix Y , X, tY, M, Ws, U, V, R, W, Wm
		real scalar iindex, tindex, windex, iter, obs, col, idx, error
		string scalar name
		real colvector s, b1, b2


		index = st_data(first::last, (id,time))
		st_view(Y, (first::last), y)
		st_view(X, (first::last), x)
		b1 = st_matrix(bname)'



		if (strlen(w) > 0) {
			windex = st_varindex(w)
			st_view(W, (first::last), w)
			M = invsym(cross(X, W, X)) * X'* diag(W)
			/* define vector weight for each N (as sum of individual weight) */
			Ws = J(N, T, .)
			for (obs = 1; obs <= last - first + 1; obs++) {    
				Ws[|index[obs, .]|] = W[obs]
			}
			Wm = rowsum(Ws :!=.)
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
			for (obs = 1; obs <= last -first + 1; obs++) {    
				R[|index[obs, .]|]= tY[obs]
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

			for (obs = 1; obs <= last -first + 1; obs++) {    
				tY[obs] = R[|index[obs, .]|] 
			}
			
			/* estimate coefficient of (Y- PCA(RES)) on b */
			b2 = M * (Y :- tY)
			error = max(abs(b2 :- b1))
			b1 = b2
		}
		V = V[1::d,.] :* sqrt(T)
		U = U :/ sqrt(T)
		MT = I(T) - cross(V, V) / T
		MI = I(N) - cross(U', U') / N
		st_store(first::last, res, tY)
		st_numscalar("r(N)", iter)
		st_numscalar("r(error)", error)
		if (strlen(tfactorlist) > 0){
			name = tokens(tfactorlist) 
			for (col = 1; col <= d; col++){
				idx = st_addvar("float", name[col])
				for (obs = first; obs <= last ; obs++) { 
					st_store(obs, idx , V[col, index[obs-first+1 ,2]])
				} 
			}
		}

		if (strlen(id1gen) > 0){
			for (col = 1; col <= d; col++){
				name =  id1gen + strofreal(col)
				idx = st_addvar("float", name)
				for (obs = first; obs <= last ; obs++) { 
					st_store(obs, 	idx , U[_index[obs -first + 1,1], col])
				} 
			}
		}
		return(&MT, &MI, &index)
	}


	void transform(pointer vector info, string scalar var, string scalar newvar, real scalar N, real scalar T, real scalar first, real scalar last){
		pointer(matrix) scalar pT
		pointer(matrix) scalar pI
		pointer(matrix) scalar pindex
		pT = info[1]
		pI = info[2]
		pindex = info[3]
		VAR = st_data((first::last), var)
		VARm = J(N, T, 0)
		for (obs = 1; obs <= last - first + 1; obs++) {  
			VARm[|(*pindex)[obs, .]|]= VAR[obs]
		}
		VARm =  *pI * VARm * *pT
		cross(*pI, *pI)
		for (obs = 1; obs <= last - first + 1; obs++) {  
			VAR[obs] = VARm[|(*pindex)[obs, .]|]
		}
		VAR = VAR 
		idx = st_addvar("double", newvar)
		st_store(first::last, idx , VAR)
	}


end
/***************************************************************************************************
ancien code
/* 
st_view(ID, (first::last), id)	
st_view(TIME, (first::last), time)
info = panelsetup(TIME, 1)
existingid	
stata(statacmd)
A = U * invsym(cross(U, U))* N * U'
Z = J(last-first+1, cols(X), .)
st_view(X, first::last, xlist)
for (t = 1; t <= T; t++){
	panelsubview(existingid, ID, t, info)
	A = A[existingid, existingid]
	Xsub = panelsubmatrix(X, t, info)
	Xsub2 = Xsub - A * Xsub/N
	Xsub2
	Z[(info[t,1])::(info[t,2]), . ] = Xsub2
}

D= cross(Z, Z)
*/
***************************************************************************************************/