		
/***************************************************************************************************

***************************************************************************************************/
program define innerregife, eclass 
	version 12
	syntax , Dimension(int) [ /// 
	id(string) time(string) idgen(string) timegen(string) /// 
	resgen(string) ///
	y(string) x(string) xname(string) yname(string) /// 
	touse(string)  wvar(string) wtype(string)  ///
	fast ///
	Absorb(string) absorbvars(string) /// 
	bstart(string) ///
	TOLerance(real 1e-9) MAXIterations(int 10000) VERBose partial  ///
	NOConstant ///
	demean /// 
	vce(string) ///
	]

	/* tempname */
	tempvar g1 g2 res
	tempname b V


	if `tolerance' < 1e-13{
		display as error "Tolerance should be higher than 1e-13"
	}
	if "`vce'" ~= ""{
		local vceoption vce(`vce')
	}
	if ("`wtype'"!=""){
		local wt [`wtype'=`wvar']
		local wvar = "`wvar'"
		local sumwt [aw=`wvar']
	}

	/* create list of new id and time name */
	forval d = 1/`dimension'{
		tempvar idfactor`d'
		local idfactorlist `idfactorlist' `idfactor`d''
		tempvar timefactor`d'
		local timefactorlist `timefactorlist' `timefactor`d''
	}


	tempname df_a
	if "`absorb'" ~= ""{
		/* case of high dimensional fixed effects : demean variables*/
		cap which hdfe.ado
		if _rc {
			di as error "hdfe.ado required when using multiple absorb variables: {stata ssc install hdfe}"
			exit 111
		}
		tempvar sample
		tempname prefix
		cap qui hdfe `y' `x'   if `touse' `wt', a(`absorbvars') gen(`prefix') sample(`sample')
		if _rc ~= 0{
			display as error "internall call to hdfe failed (error code: `=_rc'). This may be due to implicit categorical variables / time series variables in the model (i.e. i.var or L.var)". 
			exit 0
		}
		scalar `df_a' = e(df_a)
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
		local xname2 `xname'
		local yhdfe `y'
		local xhdfe `x'
	}
	else{
		/* otherwise add constant to list of regressors */
		tempvar cons
		gen `cons' = 1
		scalar `df_a' = 0
		if "`noconstant'" == "" {
			local py `y'
			local px `cons' `x' 
			local yhdfe `y'
			local xhdfe `cons' `x' 
			local xnamehdfe  _cons `xname'
			local xnamefast _cons `xname'
		}
		else{
			local py `y'
			local px `x'
			local yhdfe `y'
			local xhdfe  `x' 
			local xnamehdfe `xname'
			local xnamefast `xname'		
		}
	}


	/* create group for i and t */
	qui bys `touse' `id' : gen double `g1' = _n == 1 if `touse'
	qui replace `g1' = sum(`g1') if `touse'
	local N = `g1'[_N]

	qui bys `touse' `time'  : gen double `g2' = _n == 1 if `touse'
	qui replace `g2' = sum(`g2') if `touse'
	local T = `g2'[_N]

	qui count if `touse' 
	local touse_first = _N - r(N) + 1
	local touse_last = _N
	local obs = `touse_last'-`touse_first' + 1


	/* some checks */
	cap assert min(`T', `N') >= `dimension'
	if _rc{
		di as error "The factor structure dimension should be lower than the number of distinct values of the time variable"
		exit 0
	}
	
	cap assert max(`N', `T') <_N
	if _rc{
		di as error "More levels of FE than observations!"
		exit 3498
	}


	* initialize b
	if "`bstart'" ==""{
		qui _regress `py' `px' `wt' in `touse_first'/`touse_last', nocons
		matrix `b' = e(b)
	}
	else{
		local b `bstart'
	}

	* iterate 
	mata: info = iteration_svd("`py'", "`px'", "`wvar'", "`g1'", "`g2'", `N', `T', `dimension', `tolerance', `maxiterations', "`b'", `touse_first', `touse_last', "`idfactorlist'", "`timefactorlist'", "`res'", "`verbose'")
	tempname bend
	matrix `bend' = r(b)
	local iter = r(N)
	tempname convergence_error
	scalar `convergence_error' = r(convergence_error)

	if `iter' == `maxiterations'{
		local converged false
	}
	else{
		local converged true
	}

	if "`converged'" == "false"{
		display as text "The algorithm did not converge : convergence error is" in ye %4.3gc `convergence_error' in text " (tolerance" in ye %4.3gc `tolerance' in text")"
		display as text "Allow for more iterations with the option maxiter"

	}

	tempvar esample
	gen `esample' = `touse'


	if "`fast'" == ""{
		/* return the result from reghdfe adding i.id#c.factors and i.time#c.loading */
		foreach factor in `idfactorlist'{
			local idfactors `idfactors'	i.`time'#c.`factor'
		}
		foreach factor in `timefactorlist'{
			local timefactors `timefactors'	i.`id'#c.`factor'
		}

		cap which reghdfe
		if _rc == 111{
			di as error  `"Fixed effects feature requires command {cmd:reghdfe}"'
		}
		qui cap reghdfe `yhdfe' `xhdfe' `wt'  if `touse',  a(`absorb' `idfactors' `timefactors')  tol(`tolerance') `vceoption' keepsingletons
		if _rc ~= 0{
			display as error "internall call to reghdfe failed (error code: `=_rc'). Returning the estimate without standard errors". 
			local fast "fast"
		}
		else{
			tempname df_r
			scalar `df_r' = e(df_r) 
			tempname b V
			matrix `b' = e(b)
			matrix `V' = e(V)	
			mat colnames `b' =`xnamehdfe'
			mat colnames `V' =`xnamehdfe'
			mat rownames `V' =`xnamehdfe'

			tempname df_m
			scalar `df_m' = e(df_m)
			tempname N
			scalar `N' = e(N)
			tempname rss
			scalar `rss' = e(rss)
			tempname F
			scalar `F' = e(F)
			tempname rank
			scalar `rank' = e(rank)
			tempname mss
			scalar `mss' = e(mss)
			tempname rmse
			scalar `rmse' = e(rmse)

			tempvar esample
			qui gen `esample' = `touse'
			ereturn post `b' `V' `wt', depname(`yname') obs(`obs') esample(`esample') dof(`=`df_r'')
			ereturn scalar N = `N'
			ereturn scalar df_r = `df_r'
			ereturn scalar df_m = `df_m'
			ereturn scalar rss = `rss'
			ereturn scalar F = `F'
			ereturn scalar mss = `mss'
			ereturn scalar rmse = `rmse'	
		}
	}
	else{
		/* return the estimate without standard errors */
		local  p `: word count `xnamefast''
		tempname b V
		matrix `b' = `bend'
		matrix `V' = J(`p', `p', 0)
		mat colnames `b' =`xnamefast' 
		mat colnames `V' =`xnamefast'
		mat rownames `V' =`xnamefast'
		tempname df_r
		scalar `df_r' = `obs' - `p'
		ereturn post `b' `V', depname(`yname') obs(`obs') esample(`esample') dof(`=`df_r'')
	}



	ereturn scalar iterations = `iter'
	ereturn scalar convergence_error = `convergence_error'

	ereturn local cmd regife
	ereturn local depvar `y'
	ereturn local indepvars `xname'
	ereturn local converged `converged'
	ereturn local id `id'
	ereturn local time `time'
	ereturn local dimension `dimension'
	ereturn local title REGIFE  
	ereturn local title2 Panel structure: `id', `time'
	ereturn local title3 Factor dimension: `dimension' 
	ereturn local title4 Converged: `converged'

	Header
	ereturn display

	/* save factors, loadings and residuals */
	if "`idgen'"~=""{
		forval d = 1/`dimension'{
			qui gen `idgen'`d' = `idfactor`d''
		}
	}

	if "`timegen'"~=""{
		forval d = 1/`dimension'{
			qui gen `timegen'`d' = `timefactor`d''
		}
	}

	if "`resgen'"~=""{
		qui gen `resgen' = `res'
	}

end

/***************************************************************************************************
helper functions
***************************************************************************************************/
set matastrict on
mata:

	void iteration_svd(string scalar y, string scalar x, string scalar w, string scalar id, string scalar time, real scalar N, real scalar T, real scalar d, real scalar tolerance, real scalar maxiterations, string scalar bname, real scalar first, real scalar last, string scalar idfactorlist, string scalar timefactorlist, string scalar resgen, string scalar verbose){
		real matrix Y , X, tY, M, Ws, U, V, R, W, Wm, factors, loadings
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
			M = invsym(quadcross(X, W, X))* X' * diag(W)		
			Ws = J(N, 1, 1)
			for (obs = 1; obs <= last - first + 1; obs++) {    
				Ws[|index[obs, 1]|] = W[obs]
			}
			Ws = sqrt(Ws)
		}
		else{
			Ws = J(N, 1, 1)
			M = invsym(quadcross(X, X)) * X'
		}
		R = J(N, T, 0)
		R2 = J(N, T, 0)
		factorsfull = J(T, T, .)
		variance = J(T, T, .)
		variance2 = J(T, T, .)
		S = J(T, 1, .)
		V = J(T, T, .)
		iter = 0
		error = 1
		while ((maxiterations == 0) | (iter < maxiterations)){
			R = R2
			iter = iter + 1
			if (strlen(verbose) > 0){
				if ((mod(iter, 100)==0) & (iter > 0)){
					stata(`"display "current estimate:"')
					b1
				}
			}
			/* construct residual matrix */
			tY = Y :-  X * b1
			for (obs = 1; obs <= last -first + 1; obs++) {    
				R[|index[obs, .]|]= tY[obs]
			}
			/* do PCA of residual */
			if (strlen(w) > 0) {
				variance = cross(R, Ws, R)
			}
			else{
				variance = cross(R, R)
			}
			_symeigensystem(variance, factorsfull, S)
			/* compute R2, a low rank appproximation of R */
			variance2 = cross(factorsfull[., 1::d]', factorsfull[., 1::d]')
			R2 = R * variance2
			for (obs = 1; obs <= last -first + 1; obs++) {    
				tY[obs] = R2[|index[obs, .]|] 
			}
			/* estimate coefficient of (Y- PCA(RES)) on b */
			b2 = M * (Y :- tY)
			error = max(abs(b2 :-b1))
			b1 = b2
			if (error < tolerance){
				break
			}
		}
		factors = factorsfull[., 1::d] :* sqrt(T)
		loadings = (R * factors) :/ T
		MT = I(T) :- cross(factors', factors') / T 
		MI = I(N) :- ((loadings :* Ws) * invsym(cross((loadings :* Ws), (loadings:* Ws))) * (loadings :* Ws)')

		st_numscalar("r(N)", iter)
		st_numscalar("r(convergence_error)", error)
		st_matrix("r(b)",b1')
		names = tokens(idfactorlist)
		for (col = 1; col <= d; col++){
			idx = st_addvar("double", names[col])
			for (obs = first; obs <= last ; obs++) { 
				st_store(obs, idx , loadings[index[obs - first + 1, 1], col])
			} 
		}
		names = tokens(timefactorlist)
		for (col = 1; col <= d; col++){
			idx = st_addvar("double", names[col])
			for (obs = first; obs <= last ; obs++) { 
				st_store(obs, idx , factors[index[obs - first + 1 ,2], col])
			} 
		}
		res = Y :- tY :- X * b2
		idx = st_addvar("double", resgen)
		for (obs = first; obs <= last ; obs++) { 
			st_store(obs, idx , res[obs-first + 1])
		}
	}



	void iteration_gs(string scalar y, string scalar x, string scalar w, string scalar id, string scalar time, real scalar N, real scalar T, real scalar d, real scalar tolerance, real scalar maxiterations, string scalar bname, real scalar first, real scalar last, string scalar idfactorlist, string scalar timefactorlist, string scalar resgen, string scalar verbose){
		real matrix Y , X, tY, M, Ws, U, V, R, W, Wm, factors, loadings
		real scalar iindex, tindex, windex, iter, obs, col, idx, error
		string scalar name
		real colvector s, b1, b2
		idindex = st_data(first::last, id)
		timeindex = st_data(first::last, time)
		st_view(Y, (first::last), y)
		st_view(X, (first::last), x)
		b1 = st_matrix(bname)'
		if (strlen(w) > 0) {
			windex = st_varindex(w)
			st_view(W, (first::last), w)
			M = invsym(quadcross(X, W, X))* X' * diag(W)		
			Ws = J(N, 1, 1)
			for (obs = 1; obs <= last - first + 1; obs++) {    
				Ws[|index[obs, 1]|] = W[obs]
			}
			Ws = sqrt(Ws)
		}
		else{
			Ws = J(N, 1, 1)
			M = invsym(quadcross(X, X)) * X'
		}

		factors = J(T, d, 0.1)
		loadings = J(N, d, 0.1)
		idstorage1 = J(N, 1, 0)
		idstorage2 = J(N, 1, 0)
		timestorage1 = J(T, 1, 0)
		timestorage2 = J(T, 1, 0)
		U = J(d, d, .)
		Dx = J(d, 1, .)
		invDx = J(d, 1, .)

		V = J(d, d, .)
		iter = 0
		error = 1
		initial_iter = 30
		while ((maxiterations == 0) | (iter < maxiterations)){
			/* construct residuals */
			iter = iter + 1
			/* construct residuals  */
			tY = Y :-  X * b1
			for (r = 1 ; r <= d; r++) {
				for (inner_iter  = 1 ; inner_iter <= initial_iter; inner_iter++){
					for (obs = 1; obs <= last - first + 1; obs++) { 
						idi = idindex[obs]
						timei = timeindex[obs]   
						factor = factors[timei, r]
						idstorage1[idi]= idstorage1[idi] + tY[obs] * factor
						idstorage2[idi]= idstorage2[idi] + factor^2
					}
					for (i = 1 ; i <= N; i++){
						loadings[i, r] = idstorage1[i] / idstorage2[i]
						idstorage1[i] = 0
						idstorage2[i] = 0
					}
					for (obs = 1; obs <= last - first + 1; obs++) {    
						idi = idindex[obs]
						timei = timeindex[obs]   
						loading = loadings[idi, r]
						timestorage1[timei]= timestorage1[timei] + tY[obs] * loading
						timestorage2[timei]= timestorage2[timei] + loading^2
					}
					for (i = 1 ; i <= T; i++){
						factors[i, r] = timestorage1[i] / timestorage2[i]
						timestorage1[i] = 0
						timestorage2[i] = 0
					}
				}
				for (obs = 1; obs <= last - first + 1 ; obs++) {    
					idi = idindex[obs]
					timei = timeindex[obs]   
					tY[i] = tY[i] - loadings[idi, r] * factors[timei, r]
				}
			}
			initial_iter = 1
			for (obs = 1; obs <= last - first + 1  ; obs++) {  
				idi = idindex[obs]
				timei = timeindex[obs]   
				sum = 0
				for (r= 1 ; r <= d; r++) {
					sum = sum + loadings[idi, r] * factors[timei, r]
				}
				tY[obs] = Y[obs] - sum
			}
			/* estimate coefficient of (Y- PCA(RES)) on b */
			b2 = M * tY
			error = max(abs(b2 :-b1))
			b1 = b2
			if (error < tolerance){
				break
			}
		}

		/* scale factors and loadings */
		_symeigensystem(cross(factors, factors), U, Dx)

		for (i = 1; i <= d; i++){
			Dx[i] = sqrt(abs(Dx[i]))
		}
		for (i = 1; i <= d; i++){
			invDx[i] = 1/Dx[i]
		}
		scaledloadings = loadings * U * diag(Dx)
		_symeigensystem(cross(scaledloadings, scaledloadings), V, Dx2)
		loadings = loadings * U * diag(Dx) * V
		factors = factors *  U * diag(invDx) * V



		st_numscalar("r(N)", iter)
		st_numscalar("r(convergence_error)", error)
		st_matrix("r(b)",b1')
		names = tokens(idfactorlist)
		for (r = 1; r <= d; r++){
			idx = st_addvar("double", names[r])
			for (obs = 1; obs <= last - first + 1 ; obs++) { 
				st_store(obs, idx , loadings[idindex[obs - first + 1], r])
			} 
		}
		names = tokens(timefactorlist)
		for (r = 1; r <= d; r++){
			idx = st_addvar("double", names[r])
			for (obs = first; obs <= last ; obs++) { 
				st_store(obs, idx , factors[timeindex[obs - first + 1], r])
			} 
		}
		res = Y :- tY :- X * b2
		idx = st_addvar("double", resgen)
		for (obs = first; obs <= last ; obs++) { 
			st_store(obs, idx , res[obs-first + 1])
		}
	}



end
/***************************************************************************************************
slightly modified version from reghdfe.ado
***************************************************************************************************/

program define Header
	if !c(noisily) exit

	tempname left right
	.`left' = {}
	.`right' = {}

	local width 78
	local colwidths 1 30 51 67
	local i 0
	foreach c of local colwidths {
		local ++i
		local c`i' `c'
		local C`i' _col(`c')
	}

	local c2wfmt 10
	local c4wfmt 10
	local max_len_title = `c3' - 2
	local c4wfmt1 = `c4wfmt' + 1
	local title  `"`e(title)'"'
	local title2 `"`e(title2)'"'
	local title3 `"`e(title3)'"'
	local title4 `"`e(title4)'"'
	local title5 `"`e(title5)'"'

	// Right hand header ************************************************

	*N obs
	.`right'.Arrpush `C3' "Number of obs" `C4' "= " as res %`c4wfmt'.0f e(N)


	* Ftest
	if `"`e(chi2)'"' != "" | "`e(df_r)'" == "" {
		Chi2test `right' `C3' `C4' `c4wfmt'
	}
	else {
		Ftest `right' `C3' `C4' `c4wfmt'
	}

	* display R-squared

	if !missing(e(r2_a)) {
		.`right'.Arrpush `C3' "Adj R-squared" `C4' "= " as res %`c4wfmt'.4f e(r2_a)
	}
	if !missing(e(r2_within)) {
		.`right'.Arrpush `C3' "Within R-sq." `C4' "= " as res %`c4wfmt'.4f e(r2_within)
	}
	if !missing(e(rmse)) {
		.`right'.Arrpush `C3' "Root MSE" `C4' "= " as res %`c4wfmt'.4f e(rmse)
	}
	if !missing(e(r2_p)) {
		.`right'.Arrpush `C3' "Pseudo R2" `C4' "= " as res %`c4wfmt'.4f e(r2_p)
	}

	* iteration
	if !missing(e(iterations)) {
		.`right'.Arrpush `C3' "Iterations" `C4' "= " as res %`c4wfmt'.0f e(iterations)
	}



	// Left hand header *************************************************

	* make title line part of the header if it fits
	local len_title : length local title
	forv i=2/5 {
		if (`"`title`i''"'!="") {
			local len_title = max(`len_title',`:length local title`i'')
		}
	}

	if `len_title' < `max_len_title' {
		.`left'.Arrpush `"`"`title'"'"'
		local title
		forv i=2/5 {
			if `"`title`i''"' != "" {
				.`left'.Arrpush `"`"`title`i''"'"'
				local title`i'
			}
		}
		.`left'.Arrpush "" // Empty
	}

	* Clusters
	local kr = `.`right'.arrnels' // number of elements in the right header
	local kl = `.`left'.arrnels' // number of elements in the left header
	local N_clustervars = e(N_clustervars)
	if (`N_clustervars'==.) local N_clustervars 0
	local space = `kr' - `kl' - `N_clustervars'
	local clustvar = e(clustvar)
	forv i=1/`space' {
		.`left'.Arrpush ""
	}
	forval i = 1/`N_clustervars' {
		gettoken cluster clustvar : clustvar
		local num = e(N_clust`i')
		.`left'.Arrpush `C1' "Number of clusters (" as res "`cluster'" as text  ") " `C2' as text "= " as res %`c2wfmt'.0f `num'
	}

	HeaderDisplay `left' `right' `"`title'"' `"`title2'"' `"`title3'"' `"`title4'"' `"`title5'"'
end

program define HeaderDisplay
	args left right title1 title2 title3 title4 title5

	local nl = `.`left'.arrnels'
	local nr = `.`right'.arrnels'
	local K = max(`nl',`nr')

	di
	if `"`title1'"' != "" {
		di as txt `"`title'"'
		forval i = 2/5 {
			if `"`title`i''"' != "" {
				di as txt `"`title`i''"'
			}
		}
		if `K' {
			di
		}
	}

	local c _c
	forval i = 1/`K' {
		di as txt `.`left'[`i']' as txt `.`right'[`i']'
	}
end

program define Ftest
	args right C3 C4 c4wfmt is_svy

	local df = e(df_r)
	if !missing(e(F)) {
		.`right'.Arrpush                                ///
		`C3' "F("                              ///
		as res %4.0f e(df_m)                         ///
		as txt ","                                   ///
		as res %7.0f `df' as txt ")" `C4' "= "       ///
		as res %`c4wfmt'.2f e(F)
		.`right'.Arrpush                                ///
		`C3' "Prob > F" `C4' "= "              ///
		as res %`c4wfmt'.4f Ftail(e(df_m),`df',e(F))
	}
	else {
		local dfm_l : di %4.0f e(df_m)
		local dfm_l2: di %7.0f `df'
		local j_robust "{help j_robustsingular##|_new:F(`dfm_l',`dfm_l2')}"
		.`right'.Arrpush                                ///
		`C3' "`j_robust'"                     ///
		as txt `C4' "= " as result %`c4wfmt's "."
		.`right'.Arrpush                                ///
		`C3' "Prob > F" `C4' "= " as res %`c4wfmt's "."
	}
end

program define Chi2test

	args right C3 C4 c4wfmt

	local type `e(chi2type)'
	if `"`type'"' == "" {
		local type Wald
	}
	if !missing(e(chi2)) {
		.`right'.Arrpush                                ///
		`C3' "`type' chi2("                   ///
		as res e(df_m)                               ///
		as txt ")" `C4' "= "                         ///
		as res %`c4wfmt'.2f e(chi2)
		.`right'.Arrpush                                ///
		`C3' "Prob > chi2" `C4' "= "          ///
		as res %`c4wfmt'.4f chi2tail(e(df_m),e(chi2))
	}
	else {
		local j_robust                                  ///
		"{help j_robustsingular##|_new:`type' chi2(`e(df_m)')}"
		.`right'.Arrpush                                ///
		`C3' "`j_robust'"                     ///
		as txt `C4' "= " as res %`c4wfmt's "."
		.`right'.Arrpush                                ///
		`C3' "Prob > chi2" `C4' "= "          ///
		as res %`c4wfmt's "."
	}
end



