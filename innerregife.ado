
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


	* until reghdfe 3.0 is on ssc
	if "`fast'" == ""{
		local partial partial
	}


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
		tempvar cons
		qui gen `cons' = 1
		local xname `xname' _cons
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



	/* algorithim */
	* first reg  
	qui _regress `py' `px' `cons' `wt' in `touse_first'/`touse_last', nocons
	*qui ccemg `py' `pcx' `wt' if `touse', f(`id1' `id2')
	matrix `b' = e(b)
	qui predict `res' in `touse_first'/`touse_last'



	* iterate 
	mata: info = iteration("`py'","`res'", "`px' `cons'", "`wvar'", "`g1'", "`g2'", `N', `T', `dimension', `tolerance', `maxiterations', "`b'", `touse_first', `touse_last', "`id1gen'", "`tfactorlist'", "`verbose'")
	local iter = r(N)
	tempname convergence_error
	scalar `convergence_error' = r(convergence_error)
	tempname b
	matrix `b' = r(b)



	if `iter' == `maxiterations'{
		display as text "The algorithm did not converge : convergence error is" in ye %4.3gc `convergence_error' in text " (tolerance" in ye %4.3gc `tolerance' in text")"
		display as text "Use the maxiterations options to increase the amount of iterations"
	}

	tempvar esample
	gen `esample' = `touse'
	if "`fast'" == ""{
		/*  Errors are computed using the fact B_{Asy} has the same error distribution than the real b. beta_{ASI} can be estimated by reghdfe or by partialing dependent and regressor first. reghdfe makes it simpler to me and give same result*/ 

		if "`partial'" == ""{
			qui reghdfe `y' `x'  `cons' `wt'  in `touse_first'/`touse_last', a( `absorb' `id1'#c.(`tfactorlist')) `options'
			tempname df_a
			scalar `df_a' = e(df_a) + `T' * `dimension'
			tempname df_r
			scalar `df_r' = e(df_r) - `T' * `dimension'
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
			tempname N_hdfe
			scalar `N_hdfe' = e(N_hdfe)
			tempname N_hdfe_extended
			scalar `N_hdfe_extended' = e(N_hdfe_extended)
			tempname mobility
			scalar `mobility' = e(mobility)
			tempname M_due_to_nested
			scalar `M_due_to_nested' = e(M_due_to_nested)
			tempname df_a
			scalar `df_a' = e(df_a)
			tempname M1
			scalar `M1' = e(M1)
			tempname K1
			scalar `K1' = e(K1)
			tempname M2
			scalar `M2' = e(M2)
			tempname K2
			scalar `K2' = e(K2)
			tempname tss
			scalar `tss' = e(tss)
			tempname mss
			scalar `mss' = e(mss)
			tempname rmse
			scalar `rmse' = e(rmse)
			tempname F_absorb
			scalar `F_absorb' = e(F_absorb)
			tempname r2_a_within
			scalar `r2_a_within' = e(r2_a_within)
			tempname r2_a
			scalar `r2_a' = e(r2_a)
			tempname r2_within
			scalar `r2_within' = e(r2_within)
			tempname r2
			scalar `r2' = e(r2)
			tempname ll_0
			scalar `ll_0' = e(ll_0)
			tempname ll
			scalar `ll' = e(ll)
			tempname b V
			matrix `b' = e(b)
			matrix `V' = e(V) * e(df_r) / `df_r'
			mat colnames `b' =`xname'
			mat colnames `V' =`xname'
			mat rownames `V' =`xname'
			tempvar esample
			gen `esample' = `touse'

			ereturn post `b' `V' `wt', depname(`yname') obs(`obs') esample(`esample') dof(`=`df_r'')
			ereturn scalar df_m = `df_m'
			ereturn scalar N = `N'
			ereturn scalar rss = `rss'
			ereturn scalar F = `F'
			ereturn scalar rank = `rank'
			ereturn scalar N_hdfe = `N_hdfe'
			ereturn scalar N_hdfe_extended = `N_hdfe_extended'
			ereturn scalar mobility = `mobility' 
			ereturn scalar M_due_to_nested = `M_due_to_nested'
			ereturn scalar df_a = `df_a'
			ereturn scalar M1 = `M1'
			ereturn scalar K1 = `K1'
			ereturn scalar M2 = `M2'
			ereturn scalar K2 = `K2'
			ereturn scalar tss = `tss'
			ereturn scalar mss = `mss'
			ereturn scalar rmse = `rmse'
			ereturn scalar F_absorb = `F_absorb'
			ereturn scalar r2_a_within = `r2_a_within'
			ereturn scalar r2_a = `r2_a'
			ereturn scalar r2_within = `r2_within'
			ereturn scalar r2 = `r2'
			ereturn scalar ll_0 = `ll_0'
			ereturn scalar ll = `ll'
			ereturn scalar df_m = `df_m'
			ereturn scalar df_r = `df_r'
			ereturn local predict "regife_p"
		}
		else{
			* alternatively, I can residualize everything with respect to factors. Interestingly, I get same option
			local newlist
			foreach var in `py' `px' `cons'{
				tempvar new`var'
				mata: transform(*info[1],*info[2], *info[3], *info[4], "`wvar'", "`var'", "`new`var''", `N', `T', `touse_first', `touse_last')
				local newlist  `newlist' `new`var''
			}
			qui reg `newlist' `wt' in `touse_first'/`touse_last', nocons `options'
			scalar `df_a' =`df_a' + (`N'+`T')*`dimension'
			tempname df_r
			scalar `df_r' = e(df_r) - `df_a'
			tempname df_m
			scalar `df_m' = e(df_m)
			tempname r2
			scalar `r2' = e(r2)
			tempname b V
			matrix `b' = e(b)
			matrix `V' = e(V)* e(df_r) / `df_r'
			mat colnames `b' =`xname'
			mat colnames `V' =`xname'
			mat rownames `V' =`xname'
			tempvar esample
			gen `esample' = `touse'

			ereturn post `b' `V', depname(`yname') obs(`obs') esample(`esample') dof(`=`df_r'')
			ereturn scalar r2 = `r2'
			ereturn scalar df_a = `df_m'
			ereturn scalar df_m = `df_m'			
		}
	}
	else{
		/* don't compute error, just returns estimate */
		tempname V
		local  p `: word count `xname''
		matrix `V' = J(`p', `p', 0)
		mat colnames `b' =`xname' 
		mat colnames `V' =`xname'
		mat rownames `V' =`xname'
		tempname df_r
		scalar `df_r' = `obs' - `p'
		ereturn post `b' `V', depname(`yname') obs(`obs') esample(`esample') dof(`=`df_r'')
	}
	ereturn scalar iter = `iter'
	ereturn scalar convergence_error = `convergence_error'
	ereturn local id1 `id1'
	ereturn local id2 `id2'
	ereturn local f1 `id1gen'
	ereturn local f2 `id2gen'
	ereturn local d `dimension'
	ereturn local title REGIFE  
	ereturn local title2 Panel structure: `id1', `id2'
	ereturn local title3 Factor dimension: `dimension' 

	Header
	ereturn display

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
			Ws = J(N, 1, 1)
			for (obs = 1; obs <= last - first + 1; obs++) {    
				Ws[|index[obs, 1]|] = W[obs]
			}
			Ws = sqrt(Ws)
		}
		else{
			Ws = J(N, 1, 1)
			M = invsym(cross(X, X)) * X'
		}
		R = J(N, T, 0)
		V = J(T, T, .)
		iter = 0
		error = 1
		while ((maxiterations == 0) | (iter < maxiterations)){
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
				R = R :* Ws
			}
			/* do PCA of residual */
			_svd(R, s, V)
			if (strlen(w) > 0) {
				R = R :/ Ws
			}
			U = R[.,(1::d)] * diag(s[1::d]) 
			R = U *  V[(1::d),.] 

			for (obs = 1; obs <= last -first + 1; obs++) {    
				tY[obs] = R[|index[obs, .]|] 
			}
			/* estimate coefficient of (Y- PCA(RES)) on b */
			b2 = M * (Y :- tY)
			error = max(abs(b2 :- b1))
			b1 = b2
			if (error < tolerance){
				break
			}
		}
		
		V = V[1::d,.] :* sqrt(T)
		U = U :/ sqrt(T)
		MT = I(T) :- cross(V, V) / T
		MI = I(N) :- cross((U :* Ws)', (U:* Ws)') / (N*sum(Ws))

		st_store(first::last, res, tY)
		st_numscalar("r(N)", iter)
		st_numscalar("r(convergence_error)", error)
		st_matrix("r(b)",b1')

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
		return(&MT, &MI, &index, &Ws)
	}

	void transform(real matrix MT, real matrix MI, real matrix index, real matrix Ws, string scalar w, string scalar var, string scalar newvar, real scalar N, real scalar T, real scalar first, real scalar last){
		real matrix VARm
		real vector VAR
		real scalar obs, idx
		VAR = st_data((first::last), var)
		VARm = J(N, T, 0)
		for (obs = 1; obs <= last - first + 1; obs++) {  
			VARm[|(index)[obs, .]|]= VAR[obs]
		}
		if (strlen(w) > 0) {
			VARm =  ((MI * (VARm :* Ws) )* MT) :/Ws
		}
		else{
			VARm =  MI * VARm * MT 
		}

		for (obs = 1; obs <= last - first + 1; obs++) {  
			VAR[obs] = VARm[|index[obs, .]|]
		}
		idx = st_addvar("double", newvar)
		st_store(first::last, idx , VAR)
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
	if !missing(e(iter)) {
		.`right'.Arrpush `C3' "Number of iter" `C4' "= " as res %`c4wfmt'.0f e(iter)
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


