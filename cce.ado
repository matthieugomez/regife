program define cce, eclass sortpreserve

	version 12
	syntax anything [if] [in] [aweight fweight pweight iweight], Factors(string) [vce(string) ]


	/* syntax */
	if ("`weight'"!=""){
		local wt [`weight'`exp']
		local wtv = subinstr("`exp'","=","", .)
		local wtv2 "*`wtv'"
	}
	if regexm("`anything'", "^ *p (.*)"){
		local ccep ccep
		local varlist = regexs(1)
	}
	else if regexm("`anything'", "^ *mg (.*)"){
		local ccemg ccemg
		local varlist = regexs(1)
	}



	if "`ccemg'" ~= "" & "`weight'" == "fweight"{
		dis as error "fweight cannot be used for the ccemg estimate"
		exit 0 
	}
	if "`ccemg'" ~= "" & "`vce'" ~= ""{
		dis as error "invalid option: vce"
		exit 0 
	}
	

	/* syntax factors */
	local factors = trim("`factors'")
	if regexm("`factors'", "(.+) (.+)"){
		local id = regexs(1)
		local time = regexs(2)
	}


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



	qui count if `touse' 
	local touse_first = _N - r(N) + 1
	local touse_last = _N
	local obs = `touse_last'-`touse_first' + 1


	/* create weighted mean by time */
	sort `touse' `time'
	foreach v of varlist `y' `x' {
		tempvar `v'_p
		qui by `touse' `time': gen ``v'_p' = sum((`v')`wtv2')/sum(((`v')!=.)`wtv2') if `touse'
		qui by `touse' `time': replace ``v'_p' = ``v'_p'[_N]
		local felist `felist' ``v'_p'
	}
	if "`ccep'" ~= ""{
		if "`vce'"~=""{
			local vce vce(`vce')
		}
		qui reghdfe `y' `x' `wt' if `touse' == 1, a(`id1' `id2' `id'##c.(`felist')) `vce'
		tempname df_a
		scalar `df_a' = e(df_a) 
		tempname df_r
		scalar `df_r' = e(df_r)
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

		
		tempname b V
		mat `b' = e(b)
		mat colnames `b' =`xname'
		mat `V' = e(V)
		mat colnames `V' =`xname'
		mat rownames `V' =`xname'
		ereturn post `b' `V', depname(`yname') obs(`obs') esample(`touse') 
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
	else if "`ccemg'" ~= ""{
		* remove variables colinear for all groups (this allows to get some means that are at least meaningul)


		_rmcoll `x' `felist' `wt' in `touse_first'/`touse_last',  forcedrop 
		local vlist = r(varlist)
		local x: list vlist - felist
		local felist: list vlist - x
		* add fixed effect interacted with
		tempvar bylength
		sort `touse' `id'
		local type = cond(c(N)>c(maxlong), "double", "long")
		qui bys `touse' `id' : gen `type' `bylength' = _N 
		tempvar t
		qui by `touse' `id': gen `t' = _n == 1 if `touse'
		qui replace `t' = sum(`t')
		local N = `t'[_N]
		tempname b
		local p `: word count `x''
		mata: `b' = J(`N', `p', .)
		local iter = 0
		tempname w b1 V1 ng
		if "`wt'" ~= ""{
			mata: `w' = J(`N', 1, 0)
		}
		else{
			mata: `w' = J(`N', 1, 1)
		}
		tempvar cons



		sort `touse' `id'
		local start = `touse_first'
		while `start' <= `touse_last'{
			local base2 ""
			local iter = `iter' + 1
			local end = `start' + `bylength'[`start'] - 1
			if "`wt'" ~= ""{
				qui sum  `wt' in `start'/`end', meanonly
				mata: `w'[`iter', 1]  = st_numscalar("r(sum)")
			}

			/*  flag coefficients that are meaningless (sensible to which variable is removed) */
			local all `x' `felist'
			qui _rmcoll `all'  `wt' in `start'/`end', forcedrop 
			local p `: word count `x''
			if `=r(k_omitted)' {
				* I want to remove all variables that are jointly colinear.
				local base = r(varlist)
				local omitted: list all - base
				qui _rmcoll `omitted' `wt' in `start'/`end', forcedrop 
				local omitted = r(varlist)
				if "`omitted'" ~= "."{
					foreach v in `base'{
						qui _rmcoll `v' `omitted' `wt' in `start'/`end', forcedrop 
						if `=r(k_omitted)' == 0 {
							local base2 "`base2' `v' "
						}
					}
					cap _rmdcoll `y' `base'  `wt' in `start'/`end'
					if _rc{
						qui _regress `y' `x' `felist'  `wt' in `start'/`end'
						matrix `b1' = get(_b)
						matrix `V1'= vecdiag(get(VCE))
						local col = 0
						foreach v in `x'{
							local col = `col' + 1
							if strpos("`base2'", " `v' "){
								mata: `b'[`iter', `col']  = st_matrix("`b1'")[1,`col']
							}
						}
					}
				}
			}				
			else{
				cap _rmdcoll `y' `x' `felist'  `wt' in `start'/`end'
				if _rc == 0{
					qui _regress `y' `x' `felist'  `wt' in `start'/`end'
					matrix `b1' = get(_b)
					matrix `V1'= vecdiag(get(VCE))
					mata: `b'[`iter', .]  = st_matrix("`b1'")[1, 1::`p']
				}
			}
			local start = `end' + 1
		}
		mata: meanvar(`b', `w', "`b1'", "`V1'", "`ng'")
		matrix colnames `b1' = `xname'
		matrix colnames `ng' = `xname'
		matrix rownames `V1' = `xname'
		matrix colnames `V1' = `xname'
		ereturn post `b1' `V1', depname(`yname') obs(`obs') esample(`touse') dof(1)
		ereturn matrix Neff = `ng'

	}
	ereturn local id1 `id'
	ereturn local id2 `time'
	ereturn local predict "regife_p"
	Header
	ereturn display

end

/***************************************************************************************************
helper functions
***************************************************************************************************/

set matastrict on

mata:
	void meanvar(real matrix b, real matrix w, string scalar sb1, string scalar sV1, string scalar sng){
		b1 = J(1, cols(b), .)
		V1 =  J( cols(b), cols(b), .)
		ng = J(1, cols(b), .)
		for(i=1;i<= cols(b);++i){
			v = rowmissing(b[., i]):==0
			bnm =  select(b[., i], v)
			wnm = select(w, v)
			n = sum(v)
			wnm = wnm / sum(wnm) * n
			ng[1,i] = n
			if (n == 1){
				b1[1, i] = bnm
			}
			else if (n > 1){
				b1[1, i] = mean(bnm, wnm)
				V1[i, i] = variance(bnm, wnm)/(n-1)
			}
		}
		st_matrix(sb1, editmissing(b1,0))
		st_matrix(sV1, editmissing(V1,0))
		st_matrix(sng,  ng)
	}
end

/***************************************************************************************************

***************************************************************************************************/

/***************************************************************************************************
modified version from reghdfe.ado
***************************************************************************************************/

program define Header
	if !c(noisily) exit

	tempname left right
	.`left' = {}
	.`right' = {}

	local width 78
	local colwidths 1 30 40 67
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
	* effective

	cap local colnames: colnames e(Neff)
	if !_rc{
		local iter = 0
		tempname V
		matrix `V' = e(Neff)
		foreach name in `colnames'{
			local ++iter 
			local s = `V'[1,`iter']
			.`right'.Arrpush `C3' "Eff obs for `name'" `C4' "= " as res %`c4wfmt'.0f `s'
		}
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

