program define cce, eclass sortpreserve

	version 12
	syntax anything [if] [in] [aweight fweight pweight iweight], Factors(string) [Absorb(string) *]


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


	if "`ccemg'" ~= "" & "`absorb'" ~= ""{
		dis as error "The option absorb cannot be used for the ccemg estimate"
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
		qui reghdfe `y' `x' `wt' if `touse' == 1, a(`absorb' `id'##c.(`felist')) `options'
		ereturn local id1 `id'
		ereturn local id2 `time'
		ereturn local predict "regife_p"

		tempname b V
		mat `b' = e(b)
		mat colnames `b' =`xname'
		mat `V' = e(V)
		mat colnames `V' =`xname'
		mat rownames `V' =`xname'
		ereturn post `b' `V', depname(`yname') obs(`obs') esample(`touse') 
		ereturn display
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
		matrix rownames `V1' = `xname'
		matrix colnames `V1' = `xname'
		ereturn post `b1' `V1', depname(`yname') obs(`obs') esample(`touse') 
		ereturn display
		display as text "{lalign 26:Number of obs = }" in ye %10.0fc `obs'
		local i = 0
		display as text "{lalign 26:Number of groups = }" in ye %10.0fc `iter'
		foreach x in `xname'{
			local ++i
			display as text "{lalign 26:Number for `x' = }" in ye %10.0fc `ng'[1, `i']
			ereturn local `x'_n = `ng'[1, `i']
		}
	}

end

/***************************************************************************************************
helper functions
***************************************************************************************************/

set matastrict on

mata:
	void meanvar(real matrix b, real matrix w, string scalar sb1, string scalar sV1, string scalar sng){
		b1 = J(1, cols(b), .)
		V1 =  J(1, cols(b), .)
		ng = J(1, cols(b), .)

		for(i=1;i<= cols(b);++i){
			v = rowmissing(b[., i]):==0
			bnm =  select(b[., i], v)
			wnm = select(w, v)
			n = sum(v)
			ng[1,i] = n
			if (n == 1){
				b1[1, i] = bnm
			}
			else if (n > 1){
				b1[1, i] = mean(bnm, wnm)
				V1[1, i] = variance(bnm, wnm)/n
			}
		}
		st_matrix(sb1, editmissing(b1,0))
		st_matrix(sV1, editmissing(diag(V1),0))
		st_matrix(sng,  ng)
	}
end




