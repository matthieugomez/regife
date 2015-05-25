cap program drop regife
program define regife, eclass sortpreserve
	version 13
	syntax anything [if] [in] [aweight/], ife(varlist) Dimension(integer)  [Absorb(varlist) noCONS convergence(real 0.000001) MAXiteration(int 10000) gen(string)]


	qui ExpandFactorVariables `anything'
	local anything = r(varlist)
	tokenize `anything'
	local y `1'
	macro shift 
	local x `*'



	if "`gen'" ~= ""{
		confirm new variable `gen'
	}

	local ife = trim("`ife'")
	if regexm("`ife'", "(.+[^=]*) ([^=]*.+)"){
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
	if "`gen'" ~= ""{
		rename `res' `gen'
	}
	qui reg `res2' `xc', nocons



	/* results */
	tempname df_m
	scalar `df_m' = e(df_m) + (`N'+`T')* `dimension' + `df_a'
	tempname df_r
	scalar `df_r' = `obs' - `df_m'

	mat `b' = e(b)

	FixVarnames `oldx'
	local xnames =  r(newnames)
	mat colnames `b' =`xnames'
	mat `V' = e(V) * e(df_r)/ `df_r'
	mat colnames `V' =`xnames'
	mat rownames `V' =`xnames'

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

/***************************************************************************************************
helpfer program to handle factor variables from reghdfe
***************************************************************************************************/
cap program drop Assert
cap program drop ExpandFactorVariables
cap program drop LabelRenameVariable
cap program drop Debug
cap program drop FixVarnames

program define Assert
	syntax anything(everything equalok) [, MSG(string asis) RC(integer 198)]
	if !(`anything') {
		di as error `msg'
		exit `rc'
	}
end



program define Debug

	syntax, [MSG(string asis) Level(integer 1) NEWline COLOR(string)] [tic(integer 0) toc(integer 0)]
	
	cap mata: st_local("VERBOSE",strofreal(VERBOSE)) // Ugly hack to avoid using a global
	if ("`VERBOSE'"=="") {
		di as result "Mata scalar -VERBOSE- not found, setting VERBOSE=3"
		local VERBOSE 3
		mata: VERBOSE = `VERBOSE'
	}


	assert "`VERBOSE'"!=""
	assert inrange(`level',0, 4)
	assert (`tic'>0) + (`toc'>0)<=1

	if ("`color'"=="") local color text
	assert inlist("`color'", "text", "res", "result", "error", "input")

	if (`VERBOSE'>=`level') {

		if (`tic'>0) {
			timer clear `tic'
			timer on `tic'
		}
		if (`toc'>0) {
			timer off `toc'
			qui timer list `toc'
			local time = r(t`toc')
			if (`time'<10) local time = string(`time'*1000, "%tcss.ss!s")
			else if (`time'<60) local time = string(`time'*1000, "%tcss!s")
			else if (`time'<3600) local time = string(`time'*1000, "%tc+mm!m! SS!s")
			else if (`time'<24*3600) local time = string(`time'*1000, "%tc+hH!h! mm!m! SS!s")
			timer clear `toc'
			local time `" as result " `time'""'
		}

		if (`"`msg'"'!="") di as `color' `msg'`time'
		if ("`newline'"!="") di
	}
end


program define ExpandFactorVariables, rclass
	syntax varlist(min=1 numeric fv ts) [if] [,setname(string)] [CACHE]
	
	* If saving the data for later regressions -savecache(..)- we will need to match each expansion to its newvars
	* This mata array is used for that
	* Note: This explains why we need to wrap -fvrevar- in a loop
	if ("`cache'"!="") mata: varlist_cache = asarray_create()

	* Building the debug message may be slow, only do it if requested with verbose
	cap mata: st_local("VERBOSE",strofreal(VERBOSE))
	if ("`VERBOSE'"=="") local VERBOSE 3

	local expanded_msg `"" - variable expansion for `setname': " as result "`varlist'" as text " ->""'
	while (1) {
		gettoken factorvar varlist : varlist, bind
		if ("`factorvar'"=="") continue, break

		fvrevar `factorvar' `if' // , stub(__V__) // stub doesn't work in Stata 11.2
		local contents
		foreach var of varlist `r(varlist)' {
			LabelRenameVariable `var' // Tempvars not renamed will be dropped automatically
			if !r(is_dropped) local contents `contents' `r(varname)'
			* Yellow=Already existed, White=Created, Red=NotCreated (omitted or base)
			local color = cond(r(is_dropped), "error", cond(r(is_newvar), "input", "result"))
			if (`VERBOSE'>3) {
				local expanded_msg `"`expanded_msg' as `color' " `r(name)'" as text " (`r(varname)')""'
			}
		}
		Assert "`contents'"!="", msg("error: variable -`fvvar'- in varlist -`varlist'- in category -`setname'- is  empty after factor/time expansion")
		if ("`cache'"!="") mata: asarray(varlist_cache, "`fvvar'", "`contents'")
		local newvarlist `newvarlist' `contents'
	}

	Debug, level(3) msg(`expanded_msg')
	return clear
	return local varlist "`newvarlist'"
end

program define LabelRenameVariable, rclass
	syntax varname
	local var `varlist'
	local fvchar : char `var'[fvrevar]
	local tschar : char `var'[tsrevar]
	local is_newvar = ("`fvchar'`tschar'"!="") & substr("`var'", 1, 2)=="__"
	local name "`var'"
	local will_drop 0

	if (`is_newvar') {
		local name "`fvchar'`tschar'"
		local parts : subinstr local fvchar "#" " ", all
		local has_cont_interaction = strpos("`fvchar'", "c.")>0
		local is_omitted 0
		local is_base 0
		foreach part of local parts {
			if (regexm("`part'", "b.*\.")) local is_base 1
			if (regexm("`part'", "o.*\.")) local is_omitted 1
		}

		local will_drop = (`is_omitted') | (`is_base' & !`has_cont_interaction')
		if (!`will_drop') {
			char `var'[name] `name'
			la var `var' "[Tempvar] `name'"
			local newvar : subinstr local name "." "__", all
			local newvar : subinstr local newvar "#" "_X_", all
			* -permname- selects newname# if newname is taken (# is the first number available)
			local newvar : permname __`newvar', length(30)
			rename `var' `newvar'
			local var `newvar'
		}
	}

	return scalar is_newvar = `is_newvar'
	return scalar is_dropped = `will_drop'
	return local varname "`var'"
	return local name "`name'"
end


program define FixVarnames, rclass
	local vars `0'

	foreach var of local vars {
		local newname
		local pretyname

		* -var- can be <o.__W1__>
		if ("`var'"=="_cons") {
			local newname `var'
			local prettyname `var'
		}
		else {
			fvrevar `var', list
			local basevar "`r(varlist)'"
			local label : var label `basevar'
			local is_avge = regexm("`basevar'", "^__W[0-9]+__$")
			local is_temp = substr("`basevar'",1,2)=="__"
			local is_omitted = strpos("`var'", "o.")
			local prefix = cond(`is_omitted'>0, "o.", "")
			local name : char `basevar'[name]

			if (`is_avge') {
				local avge_str : char `basevar'[avge_equation]
				local name : char `basevar'[name]
				local prettyname `avge_str':`prefix'`name'

				local newname : char `basevar'[target]
				if ("`newname'"=="") local newname `var'
			}
			else if (`is_temp' & "`name'"!="") {
				local newname `prefix'`name'
				
				* Fix bug when the var is omitted:
				local bugmatch = regexm("`newname'", "^o\.([0-9]+)b?\.(.+)$")
				if (`bugmatch') {
					local newname = regexs(1) + "o." + regexs(2) // EG: 1o.var
				}

				local prettyname `newname'
			}
			else {
				local newname `var'
				local prettyname `newname'
			}
		}
		
		*di in red " var=<`var'> --> new=<`newname'> pretty=<`prettyname'>"
		Assert ("`newname'"!="" & "`prettyname'"!=""), ///
		msg("var=<`var'> --> new=<`newname'> pretty=<`prettyname'>")
		local newnames `newnames' `newname'
		local prettynames `prettynames' `prettyname'
	}

	local A : word count `vars'
	local B : word count `newnames'
	local C : word count `prettynames'
	Assert `A'==`B', msg("`A' vars but `B' newnames")
	Assert `A'==`C', msg("`A' vars but `C' newnames")
	
	***di as error "newnames=`newnames'"
	***di as error "prettynames=`prettynames'"

	return local newnames "`newnames'"
	return local prettynames "`prettynames'"
end

