/***************************************************************************************************

***************************************************************************************************/
program define regife, sortpreserve
	version 12
	syntax [varlist(min=1 numeric fv ts)] [if] [in] [aweight fweight pweight iweight] , Factors(string) [vce(string) Absorb(string) RESiduals(string) * ]



	local optionlist `options'

	/* syntax absorb */
	while (regexm("`absorb'", "[ ][ ]+")) {
		local absorb : subinstr local absorb "  " " ", all
	}
	local absorb : subinstr local absorb " =" "=", all
	local absorb : subinstr local absorb "= " "=", all
	

	if "`residuals'" ~= ""{
		confirm new variable `residuals'
	}


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

			local absorbgroup `absorbgroup' `=regexs(2)'
		}
		else{
			local absorbgroup `absorbgroup' `1'
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
				confirm new variable `id`i'gen'`d'
			}
		}
		else{
			local id`i' `f'
		}
		confirm var `id`i''
	}

	if ("`weight'"!=""){
		local wtype `weight'
		local wvar  `=subinstr("`exp'","=","", .)'
	}


	/* touse */
	marksample touse
	markout `touse' `id1' `id2' `wvar', strok

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


	* syntax errors
	if "`vce'" ~= ""{
		local 0 `vce'
		syntax anything [, reps(int 50) cluster(varname)]
		if "`anything'" == "bootstrap"{
			local bootstrap "yes"
		}
	}


	* mean weight
	if "`wvar'" ~= ""{
		tempvar wvar2
		qui bys `touse' `id1': gen double `wvar2' = sum(`wvar') if `touse'
		qui by `touse' `id1': replace `wvar2' = `wvar2'[_N]/_N if `touse'
		local wvar = "`wvar2'"
	}


	if "`bootstrap'" ~= "yes" {
		innerregife, dimension(`dimension') id1(`id1') id2(`id2') id1gen(`id1gen') id2gen(`id2gen') resgen(`residuals') y(`y') x(`x') yname(`yname') xname(`xname') touse(`touse') wtype(`wtype') wvar(`wvar')  absorbgroup(`absorbgroup') absorb(`absorb') vce(`vce')  `optionlist'
	}
	else{
		/* get bstart */
		qui innerregife, dimension(`dimension') id1(`id1') id2(`id2') id1gen(`id1gen') id2gen(`id2gen') resgen(`resgen') y(`y') x(`x') yname(`yname') xname(`xname') touse(`touse') wtype(`wtype') wvar(`wvar')  absorbgroup(`absorbgroup') absorb(`absorb')  `optionlist' 
		tempname bstart 
		matrix `bstart' = e(b)
		if "`cluster'" == ""{
			bootstrap,  reps(`reps') : ///
			innerregife, dimension(`dimension') id1(`id1') id2(`id2') y(`y') x(`x') yname(`yname') xname(`xname') touse(`touse') wtype(`wtype') wvar(`wvar') absorbgroup(`absorbgroup') absorb(`absorb') fast bstart(`bstart') `optionlist'
		}
		else{
			tempvar clusterid
			tsset, clear
			if "`id1'" == "`cluster'"{
				bootstrap, reps(`reps') cluster(`cluster') idcluster(`clusterid'): ///
				innerregife, dimension(`dimension') id1(`clusterid') id2(`id2') y(`y') x(`x') yname(`yname') xname(`xname') touse(`touse') wtype(`wtype') wvar(`wvar') absorbgroup(`absorbgroup') absorb(`absorb') fast bstart(`bstart') `optionlist'
			}
			else if "`id2'" == "`cluster'"{
				bootstrap, reps(`reps') cluster(`cluster') idcluster(`clusterid'): ///
				innerregife, dimension(`dimension') id1(`id1') id2(`clusterid') y(`y') x(`x') yname(`yname') xname(`xname') touse(`touse') wtype(`wtype') wvar(`wvar') absorbgroup(`absorbgroup')  fast `optionlist'
			}
			else{
				bootstrap, reps(`reps') cluster(`cluster') idcluster(`clusterid'): ///
				innerregife, dimension(`dimension') id1(`id1') id2(`id2') y(`y') x(`x') yname(`yname') xname(`xname') touse(`touse') wtype(`wtype') wvar(`wvar') absorbgroup(`absorbgroup')  fast  bstart(`bstart')`optionlist'
			}

		}
	}

end


