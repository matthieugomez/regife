/***************************************************************************************************
v0.1 07/08/2015: first release
v0.2 07/09/2015: correct normalization for loadings
v0.3 04/12/2017: correct weight
v0.4 09/01/2021: remove error when N < T + preserve tsset
***************************************************************************************************/
program define regife, sortpreserve
	version 12.0
	syntax [varlist(min=1 numeric fv ts)] [if] [in] [aweight fweight pweight iweight] , [Factors(string) ife(string) vce(string) Absorb(string) RESiduals(string) * ]

	if "`factors'" != ""{
		di as  txt "The option factors() was renamed to ife(). In the future, please use the syntax ife(`factors') to specify the factor model."
		local ife `factors'
	}
	if "`ife'" == ""{
		di as error "option ife() required"
		exit 111
	}


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

	cap which reghdfe.ado
	if _rc {
		di as error "reghdfe.ado required: {stata ssc install reghdfe}"
		exit 111
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

			local absorbvars `absorbvars' `=regexs(2)'
		}
		else{
			local absorbvars `absorbvars' `1'
		}
		macro shift
	}

	/* syntax factors */
	if regexm("`ife'", "(.*),(.*)"){
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

	local id `id1'
	local time `id2'



	if ("`weight'"!=""){
		local wtype `weight'
		local wvar  `=subinstr("`exp'","=","", .)'
	}


	/* touse */
	marksample touse
	markout `touse' `id' `time' `wvar', strok

	/*syntax varlist */  
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
		syntax anything [, CLuster(varname) *]
		if "`anything'" == "bootstrap"{
			local bootstrap "yes"
			local bootstrapoptions `options'
		}
	}


	* mean weight
	if "`wvar'" ~= ""{
		tempvar wvar2
		qui bys `touse' `id': gen double `wvar2' = sum(`wvar') if `touse'
		qui by `touse' `id': replace `wvar2' = `wvar2'[_N]/_N if `touse'
		local wvar = "`wvar2'"
	}


	if "`bootstrap'" ~= "yes" {
		innerregife, dimension(`dimension') id(`id') time(`time') idgen(`id1gen') timegen(`id2gen') resgen(`residuals') y(`y') x(`x') yname(`yname') xname(`xname') touse(`touse') wtype(`wtype') wvar(`wvar')  absorbvars(`absorbvars') absorb(`absorb') vce(`vce')  `optionlist'
	}
	else{
		/* get bstart */
		qui innerregife, dimension(`dimension') id(`id') time(`time') idgen(`id1gen') timegen(`id2gen') resgen(`resgen') y(`y') x(`x') yname(`yname') xname(`xname') touse(`touse') wtype(`wtype') wvar(`wvar')  absorbvars(`absorbvars') absorb(`absorb')  `optionlist' 
		tempname bstart 
		matrix `bstart' = e(b)
		if "`cluster'" == ""{
			bootstrap,  `bootstrapoptions' : ///
			innerregife, dimension(`dimension') id(`id') time(`time') y(`y') x(`x') yname(`yname') xname(`xname') touse(`touse') wtype(`wtype') wvar(`wvar') absorb(`absorb') absorbvars(`absorbvars') fast bstart(`bstart') `optionlist'
		}
		else{
			tempvar clusterid
			cap tsset
			cap local tssettimevar = r(timevar)
			cap local tssetpanelvar = r(panelvar)
			tsset, clear
			if "`id'" == "`cluster'"{
				local absorbvars = substr("`absorbvars'", "`id'", "`clusterid'")
				bootstrap, cluster(`cluster') idcluster(`clusterid') `bootstrapoptions': ///
				innerregife, dimension(`dimension') id(`clusterid') time(`time') y(`y') x(`x') yname(`yname') xname(`xname') touse(`touse') wtype(`wtype') wvar(`wvar') absorb(`absorb') absorbvars(`absorbvars')  fast bstart(`bstart') `optionlist'
			}
			else if "`time'" == "`cluster'"{
				local absorb = substr("`absorb'", "`time'", "`clusterid'")
				bootstrap, cluster(`cluster') idcluster(`clusterid') `bootstrapoptions': ///
				innerregife, dimension(`dimension') id(`id') time(`clusterid') y(`y') x(`x') yname(`yname') xname(`xname') touse(`touse') wtype(`wtype') wvar(`wvar') absorb(`absorb') absorbvars(`absorbvars')  fast bstart(`bstart') `optionlist'
			}
			else{
				bootstrap, cluster(`cluster') idcluster(`clusterid') `bootstrapoptions': ///
				innerregife, dimension(`dimension') id(`id') time(`time') y(`y') x(`x') yname(`yname') xname(`xname') touse(`touse') wtype(`wtype') wvar(`wvar') absorb(`absorb') absorbvars(`absorbvars')  fast  bstart(`bstart')`optionlist'
			}
			if ("`tssettimevar'" != "") {
				tsset `tssetpanelvar' `tssettimevar' 
			}
		}
	}

end




