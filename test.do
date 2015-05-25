webuse set https://github.com/matthieugomez/stata-regife/raw/master/data/
webuse  Divorce-Wolfers-AER, clear
keep if inrange(year, 1968, 1988)
egen state = group(st), label
tsset state year
egen tag = tag(state)


regife div_rate unilateral divx*   [w=stpop], f(state year) d(2)
regife div_rate L.unilateral i.year   [aw=stpop], f(state year) d(2)

tab year, gen(tempyear)
regife div_rate L.unilateral tempyear*   [aw=stpop], f(state year) d(2)


reghdfe div_rate unilateral divx* [aw=stpop], a(state year)
regife div_rate unilateral divx* [aw=stpop], f(state year) d(2)
regife div_rate unilateral divx* [aw=stpop], f(state year) d(2) a(state year)
regife div_rate unilateral divx*  i.state i.year [aw=stpop], f(state year) d(2)




xtmg div_rate unilateral, cce

tempvar touse
gen byte `touse' = 1
markout `touse' div_rate unilateral
cap drop b
gen b = .
local vlist ""
foreach v of varlist div_rate unilateral{
	cap drop `v'_p
	egen `v'_p = wtmean(`v') if `touse', by(year)
	local yearlist `vlist' `v'_p

	egen `v'_s = wtmean(`v') if `touse', by(state)
	local statelist `vlist' `v'_s

}
qui levelsof state, clean
foreach s in `=r(levels)'{
	cap reg div_rate unilateral `yearlist'    if state == `s' 
	if _rc == 0{
		qui replace b = _b[unilateral] if e(sample) == 1 & _b[unilateral] ~= 0 
	}
}
drop `yearlist' `statelist'
ttest b =0 if tag == 1



regife div_rate unilateral i.year i.state  [aw=stpop], f(state year) d(5)
regife div_rate unilateral   [aw=stpop], f(state year) d(5) a(year state)



