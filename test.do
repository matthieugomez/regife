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




use data/income-deregulation, clear
regife p20  intra_dummy, f(state year) d(2)
cce p p20  intra_dummy, f(state year)
cce mg p20  intra_dummy, f(state year)


regife p20  L.unilateral i.year   [aw=stpop], f(state year) d(2)

tab year, gen(tempyear)
regife p20  L.unilateral tempyear*   [aw=stpop], f(state year) d(2)


reghdfe p20  unilateral divx* [aw=stpop], a(state year)
regife p20  unilateral divx* [aw=stpop], f(state year) d(2)
regife p20  unilateral divx* [aw=stpop], f(state year) d(2) a(state year)
regife p20  unilateral divx*  i.state i.year [aw=stpop], f(state year) d(2)







xtmg div_rate unilateral, cce


cce p div_rate unilateral divx* [w=stpop], f(state year)
cce mg div_rate unilateral divx* [w=stpop], f(state year) 


tempvar touse
gen byte `touse' = 1
markout `touse' div_rate unilateral
cap drop b
gen b = .
local vlist ""
foreach v of varlist div_rate unilateral{
	cap drop `v'_y
	egen `v'_y = wtmean(`v') if `touse', by(year)
	local yearlist `yearlist' `v'_y
}
qui levelsof state, clean
foreach s in `=r(levels)'{
	cap reg div_rate unilateral `yearlist'    if state == `s' 
	if _rc == 0{
		qui replace b = _b[unilateral] if e(sample) == 1 & _b[unilateral] ~= 0 
	}
}
ttest b =0 if tag == 1 

ds `yearlist'
reghdfe div_rate unilateral, a(state##c.(`yearlist'))
drop `yearlist' 


regife div_rate unilateral i.year i.state  [aw=stpop], f(state year) d(5)
regife div_rate unilateral   [aw=stpop], f(state year) d(5) a(year state)



