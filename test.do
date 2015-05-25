webuse set https://github.com/matthieugomez/stata-regife/raw/master/data/
webuse  Divorce-Wolfers-AER, clear
keep if inrange(year, 1968, 1988)
egen state = group(st), label
tsset state year


regife div_rate unilateral divx*   [w=stpop], f(state year) d(2)
regife div_rate L.unilateral i.year   [aw=stpop], f(state year) d(2)

tab year, gen(tempyear)
regife div_rate L.unilateral tempyear*   [aw=stpop], f(state year) d(2)


reghdfe div_rate unilateral divx* [aw=stpop], a(state year)
regife div_rate unilateral divx* [aw=stpop], f(state year) d(2)
regife div_rate unilateral divx* [aw=stpop], f(state year) d(2) a(state year)
regife div_rate unilateral divx*  i.state i.year [aw=stpop], f(state year) d(2)




cap drop b
gen b = .
local vlist ""
foreach v of varlist div_rate unilateral divx*{
	cap drop `v'_p
	egen `v'_p = wtmean(`v'), by(year)
	local vlist `vlist' `v'_p
}
levelsof state, clean
foreach s in `=r(levels)'{
	cap reg div_rate unilateral divx* `vlist'  if year>1967 & year<1989  & state == `s' [aw=stpop]
	if _rc == 0{
		replace b = _b[unilateral] if e(sample) == 1 
	}
}
drop `vlist'
