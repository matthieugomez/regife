use data/Divorce-Wolfers-AER, clear
keep if inrange(year, 1968, 1988)
egen state = group(st), label
tsset state year


regife div_rate unilateral divx*   [w=stpop], f(state year) d(2)
reghdfe div_rate unilateral divx* [aw=stpop], a(state year)
regife div_rate unilateral divx* [aw=stpop], f(state year) d(2)
regife div_rate unilateral divx* [aw=stpop], f(state year) d(2) a(state year)
regife div_rate unilateral divx*  i.state i.year [aw=stpop], f(state year) d(2)






* check lagged and factor display correctly
regife div_rate L.unilateral i.year   [aw=stpop], f(state year) d(2)
tab year, gen(tempyear)
regife div_rate L.unilateral tempyear*   [aw=stpop], f(state year) d(2)



* check factors with time dimension gives the original data
distinct year
local nyear = r(ndistinct)
ife div_rate, f(state year) d(`nyear') gen(new) maxi(1)
assert float(div_rate) == float(new)
drop new



* compare different estimators
tsset state year
xtmg  div_rate unilateral, cce
cce mg div_rate unilateral, f(state year)
* add year trend
cce mg div_rate unilateral year, f(state year)

cce p div_rate unilateral, f(state year) vce(cluster state)
regife div_rate unilateral, f(state year) a(state year) d(2)





/* weight */
gen weight2 = mod(state, 4) + 1
expand weight2 
egen t = tag(state year)
replace div_rate = 0 if missing(div_rate)
regife div_rate unilateral   [pw=weight2] if t , f(state year) d(2)
regife div_rate unilateral , f(state year) d(2)
