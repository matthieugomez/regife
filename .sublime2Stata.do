use data/Divorce-Wolfers-AER, clear
keep if inrange(year, 1968, 1988)
egen state = group(st), label
tsset state year

regife div_rate unilateral divx*   [w=stpop], f(state year) d(2)



reghdfe div_rate unilateral divx* [aw=stpop], a(state year)
regife div_rate unilateral divx* [aw=stpop], f(state year) d(2) reps(1)
regife div_rate unilateral divx* [aw=stpop], f(state year) d(2) a(state year) reps(1)
regife div_rate unilateral divx*  i.state i.year [aw=stpop], f(state year) d(2) reps(1)






* check lagged and factor display correctly
regife div_rate L.unilateral i.year   [aw=stpop], f(state year) d(2) reps(1)
tab year, gen(tempyear)
regife div_rate L.unilateral tempyear*   [aw=stpop], f(state year) d(2) reps(1)



* check factors with time dimension gives the original data
distinct year
local nyear = r(ndistinct)
ife div_rate, f(state year) d(`nyear') gen(new) maxi(1)
assert float(div_rate) == float(new)
drop new



* compare different estimators
ccemg div_rate unilateral, f(state year)
ccep div_rate unilateral, f(state year) vce(cluster state)
regife div_rate unilateral, f(state year) a(state year) d(2)

