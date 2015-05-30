use data/Divorce-Wolfers-AER, clear
keep if inrange(year, 1968, 1988)
egen state = group(st), label
tsset state year

regife div_rate unilateral   [w=stpop], f(state year) a(state year) d(2)



reghdfe div_rate unilateral [aw=stpop], a(state year)
regife div_rate unilateral [aw=stpop], f(state year) d(2) reps(1)
regife div_rate unilateral [aw=stpop], f(state year) d(2) a(state year) reps(1)
regife div_rate unilateral  i.state i.year [aw=stpop], f(state year) d(2) reps(1)






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


* ife 
ife div_rate [aw=stpop], f(fe1 = state year) d(2)



/* rough weight example */
gen weight2 = mod(state, 4) + 1
expand weight2 
bys state year: gen temp = _n
egen gs = group(state temp)
egen gs = group(state temp)
egen t = temp == 1
replace div_rate = 0 if missing(div_rate)


regife div_rate unilateral   [fw=weight2] if t , f(state year) d(2) reps(1) maxi(0)
regife div_rate unilateral , f(state year) d(2) reps(1) maxi(0)
