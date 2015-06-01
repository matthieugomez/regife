use income-deregulation, clear


reghdfe p20 L2.intra_dummy, a(state year)
regife p20 L2.intra_dummy, f(state year) d(3) a(state year)
ccemg p20 L2.intra_dummy, f(state year) 
ccep p20 L2.intra_dummy, f(state year)



regife p20 L2.intra_dummy, f(state year) d(3) a(state year) reps(10)


hdfe div_rate2 unilateral cons [aw=stpop], a(cons state#c.year) gen(ok)


regife div_rate2 unilateral,  f(state year) d(2) partial


/* check fast and non fast give (i) good cons or not (ii) same result */
regife div_rate2 unilateral,  f(state year) d(2) 
regife div_rate2 unilateral,  f(state year) d(2)  fast

regife div_rate2 unilateral,  f(state year) a(state) d(2) 
regife div_rate2 unilateral,  f(state year) a(state) d(2)  fast






regife div_rate unilateral [aw=stpop], f(state year) d(2) a(state year) 
regife div_rate unilateral  i.state i.year [aw=stpop], f(state year) d(2) 



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
regife div_rate unilateral, f(state year) a(state year) d(2) reps(1)


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
