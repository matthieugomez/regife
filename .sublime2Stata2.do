use data/Divorce-Wolfers-AER, clear
keep if inrange(year, 1968, 1988)
egen state = group(st), label
tsset state year
gen div_rate2 = div_rate
replace div_rate2 = 0 if missing(div_rate2)

timer clear
timer on 1
regife div_rate2 unilateral,  f(state year) d(2) partial
timer off 1
timer list
