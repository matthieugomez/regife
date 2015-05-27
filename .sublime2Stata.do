use data/Divorce-Wolfers-AER, clear
keep if inrange(year, 1968, 1988)
egen state = group(st), label
tsset state year

