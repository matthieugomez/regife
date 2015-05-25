cd "/Users/Matthieu/Dropbox/Github/stata-regife"
insheet using "data/cigar.csv", clear

tsset state year
reg D.sales D.price 


gen dsales = D.sales
gen dprice = D.price

regife dsales dprice, id(state) t(year) d(2)
regife dsales dprice, id(state) t(year) d(2)


egen tempsales = mean(dsales), by(year)
egen tempprice = mean(dprice), by(year)
gen b = .
levelsof state, clean
foreach s in `r(levels)'{
	reg dsales dprice tempsales tempprice if state == `s'
	replace b = _b[dprice] if state == `s'
}
sum b
