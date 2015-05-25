cd "/Users/Matthieu/Dropbox/Github/stata-regife"
insheet using "data/cigar.csv", clear

tsset state year
regife D.sales D.price, f(state year) d(2)

gen code = mod(_n, 10)

ife sales, f(state year) d(2) gen(res)

