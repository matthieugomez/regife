/***************************************************************************************************
This simulation reproduces Table 1 in Bai (2009), adding the case of unbalanced sample
***************************************************************************************************/
tempfile filename
tempname postname
postfile `postname' N T i str20 sample b1 b2 x w using `filename', replace

local beta1 1
local beta2 3
local mu 1
local mu1 1
local mu2 5
local gamma 2
local delta 4




local N 100
local T 50

forval i = 1/100{
	drop _all
	set obs `N'
	gen N = _n
	gen lambda1 = rnormal()
	gen lambda2 = rnormal()
	gen x = lambda1 + lambda2 + rnormal()
	tempfile temp
	save `temp'

	drop _all
	set obs `T'
	gen T = _n
	gen F1 = rnormal()
	gen F2 = rnormal()
	gen w = F1 + F2 + rnormal()
	cross using `temp'

	gen eta_1 = rnormal()
	gen eta_2 = rnormal()
	gen X1 = `mu1' + lambda1 * F1 + lambda2 * F2 + lambda1 + lambda2 + F1 + F2 + rnormal()
	gen X2 = `mu2' +  lambda1 * F1 +  lambda2 * F2 + lambda1 + lambda2 + F1 + F2 + rnormal()

	gen Y = `beta1' * X1 + `beta2' * X2 + `mu' + `gamma' * x + `delta' * w +lambda1 * F1 + lambda2 * F2 + rnormal(0, 2)


	regife Y X1 X2 x w, f(N T) d(2)  fast 
	post `postname' (`N') (`T') (`i') ("balanced") (_b[X1]) (_b[X2]) (_b[x]) (_b[w])

	regife Y X1 X2 x w if mod(_n, 3) != 0, f(N T) d(2)  fast 
	post `postname' (`N') (`T') (`i') ("unbalanced") (_b[X1]) (_b[X2]) (_b[x]) (_b[w])

	display `i'

}
postclose `postname'
use `filename', clear
