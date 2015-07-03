/***************************************************************************************************

***************************************************************************************************/
tempfile filename
tempname postname
postfile `postname' N T i b1 b2 using `filename', replace
forval i = 1/500{
	local N 100
	local T 10
	local beta1 1
	local beta2 3
	local mu 5
	local gamma 2
	local delta 4

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
	gen X1 = 1 +  lambda1 * F1 +   lambda2 * F2 + lambda1 + lambda2 + F1 + F2 + rnormal()
	gen X2 = 1 +   lambda1 * F1 +  lambda2 * F2 + lambda1 + lambda2 + F1 + F2 + rnormal()

	gen Y = `beta1' * X1 + `beta2' * X2 + `mu' + `gamma' * x + `delta' * w +lambda1 * F1 + lambda2 * F2 + rnormal(9, 4)

	regife Y X1 X2 x w, f(N T) d(2) 
	post `postname' (`N') (`T') (`i') (_b[X1]) (_b[X2])
}
postclose `postname'
use `filename', clear