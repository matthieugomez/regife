/***************************************************************************************************
Simple case
Unbalanced case
Absorb case
***************************************************************************************************/
tempfile filename
tempname postname
postfile `postname' N T i str20 sample b1 b2 gamma delta sb1 sb2 sgamma sdelta using `filename', replace

local beta1 1
local beta2 3
local mu 5
local mu1 1
local mu2 1
local gamma 2
local delta 4




local N 10000
local T 10


forval i = 1/1000{
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


	regife Y X1 X2 x w, f(N T, 2)   
	if e(converged) == "true"{
		post `postname' (`N') (`T') (`i') ("balanced") (_b[X1]) (_b[X2]) (_b[x]) (_b[w]) (_se[X1]) (_se[X2]) (_se[x]) (_se[w])
	}


	regife Y X1 X2 x w, f(N T, 2) a(N)   
	if e(converged) == "true"{
		post `postname' (`N') (`T') (`i') ("absorb") (_b[X1]) (_b[X2]) (_b[x]) (_b[w]) (_se[X1]) (_se[X2]) (_se[x]) (_se[w])
	}


	regife Y X1 X2 x w if mod(_n, 3) != 0, f(N T, 2)  
	if e(converged) == "true"{
		post `postname' (`N') (`T') (`i') ("unbalanced") (_b[X1]) (_b[X2]) (_b[x]) (_b[w]) (_se[X1]) (_se[X2]) (_se[x]) (_se[w])
	}
	display `i'

}
postclose `postname'
use `filename', clear


/***************************************************************************************************
clustered errors
***************************************************************************************************/
tempfile filename
tempname postname
postfile `postname' N T i str20 sample b1 b2 gamma delta sb1 sb2 sgamma sdelta using `filename', replace

local beta1 1
local beta2 3
local mu 5
local mu1 1
local mu2 1
local gamma 2
local delta 4




local N 10000
local T 1000


forval i = 1/1000{
	drop _all
	set obs `N'
	gen N = _n
	gen lambda1 = rnormal()
	gen lambda2 = rnormal()
	gen x = lambda1 + lambda2 + rnormal()
	gen error = rnormal()
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

	gen Y = `beta1' * X1 + `beta2' * X2 + `mu' + `gamma' * x + `delta' * w +lambda1 * F1 + lambda2 * F2 + runiform() *error + rnormal(0, 2)


	regife Y X1 X2 x w, f(N T, 2)   cl(N) a(N)
	if e(converged) == "true"{
		post `postname' (`N') (`T') (`i') ("balanced") (_b[X1]) (_b[X2]) (_b[x]) (_b[w]) (_se[X1]) (_se[X2]) (_se[x]) (_se[w])
	}

}
postclose `postname'
use `filename', clear


