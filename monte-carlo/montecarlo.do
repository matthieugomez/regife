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




local Ns 100 100 100
local Ts 10 20 50

forval j = 1/3{
	local N `:word `j' of `Ns''
	local T `:word `j' of `Ts''

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


		/* simple case */
		matrix bstart = [`mu', `beta1' , `beta2', `gamma', `delta']
		regife Y X1 X2 x w, f(N T, 2)   bstart(`bstart')

		if e(converged) == "true"{
			post `postname' (`N') (`T') (`i') ("balanced") (_b[X1]) (_b[X2]) (_b[x]) (_b[w]) (_se[X1]) (_se[X2]) (_se[x]) (_se[w])
		}

		/* absorb case */
		matrix bstart = [`beta1' , `beta2', `delta']
		regife Y X1 X2 w, f(N T, 2) a(N) bstart(`bstart')  
		if e(converged) == "true"{
			post `postname' (`N') (`T') (`i') ("absorb") (_b[X1]) (_b[X2]) (_b[x]) (_se[X1]) (_se[X2]) (_se[w])
		}

		/* unbalanced case */
		matrix bstart = [`mu', `beta1' , `beta2', `gamma', `delta']
		regife Y X1 X2 x w if mod(_n, 3) != 0, f(N T, 2)  bstart(`bstart')
		if e(converged) == "true"{
			post `postname' (`N') (`T') (`i') ("unbalanced") (_b[X1]) (_b[X2]) (_b[x]) (_b[w]) (_se[X1]) (_se[X2]) (_se[x]) (_se[w])
		}
		display `i'
	}
}
postclose `postname'
use `filename', clear

drop _all
set obs 100
gen N = _n
gen lambda1 = rnormal()
gen lambda2 = rnormal()
gen x = lambda1 + lambda2 + rnormal()
tempfile temp
save `temp'


drop _all
set obs 20
gen T = _n
gen F1 = rnormal()
gen F2 = rnormal()
gen w = F1 + F2 + rnormal()
cross using `temp'

gen eta_1 = rnormal()
gen eta_2 = rnormal()
gen X1 = 3 + 4 * lambda1  + lambda2 * F2 
gen Y = 3 + 3 * X1 +  lambda1 * F1 +  lambda2 * F2 + lambda1 + lambda2 + F1 + F2 + rnormal()

ife X1, a(N) f(N T, 2) res(resx)
ivreg Y (X1 = res)
