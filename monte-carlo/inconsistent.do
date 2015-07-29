
tempfile filename
tempname postname
postfile `postname' str10 estimator beta using `filename', replace

forval i = 1/1000{
	drop _all
	set obs 100
	gen N = _n
	gen lambda1 = rnormal()^2 + 3
	gen lambda2 = 0.2 * lambda1^2 + rnormal()
	tempfile temp
	save `temp'

	drop _all
	set obs 20
	gen T = _n
	gen F1 = rnormal() + 1
	gen F2 = 0.8 * F1 + rnormal() + 2
	cross using `temp'



	gen x = rnormal()
	gen X1 = lambda1 + lambda1 * cos(F1) + lambda1 * (F1 * rnormal()^2) + rnormal()
	gen Y = 3 + 3 * X1 + 10 * (lambda1 * F1) + rnormal()

	ife X1, a(N) f(N T, 1) res(resX1)
	reghdfe Y resX1, a(N)
	post `postname' ("3step") (_b[resX1])
	regife Y X1 , a(N) f(N T, 1)
	post `postname' ("regife") (_b[X1])
}

postclose `postname'
use `filename', clear
sumup b, by(estimator)




tempfile filename
tempname postname
postfile `postname' str10 estimator beta using `filename', replace

forval i = 1/10{
	drop _all
	set obs 100
	gen N = _n
	gen lambda1 = rnormal()^2 + 3
	gen lambda2 = 0.2 * lambda1^2 + rnormal()
	tempfile temp
	save `temp'

	drop _all
	set obs 20
	gen T = _n
	gen F1 = rnormal() + 1
	gen F2 = 0.8 * F1 + rnormal() + 2
	cross using `temp'

	gen X1 = 1 + lambda1 * F1 + rnormal()
	gen Y = 3 + 3 * X1 - 10 * lambda1 * F1 + rnormal() * lambda1 * F1 + rnormal()


	ife X1, a(N) f(N T, 1) res(resX1)
	reghdfe Y resX1, a(N)
	post `postname' ("3step") (_b[resX1])
	regife Y X1 , a(N) f(N T, 1)
	post `postname' ("regife") (_b[X1])
}

postclose `postname'
use `filename', clear
sumup b, by(estimator)




tempfile filename
tempname postname
postfile `postname' str10 estimator beta using `filename', replace

forval i = 1/100{
	drop _all
	set obs 100
	gen N = _n
	gen lambda1 = rnormal()^2 + 3
	gen lambda2 = 0.2 * lambda1^2 + rnormal()
	tempfile temp
	save `temp'

	drop _all
	set obs 20
	gen T = _n
	gen F1 = rnormal() + 1
	gen F2 = 0.8 * F1 + rnormal() + 2
	cross using `temp'

	gen X1 = 1 + lambda1 * F1 + rnormal()
	gen Y = 3 + 3 * X1 - 10 * lambda1 * F1 + 10 * rnormal()^2 * lambda1 * F1 + rnormal()

	ife X1, a(N) f(N T, 1) res(resX1)
	reghdfe Y resX1, a(N)
	post `postname' ("3step") (_b[resX1])
	regife Y X1 , a(N) f(N T, 2)
	post `postname' ("regife") (_b[X1])
}

postclose `postname'
use `filename', clear
sumup b, by(estimator)



tempfile filename
tempname postname
postfile `postname' str10 estimator beta using `filename', replace

forval i = 1/100{
	drop _all
	set obs 100
	gen z = rnormal()
	gen hidden = rnormal()
	gen X1 =  z  + rnormal()
	gen Y = 3 * X1 + 2 * z + 10 * hidden^2 * z + rnormal()
	reg Y X1 z
	post `postname' ("reg") (_b[X1])
	reg X1 z
	predict res, res
	reg Y res
	post `postname' ("ivreg") (_b[res])
}
postclose `postname'
use `filename', clear
sumup b, by(estimator)






