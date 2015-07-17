webuse nlswork, clear
keep if id <= 100
regife ln_w tenure, f(id year, 1)
regife ln_w tenure, a(id) f(id year, 1)
regife ln_w tenure, a(id year) f(id year, 1)
regife ln_w tenure, f(fid = id fyear = year, 1)
drop fid* fyear*
regife ln_w tenure, a(feid = id feyear = year) f(fid = id fyear = year, 1)
regife ln_w tenure, f(id year, 1) residuals(newvar)
regife ln_w tenure, f(id year, 1) vce(bootstrap)
regife ln_w tenure, f(id year, 1) vce(bootstrap, cluster(id))
regife ln_w tenure, a(id) f(id year, 1) vce(bootstrap, cluster(id))




webuse nlswork, clear
keep if id <= 100
ife ln_w , f(fei = id fey = year, 1) res(new)
gen temp = abs(ln_w - fei1 * fey1 - new)
sum temp, meanonly
assert r(mean) < 0.001
drop temp new fei* fey*
ife ln_w , f(fei = id fey = year, 2) res(new)
gen temp = abs(ln_w - (fei1 * fey1 + fei2 * fey2 + new))
sum temp, meanonly
assert r(mean) < 0.001
drop temp new fei* fey*
ife ln_w , a(fe = id) f(fei = id fey = year, 2) res(new)
gen temp = abs(ln_w - (fe + fei1 * fey1 + fei2 * fey2 + new))
sum temp
assert r(sd) < 0.001
drop temp new fe* fei* fey*
ife ln_w , a(fe = year) f(fei = id fey = year, 2) res(new)
gen temp = abs(ln_w - (fe + fei1 * fey1 + fei2 * fey2 + new))
sum temp
assert r(sd) < 0.001
drop temp new fe* fei* fey*

ife ln_w , a(fe1 = year fe2 = id) f(fei = id fey = year, 2) res(new)
gen temp = abs(ln_w - (fe1 + fe2 + fei1 * fey1 + fei2 * fey2 + new))
sum temp
assert r(sd) < 0.001
drop temp new fe* fei* fey*




