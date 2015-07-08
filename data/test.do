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

