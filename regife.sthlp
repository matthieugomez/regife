{smcl}
{* *! version 0.2 08jul2015}{...}
{vieweralsosee "tabstat" "help tabstat"}{...}
{viewerjumpto "Syntax" "sumup##syntax"}{...}
{viewerjumpto "Description" "sumup##description"}{...}
{viewerjumpto "Options" "sumup##options"}{...}
{viewerjumpto "Examples" "sumup##examples"}{...}
{viewerjumpto "Stored results" "sumup##results"}{...}
{viewerjumpto "References" "sumup##references"}{...}
{viewerjumpto "Author" "sumup##contact"}{...}

{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:reghdfe} {hline 2}}Linear models with interactive fixed effects{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}


{p 8 15 2} {cmd:regife}
{depvar} [{indepvars}] 
{ifin} {it:{weight}}
{cmd:,} 
{opth f:actors(idvar timevar, dimensionint)} 
[{help regife##options:options}] {p_end}

{marker description}{...}
{title:Description}

{pstd}
{cmd:regife} fits a model with interactive fixed effects following Bai (2009). Optionally, it saves the estimated factors. Errors are computed following the regressions indicated in Section 6, but Monte Carlo evidence suggest bootstraps performs better in finite sample. The program requires {help reghdfe} and {help hdfe} to be installed (both are available on SSC).


{marker options}{...}
{title:Options}
{synoptset 30 tabbed}{...}
{synopthdr}
{synoptline}
{synopt :{opt f:actors(idvar timevar, dimensionint)}} id variable, time variable, and factor dimension. To save the estimates for the factor model, write {newvar}{inp:={it:factorvar}}.{p_end}
{synopt :{opt a:bsorb}{cmd:(}{help reghdfe##absvar:absvar}[...]{cmd:)}} identifiers of the fixed effects that will be absorbed. To save the estimates for the fixed effect, write {newvar}{inp:={it:absvarvar}}.{p_end}
{synopt:{opt vce}{cmd:(}{help reghdfe##vcetype:vcetype}[, {it:opt}]{cmd:)}}{it:vcetype}}
is {opt un:adjusted}/{opt ols} (default), {opt r:obust}, {opt bootrap} or {opt cl:uster} {it:clustervars}. Monte carlo evidence suggests that bootstrap performs better in finite sample{p_end}
{synopt:{opt tol:erance(#)}} specifies the tolerance criterion for convergence; default is {cmd:tolerance(1e-9)}{p_end}
{synopt:{opt max:iterations(#)}} specifies the maximum number of iterations; default is {cmd:maxiterations(5000)}. 0 corresponds to an illimited number of iterations{p_end}
{synopt :{opt res:iduals(newvar)}} save residuals {p_end}
{synopt :{opt bstart(matrix)}} start the iteration algorithm at a given value  for b{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}{cmd:fweight}s, {cmd:aweight}s and {cmd:pweight}s are allowed but should be constant within idvar; see {help weight}.{p_end}




{marker examples}{...}
{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. webuse nlswork}{p_end}
{phang2}{cmd:. keep if id <= 100}{p_end}

{pstd}Factor model in id, year of dimension 1{p_end}
{phang2}{cmd:. regife  ln_w tenure, f(id year, 1)}{p_end}
{pstd}Model including id fixed effect, and a factor model in id, year of dimension 2{p_end}
{phang2}{cmd:. regife  ln_w tenure, a(id) f(id year, 1)}{p_end}
{pstd}Model including id fixed effect, year fixed effect,  and a factor model in id, year of dimension 1{p_end}
{phang2}{cmd:. regife  ln_w tenure, a(id year) f(id year, 1)}{p_end}
{pstd}Save interactive fixed effects{p_end}
{phang2}{cmd:. regife  ln_w tenure, f(fid = id fyear = year, 1)}{p_end}
{pstd}Save fixed effects and interactive fixed effects{p_end}
{phang2}{cmd:. regife  ln_w tenure, a(feid = id fe year = year) f(fid = id fyear = year, 1)}{p_end}
{pstd}Generate residuals{p_end}
{phang2}{cmd:. regife  ln_w tenure, f(id year, 1) residuals(newvar)}{p_end}
{pstd}Bootstrap standard errros{p_end}
{phang2}{cmd:. regife  ln_w tenure, f(id year, 1) vce(bootstrap)}{p_end}
{pstd}Block bootstrap with respect to id{p_end}
{phang2}{cmd:. regife  ln_w tenure, f(id year, 1) vce(bootstrap, cluster(id))}{p_end}




{marker results}{...}
{title:Stored results}

{pstd}
{cmd:regife} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(df_r)}}residual degrees of freedom{p_end}
{synopt:{cmd:e(df_m)}}model degrees of freedom{p_end}
{synopt:{cmd:e(F)}}residual sum of squares{p_end}
{synopt:{cmd:e(mss)}}total sum of squares{p_end}
{synopt:{cmd:e(rss)}}residual sum of squares{p_end}
{synopt:{cmd:e(iterations)}}number of iterations{p_end}
{synopt:{cmd:e(error)}}convergence error{p_end} 

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:regife}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(indepvars)}}names of independent variables{p_end}
{synopt:{cmd:e(id)}}id variable {p_end}
{synopt:{cmd:e(time)}}time variable {p_end}
{synopt:{cmd:e(dimension)}}dimension{p_end}
{synopt:{cmd:e(congerged)}}did the algorithm converge?{p_end}




{marker references}{...}
{title:References}

{pstd}
{cmd:regife} implements the estimate proposed by:

{phang}
Jushan Bai. "Panel Data Models with Interactive Fixed Effects".
{it:Econometrica, 77.4 (2009): 1229-1279.}
{p_end}


{marker contact}{...}
{title:Author}

{phang}
Matthieu Gomez

{phang}
Department of Economics, Princeton University

{phang}
Please report issues on Github
{browse "https://github.com/matthieugomez/stata-regife":https://github.com/matthieugomez/stata-regife}
{p_end}


