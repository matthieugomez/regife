

# regife

The command `regife` estimates models with interactive fixed effects (Bai 2009)


```
use "data/income-deregulation", clear
```

Model with interactive fixed effects (factor model of 3 dimension)
```
regife p30 intra_dummy, f(state year) d(3) reps(100)
```

Model with interactive fixed effects (state + year fe + factor model of 2 dimensions)

```
regife p30 intra_dummy,  a(state year)  f(state year) reps(100)
```


- The command handles unbalanced panels (ie missing observation for a given id x time) by imputing observations as described in the appendix of Bai 2009.
- The command `regife` is estimated on the residuals after removing the fixed effect specified in `absorb`. This is correct as long as the fixed effects are compatible with a factor model (in particular id fe and/or time fe)




## Standard errors


You should estimate standard errors by bootstrap, by specifying the option `reps`. If both `reps` and `cluster` are specified, standard errors are computed by block bootstrap.

```
regife p30 intra_dummy, f(state year) d(2) a(state year) cl(state) reps(50)
```


## Predict

You can save the loadings and factors using the symbol `=` in the option `factors`

```
regife p30 intra_dummy, f(loading_state=state factor_year=year) d(2) a(state year)  
```


After saving them,  you can generate the estimated factor using `predict`:

```
predict factors, f
```
- `f` returns the factor term
- `res` returns residuals without the factor term
- `xb` returns prediction without the factor term
- `resf` returns residuals augmented with the factor term
- `xbf` returns prediction augmented with the factor term

To use the option `f`, `xb` and `resf`, you need to save the interactive fixed effects first








## Syntax
The syntax is

```
regife depvar [indepvars]  [weight] [if] [in], 
	Factors(idvar timevar) Dimension(integer)  [
	Absorb(string) noCONS 
	reps(int 0) cluster(clustervars)
	TOLerance(real 1e-6) MAXIterations(int 10000) 
	]
```






# cce (Pesaran 2006)

The command `ccemg` and `ccep` correspond respectively to Pesaran (2006) Common Correlated Effects Mean Group estimator (CCEMG) and Common Correlated Effects Pooled estimator (CCEP). 

Like with year fixed effect, these commands generate the mean value of regressors at each time accross groups. and add them as regressor. After this step,
- `ccemg` runs the new model within each group and takes the average of beta accross all groups. Errors can be computed with the option `vce`
- `ccep` runs the new model on the pooled sample, interacting the mean regressors with group dummies. Errors are obtained by ttest for the betas.

These estimates are easy to compute manually compared to Bai (2009) estimate. I include them in this package for easier comparability. The ccemg estimate is also available in the Stata package [`xtmg`](https://ideas.repec.org/c/boc/bocode/s457238.html). 

The syntax for these estimates

```
ccemg p30 intra_dummy, f(state year)
ccep p30 intra_dummy, f(state year) vce(cluster state)
```




# ife
The command `ife` estimates a factor model for a given variable. 

Contrary to Stata usual `pca` command, 
- `ife` handles panel data (ie dataset where each row represents an id and a time) 
- `ife` handles unbalanced panel data : in the first step, missing observations are set to zero and a factor model is estimated.  In a second step, missing observations are replaced by the predicted value of the factor model, etc until convergence. This corresponds to the algorithm described in Stock and Watson (1998).




To generate the loadings and/or the factors, use the lhs of `=`
```
ife p30, f(loading_state=state factor_year=year)  d(2)
```

To directly generate the residual of the factor model, use `residuals`

```
ife p30, f(state year) residuals(p30_res)
```





By default, `ife` demeans the variable and estimates a factor model on it. If you want to estimate a `pca`, you probably want to demean with respect to id and/or time. To do so, use the option `absorb`.

```
ife p30, a(state) f(state year)  d(2) residuals(p30_res)
ife p30, a(state year) f(state year)  d(2) residuals(p30_res)
```




# Installation

If you have Stata >= 13

```
net install hdfe, from (https://raw.githubusercontent.com/sergiocorreia/reghdfe/master/package/)
net install regife, from(https://github.com/matthieugomez/stata-regife/raw/master/)
```
(`regife` requires `hdfe` if you want to estimate models with multiple fixed effects)



With Stata 12 or older, download the zipfiles of the repositories and run in Stata
```
net install hdfe, from("SomeFolderReghdfe")
net install regife, from("SomeFolderRegife")
```