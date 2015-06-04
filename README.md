

# regife

The command `regife` estimates models with interactive fixed effects (Bai 2009)


```
use "data/income-deregulation", clear
```

Model with interactive fixed effects (3 factors)
```
regife p30 intra_dummy, f(state year) d(3)
```

Model with interactive fixed effects (2 factors beyond state + year fe)

```
regife p30 intra_dummy,  a(state year)  f(state year) d(2)
```





## Standard errors

Standard errors returned by the commands are the ones obtained by the regression with time factors as regressors (interacted with id dummy).

```
regife p30 intra_dummy, f(state year) d(2) a(state year) reps(0)`
```

You should probably estimate it by bootstrap, in this case you can directly specify the option `reps`. If both `reps` and `cluster` are specified, standard errors are computed by block bootstrap.

```
regife p30 intra_dummy, f(state year) d(2) a(state year) cl(state) reps(50)
```


## Predict

You can save the loadings and factors using the symbol `=` in the option `factors`

```
regife p30 intra_dummy, f(fs=state fy=year) d(2) a(state year)  
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


The `absorb` option: the command `regife` is estimated on the residuals after removing the fixed effect specified in `absorb`. The fixed effect specified in `absorb` *must* be compatible with the interactive fixed effect model (in particular, it's fine to add id and year fixed effect). The syntax for `absorb` is the same than `reghdfe`.



# ife
The command `ife` estimates a factor model for a given variable. Contrary to Stata usual `pca` command, 
- this allows to estimate PCA on a dataset in a long form (in particular panel data)
- it handles unbalanced panels using an algorithm akin to Stock and Watson (1998). 

Missing combinations id x date in the dataset are considered to be missing, not zero.

```
ife p30 intra_dummy, f(fs=state fy=year)  d(2)
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




# Installation

`regife` requires `reghdfe` and `hdfe`:
```
net install reghdfe, from (https://raw.githubusercontent.com/sergiocorreia/reghdfe/updated_mata/package/)
net install hdfe, from (https://raw.githubusercontent.com/sergiocorreia/reghdfe/master/package/)
net install regife, from(https://github.com/matthieugomez/stata-regife/raw/master/)
```
