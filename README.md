

# regife

The command `regife` estimates models with interactive fixed effects (Bai 2009)


```
use "data/Divorce-Wolfers-AER", clear
egen state = group(st), label
keep if inrange(year, 1968, 1988) 
```

Model with interactive fixed effects (3 factors)
```
regife div_rate unilateral, f(state year) d(3)
```

Model with interactive fixed effects (3 factors) beyond state + year fe

```
regife div_rate unilateral,  f(state year) a(state year) d(3)
```



## Standard errors
Correct standard errors must be obtained by boostrap : just use the option `reps`:

```
regife div_rate unilateral, f(state year) d(2) a(state year) reps(50)
```

To block bootstrap, use the option `cluster`

```
regife div_rate unilateral, f(state year) d(2) a(state year) reps(50) cl(state)
```


## Predict

Save the interactive fixed effect using the symbol `=` in the option `ife`

```
regife div_rate unilateral, f(fs=state fy=year) d(2) a(state year) reps(50)
```


After saving these fixed effects, you can generate the estimated factor using `predict`:

```
regife div_rate unilateral, f(fs=state fy=year)  d(2)
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
regife depvar [indepvars]  [aweight pweight fweight] [if] [in], ///
	Factors(idvar timevar) Dimension(integer)  [
	Absorb(string) noCONS 
	TOLerance(real 1e-6) MAXIterations(int 10000) 
	]
```


Note on the `absorb` option: the command `regife` is estimated on the residuals after removing the fixed effect specified in `absorb`. The fixed effect specified in `absorb` *must* be compatible with the interactive fixed effect model (although currently `regife` does not check it is the case). The syntax for `absorb` is the same than `reghdfe`.



# ife
The command `ife` estimates a factor model for a given variable. Contrary to Stata usual `pca` command, 
- this allows to estimate PCA on a dataset in a long form (in particular panel data)
- it handles unbalanced panels using an algorithm akin to Stock and Watson (1998). 

Missing combinations id x date in the dataset are considered to be missing, not zero.

```
ife div_rate unilateral, f(fs=state fy=year)  d(2)
```

# cce (Pesaran 2006)

The command `ccemg` and `ccep` correspond respectively to Pesaran (2006) Common Correlated Effects Mean Group estimator (CCEMG) and Common Correlated Effects Pooled estimator (CCEP). 

Like with year fixed effect, these commands generate the mean value of regressors at each time accross groups. and add them as regressor. After this step,
- `ccemg` runs the new model within each group and takes the average of beta accross all groups.
- `ccep` runs the new model on the pooled sample, interacting the mean regressors with group dummies

These estimates are easy to compute manually compared to Bai (2009) estimate. I included in this package for easier comparability. The ccemg estimate is also available in the Stata package [`xtmg`](https://ideas.repec.org/c/boc/bocode/s457238.html). 

The syntax for these estimates

```
ccemg depvar [indepvars]  [aweight pweight fweight] [if] [in], ///
	Factors(idvar timevar)
ccep depvar [indepvars]  [aweight pweight fweight] [if] [in], ///
	Factors(idvar timevar) vce(vceoption)
```



# Installation

```
net install regife, from(https://github.com/matthieugomez/stata-regife/raw/master/)
```

If you want to use the option `absorb`, you must download the command `hdfe` 

```
net install hdfe, from (https://raw.githubusercontent.com/sergiocorreia/reghdfe/master/package/)
```