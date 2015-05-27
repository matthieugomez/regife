

# regife

The command `regife` estimates models with interactive fixed effects (Bai 2009)


```
use "data/Divorce-Wolfers-AER", clear
egen state = group(st), label
keep if inrange(year, 1968, 1988) 
```

Model with interactive fixed effects
```
regife div_rate unilateral, f(state year) d(3)
```

Model with state / year fe + interactive fixed effects:

```
regife div_rate unilateral,  f(state year) a(state year) d(3)
```



## Standard errors
Correct standard errors can be obtained by boostrap. Just use the option `reps`:

```
regife div_rate unilateral, f(state year) d(2) a(state year) reps(50)
```

To block bootstrap by state:

```
regife div_rate unilateral, f(state year) d(2) a(state year) reps(50) cl(state)
```


## Predict

Save the interactive fixed effect using the symbol `=` in the option `ife`

```
regife div_rate unilateral , f(fs=state fy=year) d(2) a(state year) reps(50)
```


Generate the estimated factor using `predict`:

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


The command `regife` is estimated on the residuals after removing the fixed effect specified in `absorb`. The fixed effect specified in `absorb` *must* be compatible with the interactive fixed effect model (although currently `regife` does not check it is the case). The syntax for `absorb` is the same than `reghdfe`.



# ife
The command `ife` estimates a factor model for a given variable. Contrary to Stata usual `pca` command, this allows to estimate PCA on a dataset in a long form (in particular panel data). Moreover, it handles unbalanced panels using an algorithm akin to Stock and Watson (1998).

```
ife div_rate unilateral, f(fs=state fy=year)  d(2)
```


# Installation


```
net install regife , from(https://github.com/matthieugomez/stata-regife/raw/master/)
```

If you want to use the option `absorb`, you should download the command `hdfe` 

```
cap ado uninstall reghdfe
net install hdfe, from (https://raw.githubusercontent.com/sergiocorreia/reghdfe/master/package/)
```