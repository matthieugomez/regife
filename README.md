

# regife

The command `regife` estimates models with interactive fixed effects (Bai 2009)



The syntax is

```
regife depvar [indepvars]  [if] [in], Factors(idvar timevar) Dimension(integer)  [
	Absorb(string) noCONS 
	TOLerance(real 1e-6) MAXIterations(int 10000) 
	GENerate(newvarname)
]
```

# example

```
use "data/Divorce-Wolfers-AER", clear
egen state = group(st), label
keep if inrange(year, 1968, 1988) 
```

Model with state / year fe
```
reghdfe div_rate unilateral divx* [aw=stpop], a(state year) vce(cluster state)
```

Model with state / year fe + interactive fixed effects:

```
regife div_rate unilateral divx* [w=stpop],  f(state year) a(state year) d(4) vce(cluster state)
```



### save

Save the interactive fixed effect using the symbol `=` in the option `ife`

```
regife depvar [indepvars],  f(f1=firm f2=year) 
```

Check that all the equations give the same coefficients:

```
regife depvar [indepvars], f(fi=id ft=timevar)  d(2)
reg y [indepvars] i.timevar#c.fi1 i.timevar#c.fi2
reg y [indepvars] i.idvar#c.ft1 i.idvar#c.ft2
```

Obtain the factor using `predict`:

```
regife depvar [indepvars],  f(f1=firm f2=year) 
predict factors, f
```
- `f` returns the factor term
- `res` returns residuals without the factor term
- `xb` returns prediction without the factor term
- `resf` returns residuals augmented with the factor term
- `xbf` returns prediction augmented with the factor term

To use the option `f`, `xb` and `resf`, you need to save the interactive fixed effects first


### absorb

The command `regife` is estimated on the residuals after removing the fixed effect specified in `absorb`. The fixed effect specified in `absorb` *must* be compatible with the interactive fixed effect model (although currently `regife` does not check it is the case). The syntax for `absorb` is the same than `reghdfe`.



# ife
The command `ife` estimates a factor model for a given variable. Contrary to Stata usual `pca` command, this allows to estimate PCA on a dataset in a long form (in particular panel data). Moreover, it handles unbalanced panels using an algorithm akin to Stock and Watson (1998).

```
ife x, f(idvar timevar) d(2) gen(res)
```


# Standard errors
Standard errors reported are the ones of the dependent variable minus estimated factor structore over the regressor (adjusted for (T+N)*d number of freedoms corresponding to interactive fixed effects).
They are incorrect.

Correct standard errors can be obtained by boostrap. Just use the option `reps` (with `cluster` if you want to compute the right errors)

```
egidfe div_rate unilateral divx* , f(state year) d(2) a(state year) reps(50)
```

To block bootstrap by state:

```
egidfe div_rate unilateral divx* , f(state year) d(2) a(state year) reps(50) cl(state)
```


# installation


```
net install regife , from(https://github.com/matthieugomez/stata-regife/raw/master/)
```

`regife` requires the command `hdfe` if you use the option `absorb`

```
cap ado uninstall reghdfe
net install hdfe, from (https://raw.githubusercontent.com/sergiocorreia/reghdfe/master/package/)
```