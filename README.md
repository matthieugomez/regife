
The command `regife` estimates models with interactive fixed effects (Bai 2009)

# Syntax
The syntax is

```
regife y xvarlist, Factors(idvar timevar)  Dimension(integer)  [
	Absorb(string) noCONS 
	convergence(real 0.000001) MAXiteration(int 500) 
	GENerate(newvarname)
]
```


### Absorb
Specify supplementary fixed effect using the option `absorb` (which allows for same syntax than `reghdfe`). Fixed effect specified in `absorbs` must be compatible with an interactive fixed effect model with respect to `id` and `time` (`regife` is estimated on the residuals after removing the fixed effect specified in `absorb`)


### Save

Save the interactive fixed effect using the symbol `=` in the option `ife`

```
regife y xvarlist,  f(fi=id ft=timevar) 
```

Check that all the equations give the same coefficients:

```
regife y xvarlist, f(fi=id ft=timevar)  d(2)
reg y xvarlist i.timevar#c.fi1 i.timevar#c.fi2
reg y xvarlist i.idvar#c.ft1 i.idvar#c.ft2
```

Directly save the factors using the option `gen`

```
regife y xvarlist, f(idvar timevar) d(2) gen(res)
```


### ife

The command `ife` estimates a factor model for a given variable. Contrary to Stata usual `pca` command, this allows to estimate PCA on a dataset in a long form (in particular panel data). Moreover, it handles unbalanced panels using an algorithm akin to Stock and Watson (1998).

```
ife x, f(idvar timevar) d(2) gen(res)
```

# Installation

`regife` requires the command `hdfe` if you use the option `absorb`

```
cap ado uninstall reghdfe
net install hdfe, from (https://raw.githubusercontent.com/sergiocorreia/reghdfe/master/package/)
net install regife , from(https://github.com/matthieugomez/stata-regife/raw/master/)
```