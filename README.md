
The command `regife` estimates models with interactive fixed effects (Bai 2009)

# regife


The syntax is

```
regife depvar [indepvars]  [if] [in], Factors(idvar timevar) Dimension(integer)  [
	Absorb(string) noCONS 
	TOLerance(real 1e-6) MAXIterations(int 10000) 
	GENerate(newvarname)
]
```


#

webuse set https://github.com/matthieugomez/stata-regife/raw/master/data/
webuse  Divorce-Wolfers-AER, clear
egen state = group(st), label

reghdfe div_rate unilateral divx* if inrange(year, 1968, 1988)  [aw=stpop], a(state year)
regife div_rate unilateral divx* if inrange(year, 1968, 1988)  [w=stpop], f(state year) d(2)

### Save

Save the interactive fixed effect using the symbol `=` in the option `ife`

```
regife depvar [indepvars],  f(fi=id ft=timevar) 
```

Check that all the equations give the same coefficients:

```
regife depvar [indepvars], f(fi=id ft=timevar)  d(2)
reg y [indepvars] i.timevar#c.fi1 i.timevar#c.fi2
reg y [indepvars] i.idvar#c.ft1 i.idvar#c.ft2
```

Directly save the factors using the option `gen`

```
regife depvar [indepvars], f(idvar timevar) d(2) gen(res)
```

### Absorb
The command `regife` is estimated on the residuals after removing the fixed effect specified in `absorb`. The fixed effect specified in `absorb` *must* be compatible with the interactive fixed effect model (although currently `regife` does not check it is the case). The syntax for `absorb` is the same than `reghdfe`.


# ife
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