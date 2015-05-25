
The command `regife` estimates models with interactive fixed effects (Bai 2009)

# Syntax
The syntax is

```
regife y xvarlist, ife(idvar timevar)  Dimension(integer)  
[ Absorb(string) noCONS 
  convergence(real 0.000001) MAXiteration(int 500) 
  gen(string)
]
```


### Absorb
Specify supplementary fixed effect using the option `absorb` (which allows for same syntax than `reghdfe`). Fixed effect specified in `absorbs` must be compatible with an interactive fixed effect model with respect to `id` and `time` (`regife` is estimated on the residuals after removing the fixed effect specified in `absorb`)


### Save

Save the interactive fixed effect using the symbol `=` in the option `ife`

```
regife y xvarlist,  ife(fi=id ft=timevar) 
```

Check that all the equations give the same coefficients:

```
regife y xvarlist, ife(fi=id ft=timevar)  d(2)
reg y xvarlist i.timevar#c.fi1 i.timevar#c.fi2
reg y xvarlist i.idvar#c.ft1 i.idvar#c.ft2
```

Directly save the factors using the option `gen`

```
regife y xvarlist, ife(idvar timevar) d(2) gen(res)
```

# Installation

```
net install regife , from(https://github.com/matthieugomez/stata-regife/raw/master/)
```