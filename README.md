# stata-regife

A Stata command to compute Interactive fixed effect (Bai 2009)


## Syntax
Syntax is

```
regife y xvarlist, Id(idvar) Time(timevar) Dimension(integer)  ///
[Absorb(string) noCONS convergence(real 0.000001) MAXiteration(int 500) gen(string)]
```


### Absorb
You can specified supplementary fixed effect using the option `command`. 

Important: fixed effect specified in `absorbs` must be compatible with an interactive fixed effect model with respect to `id` and `time`. Mathematically, `regife` is estimated on the residuals after removing the fixed effect. 


### Save
Save the factors for residuals using the option `gen`

```
regife y xvarlist, i(f1=idvar) time(f2=timevar) gen(res)
```
Save the interactive fixed effect using the symbol equal

```
regife y xvarlist, i(f1=idvar) time(f2=timevar) d(2)
```

You can check that all the equations give the same coefficients for `x`

```
regife y xvarlist, i(f1=idvar) t(f2=timevar) d(2)
reg y xvarlist i.timevar#c.f1
reg y xvarlist i.idvar#c.f2
```

