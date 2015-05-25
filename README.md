
A Stata command to compute Interactive fixed effect (Bai 2009)


The syntax is

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
regife y xvarlist, i(fi=idvar) time(ft=timevar) gen(res)
```
Save the interactive fixed effect using the option `absorb`

```
regife y xvarlist, i(fi=idvar) time(ft=timevar)
```

You can check that all the equations give the same coefficients:

```
regife y xvarlist, i(f1=idvar) t(f2=timevar) d(2)
reg y xvarlist i.timevar#c.fi1 i.timevar#c.fi2
reg y xvarlist i.idvar#c.ft1 i.idvar#c.ft2
```

