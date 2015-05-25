# stata-regife

A Stata command to compute Interactive fixed effect (Bai 2009)


## Syntax
Syntax is

```
regife yvar xvar, Id(idvar) Time(timevar) Dimension(integer)  ///
[Absorb(string) noCONS convergence(real 0.000001) MAXiteration(int 500) gen(string)]
```


### Absorb
You can specified supplementary fixed effect using the option `command`. Mathematically, `regife` is estimated on the residuals after removing the fixed effect. This mean fixed effect specified in `absorbs` must be compatible with an interactive fixed effect model with respect to `id` and `time`


### Save
Save the factors for residuals using the option `gen`

```
regife yvar xvar, Id(f1=idvar) Time(f2=timevar) gen(res)
```
Save the interactive fixed effect using the symbol equal

```
regife yvar xvar, Id(f1=idvar) Time(f2=timevar) d(2)
```

You can check that the following model gives the same coefficient:

```
reg yvar xvar i.timevar#c.f1
```

