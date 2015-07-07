

# regife (Bai 2009)

The command `regife` estimates models with interactive fixed effects (Bai 2009). 


### Syntax

`regife` requires a formula and the option `factors`, composed of an id variable, a time variable, and the dimension.

```
insheet "data/cigar.csv", clear
regife sales price, f(state year, 3)
```




### Absorb
Impose id or time fixed effect with the option `absorb`. 

```
regife sales price, f(state year, 2)  a(state year)
```
In my experience, the convergence is much faster when id or time fixed effects are specified.




### Unbalanced Panel
The command handles unbalanced panels (ie missing observation for a given id, time) as described in the appendix of Bai 2009. In this case,  standard errors should be estimated by bootstrap.



### Standard errors
Robust standard errors can be specified with the option `vce`.  

```
regife sales price, f(state year, 2) a(state year) vce(cluster state) 
```
Except for bootstraped standard errors, the `vce` option is simply passed to the command regressing y on x and covariates of the form `i.id#c.year` and `i.year#c.id` (as discussed in section 6 of of Bai 2009).



In my experience, bootstraped errors are much more performant in finite sample:
```
regife sales price, f(state year, 2)  vce(bootstrap, reps(100))
```

Obtain standard errors by block bootstrap:
```
regife sales price, f(state year, 2)  vce(bootstrap, cluster(state))
```


### Weights
Weights are supported but should be constant within id


### Convergence
The iteration algorithm can be modified using the option `tolerance` (default to 1e-9) or `maxiteration` (default to 10000).



### Save factors
Save loadings and/or factors by specifying new variable names using `=`
```
regife sales price, f(loading_state=state factor_year=year, 2) 
```






# ife
The command `ife` estimates a factor model for some variable

- Contrary to Stata `pca` command, 
 - `ife` handles dataset in long forms (ie `ife` works on datasets where each row represents an id and a time, rather than datasets where each row represents an id and each variable a time).
 - `ife` handles unbalanced panel data : in the first step, missing observations are set to zero and a factor model is estimated.  In a second step, missing observations are replaced by the predicted value of the factor model, etc until convergence. This corresponds to the algorithm described in Stock and Watson (1998).


- To save loadings and/or the factors, use the lhs of `=`
 ```
 ife sale, f(loading_state=state factor_year=year)  d(2)
 ```

- By default, `ife` does not demean the variable. If you want to estimate a PCA, you probably want to demean the variable with respect to id and/or time. To do so, use the option `absorb`. 


 ```
 ife sale, a(fe_state = state) f(factors = state loading = year, 2)  
 ife sale, a(fe_state = state fe_year = year) f(factors = state loading = year, 2) 
 ```

- The residual from the overall factor model can then be obtained in the following way:

 ```
 gen predict = p30 - (fe_state + fe_year + factors_1 * loading_1 + factors_2 * loading_2)
 ```

 Instead of saving each part of the factor model, obtain directly the residuals using the `residuals` option

 ```
 ife sale, f(state year, 2) residuals(p30_res)
 ife sale, a(state year) f(state year, 2) residuals(p30_res)
 ```




# cce (Pesaran 2006)

The command `ccemg` and `ccep` correspond respectively to Pesaran (2006) Common Correlated Effects Mean Group estimator (CCEMG) and Common Correlated Effects Pooled estimator (CCEP). 

Like with year fixed effect, these commands generate the mean value of regressors at each time accross groups. and add them as regressor. After this step,
- `ccemg` runs the new model within each id and averages the betas accross all ids. 
- `ccep` runs the new model on the pooled sample, interacting the newly created variables with group dummies. 

`ccep` relies on `reghdfe` and is generally faster than `ccemg`.




# Installation
`regife` requires [`reghdfe` and `hdfe`](https://github.com/sergiocorreia/reghdfe) with version 3.0+

If you have Stata 13+

```
net install hdfe, from (https://raw.githubusercontent.com/sergiocorreia/reghdfe/master/package/)
net install reghdfe, from (https://raw.githubusercontent.com/sergiocorreia/reghdfe/master/package/)
net install regife, from(https://github.com/matthieugomez/stata-regife/raw/master/)
```



With Stata 12 or older, download the zipfiles of the repositories and run in Stata the following commands:
```
net install hdfe, from("SomeFolderHdfe")
net install reghdfe, from("SomeFolderReghdfe")
net install regife, from("SomeFolderRegife")
```