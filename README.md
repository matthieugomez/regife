

# regife

The command `regife` estimates models with interactive fixed effects (Bai 2009)



## Syntax

Estimate models with interactive fixed effects by specifying the option `factors`, containing an id variable, a time variable, and the dimension.

```
use "data/income-deregulation", clear
regife p30 intra_dummy, f(state year, 3)
```


## Save factors
Save the loadings and factors by specifying new variable names using `=`
```
regife p30 intra_dummy, f(loading_state=state factor_year=year, 2) 
```


## Absorb
Impose id or time fixed effect with the option `absorb`. The convergence is generally *much* faster when id or time fixed effects are specified.

```
regife p30 intra_dummy,  a(state year)  f(state year, 2)
```
estimates a model with state fe + year fe + factor model of 2 dimensions





## Unbalanced Panel
The command handles unbalanced panels (ie missing observation for a given id x time) as described in the appendix of Bai 2009. In this case,  *standard errors should be estimated by bootstrap* 

## Weights
Weights are supported but should be constant within id


## Standard errors
Robust standard errors can be specified with the option `vce`. 

In all cases except boostrap, the option is simply passed to the regression of y on x and covariates of the form `i.id#c.year` and `i.year#c.id` (as discussed in section 6 of of Bai 2009).


```
regife p30 intra_dummy, f(state year, 2) a(state year) vce(cluster state) 
```

To compute boostraped standard errors
```
regife p30 intra_dummy, f(state year, 2)  vce(bootstrap, reps(100) cluster(state))
```




## Convergence
Modify when the iteration stops by using the option `tolerance` (default to 1e-9) or `maxiteration` (default to 10000)





# cce (Pesaran 2006)

The command `ccemg` and `ccep` correspond respectively to Pesaran (2006) Common Correlated Effects Mean Group estimator (CCEMG) and Common Correlated Effects Pooled estimator (CCEP). 

Like with year fixed effect, these commands generate the mean value of regressors at each time accross groups. and add them as regressor. After this step,
- `ccemg` runs the new model within each group and takes the average of beta accross all groups. Errors can be computed with the option `vce`
- `ccep` runs the new model on the pooled sample, interacting the mean regressors with group dummies. 

These estimates can be helpful to start the iteration for the `regife` estimator using the option `bstart'.





# ife
The command `ife` estimates a factor model for a variable

- Contrary to Stata usual `pca` command, 
 - `ife` handles panel data (ie dataset where each row represents an id and a time) 
 - `ife` handles unbalanced panel data : in the first step, missing observations are set to zero and a factor model is estimated.  In a second step, missing observations are replaced by the predicted value of the factor model, etc until convergence. This corresponds to the algorithm described in Stock and Watson (1998).


- To generate the loadings and/or the factors, use the lhs of `=`
 ```
 ife p30, f(loading_state=state factor_year=year)  d(2)
 ```

- By default, `ife` does not demean the variable. If you want to estimate a PCA, you probably want to demean the variable with respect to id and/or time. To do so, use the option `absorb`. 


 ```
 ife p30, a(fe_state = state) f(factors = state loading = year, 2)  
 ife p30, a(fe_state = state fe_year = year) f(factors = state loading = year, 2) 
 ```

- The residual from the overall factor model can then be obtained in the following way:

 ```
 gen predict = p30 - (fe_state + fe_year + factors_1 * loading_1 + factors_2 * loading_2)
 ```

 Instead of saving each part of the factor model, obtain directly the residuals using the `residuals` option

 ```
 ife p30, f(state year, 2) residuals(p30_res)
 ife p30, a(state year) f(state year, 2) residuals(p30_res)
 ```



## Installation

If you have Stata >= 13

```
net install hdfe, from (https://raw.githubusercontent.com/sergiocorreia/reghdfe/master/package/)
net install reghdfe, from (https://raw.githubusercontent.com/sergiocorreia/reghdfe/master/package/)
net install regife, from(https://github.com/matthieugomez/stata-regife/raw/master/)
```
(`regife` requires `hdfe` if you want to estimate models with multiple fixed effects)



With Stata 12 or older, download the zipfiles of the repositories and run in Stata
```
net install hdfe, from("SomeFolderHdfe")
net install reghdfe, from("SomeFolderReghdfe")
net install regife, from("SomeFolderRegife")
```