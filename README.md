

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




### Unbalanced Panel
The command handles unbalanced panels (ie missing observation for a given id, time) as described in the appendix of Bai 2009. 



### Standard errors
The `vce` option allows to compute robust standard errors 

```
regife sales price, f(state year, 2) a(state year) vce(cluster state) 
```

Except for bootstrap, the `vce` option is simply passed to the command regressing y on x and covariates of the form `i.id#c.year` and `i.year#c.id` (as discussed in section 6 of of Bai 2009).

In presence of correlation, the estimate for beta isbiased (See Theorem 3 in Bai 2009). Instead of computing robust standard errors, you may want to add enough factors until residuals are i.i.d.


For small T, it seems wiser to bootstrap the standard errors, (see the following [monte carlo results](monte-carlo/result.png).
```
regife sales price, f(state year, 2)  vce(bootstrap, reps(100))
```

To obtain standard errors by block bootstrap:
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
regife sales price, a(fe = state) f(loading_state=state factor_year=year, 2) 
```

To obtain residuals, directly use the option `residuals`


```
regife sales price, f(state year, 2) residuals(newres)
```


### Speed
`regife` can be quite slow, and typically, a high number of iterations is required until convergence. 

- You can start the convergence at a given `beta` using `bstart`
- The more correlated X and the factor model are, the slower the Gauss-Seidel method. A first consequence is that regife is slow in these cases when estimates are far from the OLS estimates. A second consequence is that adding id or time fixed effects makes the convergence faster.
- I've written a [similar command](https://github.com/matthieugomez/PanelFactorModels.jl) in Julia, which is more than 100x faster



# ife
The command `ife` estimates a factor model for some variable

- Contrary to Stata `pca` command, 
 - `ife` handles dataset in long forms (ie `ife` works on datasets where each row represents an id and a time, rather than datasets where each row represents an id and each variable a time).
 - `ife` handles unbalanced panel data : in the first step, missing observations are set to zero and a factor model is estimated.  In a second step, missing observations are replaced by the predicted value of the factor model, etc until convergence. This corresponds to the algorithm described in Stock and Watson (1998).


- To save loadings and/or the factors, use the lhs of `=`
 ```
 ife sale, f(loading_state=state factor_year=year)  d(2)
 ```

 If you want to obtain residuals, you can directly use the `residuals` option

 ```
 ife sale, f(state year, 2) residuals(p30_res)
 ```

- By default, `ife` does not demean the variable. If you want to estimate a PCA, you probably want to demean the variable with respect to id and/or time. To do so, use the option `absorb`. 


 ```
 ife sale, a(fe_state = state fe_year = year) f(factors = state loading = year, 2) 
 ife sale, a(state year) f(state year, 2) residuals(p30_res)
 ```








# Installation

### regife
`regife` is now available on ssc
```
ssc install regife
```


To install the latest version  on Github 
- with Stata13+
	```
	net install regife, from(https://github.com/matthieugomez/stata-regife/raw/master/)
	```

- with Stata 12 or older, download the zipfiles of the repositories and run in Stata the following commands:
	```
	net install regife, from("SomeFolderRegife")
	```
