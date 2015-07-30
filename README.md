

# regife (Bai 2009)

The command `regife` estimates models with interactive fixed effects using a least square estimate (Bai 2009).

For an observation `i`, denote `(jλ(i), jf(i))` the associated pair (id x time).  This package estimates the set of coefficients `beta`, of factors `(f1, .., fr)` and of loadings `(λ1, ..., λr)` that solve

![minimization](img/minimization.png)







### Syntax

`regife` requires a formula and the option `factors`, composed of an id variable, a time variable, and the dimension.

```
insheet "data/cigar.csv", clear
regife sales price, f(state year, 3)
```


### Options

#### Standard errors
The `vce` option allows to compute robust standard errors 

```
regife sales price, f(state year, 2) a(state year) vce(cluster state) 
```


#### Absorb
Impose id or time fixed effect with the option `absorb`. 

```
regife sales price, f(state year, 2)  a(state year)
```
#### Unbalanced Panel
The command handles unbalanced panels (ie missing observation for a given id, time) as described in the appendix of Bai 2009. 





#### Weights
Weights are supported but should be constant within id

#### Save factors
Save loadings and/or factors by specifying new variable names using `=`

```
regife sales price, a(fe = state) f(loading_state=state factor_year=year, 2) 
```

To obtain residuals, directly use the option `residuals`


```
regife sales price, f(state year, 2) residuals(newres)
```




## FAQ
#### When should I use interactive fixed effects?
Time fixed effects assume aggregate shocks impact each individual in the same way. In contrast, interactive fixed effects allow individuals to have different exposure to aggregate shocks. 


Another intuition for the difference is to consider the model 
```
Y = X'b + g(i, t) + ϵ
```
Fixed effect correspond to a first order taylor approximation of `g`, while interactive fixed effects correspond to a second order expansion of `g`.

#### Can't I just remove the endogeneity by replacing X with the residuals of X on a factor model?
For models with fixed effect, one can obtain consistent estimate by (i) demeaning regressors and (ii) use the residuals in the original regression.
In contrast, this method does not work i n models with interactive fixed effects. The intuition is that this kind of method (based on the FWL theorem) relies on linear projections, but factor models are non linears.


#### How are standard errors computed?
The standard errors are the ones obtained by a regression of y on x and covariates of the form `i.id#c.year` and `i.year#c.id`. This method is hinted in section 6 of of Bai (2009).

[Monte carlo evidence](monte-carlo/result.png) suggest to bootstrap the standard errors for small T.
```
regife sales price, f(state year, 2)  vce(bootstrap, reps(100))
```

To obtain standard errors by block bootstrap:
```
regife sales price, f(state year, 2)  vce(bootstrap, cluster(state))
```

#### What if I don't know the number of factors?
As proven in Moon Weidner (2015), overestimating the number of factors does not threaten the consistency of the interactive fixed effect estimates. The intuition is that irrelevant factors behave similarly to irrelevant covariates in a traditional OLS. A rule of thumb is to add factors until the result is not sensible to the number of factors.

#### Does regife implement the bias correction term in Bai (2009)?
In presence of correlation, the estimate for beta is biased (See Theorem 3 in Bai 2009 that derives the correction term). However, `regife` does not implement any correction. You may want to add enough factors until residuals are approximately i.i.d.


#### How can I speedup regife?

`regife` can be quite slow: tl convergence. 

- You can start the convergence at a given `beta` using `bstart`
- You can decrease the `tolerance` (default to 1e-9) or `maxiteration` (default to 10000).
- The Gauss-Seidel method requires a high number of iteration especially when X and the factor model are very correlated : this means `regife` is slow in those cases where the interactive fixed effect estimates substantially differ from the OLS estimates. For the same reason, adding id or time fixed effects generally makes the convergence much faster.
- I've written a [similar command](https://github.com/matthieugomez/PanelFactorModels.jl) in Julia, which is more than 100x faster


#### Why does regife return different estimates than the phht package in R?
The phht package in R also allows to compute the interactive fixed effect estimate in the case of balanced panels. This package returns wrong estimates in the case without fixed effects. 




# ife
The command `ife` estimates a factor model for a unique variable

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
`regife` is now available on ssc
```
ssc install regife
```

To install the latest version  on Github (or the one including ife)
- with Stata13+
	```
	net install regife, from(https://github.com/matthieugomez/stata-regife/raw/master/)
	```

- with Stata 12 or older, download the zipfiles of the repositories and run in Stata the following commands:
	```
	net install regife, from("SomeFolderRegife")
	```
