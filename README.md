

# regife (Bai 2009)

The command `regife` estimates models with interactive fixed effects following Bai (2009).

For an observation `i`, denote (`jλ(i)`, `jf(i)`) the associated pair (`id` x `time`).  The command estimates models of the form 

![model](img/model.png)


This algorithm returns the set of coefficients `β`, of factors `(f1, .., fr)` and of loadings `(λ1, ..., λr)` that minimize

![minimization](img/minimization.png)







## Syntax

`regife` requires a formula and the option `factors`, that specifies the id variable, the time variable, and the dimension:

```
insheet "data/cigar.csv", clear
regife sales price, f(state year, 3)
```




#### Absorb
You can impose id and/or time fixed effect by estimating models of the form

![model](img/femodel.png)

Just use the option `absorb` along with the option `factors`:

```
regife sales price, f(state year, 2)  a(state year)
```





#### Unbalanced Panel
The command handles unbalanced panels (ie missing observation for a given id, time) as described in the appendix of Bai 2009. 





#### Weights
Weights are supported but should be constant within id

#### Save factors
Save loadings and/or factors by specifying new variable names at the left hand side of `=`

```
regife sales price, a(fe_state=state) f(ife_state=state ife_year=year, 2) 
```

To save residuals, use the option `residuals`


```
regife sales price, f(state year, 2) residuals(newres)
```




## FAQ
#### When should I use interactive fixed effects?
Time fixed effects assume aggregate shocks impact each individual in the same way. In contrast, interactive fixed effects allow individuals to have different exposure to aggregate shocks. 

You can find such models in the following articles:

- Eberhardt, Helmers, Strauss (2013) *Do spillovers matter when estimating private returns to R&D?*
- Hagedorn, Karahan, Movskii (2015) *Unemployment Benefits and Unemployment in the Great Recession: The Role of Macro Effects*
- Hagedorn, Karahan, Movskii (2015) *The impact of unemployment benefit extensions on employment: the 2014 employment miracle?* 
- Totty (2015) *The Effect of Minimum Wages on Employment: A Factor Model Approach*

#### How are standard errors computed?
The `vce` option is passed to a regression of y on x and covariates of the form `i.id#c.year` and `i.year#c.id`. You can use any `vce` option availbe in `reghdfe`. This way of computing standard errors is hinted in section 6 of of Bai (2009).


```
regife sales price, f(state year, 2) a(state year) vce(cluster state) 
```


That being said, personal [Monte carlo evidence](monte-carlo/montecarlo.do) suggest to bootstrap the standard errors for small T.
```
regife sales price, f(state year, 2)  vce(bootstrap, reps(100))
regife sales price, f(state year, 2)  vce(bootstrap, cluster(state))
```

#### What if I don't know the number of factors?
As proven in Moon Weidner (2015), overestimating the number of factors still returns consistent estimates: irrelevant factors behave similarly to irrelevant covariates in a traditional OLS. A rule of thumb is to check that your estimate stays constant when you add more factors.

#### Does regife implement the bias correction term in Bai (2009)?
In presence of cross or time correlation beyond the factor structure, the estimate for beta is biased (but still consistent): see Theorem 3 in Bai 2009, which derives the correction term in special cases. However, `regife` does not implement any correction. You may want to add enough factors until residuals are approximately i.i.d.


#### How can I speedup the convergence?

- Start the convergence at a given `beta` using `bstart`.
- Decrease the `tolerance` (default to 1e-9) or `maxiteration` (default to 10000).
- The algorithm used in `regife` requires a lot of iterations when interactive fixed effects are correlated with the RHS variable. This means `regife` is slow exactly in those cases where the interactive fixed effect estimates substantially differ from the OLS estimates. For the same reason, adding id or time fixed effects generally makes the convergence much faster.
- I've written a [similar command](https://github.com/matthieugomez/PanelFactorModels.jl) in Julia, which is more than 100x faster


#### Can't `β` be simply estimated by replacing X with the residuals of X on a factor model?
For models with fixed effect, an equivalent way to obtain β is to first demean regressors within groups and then use their residuals in the original regression.
In contrast, this method does not work with models with interactive fixed effects. While fixed effects are linear projections (so that the Frisch-Waugh-Lovell theorem holds), factor models are non linear projections.


# Installation
`regife` is now available on ssc. It requires reghdfe

```
ssc install reghdfe
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
