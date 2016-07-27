

# regife (Bai 2009)

The command `regife` estimates models with interactive fixed effects following Bai (2009).

For an observation `i`, denote (`jλ(i)`, `jf(i)`) the associated pair (`id` x `time`).  The command estimates models of the form 

![model](img/model.png)


The model is estimated by least square, i.e. by finding the coefficients `β`, of factors `(f1, .., fr)` and of loadings `(λ1, ..., λr)` that minimize

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
#### When should one use interactive fixed effects models?
Some litterature using this estimation procedure::

- Eberhardt, Helmers, Strauss (2013) *Do spillovers matter when estimating private returns to R&D?*
- Hagedorn, Karahan, Movskii (2015) *Unemployment Benefits and Unemployment in the Great Recession: The Role of Macro Effects*
- Hagedorn, Karahan, Movskii (2015) *The impact of unemployment benefit extensions on employment: the 2014 employment miracle?* 
- Totty (2015) *The Effect of Minimum Wages on Employment: A Factor Model Approach*

#### How are standard errors computed?
Errors are obtained by regressing y on x and covariates of the form `i.id#c.year` and `i.year#c.id`. This way of computing standard errors is hinted in section 6 of of Bai (2009).


#### Does this command implement the bias correction term in Bai (2009)?
In presence of cross or time correlation beyond the factor structure, the estimate for beta is consistent but biased (see Theorem 3 in Bai 2009, which derives the correction term in special cases). However, this package does not implement any correction. You may want to check that your residuals are approximately i.i.d.


## References
- Bai, Jushan. *Panel data models with interactive fixed effects.* (2009) Econometrica 
- Ilin, Alexander, and Tapani Raiko. *Practical approaches to principal component analysis in the presence of missing values.* (2010) The Journal of Machine Learning Research 11 
-  Koren, Yehuda. *Factorization meets the neighborhood: a multifaceted collaborative filtering model.* (2008) Proceedings of the 14th ACM SIGKDD international conference on Knowledge discovery and data mining. 
- Raiko, Tapani, Alexander Ilin, and Juha Karhunen. *Principal component analysis for sparse high-dimensional data.* (2008) Neural Information Processing.
- Srebro, Nathan, and Tommi Jaakkola. *Weighted low-rank approximations* (2010) The Journal of Machine Learning Research 11 
- Nocedal, Jorge and Stephen Wright *An Inexact Levenberg-Marquardt method for Large Sparse Nonlinear Least Squares*  (1985) The Journal of the Australian Mathematical Society

#### How can I speedup the convergence?

- Start the convergence at a given `beta` using `bstart`.
- Decrease the `tolerance` (default to 1e-9) or `maxiteration` (default to 10000).
- Save your dataset in `.csv` and use a similar [command](https://github.com/matthieugomez/PanelFactorModels.jl) in Julia, which is much more faster.


#### Can `β` be estimated by replacing X with the residuals of X on a factor model?
For models with fixed effect, an equivalent way to obtain β is to first demean regressors within groups and then regress `y` on these residuals instead of the original regressors.
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
