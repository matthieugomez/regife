
# regife

The command `regife` estimates linear models with interactive fixed effects following [Bai (2009)](https://onlinelibrary.wiley.com/doi/pdf/10.3982/ECTA6135).

**Note:** Interactive fixed effects are factor models (loadings x factors). If you want *interacted* fixed effects (e.g., `i.state#i.year`), use [reghdfe](https://github.com/sergiocorreia/reghdfe) instead.

For an observation `i`, denote (`j_id(i)`, `j_time(i)`) the associated pair (`id` x `time`). The command estimates models of the form

![model](img/model.png)

by finding the coefficients `beta`, factors `(f1, .., fd)`, and loadings `(lambda1, ..., lambdad)` that minimize

![minimization](img/minimization.png)

The command is slow. A much faster implementation is available in [Julia](https://github.com/matthieugomez/InteractiveFixedEffectModels.jl).

## Installation

`regife` is available on SSC. It requires [reghdfe](https://github.com/sergiocorreia/reghdfe). Using `absorb()` with multiple variables additionally requires [hdfe](https://ideas.repec.org/c/boc/bocode/s457928.html).

```stata
ssc install reghdfe
ssc install regife
```

To install the latest version from GitHub:

```stata
net install regife, from("https://raw.githubusercontent.com/matthieugomez/regife/main/")
```

## Syntax

### `regife`

```
regife depvar [indepvars] [if] [in] [weight], ife(idvar timevar, d) [options]
```

#### Required option

- **`ife(idvar timevar, d)`** specifies the id variable, time variable, and the dimension `d` of the factor model.

#### Optional

- **`absorb(absvar [...])`** absorbs standard fixed effects (passed to `reghdfe`). Example: `absorb(state year)`.
- **`vce(vcetype)`** specifies the variance-covariance estimator: `unadjusted` (default), `robust`, `bootstrap`, or `cluster clustvar`. Monte Carlo evidence suggests bootstrap performs better in finite samples.
- **`residuals(newvar)`** saves residuals to a new variable.
- **`tolerance(#)`** convergence tolerance; default `1e-9`.
- **`maxiterations(#)`** maximum iterations; default `10000`.
- **`bstart(matrix)`** provides a starting value for the coefficient vector.

Weights (`fweight`, `aweight`, `pweight`) are supported but must be constant within `idvar`.

### `ife`

`ife` extracts factors and loadings from a single variable (without regression). This is useful for 3-step estimation: extract factors from a variable, then use them as controls in a second regression.

```
ife varname [if] [in] [weight], factors(idvar timevar, d) [absorb(absvar) residuals(newvar)]
```

## Examples

### Basic usage

```stata
webuse nlswork
keep if id <= 100
regife ln_w tenure, ife(id year, 1)
```

### With absorbed fixed effects

Impose id and/or time fixed effects alongside the interactive fixed effects:

```stata
regife ln_w tenure, ife(id year, 2) absorb(id year)
```

### Save factors and loadings

Specify new variable names at the left-hand side of `=` inside `ife()`:

```stata
regife ln_w tenure, ife(ife_id=id ife_year=year, 2)
```

This creates variables `ife_id1`, `ife_id2` (loadings) and `ife_year1`, `ife_year2` (factors). Similarly for `absorb()`:

```stata
regife ln_w tenure, absorb(fe_id=id) ife(ife_id=id ife_year=year, 2)
```

### Save residuals

```stata
regife ln_w tenure, ife(id year, 2) residuals(newres)
```

### Bootstrap standard errors

```stata
regife ln_w tenure, ife(id year, 1) vce(bootstrap)
regife ln_w tenure, ife(id year, 1) vce(bootstrap, cluster(id))
```

### Extract factors with `ife`

```stata
ife ln_w, factors(fei=id fey=year, 2) absorb(id) residuals(res)
```

This decomposes `ln_w` (after absorbing the `id` fixed effect) into loadings (`fei1`, `fei2`), factors (`fey1`, `fey2`), and residuals (`res`).

### Unbalanced panels

The command handles unbalanced panels (i.e., missing observations for a given id-time pair) as described in the appendix of Bai (2009).


## FAQ

### When should one use interactive fixed effects models?

Interactive fixed effects are useful when unobserved heterogeneity has a factor structure -- for instance, when units respond differently to common shocks. Some applications:

- Eberhardt, Helmers, Strauss (2013) *Do Spillovers Matter When Estimating Private Returns to R&D?*
- Hagedorn, Karahan, Manovskii (2015) *Unemployment Benefits and Unemployment in the Great Recession: The Role of Macro Effects*
- Hagedorn, Karahan, Manovskii (2015) *The Impact of Unemployment Benefit Extensions on Employment: The 2014 Employment Miracle?*
- Totty (2015) *The Effect of Minimum Wages on Employment: A Factor Model Approach*

### How are standard errors computed?

Standard errors are obtained by regressing `y` on `x` and covariates of the form `i.id#c.factor` and `i.time#c.loading`, as suggested in Section 6 of Bai (2009). Monte Carlo evidence suggests bootstrap standard errors (`vce(bootstrap)`) perform better in finite samples.

### Does this command implement the bias correction in Bai (2009)?

No. In the presence of cross-sectional or time-series correlation beyond the factor structure, the estimator for beta is consistent but biased (see Theorem 3 in Bai 2009, which derives the correction term in special cases). This package does not implement any bias correction. You may want to check that your residuals are approximately i.i.d.


## Changes

### v0.6 (2026-03-29)
Fix cluster bootstrap. In previous versions, `vce(bootstrap, cluster())` silently ignored the `absorb()` option and used the wrong panel id for resampling when variable names were abbreviated. Point estimates from `regife` without bootstrap were unaffected. Standard errors from `vce(bootstrap)` without `cluster()` were also unaffected.

### v0.4 (2021-09-01)
Remove error when N < T. Preserve tsset.

### v0.3 (2017-04-12)
Correct weight handling.

### v0.2 (2015-07-09)
Correct normalization for loadings.

### v0.1 (2015-07-08)
First release.


## References

Bai, Jushan. *Panel Data Models with Interactive Fixed Effects.* Econometrica, 2009.


## Citation

```
Matthieu Gomez, 2015. "REGIFE: Stata module to estimate linear models with interactive fixed effects."
Statistical Software Components, Boston College Department of Economics.
https://ideas.repec.org/c/boc/bocode/s458042.html
```
