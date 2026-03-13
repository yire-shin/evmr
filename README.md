evmr: Extreme Value Modeling for r-Largest Order Statistics
================

# evmr

`evmr` is an R package for **extreme value modeling using the r-largest
order statistics framework**.

The package provides tools for fitting, analyzing, and comparing extreme
value models that use the **largest r observations within each block**.

The package is designed for applications in

- hydrology
- climatology
- environmental science
- extreme risk analysis

where multiple extreme events may occur within the same block.

------------------------------------------------------------------------

# Installation

``` r
# install.packages("remotes")
remotes::install_github("yire-shin/evmr")
```

------------------------------------------------------------------------

# Workflow

A typical workflow in `evmr` is

``` text
random generation → model fitting → return level estimation → profile likelihood → r selection
```

For each supported model, the package provides a consistent set of
functions:

| function | description |
|----|----|
| random generator | simulate r-largest order statistics |
| `.fit()` | fit the model by maximum likelihood |
| `.rl()` | estimate return levels |
| `.prof()` | obtain profile likelihood confidence intervals |
| `Edtest()` | perform entropy difference based sequential testing for selecting r |

------------------------------------------------------------------------

# rK4D Model

The **rK4D model** is based on the **four-parameter Kappa
distribution**.

### Random generation

``` r
x <- rk4dr(
  n = 50, r = 3,
  loc = 10, scale = 2,
  shape1 = 0.1, shape2 = 0.1
)

head(x$rmat)
```

### Model fitting

``` r
fit <- rk4d.fit(x$rmat, num_inits = 5)
fit$mle
```

### Return levels

``` r
rk4d.rl(fit)
```

### Profile likelihood confidence intervals

``` r
rk4d.prof(fit, m = 100, xlow = 12, xup = 25)
```

### Entropy difference test

``` r
rk4dEdtest(x$rmat)
```

### Reference

Shin, Y., & Park, J.-S. (2023).  
Modeling climate extremes using the four-parameter kappa distribution
for r-largest order statistics.  
*Weather and Climate Extremes*.  
<https://doi.org/10.1016/j.wace.2022.100533>

Hosking, J. R. M. (1994).  
*The four-parameter Kappa distribution*. Cambridge University Press.

Martins, E. S., & Stedinger, J. R. (2000).  
Generalized maximum-likelihood generalized extreme-value quantile
estimators for hydrologic data.  
*Water Resources Research*.  
<https://doi.org/10.1029/1999WR900330>

Coles, S., & Dixon, M. (1999).  
Likelihood-based inference for extreme value models.  
*Extremes*.  
<https://doi.org/10.1023/A:1009905222644>

------------------------------------------------------------------------

# rGLO Model

The **rGLO model** is based on the **generalized logistic
distribution**.

### Random generation

``` r
x <- rglor(
  n = 50, r = 3,
  loc = 10, scale = 2,
  shape = 0.1
)

head(x$rmat)
```

### Model fitting

``` r
fit <- rglo.fit(x$rmat, num_inits = 5)
fit$mle
```

### Return levels

``` r
rglo.rl(fit)
```

### Profile likelihood confidence intervals

``` r
rglo.prof(fit, m = 100, xlow = 12, xup = 25)
```

### Entropy difference test

``` r
rgloEdtest(x$rmat)
```

### Reference

Shin, Y., & Park, J.-S. (2024).  
Generalized logistic model for r-largest order statistics with
hydrological application.  
*Stochastic Environmental Research and Risk Assessment*.  
<https://doi.org/10.1007/s00477-023-02642-7>

------------------------------------------------------------------------

# rGGD Model

The **rGGD model** is based on the **generalized Gumbel distribution**.

### Random generation

``` r
x <- rggdr(
  n = 50, r = 3,
  loc = 10, scale = 2,
  shape = 0.1
)

head(x$rmat)
```

### Model fitting

``` r
fit <- rggd.fit(x$rmat, num_inits = 5)
fit$mle
```

### Return levels

``` r
rggd.rl(fit)
```

### Profile likelihood confidence intervals

``` r
rggd.prof(fit, m = 100, xlow = 12, xup = 25)
```

### Entropy difference test

``` r
rggdEdtest(x$rmat)
```

### Reference

Shin, Y., & Park, J.-S. (2025).  
Generalized Gumbel model for r-largest order statistics with application
to peak streamflow.  
*Scientific Reports*.  
<https://doi.org/10.1038/s41598-024-83273-y>

Jeong et al. (2014).  
A three-parameter kappa distribution with hydrologic application: A
generalized Gumbel distribution.  
<https://doi.org/10.1007/s00477-014-0865-8>

------------------------------------------------------------------------

# rGD Model

The **rGD model** is based on the **Gumbel distribution**.

### Random generation

``` r
x <- rgdr(
  n = 50, r = 3,
  loc = 10, scale = 2
)

head(x$rmat)
```

### Model fitting

``` r
fit <- rgd.fit(x$rmat, num_inits = 5)
fit$mle
```

### Return levels

``` r
rgd.rl(fit)
```

### Profile likelihood confidence intervals

``` r
rgd.prof(fit, m = 100, xlow = 12, xup = 25)
```

### Entropy difference test

``` r
rgdEdtest(x$rmat)
```

### Reference

Coles, S., & Dixon, M. (1999).  
Likelihood-based inference for extreme value models.  
*Extremes*.  
<https://doi.org/10.1023/A:1009905222644>

------------------------------------------------------------------------

# rLD Model

The **rLD model** is based on the **logistic distribution**.

### Random generation

``` r
x <- rldr(
  n = 50, r = 3,
  loc = 10, scale = 2
)

head(x$rmat)
```

### Model fitting

``` r
fit <- rld.fit(x$rmat, num_inits = 5)
fit$mle
```

### Return levels

``` r
rld.rl(fit)
```

### Profile likelihood confidence intervals

``` r
rld.prof(fit, m = 100, xlow = 12, xup = 25)
```

### Entropy difference test

``` r
rldEdtest(x$rmat)
```

### Reference

Ahmad et al. (1988).  
Log-logistic flood frequency analysis.  
*Journal of Hydrology*.  
<https://doi.org/10.1016/0022-1694(88)90015-7>

------------------------------------------------------------------------

# Unified Model Comparison Using Real Data

The package also provides a unified wrapper function, `evmr()`, for
fitting multiple models simultaneously.

The following example uses the built-in `bangkok` dataset.

``` r
library(evmr)

data(bangkok)

evmr(
  bangkok,
  models = c("rk4d", "rglo", "rggd", "rgd", "rld"),
  num_inits = 5
)
```

This returns a combined summary table containing

- model name
- parameter estimates
- likelihood values
- standard errors
- return level estimates

------------------------------------------------------------------------

# Choosing r

In practice, the value of `r` is usually unknown.

Small values of `r` may waste information, while overly large values may
introduce bias.

The package provides entropy difference based sequential testing
procedures through the `.edtest()` functions to help determine a
reasonable value of `r`.

Reference:

Bader, B., Yan, J., & Zhang, X. (2017).  
Automated selection of r for the r-largest order statistics approach.  
*Statistics and Computing*.  
<https://doi.org/10.1007/s11222-016-9697-3>

------------------------------------------------------------------------

# References

Bader, B., Yan, J., & Zhang, X. (2017).  
Automated selection of r for the r-largest order statistics approach.  
*Statistics and Computing*.  
<https://doi.org/10.1007/s11222-016-9697-3>

Coles, S., & Dixon, M. (1999).  
Likelihood-based inference for extreme value models.  
*Extremes*.  
<https://doi.org/10.1023/A:1009905222644>

Hosking, J. R. M. (1994).  
*The four-parameter Kappa distribution*. Cambridge University Press.

Martins, E. S., & Stedinger, J. R. (2000).  
Generalized maximum-likelihood generalized extreme-value quantile
estimators for hydrologic data.  
*Water Resources Research*.  
<https://doi.org/10.1029/1999WR900330>

Ahmad et al. (1988).  
Log-logistic flood frequency analysis.  
*Journal of Hydrology*.  
<https://doi.org/10.1016/0022-1694(88)90015-7>

Jeong et al. (2014).  
A three-parameter kappa distribution with hydrologic application: A
generalized Gumbel distribution.  
<https://doi.org/10.1007/s00477-014-0865-8>

Shin, Y., & Park, J.-S. (2023).  
Modeling climate extremes using the four-parameter kappa distribution
for r-largest order statistics.  
*Weather and Climate Extremes*.  
<https://doi.org/10.1016/j.wace.2022.100533>

Shin, Y., & Park, J.-S. (2024).  
Generalized logistic model for r-largest order statistics with
hydrological application.  
*Stochastic Environmental Research and Risk Assessment*.  
<https://doi.org/10.1007/s00477-023-02642-7>

Shin, Y., & Park, J.-S. (2025).  
Generalized Gumbel model for r-largest order statistics with application
to peak streamflow.  
*Scientific Reports*.  
<https://doi.org/10.1038/s41598-024-83273-y>
