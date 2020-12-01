# copula.est

<!-- badges: start -->

<!-- badges: end -->

The goal of R function copula.est is to estimate copula with or without PQD constraint using Bernstein polynomials.

## Description 

This function implements copula estimation methods developed in the paper:

Lu, Lu, and Sujit K. Ghosh. *“Nonparametric Estimation and Testing for Positive Quadrant Dependent Bivariate Copula.”* Journal of Business & Economic Statistics
https://doi.org/10.1080/07350015.2020.1855186

## Usage

``` r
copula.est(X, m1 = NULL, m2 = NULL, is.pqd = F)
```

## Arguments

**X**: a n*2 matrix of data.

**m1, m2**: degrees of Bernstein polynomial estimated by grid search if m1 = NULL and m2 = NULL.

**is.pqd**: is.pqd = F for unconstrained copula; is.pqd = T for PQD-constrained copula. The default is is.pqd = F.


## Value

**theta.matrix**: A matrix of estimated theta_(k1, k2), k1 = 1, ..., m1-1, k2 = 1, ..., m2-1.

**m1**: estimated degree m1.

**m2**: estimated degree m2.

## Note

This function also outputs a contour plot of the estimated coupla.

## Examples

``` r
# install.packages("copula")
library(copula)
set.seed(123)
# sample n=100 iid data from bivarate Frank copula with theta.true equal to -2
n=100
theta.true <- -2
fr <- frankCopula(theta.true, 2)
X <- rCopula(n, fr)

copula.est(X, m1 = NULL, m2 = NULL, is.pqd = F)

$theta.matrix
              [,1]         [,2]          [,3]
 [1,] 4.741365e-19 4.992646e-18 -7.086946e-18
 [2,] 1.397889e-03 1.397889e-03  3.759235e-02
 [3,] 1.397889e-03 1.397889e-03  1.145154e-01
 [4,] 1.397889e-03 1.397889e-03  1.914385e-01
 [5,] 1.397889e-03 1.397889e-03  2.513979e-01
 [6,] 1.397889e-03 4.814543e-02  2.981454e-01
 [7,] 4.590606e-02 1.250685e-01  3.750685e-01
 [8,] 1.228291e-01 2.019916e-01  4.519916e-01
 [9,] 1.755439e-01 2.789147e-01  5.289147e-01
[10,] 1.755439e-01 3.206326e-01  5.706326e-01
[11,] 1.755439e-01 3.975557e-01  6.475557e-01
[12,] 1.944408e-01 4.444408e-01  6.944408e-01

$m1
[1] 13

$m2
[1] 4
```
