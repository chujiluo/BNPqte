# BNPqte

This **R** package implements the paper, Xu D, Daniels MJ, Winterstein AG.
*A Bayesian nonparametric approach to causal inference on quantiles.* Biometrics. 2018;74(3):986-996. doi:10.1111/biom.12863,
which proposes a Bayesian nonparametric approach to estimate multiple quantile casual effects simultaneously. In specific, 
the p<sup>th</sup> quantile causal effect is defined as
<p align="center">
  <img src="https://render.githubusercontent.com/render/math?math=F_1^{-1}(p) - F_0^{-1}(p),">
</p>
and p = F<sub>t</sub>(y) can be identified as
<p align="center">
  <img src="https://render.githubusercontent.com/render/math?math=F_t(y) = \int P(Y<y|e(X)=e(x),T=t)dG(e(x)),">
</p>
where Y is the continuous outcome, X are the covariates, T is the binary treatment and e(X) is the propensity score (PS).
They estimate the PS e(X) by fitting a probit Bayesian Additive Regression Trees (BART) model on T ~ X. 
For each posterior sample of the PS e<sup>{k}</sup>, they estimate the conditional distribution of Y | e<sup>{k}</sup>, T = t
by fitting a weight dependent Dirichlet process mixture model (DPMM) of multivariate Normals on (Y, e<sup>{k}</sup>)
for the treatment group (T = 1) and the control group (T = 0) separately.
Finally, the unknown distribution of the PS is approximated using Bayesian bootstrap (Rubin, 1981).
Hence, the cumulative distribution function F<sub>t</sub>(y) is estimated by
<p align="center">
  <img src="https://render.githubusercontent.com/render/math?math=\hat{F}_t(y) = \frac{1}{KL} \sum_{k=1}^K \left\{ \sum_{i=1}^n u_i^k \left[ \sum_{l=1}^L \hat{F}^{\{kl\}}(y | e^{\{k\}} (x_i), T=t) \right] \right\}">
</p>
where k is the indices of the BART posterior samples, l is the indices of the DPMM posterior samples and (u<sub>1</sub>, ..., u<sub>n</sub>) are samples from
Dirichlet(1, ..., 1) distribution.

Based on the context of the paper, this package has the following functions:

1. `wbart`, `pbart` and `lbart` for regular BART, probit BART and logit BART, which are inheritted from the CRAN **R** package **BART** but have two modifications

    1.1. add an alternative splitting probability <img src="https://render.githubusercontent.com/render/math?math=p(d)=\alpha^d"> where *d* is the depth of the node;
  this exponential probability is from Rockova and Saha (2019), which ensures the optimal posterior concentration up to a log factor;
  
    1.2. adjust the variable inclusion proportions to mixed-type covariates and add Metropolis importance, according to Luo and Daniels (2020);
  
2. `DPMdensity` and `DPMcdensity` for joint density and conditional density estimation using DPMM of multivariate Normals,
and the posterior distribution of the DPMM is sampled using blocked Gibbs sampling introduced in Ishwaran (2000) and Ishwaran and James (2001)；

3. `qte` for quantile casual effects estimation; parallel computing for fitting multiple DPMMs is available in `qte` by setting `mc.cores` greater than 1;

4. `predict` and `plot` S3 Methods are available for functions `DPMdensity`, `DPMcdensity` and `qte`.

To increase the computational efficiency, the actual computation in **BNPqte** package is carried out in C++ and **RcppArmadillo**. In addition,
the **BNPqte** package, inheritted from **BART** package, takes advantage of multi-threading via forking as provided by the **parallel** package and OpenMP
when available and supported by the platform.
