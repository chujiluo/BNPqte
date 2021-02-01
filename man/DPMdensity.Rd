\name{DPMdensity}
\alias{DPMdensity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Semiparametric Bayesian density estimation using a DPM of multivariate Normal
}
\description{
Create a (hierarchical) Dirichlet process mixture of conjugate multivariate Normal distribtions. Posterior is sampled using blocked Gibbs sampling (Ishwaran and James, \emph{JASA}, 2001).
}
\usage{
DPMdensity(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
  a matrix giving the data from which the density estimate is to be computed.
  }
  \item{ngrid}{
  number of grid points where the density estimate is evaluated. This is only used if the dimension of \eqn{y} is equal to \eqn{2}. By default, \eqn{ngrid=1000L}.
  }
  \item{grid}{
  matrix of dimension \eqn{ngrid*d} of grid points where the density estimate is evaluated, where \eqn{d=ncol(y)}. This is only used if the dimension of \eqn{y} is equal to \eqn{2}. By default, \eqn{grid=NULL} and the \eqn{grid} is chosen according to the range of the data.
  }
  \item{updateAlhpa}{
  a logical variable indicating whether the concentration parameter \eqn{alpha} is updated (\eqn{TRUE}) or not (\eqn{FALSE}). By default, \eqn{updateAlpha=TRUE}.
  }
  \item{useHyperpriors}{
  a logical variable indicating whether hyperpriors are assigned to the parameters of the base distribution (Normal-Inverse-Wishart) (\eqn{TRUE}) or not (\eqn{FALSE}). By default, \eqn{useHyperpriors=TRUE}.
  }
  \item{nclusters}{
  number of clusters used in the blocked Gibbs sampling. By default, \eqn{nclusters=50L}.
  }
  \item{nskip}{
  number of initial iterations of the MCMC to be discarded. By default, \eqn{nskip=1000L}.
  }
  \item{ndpost}{
  number of posterior samples returned to the user. \eqn{ndpost*keepevery} iterations of the MCMC after the burn-in period. By default, \eqn{ndpost=1000L}.
  }
  \item{keepevery}{
  every \eqn{keepevery} posterior sample is kept to be returned to the user. By default, \eqn{keepevery=1L}.
  }
  \item{state}{
  a list giving the current value of the parameters (and hyperparameters). This list is obtained from the returns of this function and used if the current analysis is the continuation of a previous analysis. By default, \eqn{state=NULL}.
  }
  \item{status}{
  a logical variable indicating whether this run is new (\eqn{status=TRUE}) or the continuation of a previous analysis (\eqn{status=FALSE}). In the latter case, the current value of the parameters (and hyperparameters) must be specified in the object \eqn{state}.
  }
  \item{alpha}{
  a positive real number giving the concentration parameter of the Dirichlet process \eqn{DP(alpha, G0)}. This variable is fixed if \eqn{updateAlpha=FALSE}, and random with the prior \eqn{Gamma(alpha | a0, b0)} otherwise. By default, \eqn{alpha=2.0}.
  }
  \item{a0}{
  a positive real number giving the shape parameter of the Gamma prior for the concentration parameter \eqn{alpha}, i.e. \eqn{alpha ~ Gamma(a0, b0)}. This is only used if \eqn{updateAlpha=TRUE}. By default, \eqn{alpha=2.0}.
  }
  \item{b0}{
  a positive real number giving the rate parameter of the Gamma prior for the concentration parameter \eqn{alpha}, i.e. \eqn{alpha ~ Gamma(a0, b0)}. This is only used if \eqn{updateAlpha=TRUE}. By default, \eqn{b0=1.0}.
  }
  \item{m}{
  a vector giving the mean parameter of the base distribution, i.e. \eqn{(zeta, Omega) ~ Normal(zeta | m, Omega/lambda) x Inverse-Wishart(Omega | nu, Psi)}. This parameter is fixed if \eqn{useHyperpriors=FALSE} and random with a hyperprior \eqn{Normal(m | m0, S0)} otherwise. By default, \eqn{m=colMeans(y)}.
  }
  \item{m0}{
  a vector giving mean hyperparameter of the Normal hyperprior for the mean parameter \eqn{m}, i.e. \eqn{m ~ Normal(m0, S0)}. This is only used if \eqn{useHyperpriors=TRUE}. By default, \eqn{m0=colMeans(y)}.
  }
  \item{S0}{
  a positive definite matrix giving the covariance hyperparameter of the Normal hyperprior for the mean parameter \eqn{m}, i.e. \eqn{m ~ Normal(m0, S0)}. This is only used if \eqn{useHyperpriors=TRUE}. If nothing is specified for \eqn{S0}, then \eqn{S0} is a diagonal matrix with diagonal elements \eqn{(r_{y_j}/4)^2} (see Taddy 2008 for details).
  }
  \item{lambda}{
  a positive real number giving the scailing parameter of the base distribution, i.e. \eqn{(zeta, Omega) ~ Normal(zeta | m, Omega/lambda) x Inverse-Wishart(Omega | nu, Psi)}. This parameter is fixed if \eqn{useHyperpriors=FALSE} and random with a hyperprior \eqn{Gamma(lambda | gamma1, gamma2)} otherwise. By default, \eqn{lambda=0.5}.
  }
  \item{gamma1}{
  a positive real number giving the shape hyperparameter of the Gamma hyperprior for \eqn{lambda}, i.e. \eqn{lambda ~ Gamma(gamma1, gamma2)}. This is only used if \eqn{useHyperpriors=TRUE}. By default, \eqn{gamma1=3.0}.
  }
  \item{gamma2}{
  a positive real number giving the rate hyperparameter of the Gamma hyperprior for \eqn{lambda}, i.e. \eqn{lambda ~ Gamma(gamma1, gamma2)}. This is only used if \eqn{useHyperpriors=TRUE}. By default, \eqn{gamma2=2.0}.
  }
  \item{nu}{
  degrees of freedom \eqn{(> ncol(y)-1)} for the base \eqn{Inverse-Wishart(Omega | nu, Psi)} distribution. By default, \eqn{nu=ncol(y)+2}.
  }
  \item{Psi}{
  a positive definite matrix giving the scale matrix parameter for the base \eqn{Inverse-Wishart(Omega | nu, Psi)} distribution. This parameter is fixed if \eqn{useHyperpriors=FALSE} and random with a hyperprior \eqn{Wishart(Psi | nu0, Psi0)}. If nothing is specified for \eqn{Psi}, then \eqn{Psi} is a diagonal matrix with diagonal elements \eqn{(r_{y_j}/4)^2} (see Taddy 2008 for details).
  }
  \item{nu0}{
  degrees of freedom \eqn{(> ncol(y)-1)} of the hyperprior distribution for \eqn{Psi}, i.e. \eqn{Psi ~ Wishart(nu0, Psi0)}. This is only used if \eqn{useHyperpriors=TRUE}. By default, \eqn{nu0=ncol(y)+2}.
  }
  \item{Psi0}{
  a positive definite matrix giving the scale matrix hyperparameter for \eqn{Psi}, i.e. \eqn{Psi ~ Wishart(nu0, Psi0)}. This is only used if \eqn{useHyperpriors=TRUE}. If nothing is specified for \eqn{Psi0}, then \eqn{Psi0=S0/nu0}. (see Taddy 2008 for details).
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
