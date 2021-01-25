\name{glme}
\alias{glme}
\title{Generalized Linear Mixed Effects Models}
\usage{
glme(fixed, data, random, correlation, weights, subset,
     method, na.action, control, contrasts, keep.data)
}
\description{
This function fits a linear mixed effect model with generalized inference.
}
\arguments{
  \item{fixed}{a linear model formula, with the response on the left of a  operator and an expression involving parameters and covariates on the right.}
  \item{data}{an optional data frame containing the variables named in model, fixed, random, correlation, weights, subset, and naPattern. By default the variables are taken from the environment from which glme is called.}
  \item{random}{a two-sided linear formula of the form \code{f1+...+fn ~ x1+...+xm}, or a list of two-sided formulas of the form \code{f1 ~ x1+...+xm}, with possibly different models for different parameters. The \code{f1,...,fn} are the names of parameters included on the right hand side of model and the \code{x1+...+xm} expressions define linear models for these parameters. On the right hand side of the formula(s) indicates a single fixed effects for the corresponding parameter(s).}
  \item{correlation}{an optional corStruct object describing the within-group correlation structure}
  \item{weights}{an optional varFunc object or one-sided formula describing the within-group heteroscedasticity structure.}
  \item{subset}{an optional expression indicating the subset of the rows of data that should be used in the fit. This can be a logical vector, or a numeric vector indicating which observation numbers are to be included, or a character vector of the row names to be included. All observations are included by default.}
  \item{method}{a character string. If "GM" the model is fit by generalized inference. If "REML" the model is fit by maximizing the restricted log-likelihood. If "ML" the log-likelihood is maximized. Defaults to "GM".}
  \item{na.action}{a function that indicates what should happen when the data contain NAs.}
  \item{control}{a list of control values for the estimation algorithm to replace the default values returned.}
  \item{contrasts}{an optional list. See the contrasts.arg of model.matrix.default.}
  \item{keep.data}{logical: should the \code{data} argument (if supplied
   and a data frame) be saved as part of the model object.}
}

\value{
  \item{fixed}{returns the coefficient estimations and model summary of the fixed part.}
  \item{sd}{returns the standard deviation of random effects.}
  \item{coefficients}{returns the coefficient estimations of the fixed and random part of the mixed model.}
}

\author{
Sam Weerahandi, Berna Yazici, Ching-Ray Yu, Mustafa Cavus
}


\references{
Yu, C.R., Kelly H.Z., Carlsson, M.O., and Weerahandi, S. (2015) Generalized Estimation of the BLUP in Mixed-Effects Models: A Comparison with ML and REML, Communications in Statistics - Simulation and Computation, 44:3, 694-704, https://doi.org/10.1080/03610918.2013.790445

Weerahandi, S. and Yu, CR. (2020) Exact distributions of statistics for making inferences on mixed models under the default covariance structure. Journal of Statistical Distributions and Applications, 7:4, https://doi.org/10.1186/s40488-020-00105-w

Gamage, J., Mathew, T., and Weerahandi, S. (2013) Generalized prediction intervals for BLUPs in mixed models, Journal of Multivariate Analysis, 120, 226 - 233, https://doi.org/10.1016/j.jmva.2013.05.011.
}

\examples{
library(nlme)
library(glme)
glme(distance ~ age + Sex, data = Orthodont, random = ~ age|Subject, method = "GM")
}
