#########################################
# ---- Functions for setting classes ----
#########################################
# Content:
#
# Set the following classes:
#  - "penfaData",
#  - "penfaSampleStats"
#  - "penfaModel"
#  - "penfaPenalty"
#  - "penfa"
#
# ------------------------------------------------------------------------------

#' S4 Class for describing the input data
#'
#' @description The \code{penfaData} class gives information on the data set
#'   provided in input for analysis. This class is an adaptation of the
#'   \code{lavData} class from the
#'   [lavaan](https://CRAN.R-project.org/package=lavaan) package.
#'
#' @slot ngroups Integer. The number of groups.
#' @slot group Character. The observed variables defining the groups.
#' @slot group.label Character. The group labels, that is, the values of the
#'   \code{group} variable, if any.
#' @slot std.ov Logical indicating whether the observed variables should be
#'   standardized.
#' @slot nobs List of the effective number of observations in each group.
#' @slot norig List of the original number of observations in each group.
#' @slot ov.names List of the observed variable names in each group.
#' @slot ov List of details at the observed variable level.
#' @slot case.idx List of the case (i.e., observation) indices in each group.
#' @slot X List. Local copy of the input data set split into groups.
#'
#' @seealso \code{\link{penfa}}
#'
#'
#'
#' @examples
#'
#' data(ccdata)
#'
#' syntax = 'help  =~   h1 + h2 + h3 + h4 + h5 + h6 + h7 + 0*v1 + v2 + v3 + v4 + v5
#'           voice =~ 0*h1 + h2 + h3 + h4 + h5 + h6 + h7 +   v1 + v2 + v3 + v4 + v5'
#'
#' alasso_fit <- penfa(## factor model
#'                     model  = syntax,
#'                     data   = ccdata,
#'                     std.lv = TRUE,
#'                     ## penalization
#'                     pen.shrink = "alasso",
#'                     eta = list(shrink = c("lambda" = 0.01), diff = c("none" = 0)),
#'                     ## automatic procedure
#'                     strategy = "auto")
#'
#' alasso_fit@Data
#' str(alasso_fit@Data)
#'
#'
setClass("penfaData",
         slots = c(
           ngroups     = "integer",
           group       = "character",
           group.label = "character",
           std.ov      = "logical",
           nobs        = "list",
           norig       = "list",
           ov.names    = "list",
           ov          = "list",
           case.idx    = "list",
           X           = "list"))


#' S4 Class for describing the sample moments
#'
#' @description The \code{penfaSampleStats} class provides information on the
#'   sample moments of the factor analysis model. This
#'   class is an adaptation of the \code{lavSampleStats} class from the
#'   [lavaan](https://CRAN.R-project.org/package=lavaan) package.
#'
#' @slot var List of the variances of the observed variables in every group.
#' @slot cov List of the covariance matrices of the observed variables
#'   in every group.
#' @slot mean List of the means of the observed variables in every group.
#' @slot group.w List of group weights.
#' @slot nobs List of the effective number of observations for every group.
#' @slot ntotal Integer. Total number of observations across all groups.
#' @slot ngroups Integer. Number of groups.
#' @slot icov List of the inverse matrices of the covariance matrices of the
#'   observed variables in every group.
#' @slot cov.log.det List of the logarithms of the determinants of the
#'   covariance matrices of the observed variables for every group.
#'
#' @seealso \code{\link{penfa}}
#'
#'
#'
#'
#' @examples
#'
#' data(ccdata)
#'
#' syntax = 'help  =~   h1 + h2 + h3 + h4 + h5 + h6 + h7 + 0*v1 + v2 + v3 + v4 + v5
#'           voice =~ 0*h1 + h2 + h3 + h4 + h5 + h6 + h7 +   v1 + v2 + v3 + v4 + v5'
#'
#' alasso_fit <- penfa(## factor model
#'                     model  = syntax,
#'                     data   = ccdata,
#'                     std.lv = TRUE,
#'                     ## penalization
#'                     pen.shrink = "alasso",
#'                     eta = list(shrink = c("lambda" = 0.01), diff = c("none" = 0)),
#'                     ## automatic procedure
#'                     strategy = "auto")
#'
#' alasso_fit@SampleStats
#'
setClass("penfaSampleStats",
         slots = c(
           var         = "list",
           cov         = "list",
           mean        = "list",
           group.w     = "list",
           nobs        = "list",
           ntotal      = "integer",
           ngroups     = "integer",
           icov        = "list",
           cov.log.det = "list"))


#' S4 Class for internal representation of a factor model
#'
#' @description The \code{penfaModel} class gives the internal matrix
#'   representation of a factor analysis model. Note that this representation
#'   summarizes the characteristics of the model itself (e.g., number of items,
#'   number of factors, parameter indices, etc), without information on the
#'   penalization process (see \code{\linkS4class{penfaPenalty}} for that
#'   aspect). This class is an adaptation of the \code{lavModel} class from the
#'   [lavaan](https://CRAN.R-project.org/package=lavaan) package.
#'
#'
#' @slot GLIST List. The model matrices and vectors: "lambda" for the factor
#'   loading matrix, "psi" for the covariance matrix of the unique factors,
#'   "phi" for the covariance matrix of the common factors, "tau" for the
#'   intercept vector, and "kappa" for the vector of factor means. In case of a
#'   multiple-group analysis, the elements of each group are presented
#'   sequentially.
#' @slot dimNames List. Dimension names (row names and column names) of every
#'   model matrix and vector.
#' @slot isSymmetric Logical vector declaring whether each model matrix/vector
#'   is symmetric.
#' @slot mmSize Integer vector specifying the size (unique elements only) of
#'   each model matrix/vector.
#' @slot meanstructure Logical. It declares whether the model includes a
#'   meanstructure.
#' @slot ngroups Integer. The number of groups.
#' @slot nmat Integer vector specifying the number of model matrices/vectors for
#'   each group.
#' @slot nvar Integer vector specifying the number of observed variables in each
#'   group.
#' @slot num.idx List of the indices of the observed variables in each group.
#' @slot nx.free Integer. The number of parameters of the factor model. This
#'   count does not include the fixed parameters, but it does include the
#'   parameters that will be penalized (if any) during optimization. (see
#'   \code{\linkS4class{penfaPenalty}} for additional details in this respect).
#' @slot nx.user Integer. The total count of the parameters that are being
#'   estimated and the ones that have been fixed.
#' @slot m.free.idx List. For each model matrix, the indices of the elements to
#'   be estimated (i.e., non-fixed). The counter starts at 1 for every model
#'   matrix.
#' @slot x.free.idx List. For each model matrix, the indices of the elements to
#'   be estimated (i.e., non-fixed). The counter continues from the previous
#'   model matrix.
#' @slot m.user.idx List. Much like \code{m.free.idx}, but it also contains the
#'   indices of the parameters that have been fixed by the user.
#' @slot x.user.idx List. Much like \code{x.free.idx}, but it also contains the
#'   indices of the parameters that have been fixed by the user.
#' @slot x.free.var.idx Vector of integers denoting the indices corresponding to
#'   the unique variances.
#'
#' @seealso \code{\link{penfa}}
#'
#'
#'
#' @examples
#'
#' data(ccdata)
#'
#' syntax = 'help  =~   h1 + h2 + h3 + h4 + h5 + h6 + h7 + 0*v1 + v2 + v3 + v4 + v5
#'           voice =~ 0*h1 + h2 + h3 + h4 + h5 + h6 + h7 +   v1 + v2 + v3 + v4 + v5'
#'
#' alasso_fit <- penfa(## factor model
#'                     model  = syntax,
#'                     data   = ccdata,
#'                     std.lv = TRUE,
#'                     ## penalization
#'                     pen.shrink = "alasso",
#'                     eta = list(shrink = c("lambda" = 0.01), diff = c("none" = 0)),
#'                     ## automatic procedure
#'                     strategy = "auto")
#'
#' alasso_fit@Model
#'
setClass("penfaModel",
         slots = c(
           GLIST          = "list",
           dimNames       = "list",
           isSymmetric    = "logical",
           mmSize         = "integer",
           meanstructure  = "logical",
           ngroups        = "integer",
           nmat           = "integer",
           nvar           = "integer",
           num.idx        = "list",
           nx.free        = "integer",
           nx.user        = "integer",
           m.free.idx     = "list",
           x.free.idx     = "list",
           m.user.idx     = "list",
           x.user.idx     = "list",
           x.free.var.idx = "integer"))


#' S4 Class for describing the penalization process
#'
#' @description The \code{penfaPenalty} class provides information on the
#'   penalization process, such as the user-specified penalty functions, the
#'   optimal values of the tuning parameters, and the penalty matrices at
#'   convergence.
#'
#' @slot strategy Character. The strategy used for the selection of the tuning
#'   parameter(s). If \code{strategy = "auto"}, the optimal values of the tuning
#'   parameters are determined via the automatic tuning parameter procedure; if
#'   \code{strategy = "fixed"}, a penalized factor model with the values of the
#'   tuning parameters stored in the option \code{eta} is estimated.
#' @slot penalty List. A list of the user-specified penalty functions for
#'   sparsity ("shrink") and parameter equivalence ("diff").
#' @slot tuning List. A named list containing the optimal values of the tuning
#'   parameter(s) if \code{strategy = "auto"} or the user-specified fixed
#'   values of the tuning parameter(s) if \code{strategy = "fixed"}. The list
#'   has two components with names "shrink" and "diff", and refers to the tuning
#'   parameters used for shrinkage and group equivalence, respectively. The
#'   components of the list are, in turn, the named vectors specifying the type
#'   of parameter matrices or vectors that were penalized.
#' @slot pmat List. A named list containing the names of the parameter matrices
#'   and vectors that were penalized for sparsity ("shrink") and/or group
#'   equivalence ("diff").
#' @slot pen.idx List. A named list with the indices of the parameters that were
#'   penalized for sparsity ("shrink") and/or group equivalence ("diff").
#' @slot Sh.info List. A list of the penalization terms, vectors and matrices
#'   evaluated at the optimal values of the tuning parameters. In particular,
#'   its argument \code{S.h} returns the estimated penalty matrix. If the factor
#'   model is penalized only through a shrinkage penalty (i.e.,
#'   \code{pen.shrink} is not \code{'none'}), and there is no penalization on
#'   the differences (i.e., \code{pen.diff = 'none'}), then \code{S.h} is a
#'   diagonal matrix whose elements precisely quantify the extent to which each
#'   model parameter has been penalized.
#' @slot extra List. A list possibly containing additional information on the
#'   penalization process, such as the hyperparameter values for some penalty
#'   functions (e.g., for the alasso, the value of the exponent and the adaptive
#'   weights.)
#' @slot automatic List. If \code{strategy = "auto"}, it contains information on
#'   the automatic multiple tuning parameter procedure, such as the optimal
#'   values of the tuning parameters, the convergence status, the specified
#'   value of the influence factor, the number of necessary iterations, and the
#'   tolerance level.
#'
#' @seealso \code{\link{penfa}}
#'
#'
#'
#' @examples
#'
#' data(ccdata)
#'
#' syntax = 'help  =~   h1 + h2 + h3 + h4 + h5 + h6 + h7 + 0*v1 + v2 + v3 + v4 + v5
#'           voice =~ 0*h1 + h2 + h3 + h4 + h5 + h6 + h7 +   v1 + v2 + v3 + v4 + v5'
#'
#' alasso_fit <- penfa(## factor model
#'                     model  = syntax,
#'                     data   = ccdata,
#'                     std.lv = TRUE,
#'                     ## penalization
#'                     pen.shrink = "alasso",
#'                     eta = list(shrink = c("lambda" = 0.01), diff = c("none" = 0)),
#'                     ## automatic procedure
#'                     strategy = "auto")
#'
#' alasso_fit@Penalize
#'
#' str(alasso_fit@Penalize)
#'
setClass("penfaPenalty",
         slots = c(
           strategy  = "character",
           penalty   = "list",
           tuning    = "list",
           pmat      = "list",
           pen.idx   = "list",
           Sh.info   = "list",
           extra     = "list",
           automatic = "list"))


#' S4 Class for describing a \code{penfa} model
#'
#' @description The \code{penfa} class represents a (fitted) penalized factor
#'   analysis model. It contains a description of the model as specified by the
#'   user, a summary of the data, an internal matrix representation, the
#'   fitting results, and the penalized quantities.
#'
#' @section Objects from the Class:
#'
#' Objects can be created via the \code{\link{penfa}} function.
#'
#' @section Slots:
#'
#' \describe{
#'  \item{\code{version}:}{The \code{penfa} package version used to create this object.}
#'  \item{\code{call}:}{The function call as returned by \code{match.call()}.}
#'  \item{\code{timing}:}{The elapsed time (user + system) for various parts of
#'  the program as a list, including the total time.}
#'  \item{\code{Options}:}{Named list of options that were provided by the user
#'  or filled-in automatically. See \code{\link{penfaOptions}} for additional
#'  details.}
#'  \item{\code{ParTable}:}{Named list describing the model parameters.
#'  Can be coerced to a data.frame. This is also called "parameter table".
#'  It includes information on the fixed, free and penalized parameters, their
#'  indices, the active penalization strategies ("none", "shrink", "diff", or
#'  "shrink + diff"), the starting values, the estimated parameters and the
#'  associated standard errors.}
#'  \item{\code{pta}:}{Named list containing parameter table attributes, like
#'  observed and latent variable names, their indices, and the number of groups.}
#'  \item{\code{Data}:}{Object of internal class \code{"penfaData"}; contains
#'  information about the data set. See the \code{\linkS4class{penfaData}} class
#'  for additional details.}
#'  \item{\code{SampleStats}:}{Object of internal class \code{"penfaSampleStats"};
#'  contains the sample statistics. See the \code{\linkS4class{penfaSampleStats}}
#'  class for additional details.}
#'  \item{\code{Model}:}{Object of internal class \code{"penfaModel"}: the internal
#'  (matrix) representation of the model. See the \code{\linkS4class{penfaModel}}
#'  class for additional details.}
#'  \item{\code{Optim}:}{List. Information about the optimization process. This
#'  includes the estimated parameters (\code{x}), the number of estimated
#'  parameters (\code{npar}), the number of trust-region iterations
#'  (\code{iterations}), the value of the penalized objective function
#'  (\code{fx.pen}), the value of the unpenalized objective
#'  function (\code{fx.unpen}), the penalized log-likelihood (\code{logl.pen};
#'  this is equal to \code{fx.pen} multiplied by (-1)), the unpenalized
#'  log-likelihood (\code{logl.unpen}; this is equal to \code{fx.unpen} multiplied
#'  by (-1)), the penalized gradient (\code{dx.pen}), the penalized Hessian/Fisher
#'  information matrix (\code{hessian.pen}), the list of control arguments for
#'  the trust-region algorithm (\code{control}), and how many times the
#'  objective function became non-positive definite during the estimation
#'  process (\code{npd}). If penfa was called with the option \code{verbose =
#'  TRUE}, the following additional arguments coming from the trust-region
#'  function \code{trust} are reported in the \code{Optim} slot: \code{argpath},
#'  \code{argtry}, \code{type}, \code{accept}, \code{radii}, \code{rho},
#'  \code{fx.val}, \code{fx.valtry}, \code{change}, \code{stepnorm}. See the
#'  manual page of \code{trust} from the \code{trust} package for an overview of
#'  these quantities.}
#'  \item{\code{Penalize}:}{Object of internal class \code{"penfaPenalty"}; contains
#'  information about the penalization. See the \code{\linkS4class{penfaPenalty}} for
#'  additional details.}
#'  \item{\code{Implied}:}{List. Model-implied moments (covariance matrix and
#'  mean vector).}
#'  \item{\code{Vcov}:}{List. Information about the covariance matrix (vcov) of
#'  the model parameters. This slot includes the following quantities:  the type
#'  of penalized information matrix used in the model (either Hessian or Fisher;
#'  \code{information}), the vcov matrix of parameters (\code{vcov}), whether
#'  the convergence checks on the penalized gradient and the penalized
#'  information matrix were satisfied (\code{solution}), whether the employed
#'  information matrix was positive-definite (\code{pdef}),  whether the
#'  estimated factor solution was admissible (\code{admissibility}), the
#'  standard errors computed according to the Bayesian result from the
#'  information matrix reported in \code{information} (\code{se}), and the 95%
#'  confidence intervals (\code{ci}).}
#'  \item{\code{Inference}:}{List. Information on effective degrees of the model
#'  and information criteria for model selection. This slot reports the
#'  following quantities: effective degree of freedom for each parameter
#'  (\code{edf.single}), total edf (\code{edf}), influence matrix
#'  (\code{influence.mat}), generalized information criteria (\code{IC}), such
#'  as AIC and BIC.}
#'  \item{\code{external}:}{List. Empty slot.}
#'}
#'
#' @section  Methods:
#'
#' The following methods are available for an object of class
#' \code{\linkS4class{penfa}}:
#'
#' \describe{
#' \item{show}{\code{signature(object = "penfa")}: Prints a short summary of the
#' estimation process, including the optimization method, the specified penalty
#' functions, the convergence status, the number of iterations, the tuning
#' selection strategy, and the effective degrees of freedom. See the manual page of
#' \code{show,penfa-method} for details. } \item{summary}{\code{signature(object
#' = "penfa", header = TRUE, estimates = TRUE, ci = TRUE, level =} \code{ 0.95,
#' nd = 3L, cutoff = 0.05, extra = TRUE)}: Prints a summary of the model
#' parameter estimates, and the optimization process. See the manual page of
#' \code{summary,penfa-method} for details.} \item{coef}{\code{signature(object
#' = "penfa", type = "free", labels = TRUE)}: Returns the estimates of the
#' parameters in the model as a named numeric vector. See the manual page of
#' \code{coef,penfa-method} for details.} \item{fitted}{\code{signature(object =
#' "penfa", labels = TRUE)}: Returns a list of the model-implied moments (per
#' group). See the manual page of \code{fitted,penfa-method} for details.}
#' }
#'
#' @seealso \code{\link{penfa}}, \code{\link{penfaParEstim}}
#'
#'
#' @export
#'
#'
#' @references Geminiani, E., Marra, G., & Moustaki, I. (2021). "Single- and
#'   Multiple-Group Penalized Factor Analysis: A Trust-Region Algorithm Approach
#'   with Integrated Automatic Multiple Tuning Parameter Selection."
#'   Psychometrika, 86(1), 65-95. \doi{10.1007/s11336-021-09751-8}
#'
#'
#'
#'
#'
setClass("penfa",
         slots = c(
           version     = "character",        # package version
           call        = "call",             # matched call
           timing      = "list",             # timing information
           Options     = "list",             # options
           ParTable    = "list",             # parameter table user-specified model
           pta         = "list",             # parameter table attributes
           Data        = "penfaData",        # full data
           SampleStats = "penfaSampleStats", # sample statistics
           Model       = "penfaModel",       # internal model representation
           Optim       = "list",             # optimizer results
           Penalize    = "penfaPenalty",     # information on penalization
           Implied     = "list",             # model-implied moments
           Vcov        = "list",             # vcov, std.err, conf. intervals
           Inference   = "list",             # effective dof, IC and more
           external    = "list"))            # optional empty slot
