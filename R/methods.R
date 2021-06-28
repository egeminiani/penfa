############################################################
# ---- Methods for printing and summarizing the output ----
############################################################
# Content:
#
# - set methods coef(), show(), summary(), fitted() for a 'penfa' class object
# - set method show() for a 'penfaData' class object
#
# - object_print_short_summary : controls how a 'penfa' object is printed/summarized
# - object_print_header : displays information on convergence status
# - penfadata_print_short : displays information on the data set structure
# - object_print_optim : displays information on the optimization process and edf (if model converged)
# - object_print_penalty : displays information on the penalization
# - penfaParEstim : displays information on the parameter estimates
# - print.penfaParEstim : controls how an object of class 'penfaParEstim" is printed
# - .makeNames : creates the names to be used when printing the parameter estimates
# - print.penfa.vector : print method for object of class "penfa.vector"
# - print.penfa.matrix.symmetric : controls how an object of class 'penfa.matrix.symmetric" is printed
#                                 (i.e., only lower triangle of a symmetric matrix)
# - plot.penfaMat : controls the plot of a penalty matrix of class 'penfaPenMat'
# - penfaOut : prints the output in the form of the estimated parameter matrices
# - penmat : extract the penalty matrix or matrices (shrinkage, differences, and full)
# - logLik.penfa : prints the loglikelihood (no penalty term) from a fitted penfa model
#
# Some of these functions are adaptations from the routines present in the lavaan package (https://CRAN.R-project.org/package=lavaan)
# -----------------------------------------------------------------------------------------------------------------------

#' Coefficients from a \code{penfa} object
#'
#' An S4 method returning the estimates of the model parameters
#'
#' @param object An object of class \code{penfa}, found as a result of a call
#'   to \code{penfa}.
#' @param type Character. If \code{type="free"}, only the estimated parameters
#'   (both penalized and unpenalized) are returned. If \code{type="user"}, all
#'   parameters listed in the parameter table are returned, including fixed
#'   parameters.
#' @param labels Logical. If \code{TRUE}, parameters are returned with their
#'   names.
#'
#' @seealso \code{\link{penfa}}, \code{\link{penfa-class}}
#'
#'
#' @examples
#'
#' ## --- Continuing the Examples from penfa manual page:
#' \dontshow{require(utils)
#' example("penfa", echo = FALSE)}
#' coef(alasso_fit)
#'
#'
#'
#' @export
#'
setMethod("coef", "penfa",
          function(object, type = "free", labels = TRUE) {
            inspect_coef(object     = object,
                         type       = type,
                         add.labels = labels,
                         add.class  = TRUE)
            })


#' Display a \code{penfa} object
#'
#' An S4 method printing a short summary of the estimation process, including
#' the optimization method, the specified penalty functions, the convergence
#' status, the number of iterations, the tuning selection strategy, and the
#' effective degrees of freedom.
#'
#' @param object An object of class \code{penfa}, found as a result of a call
#'   to \code{penfa}.
#'
#' @seealso \code{\link{penfa}}, \code{\link{penfa-class}}
#'
#'
#' @examples
#'
#' ## --- Continuing the Examples from penfa manual page:
#' \dontshow{require(utils)
#' example("penfa", echo = FALSE)}
#' alasso_fit
#' alasso_fitMG
#'
#'
#' @export
setMethod("show", "penfa",
          function(object) {
            # show only basic information
            object_print_short_summary(object)
          })


#' Summary constructor for a \code{penfa} object
#'
#' An S4 method printing a summary of the model parameter estimates for an
#' object of class \code{penfa}.
#'
#' @param object An object of class \code{penfa}, found as a result of a call
#'   to \code{penfa}.
#' @param header Logical. If \code{TRUE}, the header section is printed. The
#'   header contains relevant information about the data, the fitted model, the
#'   optimization process, and the penalization strategy, including, for
#'   instance, the employed penalties, the estimated effective degrees of
#'   freedom (\emph{edf}), the optimal values of the tuning parameter(s), the
#'   GBIC and many others.
#' @param estimates Logical. If \code{TRUE}, a section with the parameter
#'   estimates is printed out.
#' @param ci Logical. If \code{TRUE}, confidence intervals are added to the
#'   parameter estimates section.
#' @param level Logical. It denotes the significance level used for the
#'   statistical tests.
#' @param nd Integer. It determines the number of digits after the decimal point
#'   to be printed in the parameter estimates section.
#' @param cutoff Numeric. Standard errors and confidence intervals for the
#'   penalized parameter estimates falling below the \code{cutoff} value are not
#'   displayed. Confidence intervals for the parameters that have been penalized
#'   and shrunken to zero must be treated with caution.
#' @param extra Logical. If \code{TRUE}, additional information on the model
#'   are displayed.
#'
#' @seealso \code{\link{penfa}}, \code{\link{penfa-class}}
#'
#' @examples
#'
#' ## --- Continuing the Examples from penfa manual page:
#' \dontshow{require(utils)
#' example("penfa", echo = FALSE)}
#' summary(alasso_fit)
#' summary(alasso_fitMG)
#'
#'
#' @export
setMethod("summary", "penfa",
          function(object,
                   header    = TRUE,
                   estimates = TRUE,
                   ci        = TRUE,
                   level     = 0.95,
                   nd        = 3L,
                   cutoff    = 0.05,
                   extra     = TRUE){

            # return object
            res <- list()
            # print the 'long' summary (i.e., short + extra = TRUE)
            if(header) {
              object_print_short_summary(object, nd = nd, extra = extra)
            }

            if(estimates) {
              PE <- penfaParEstim(object, ci = ci, level = level, remove.nonfree = FALSE,
                                  output = "text", header = TRUE)
              print(PE, cutoff = cutoff, nd = nd)
              res$PE <- as.data.frame(PE)
            }

            invisible(res)
          })


#' Model-implied moments for a \code{penfa} object
#'
#' An S4 method returning the model-implied moments for an object of class
#' \code{penfa}. For every group, a list with the model-implied  moments is returned:
#' \code{cov} contains the implied covariance matrix, and \code{mean} the implied
#' mean vector. If just the covariance matrix is analyzed, only the \code{cov}
#' argument is returned.
#'
#' @param object An object of class \code{penfa}, found as a result of a call
#'   to \code{penfa}.
#' @param labels Logical. If \code{TRUE}, the model-implied moments are named
#'   according to the item names used in the model syntax.
#'
#' @seealso \code{\link{penfa}}, \code{\link{penfa-class}}
#'
#'
#' @examples
#'
#' ## --- Continuing the Examples from penfa manual page:
#' \dontshow{require(utils)
#' example("penfa", echo = FALSE)}
#' fitted(alasso_fit)
#' fitted(alasso_fitMG)
#'
#'
#' @export
setMethod("fitted", "penfa",
          function(object, labels = TRUE) {
            object_inspect_implied(object,
                                   add.labels = labels,
                                   add.class  = TRUE,
                                   drop.list.single.group = TRUE)
          })


#' Display details on the input data
#'
#' An S4 method showing information on the input data, including the number of
#' observations. In case of a multiple-group analysis, the sample sizes for each
#' group are displayed.
#'
#' @param object An object of class \code{penfaData}, found in the \code{Data}
#'   slot from a \code{penfa} class object.
#'
#' @seealso \code{\link{penfaData-class}}
#'
#'
#' @examples
#'
#' ## --- Continuing the Examples from penfa manual page:
#' \dontshow{require(utils)
#' example("penfa", echo = FALSE)}
#' alasso_fit@Data
#' alasso_fitMG@Data
#'
#' @export
setMethod("show", "penfaData",
          function(object) {
            # print 'penfaData' object
            penfadata_print_short(object) })


#' Display details on the penalization
#'
#' An S4 method showing information on the penalization process, including the
#' employed penalty functions and the model matrices they affect. Additionally,
#' it reports the optimal values of the tuning parameters and the tuning
#' parameter selection strategy. If the automatic procedure was used, the output
#' would also show the value of the influence factor, and the number of two-steps
#' iterations.
#'
#' @param object An object of class \code{penfaPenalty}, found in the
#'   \code{Penalize} slot from an object of \code{penfa} class.
#'
#' @seealso \code{\link{penfaPenalty-class}}
#'
#'
#' @examples
#'
#' ## --- Continuing the Examples from penfa manual page:
#' \dontshow{require(utils)
#' example("penfa", echo = FALSE)}
#' alasso_fit@Penalize
#' alasso_fitMG@Penalize
#'
#'
#' @export
setMethod("show", "penfaPenalty",
          function(object){
            # print penfaPenalty object
            object_print_penalty(object,
                                 extra = TRUE,
                                 header = TRUE)
          })



# object_print_short_summary
object_print_short_summary <- function(object, nd = 3L, extra = FALSE){

  # 1. print header
  object_print_header(object)

  # 2. print penfadata
  penfadata_print_short(object, nd = nd, extra = extra)

  # 3. print optim info
  object_print_optim(object, nd = nd, extra = extra)

  # 4. print penalization info
  object_print_penalty(object, nd = nd, extra = extra)
}

# object_print_header
object_print_header <- function(object) {

  cat(sprintf("penfa %s ", utils::packageDescription("penfa", fields="Version")))
  # Convergence or not?
  if(object@Optim$converged){
    cat(sprintf("reached convergence\n"))
    }else{
      cat(sprintf("did NOT reach convergence\n"))
    }
  cat("\n")
}

# penfadata_print_short
penfadata_print_short <- function(object, nd = 3L, extra = FALSE){

  if(inherits(object, "penfa")){
    penfadata <- object@Data
  }else if(inherits(object, "penfaData")){
    penfadata <- object
  }

  num.format  <- paste("%", max(8L, nd + 5L), ".", nd, "f", sep = "")

  # cat("Data information:\n\n")
  c1 <- c2 <- character(0L)

  # number of observations
  if(penfadata@ngroups == 1L) {
    c1 <- c(c1, "Number of observations")
    c2 <- c(c2, penfadata@nobs[[1L]])
  } else {
    c1 <- c(c1, "Number of observations per group:");
    c2 <- c(c2, "");
    for(g in 1:penfadata@ngroups) {
      c1 <- c(c1, sprintf("  %-40s", penfadata@group.label[[g]]))
      c2 <- c(c2, penfadata@nobs[[g]])
    } # g
  }

  if(extra){
    # holds for one group only;
    c1 <- c(c1, "Number of groups", "Number of observed variables", "Number of latent factors")
    c2 <- c(c2, penfadata@ngroups, object@pta$nvar[[1]], object@pta$nfac[[1]])
  }

  # empty row
  c1 <- c(c1, ""); c2 <- c(c2, "")

  # format c1/c2
  c1 <- format(c1, width = 42L)
  c2 <- format(c2, width = 18L + max(0, (nd - 3L)) * 4L, justify = "right")

  # create character matrix
  M <- cbind(c1, c2, deparse.level = 0)
  colnames(M) <- rep("",  ncol(M))
  rownames(M) <- rep(" ", nrow(M))

  # print
  utils::write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)
  invisible(M)
}

# object_print_optim
object_print_optim <- function(object, nd = 3L, extra = FALSE){

  #cat("Optimization information:\n\n")

  num.format  <- paste("%", max(8L, nd + 5L), ".", nd, "f", sep = "")

  c1 <- c("Estimator", "Optimization method", "Information", "Strategy")

  if(all(unlist(object@Penalize@tuning) == 0)){ # at least one non-zero tuning?
    estimator <- "MLE"
  }else{
    estimator <- "PMLE"
  }

  c2 <- c(estimator, "trust-region",
          object@Options$information,
          object@Penalize@strategy)

  if(object@Penalize@strategy == "fixed"){
    c1 <- c(c1, "Number of iterations")
    c2 <- c(c2, object@Optim$iterations)
  }else{ # auto
    c1 <- c(c1,  "Number of iterations (total)", "Number of two-steps (automatic)")
    c2 <- c(c2,  object@Optim$iterations  + object@Penalize@automatic$iter.inner,
            object@Penalize@automatic$iter)
  }

  if(extra){
    if(object@Penalize@strategy == "auto"){
      c1 <- c(c1, "Influence factor")
      c2 <- c(c2, object@Penalize@automatic$gamma)
    }

    c1 <- c(c1, "Number of parameters:", "  Free", "  Penalized")
    c2 <- c(c2, "", sum(object@ParTable$type=="free"), sum(object@ParTable$type=="pen") )
    # qstar <- length(Reduce(union, unlist(object@Penalize@pen.idx, recursive = FALSE))); object@Optim$npar-qstar
  }

  if(!is.null(object@Inference$edf)){
    c1 <- c(c1, "Effective degrees of freedom")
    c2 <- c(c2, sprintf(num.format, object@Inference$edf))
  }

  if(extra){
    c1 <- c(c1, "GIC", "GBIC")
    c2 <- c(c2, sprintf(num.format, object@Inference$IC$AIC),
            sprintf(num.format, object@Inference$IC$BIC))
  }

  #empty last row
  c1 <- c(c1, ""); c2 <- c(c2, "")

  # format
  c1 <- format(c1, width = 42L)
  c2 <- format(c2, width = 18L + max(0, (nd - 3L)) * 4L, justify = "right")

  # character matrix
  M <- cbind(c1, c2, deparse.level = 0)
  colnames(M) <- rep("",  ncol(M))
  rownames(M) <- rep(" ", nrow(M))

  # print
  utils::write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)
  invisible(M)
}

# object_print_penalty
object_print_penalty <- function(object, nd = 3L, extra = FALSE, header = FALSE){ # header for print penfaPenalty object

  num.format  <- paste("%", max(8L, nd + 5L), ".", nd, "f", sep = "")
  c1 <- c2 <- character(0L)

  if(inherits(object, "penfa")){
    object <- object@Penalize
  }

  stopifnot(inherits(object, "penfaPenalty"))

  if(all(unlist(object@tuning) == 0) | all(unlist(object@penalty == "none"))){
    # full MLE, display no info
    # c1 <- c(c1, "Penalty", "Optimal tuning parameter")
    # c2 <- c(c2, "none",    "-" )

    # display info only for penfaPenalty class object
    if(extra && header){
      c1 <- c(c1, "Penalty function")
      c2 <- c(c2, "none")
    }

  }else{

    # PMLE
    real.tunings <- object@tuning
    # select nonzero tunings
    which_tun_to_select <- lapply(real.tunings, function(x){ x != 0})
    real.tunings <- mapply(function(x, i) x[i],
                           real.tunings, which_tun_to_select)
    real.tunings <- real.tunings[lapply(real.tunings,length)>0]


    # Is there sparsity/invariance?
    tun.shrink <- real.tunings$shrink; lx.shrink  <- length(tun.shrink)
    tun.diff <- real.tunings$diff; lx.diff  <- length(tun.diff)
    penalties.add.tun <- c("alasso", "scad", "mcp")

    if(header){
      c1 <- c(c1, "Strategy")
      c2 <- c(c2, object@strategy)
      if(object@strategy == "auto"){
        c1 <- c(c1, "Influence factor", "Number of two-steps (automatic)", "")
        c2 <- c(c2, object@automatic$gamma, object@automatic$iter, "")
      }
    }

    # Info on penalty function
    if(sum(unlist(object@penalty) == "none") == 1){ # only 1 penalty
      c1 <- c(c1, "Penalty function:")
      c2 <- c(c2, "")

      # which one?
      purpose   <- ifelse(lx.shrink > 0, "  Sparsity", "  Invariance")
      right.pen <- ifelse(purpose == "  Sparsity", object@penalty$shrink[1], object@penalty$diff[1])
      c1        <- c(c1, purpose, "")
      c2        <- c(c2, right.pen, "")

      if(extra){
        # If alasso, scad, mcp --> additional tuning;
        if(right.pen %in% penalties.add.tun){
          c1 <- c(c1, "Additional tuning parameter", sprintf("  %s", right.pen), "")
          # in case it's alasso, discard the "weights" argument
          list.add.tun <- object@extra[names(object@extra) != "weights"]
          c2 <- c(c2,  "", unlist(list.add.tun), "")
        } #else do nothing
      }

    }else{ # 2 penalties
      c1 <- c(c1, "Penalty functions:", "  Sparsity", "  Invariance", "")
      pen4shrinkage <- object@penalty$shrink[1]
      pen4diff      <- object@penalty$diff[1]
      c2 <- c(c2, "", pen4shrinkage, pen4diff, "")
      if(extra){

        # two *different* penalties having *both* an additional tuning
        if(all(c(pen4shrinkage, pen4diff) %in% penalties.add.tun) & pen4shrinkage != pen4diff){
          txt <- "Additional tuning parameters"
          c1 <- c(c1, txt, sprintf("  %s", pen4shrinkage), sprintf("  %s", pen4diff), "")
          list.add.tun <- object@extra[names(object@extra) != "weights"]
          # remove "a." to get the penalty name
          str.pen <- unlist(lapply(names(object@extra), function(x){ substring(x, 3) }))
          c2 <- c(c2,  "", unlist(list.add.tun)[which(str.pen == pen4shrinkage)],
                  unlist(list.add.tun)[which(str.pen == pen4diff)] , "")
          }else if(all(c(pen4shrinkage, pen4diff) %in% penalties.add.tun) & pen4shrinkage == pen4diff){
            # same penalty *with additional tuning* twice
            txt  <- "Additional tuning parameter"
            c1 <- c(c1, txt, sprintf("  %s", pen4shrinkage), "")
            list.add.tun <- object@extra[names(object@extra) != "weights"]
            str.pen <- unlist(lapply(names(object@extra), function(x){ substring(x, 3) }))
            c2 <- c(c2,  "", unlist(list.add.tun)[which(str.pen == pen4shrinkage)], "")
          }else if(any(c(pen4shrinkage, pen4diff) %in% penalties.add.tun)){
            # different penalties with one not having additional tuning
            txt  <- "Additional tuning parameter"
            pen.vector <- c(pen4shrinkage, pen4diff)
            two.tuning.pen <- pen.vector[which(pen.vector %in% penalties.add.tun)]
            c1 <- c(c1, txt, sprintf("  %s", two.tuning.pen), "")
            list.add.tun <- object@extra[names(object@extra) != "weights"]
            c2 <- c(c2,  "", unlist(list.add.tun), "")
          }else{
            # none of the two penalties has an additional tuning (e.g. lasso and ridge)
            # do nothing
        }
      } # extra
    }

    if(extra){
      nm.shrink  <- names(tun.shrink)
      nm.diff    <- names(tun.diff)
      # Single or multiple tunings?
      c1 <- c(c1, c(ifelse(length(unlist(real.tunings)) == 1, "Optimal tuning parameter:", "Optimal tuning parameters:")))
      c2 <- c(c2, "")

      dt <- data.frame(ch.vec = c("Factor loadings", "Intercepts", "Factor (co)variances",
                                  "Factor means", "Unique variances"),
                       pmat.nm =  c("lambda", "tau", "phi", "kappa", "psi"))

      if(lx.shrink > 0){
        c1 <- c(c1, "  Sparsity")
        c2 <- c(c2, "")
        common <- merge(as.data.frame(nm.shrink), dt, by.x = "nm.shrink", by.y = "pmat.nm",
                        sort = FALSE) # no sort, so maintain the right order
        c1 <- c(c1, sprintf("   - %-37s", common$ch.vec))
        c2 <- c(c2, sprintf(num.format, tun.shrink))
      }

      if(lx.diff > 0){
        c1 <- c(c1, "  Invariance")
        c2 <- c(c2, "")
        common <- merge(as.data.frame(nm.diff), dt, by.x = "nm.diff", by.y = "pmat.nm", sort = FALSE)
        c1 <- c(c1, sprintf("   - %-37s", common$ch.vec))
        c2 <- c(c2, sprintf(num.format, tun.diff))
      }

    }
  } # PMLE

  # empty last row
  c1 <- c(c1, ""); c2 <- c(c2, "")

  # format
  c1 <- format(c1, width = 42L)
  c2 <- format(c2, width = 18L + max(0, (nd - 3L)) * 4L, justify = "right")

  # character matrix
  M <- cbind(c1, c2, deparse.level = 0)
  colnames(M) <- rep("",  ncol(M))
  rownames(M) <- rep(" ", nrow(M))

  # print
  utils::write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)

  invisible(M)
}

#' Print parameter estimates in table format
#'
#' @description The parameter estimates of the penalized factor analysis model
#'   in each group.
#'
#' @param object An object of class \code{\linkS4class{penfa}}.
#' @param se Logical. If \code{TRUE}, it includes a column with the standard
#'   errors.
#' @param ci Logical. If \code{TRUE}, the confidence intervals are added to the
#'   output.
#' @param level The confidence level, default is 0.95.
#' @param remove.nonfree 	Logical. If \code{TRUE}, it filters the output and
#'   removes all rows with fixed (that is, neither free, nor penalized)
#'   parameters.
#' @param output Character. If "data.frame", the parameter table is displayed as
#'   a standard formatted data.frame. If "text", the parameter table is
#'   displayed with subsections (as used by the \code{summary} function).
#' @param header 	Logical, only used if \code{output = "text"}. If \code{TRUE},
#'   it prints a header on top of the parameter list with details on the group
#'   levels and the information matrix used during optimization by the
#'   trust-region algorithm.
#'
#' @seealso \code{\link{penfa}}
#'
#'
#' @export
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
#' penfaParEstim(alasso_fit)
#'
penfaParEstim <- function(object,
                     se     = TRUE,
                     ci     = TRUE,
                     level  = 0.95,
                     remove.nonfree = FALSE,
                     output = "data.frame",
                     header = FALSE) {

  PARTABLE <- as.data.frame(object@ParTable, stringsAsFactors = FALSE)
  LIST <- PARTABLE[,c("lhs", "op", "rhs", "type", "penalty", "free")]

  if(!is.null(PARTABLE$user)) {
    LIST$user <- PARTABLE$user
  }
  if(!is.null(PARTABLE$group)) {
    LIST$group <- PARTABLE$group
  } else {
    LIST$group <- rep(1L, length(LIST$lhs))
  }
  if(!is.null(PARTABLE$est)){
    LIST$est <- PARTABLE$est
  }

  # add se
  if(se) {
    LIST$se <- PARTABLE$se
    # handle tiny SEs
    LIST$se <- ifelse(LIST$se < sqrt(.Machine$double.eps), 0,  LIST$se)
    tmp.se  <- ifelse(LIST$se < sqrt(.Machine$double.eps), NA, LIST$se)
  }

  ind.pen.p   <- which(LIST$type == "pen")
  ind.free.p  <- which(LIST$type == "free")
  ind.fixed.p <- which(LIST$type == "fixed")

  # confidence interval
  if(se & ci) {
    LIST$ci.upper <- LIST$ci.lower <- NA
    ci  <- compute_CI(param = object@Optim$x, VCOV = object@Vcov$vcov,
                      std.err = object@Vcov$se, modpenalty = object@Penalize,
                      options = object@Options, level = level)
    if(length(ind.fixed.p) > 0){
      LIST$ci.lower[-ind.fixed.p] <- ci[,1]
      LIST$ci.upper[-ind.fixed.p] <- ci[,2]
      LIST$ci.lower[ind.fixed.p] <- LIST$ci.upper[ind.fixed.p] <- LIST$est[ind.fixed.p]
    }
  }

  # if single group, remove group column
  if(object@Data@ngroups == 1L) LIST$group <- NULL

  # remove non-free paramters
  if(remove.nonfree) {
    nonfree.idx <- which( LIST$free == 0L )
    if(length(nonfree.idx) > 0L) {
      LIST <- LIST[-nonfree.idx,]
    }
  }
  # remove 'free' column and LIST$user
  LIST$free <- LIST$user <- NULL

  # Collect info on significance level by renaming the columns ci.lower and ci.upper
  alpha   <- (1 - level)/2; alpha <- c(alpha, 1 - alpha)
  names(LIST)[names(LIST) %in% c("ci.lower", "ci.upper")] <- paste0("ci_",alpha)

  if(output == "text") {
    class(LIST) <- c("penfaParEstim", "penfa.data.frame", "data.frame")
    if(header) {
      attr(LIST, "information") <- object@Options$information
      attr(LIST, "group.label") <- object@Data@group.label
      attr(LIST, "header")      <- header
    }
  } else {
    class(LIST) <- c("penfa.data.frame", "data.frame")
  }
  LIST
}

#' @export
# print.penfaParEstim
print.penfaParEstim <- function(x, ..., cutoff = 0.05, nd = 3L) {

  # format for numeric values
  num.format  <- paste("%", max(8L, nd + 5L), ".", nd, "f", sep = "")
  char.format <- paste("%", max(8L, nd + 5L), "s", sep = "")
  # output sections
  GSECTIONS <- c("Latent Variables", "Covariances", "Intercepts", "Variances")

  # header?
  header <- attr(x, "header")
  if(is.null(header)) {
    header <- FALSE
  }
  if(header) {
    cat("\nParameter Estimates:\n")
    # info about standard errors (if we have x$se only) : information matrix
    if(!is.null(x$se)){
      # container
      c1 <- c2 <- character(0L)
      # c1      <- c(c1, "Information")
      # tmp.txt <- attr(x, "information")
      # c2      <- c(c2, paste(toupper(substring(tmp.txt, 1, 1)), substring(tmp.txt, 2), sep = ""))

      # format c1/c2
      c1 <- format(c1, width = 38L)
      c2 <- format(c2, width = 13L + max(0, (nd - 3L)) * 4L, justify = "right")

      # create character matrix
      M <- cbind(c1, c2, deparse.level = 0)
      colnames(M) <- rep("",  ncol(M))
      rownames(M) <- rep(" ", nrow(M))

      # print
      utils::write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)
    }
  }

  # number of groups
  if(is.null(x$group)) {
    ngroups <- 1L
    x$group <- rep(1L, length(x$lhs))
  } else {
    ngroups <- partable_ngroups(x)
  }

  # round to 3 digits after the decimal point
  y <- as.data.frame( lapply(x, function(x) { if(is.numeric(x)) { sprintf(num.format, x) }else{ x } }),
                      stringsAsFactors = FALSE)

  # always remove op/group/rhs/penalty columns
  y$op <- y$group <- y$rhs <- y$penalty <- NULL
  alpha <- sapply(strsplit(names(x)[which(grepl("ci_", names(x)))], "ci_"), "[", 2)
  alpha <- as.numeric(alpha)

  # convert to character matrix
  m <- as.matrix(format.data.frame(y, na.encode = FALSE, justify = "right"))

  # use empty row names
  rownames(m) <- rep("", nrow(m))

  # handle se == 0.0
  if(!is.null(x$se)){
    se.idx <- which(x$se == 0)
    if(length(se.idx) > 0L) {
      m[se.idx, "se"] <- ""
      if(!is.null(x$z)) {
        m[se.idx, "z"] <- ""
      }
      if(!is.null(x$pvalue)) {
        m[se.idx, "pvalue"] <- ""
      }
    }

    # handle se == NA
    se.idx <- which(is.na(x$se))
    if(length(se.idx) > 0L) {
      if(!is.null(x$z)) {
        m[se.idx, "z"] <- ""
      }
      if(!is.null(x$pvalue)) {
        m[se.idx, "pvalue"] <- ""
      }
    }

    # hide CI for penalized parameters shrunken to zero
    # hide.idx <- which(abs(x$est) < cutoff & x$type == "pen")
    # if we want only for penalized factor loadings
    hide.idx <- which(abs(x$est) < cutoff & (x$penalty == "shrink" | x$penalty == "shrink + diff"))
    if(length(hide.idx) > 0){

      if(!is.null(x$se)){
        m[hide.idx, "se"] <- ""
        m[hide.idx, which(  grepl("ci_", attr(m, "dimnames")[2][[1]])) ] <- ""
      }
    }
    if(cutoff >= 0 & cutoff < 0.01){
      warning("penfa WARNING: Confidence intervals for penalized parameters shrunken to zero must be treated with caution.")
    }
  }

  # rename some column names
  colnames(m)[ colnames(m) ==    "lhs" ] <- ""
  colnames(m)[ colnames(m) ==     "op" ] <- ""
  colnames(m)[ colnames(m) ==    "rhs" ] <- ""
  colnames(m)[ colnames(m) ==    "type" ] <- "Type   "
  colnames(m)[ colnames(m) ==    "est" ] <- "Estimate"
  colnames(m)[ colnames(m) ==     "se" ] <- "Std.Err"
  colnames(m)[ which(grepl("ci_", colnames(m))) ] <- sprintf("%.*g%%", 3, 100 * alpha)

  # format column names
  colnames(m) <- sprintf(char.format, colnames(m))

  # group-specific sections
  for(g in 1:ngroups) {

    # group header
    if(ngroups > 1L) {
      group.label <- attr(x, "group.label")
      cat("\n\n")
      cat("Group ", g, " [", group.label[g], "]:\n", sep="")
    }

    # ov/lv names
    # use group instead of block
    ov.names <- partable_vnames(x, "ov", group = g)
    lv.names <- partable_vnames(x, "lv", group = g)

    # group-specific sections
    for(s in GSECTIONS) {
      if(s == "Latent Variables") {
        row.idx <- which( x$op == "=~" & !x$lhs %in% ov.names & x$group == g)
        if(length(row.idx) == 0L) next
        m[row.idx,1] <- .makeNames(x$rhs[row.idx], x$label[row.idx])
      } else if(s == "Covariances") {
        row.idx <- which(x$op == "~~" & x$lhs != x$rhs & x$group == g)
        if(length(row.idx) == 0L) next
        # make distinction between residual and plain
        y.names <- unique(ov.names)
        PREFIX <- rep("", length(row.idx))
        PREFIX[ x$rhs[row.idx] %in% y.names ] <- "  ."
        m[row.idx,1] <- .makeNames(x$rhs[row.idx], x$label[row.idx], PREFIX = PREFIX)
      } else if(s == "Intercepts") {
        row.idx <- which(x$op == "~1" & x$group == g)
        if(length(row.idx) == 0L) next
        # make distinction between intercepts and means
        y.names <- unique(ov.names)
        PREFIX <- rep("", length(row.idx))
        PREFIX[ x$lhs[row.idx] %in% y.names ] <- "  ."
        m[row.idx,1] <- .makeNames(x$lhs[row.idx], x$label[row.idx], PREFIX = PREFIX)
      } else if(s == "Variances") {
        row.idx <- which(x$op == "~~" & x$lhs == x$rhs & x$group == g)
        if(length(row.idx) == 0L) next
        # make distinction between residual and plain
        y.names <- unique(ov.names)
        PREFIX <- rep("", length(row.idx))
        PREFIX[ x$rhs[row.idx] %in% y.names ] <- "  ."
        m[row.idx,1] <- .makeNames(x$rhs[row.idx], x$label[row.idx], PREFIX = PREFIX)
      } else {
        row.idx <- integer(0L)
      }

      # two formatting types:
      #  - for regular (nothing to do, except row/colnames)
      #  - for Latent Variables and Covariances 'bundle' the output per lhs element

      # bundling
      if(s %in% c("Latent Variables", "Covariances")) {
        nel <- length(row.idx)
        M <- matrix("", nrow = nel*2, ncol = ncol(m))
        colnames(M) <- colnames(m)
        rownames(M) <- rep("", NROW(M))
        LHS <- paste(x$lhs[row.idx], x$op[row.idx])
        lhs.idx <- seq(1, nel*2L, 2L)
        rhs.idx <- seq(1, nel*2L, 2L) + 1L
        if(s == "Covariances") {
          # make distinction between residual and plain
          y.names <- unique( ov.names )
          PREFIX <- rep("", length(row.idx))
          PREFIX[ x$lhs[row.idx] %in% y.names ] <- "."
        } else {
          PREFIX <- rep("", length(LHS))
        }
        M[lhs.idx, 1] <- sprintf("%1s%-15s", PREFIX, LHS)
        M[rhs.idx,  ] <- m[row.idx,]
        # avoid duplicated LHS labels
        if(nel > 1L) {
          del.idx <- integer(0)
          old.lhs <- ""
          for(i in 1:nel) {
            if(LHS[i] == old.lhs) {
              del.idx <- c(del.idx, lhs.idx[i])
            }
            old.lhs <- LHS[i]
          }
          if(length(del.idx) > 0L) {
            M <- M[-del.idx,,drop=FALSE]
          }
        }
        cat("\n", s, ":\n", sep = "")
        #cat("\n")
        print(M, quote = FALSE)
      }else { # Regular
        M <- m[row.idx,,drop=FALSE]
        colnames(M) <- colnames(m)
        rownames(M) <- rep("", NROW(M))
        cat("\n", s, ":\n", sep = "")
        print(M, quote = FALSE)
      }
    } # GSECTIONS

  } # groups

  cat("\n")
  invisible(m)
}

# .makeNames
.makeNames <- function(NAMES, LABELS, PREFIX = NULL) {

  W <- 14
  if(is.null(PREFIX)) {
    PREFIX <- rep("", length(NAMES))
  }

  multiB <- FALSE
  if(any(nchar(NAMES) != nchar(NAMES, "bytes")))
    multiB <- TRUE
  if(any(nchar(LABELS) != nchar(LABELS, "bytes")))
    multiB <- TRUE
  # labels?
  l.idx <- which(nchar(LABELS) > 0L)
  if(length(l.idx) > 0L) {
    if(!multiB) {
      LABELS <- abbreviate(LABELS, 4)
      LABELS[l.idx] <- paste(" (", LABELS[l.idx], ")", sep="")
      MAX.L <- max(nchar(LABELS))
      NAMES <- abbreviate(NAMES, minlength = (W - MAX.L),
                          strict = TRUE)
    } else {
      MAX.L <- 4L
    }
    NAMES <- sprintf(paste("%-", (W - MAX.L), "s%", MAX.L, "s",
                           sep=""), NAMES, LABELS)
  } else {
    if(!multiB) {
      NAMES <- abbreviate(NAMES, minlength = W, strict = TRUE)
    } else {
      NAMES <- sprintf(paste("%-", W, "s", sep=""), NAMES)
    }
  }
  char.format <- paste("%3s%-", W, "s", sep = "")
  sprintf(char.format, PREFIX, NAMES)
}

#' @export
# print.penfa.vector
print.penfa.vector <- function(x, ..., nd = 3L) {
  y <- unclass(x)
  print( round(y, nd), ... )
  invisible(x)
}

#' @export
# print.penfa.matrix.symmetric
print.penfa.matrix.symmetric <- function(x, ..., nd = 3L) {
  y <- x; y <- unclass(y)
  ll <- lower.tri(x, diag=TRUE)
  y[ll] <- format(round(x[ll], digits=nd)); y[!ll] <- ""
  if (!is.null(colnames(x))) {
    colnames(y) <- abbreviate(colnames(x), minlength = nd + 3L)
  }
  print(y, ..., quote = FALSE)
  invisible(x)
}

#' @export
# plot.penfaMat
plot.penfaPenMat <- function(x, ...,
                             title = "Penalty matrix (absolute value and log scale)",
                             col = NULL){

  stopifnot(inherits(x, "penfaPenMat"))

  # abs to handle negative values when using difference penalties
  # log1p to handle zeroes
  logx <- log1p(abs(x))
  if(is.null(col)){
    col <- cartography::carto.pal(pal1 = "blue.pal", n1 = 20)
  }

  p <- plotly::plot_ly(x = rownames(logx),
                       y = colnames(logx),
                       z = logx,
                       type = "heatmap",
                       colors = grDevices::colorRampPalette(col)(50))
  plotly::layout(p, yaxis = list(autorange = "reversed"),
                 title = list(text = title))
}


#' Print estimated parameter matrices
#'
#' @description A utility that extracts the estimated parameter matrices
#'   and vectors of the penalized factor model for each group and rounds them to
#'   the specified number of decimal digits.
#'
#' @param object An object of class \code{\linkS4class{penfa}}, that is, a
#'   fitted penalized factor model.
#' @param which Character denoting the name of the estimated matrix or vector to
#'   display. Possible values are "lambda", "psi", "phi", "tau", and "kappa".
#'   Multiple elements can be specified. By default, all estimated matrices are
#'   shown.
#' @param ... Additional options.
#' @param nd The number of decimal digits to be used.
#'
#' @return List of the estimated parameter matrices and vectors for each group.
#'
#' @export
#'
#' @seealso \code{\link{penfa}}
#'
#'
#'
penfaOut <- function(object, which = c("lambda", "psi", "phi", "tau", "kappa"),
                     ..., nd = 3L){

  y     <- set_names2penfaOut(object)
  GLIST <- y@Model@GLIST
  mat   <- lapply(GLIST, function(x) round(x, nd, ...) )
  mat
  idx <- c()

  for(i in 1:length(which)){
    # admissible name?
    name.check <- which[i] %in% c("lambda", "psi", "phi", "tau", "kappa")
    if(name.check==FALSE){
      stop("penfa ERROR: \"which\" must be any of \"lambda\", \"psi\", \"phi\", \"tau\", or \"kappa\"")
    }
    idx <- c(idx, which(names(mat) %in% which[i]))
  }
  mat <- mat[sort(idx)]
  return(mat)
}



set_names2penfaOut <- function(object){
  for(i in 1:length(object@Model@GLIST)){
    rownames(object@Model@GLIST[[i]]) <- object@Model@dimNames[[i]][[1]]
    colnames(object@Model@GLIST[[i]]) <- object@Model@dimNames[[i]][[2]]
  }
  # # If multiple groups, add group variable name and pick the right group matrices
  # ngroups <- object@Model@ngroups
  # if(ngroups > 1){
  #   GLIST.group <- vector("list", length = ngroups)
  #   names(GLIST.group) <- object@Data@group.label
  #   for(g in 1:ngroups){
  #     nmat  <- object@Model@nmat
  #     mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
  #     GLIST.group[[g]] <- object@Model@GLIST[mm.in.group]
  #   } # g
  #   object@Model@GLIST <- GLIST.group
  # }
  return(object)
}


#' Extract estimated penalty matrix
#'
#' @description A utility that extracts the estimated penalty matrix from a
#'   fitted object of class \code{penfa}.
#'
#' @param x An object of class \code{\linkS4class{penfa}}, that is, a
#'   fitted penalized factor model.
#' @param type Character denoting the type of penalization. Type equal to
#'   \code{"full"} returns the complete penalty matrix; \code{type="shrink"}
#'   returns the penalty matrix for shrinkage; \code{type="diff"}
#'   the penalty matrix for parameter equivalence. The matrix
#'   returned by \code{type="full"} is the sum of all the \code{shrink} and
#'   \code{diff} penalty sub-matrices.
#' @param which Character prompting the extraction of the penalty matrix
#'   component corresponding to the specified model matrix. It is only valid
#'   when \code{type="shrink"} or \code{type="diff"}. Possible values are
#'   \code{"lambda"}, \code{"psi"}, \code{"phi"}, \code{"tau"}, \code{"kappa"}
#'   and \code{"none"}.  Only the model matrices penalized during model fitting
#'   (i.e., in the \code{penfa} call) can appear in the \code{which} argument.
#'
#' @return A penalty matrix of class \code{penfaPenMat}. If multiple elements
#'   are specified in the \code{which} argument, a list of penalty matrices (one
#'   for each element, and each of class \code{penfaPenMat}) is returned.
#'
#' @export
#'
#' @seealso \code{\link{penfa}}
#'
#'
penmat <- function(x, type = "full", which = NULL){

  stopifnot(inherits(x, "penfa"))
  stopifnot(type %in% c("full", "shrink", "diff"))

  if(type == "full"){
    out <- x@Penalize@Sh.info$S.h
    if(!is.null(which)){
      warning("penfa WARNING: ignored \"which\" argument when type=\"full\".\n")
    }

  }else{

    if(type == "shrink"){
      out.tmp <- x@Penalize@Sh.info$SS.shrink
    }else if(type == "diff"){
      out.tmp <- x@Penalize@Sh.info$SS.diff
    }

    # Check if the specified matrix can be extracted
    if(!is.null(which)){
      idx <- c()
      # Extract the matrices of interest
      for(i in 1:length(which)){
        # admissible name?
        name.check <- which[i] %in% c("lambda", "psi", "phi", "tau", "kappa", "none")
        if(!name.check){
          stop("penfa ERROR: \"which\" must be any of \"lambda\", \"psi\", \"phi\", \"tau\", \"kappa\" or \"none\".")
        }

        # name in the penalty?
        existence.check <- which[i] %in% names(out.tmp)
        if(!existence.check){
          stop("penfa ERROR: the model matrix \"", which[i], "\" specified in the \"which\" argument was not found in the penalty matrix.\n")
        }
        # Good to go
        idx <- c(idx, which(names(out.tmp) %in% which[i]))
      }

      out <- out.tmp[idx]
      # If just one matrix, return matrix (not list)
      if(length(idx)==1){
        out <- out[[1]]
      }

    }else{
      # 'which' was not specified
      out <- out.tmp
    }
  }
  return(out)
}

#' @export
# logLik (to be used with AIC/BIC functions)
logLik.penfa <- function(object, ...){

  # Only compute AIC/BIC when the model converged
  if(object@Model@nx.free > 0L && object@Optim$converged &&
     all(object@Vcov$solution) && object@Vcov$admissibility){

    logl <- object@Optim$logl.unpen # (unpenalized) log-likelihood
    attr(logl, "df")   <- object@Inference$edf  # the effective dof
    attr(logl, "nobs") <- object@SampleStats@ntotal
    class(logl)        <- "logLik"

  }else{
    stop("penfa ERROR: The model did not converge.")
  }
  logl
}
