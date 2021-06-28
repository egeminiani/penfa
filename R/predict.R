##########################################
# ---- Functions for predict() output ----
##########################################
# Content:
#
# - penfaPredict : computes factor scores from a fitted penfa model object.
# - predict_eta_normal : factor scores for the normal case, computed according to
#                        the classic 'regression' method.
# - computeEETA, computeEY, computeEY.LISREL : compute necessary quantities.
# - predict_eta_bartlett : factor scores for the normal case, computed according
#                          to Bartlett's method.
#
# Some of these functions are adaptations from the routines present in the lavaan package (https://CRAN.R-project.org/package=lavaan)
# -----------------------------------------------------------------------------------------------------------------------

#' Compute the factor scores from a fitted \code{penfa} model
#'
#' @description The \code{penfaPredict} function estimates the factor scores
#'   from a fitted penalized factor model. The factor scores are the estimated
#'   values ("predictions") of the common factors.
#'
#' @param object An object of class \code{\linkS4class{penfa}}.
#' @param newdata An optional data frame containing the same variables as the
#'   ones appearing in the original data frame used for fitting the model in
#'   \code{object}.
#' @param method Character indicating the method for computing the factor
#'   scores. Possible options are \code{"regression"} and \code{"bartlett"}. For
#'   the normal linear continuous case, the regression method is equivalent to
#'   the Empirical Bayes Method (EBM), whereas Bartlett's strategy is equivalent
#'   to maximum likelihood's method.
#' @param label Logical. If \code{TRUE}, the columns are labeled.
#' @param append.data Logical. If \code{TRUE}, the original data set (or the
#'   data set provided in the \code{newdata} argument) is appended to the factor
#'   scores.
#' @param assemble Logical. If \code{TRUE}, the factor scores from each group
#'   are assembled in a single data frame of the same dimensions as the original
#'   data set and with a group column defining the groups.
#'
#' @export
#'
#' @seealso \code{\link{penfa}}
#'
#'
#'
#' @references Geminiani E. (2020), "A penalized likelihood-based framework for
#'   single and multiple-group factor analysis models" (Doctoral dissertation,
#'   University of Bologna). Available at
#'   \url{http://amsdottorato.unibo.it/9355/}.
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
#'                     strategy = "auto",
#'                     gamma = 4)
#'
#' fscores <- penfaPredict(alasso_fit)
#'
penfaPredict <- function(object,
                         newdata     = NULL,
                         method      = "regression",
                         label       = TRUE,
                         append.data = FALSE,
                         assemble    = FALSE){

  stopifnot(inherits(object, "penfa"))
  model       <- object@Model
  data        <- object@Data
  samplestats <- object@SampleStats
  implied     <- object@Implied
  pta         <- object@pta

  # need full data set supplied
  if(is.null(newdata)){
    # if(data@data.type != "full"){
    #   stop("penfa ERROR: sample statistics were used for fitting and newdata is empty")
    # }else
    if(is.null(data@X[[1]])){
      stop("penfa ERROR: no local copy of data")
    }else{
      data.obs <- data@X
    }
  }else{
    OV <- data@ov
    newData <- penfaData(data = newdata, group = data@group,
                         ov.names = data@ov.names,
                         options = list(std.ov = data@std.ov,
                                        group.label=data@group.label,
                                        warn = FALSE))
    data.obs <- newData@X
  }

  ok  <- is_Admissible(object@Model, verbose = FALSE)
  out <- predict_eta(object = NULL, model = model, data = data,
                     samplestats = samplestats, implied = implied, data.obs = data.obs,
                     method = method)

  # remove dummy lv? (removes attr!)
  out <- lapply(seq_len(data@ngroups), function(g){ out[[g]] })

  # append original/new data? (also remove attr)
  if(append.data){
    out <- lapply(seq_len(data@ngroups), function(g) {
      ret <- cbind(out[[g]], data.obs[[g]])
      ret
    })
  }

  # label
  if(label){
    for(g in seq_len(data@ngroups)){
      if(append.data) {
        colnames(out[[g]]) <- c(pta$vnames$lv[[1]], data@ov.names[[1]])
      }else{
        colnames(out[[g]]) <- pta$vnames$lv[[1]]
      }
    }

    # group.labels
    if(data@ngroups > 1L) {
      names(out) <- data@group.label
    }
  } # label

  # matrix class
  out <- lapply(out, "class<-", c("matrix"))

  if(data@ngroups == 1L) {
    res <- out[[1L]]
  } else {
    res <- out
  }

  # assemble multiple groups into a single data.frame?
  if(data@ngroups > 1L && assemble) {
    if(!is.null(newdata)) {
      data <- newData
    }
    DATA <- matrix(as.numeric(NA), nrow = sum(unlist(data@norig)),
                   ncol = ncol(out[[1L]])) # assume == per g
    colnames(DATA) <- colnames(out[[1L]])
    for(g in seq_len(data@ngroups)) {
      DATA[ data@case.idx[[g]], ] <- out[[g]]
    }
    DATA <- as.data.frame(DATA, stringsAsFactors = FALSE)
    if(!is.null(newdata)) {
      DATA[, data@group] <- newdata[, data@group ]
    } else {
      # add group
      DATA[, data@group ] <- rep(as.character(NA), nrow(DATA))
      # we will loose the group label of omitted variables
      DATA[unlist(data@case.idx), data@group] <- rep(data@group.label, unlist(data@nobs))
    }
    res <- DATA
  }

  res
}

# predict_eta
predict_eta <- function(object = NULL, model = NULL, data = NULL, samplestats = NULL,
                        implied = NULL, data.obs = NULL, method = "regression"){
    # full object?
    if(inherits(object, "penfa")) {
      data <- object@Data
    }else{
      stopifnot(!is.null(data))
    }
    method <- tolower(method)

    # alias
    if(method == "regression") {
      method <- "ebm"
    }else if(method == "bartlett" || method == "bartlet") {
      method <- "ml"
    }

    # normal case
    if(all(data@ov$type == "numeric")) {
      if(method == "ebm"){
        out <- predict_eta_normal(object = object, model = model, data = data,
                                  implied = implied, samplestats = samplestats,
                                  data.obs = data.obs)
      }else if(method == "ml"){
        out <- predict_eta_bartlett(object = object, model = model, data = data,
                                    implied = implied, samplestats = samplestats,
                                    data.obs = data.obs)
      }else{
        stop("penfa ERROR: unknown method: ", method)
      }
    }else{
      stop("penfa ERROR: only numerical observed variables are allowed")
    }
    out
}

# predict_eta_normal
predict_eta_normal <- function(object = NULL, model = NULL, data = NULL,
                               samplestats = NULL, implied = NULL,
                               data.obs = NULL){
  # full object
  if(inherits(object, "penfa")) {
    model       <- object@Model
    data        <- object@Data
    samplestats <- object@SampleStats
    implied     <- object@Implied
  }else{
    stopifnot(!is.null(model), !is.null(data), !is.null(samplestats), !is.null(implied))
  }

  if(is.null(data.obs)) {
    data.obs <- data@X
    newdata.flag <- FALSE
  } else {
    newdata.flag <- TRUE
  }

  LAMBDA    <- computeLAMBDA(model = model)
  Sigma.hat <- implied$cov
  Sigma.inv <- lapply(Sigma.hat, MASS::ginv)
  VETA      <- computePHI(model = model)
  EETA      <- computeEETA(model = model, samplestats = samplestats)
  EY        <- computeEY(model = model, samplestats = samplestats)
  FS        <- vector("list", length = data@ngroups)

  for(g in 1:data@ngroups) {
    data.obs.g <- data.obs[[g]]
    VETA.g     <- VETA[[g]]
    EETA.g     <- EETA[[g]]
    LAMBDA.g   <- LAMBDA[[g]]
    EY.g       <- EY[[g]]
    Sigma.inv.g <- Sigma.inv[[g]]
    nfac <- ncol(VETA[[g]])
    if(nfac == 0L){
      FS[[g]] <- matrix(0, data@nobs[[g]], nfac)
      next
    }
    # center data
    Yc <- t(t(data.obs.g) - EY.g)
    # factor score coefficient matrix 'C'
    FSC <- VETA.g %*% t(LAMBDA.g) %*% Sigma.inv.g
    # compute factor scores
    FS.g <- t(FSC %*% t(Yc) + EETA.g)
    FS[[g]] <- FS.g
  }
  FS
}

# computeEETA
computeEETA <- function(model = NULL, GLIST = NULL, samplestats = NULL){
  if (is.null(GLIST))
    GLIST <- model@GLIST

  ngroups <- model@ngroups
  nmat    <- model@nmat
  EETA    <- vector("list", length = ngroups)
  for(g in 1:ngroups){
    mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]
    MLIST       <- GLIST[mm.in.group]
    EETA.g      <- as.vector(matrix(0, ncol(MLIST$lambda), 1L))
    EETA[[g]]   <- EETA.g
  }
  EETA
}

# computeEY
computeEY <- function(model = NULL, GLIST = NULL, samplestats = NULL){
  if(is.null(GLIST))
    GLIST <- model@GLIST

  ngroups <- model@ngroups
  nmat    <- model@nmat
  EY <- vector("list", length = ngroups)
  for(g in 1:ngroups){
    mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]
    MLIST <- GLIST[mm.in.group]
    EY[[g]] <- computeEY.LISREL(MLIST = MLIST, sample.mean = samplestats@mean[[g]])
  }
  EY
}

# computeEY.LISREL
computeEY.LISREL <- function (MLIST = NULL, sample.mean = NULL){
  LAMBDA <- MLIST$lambda
  if(!is.null(MLIST$tau)){
    TAU <- MLIST$tau
  }else{
    TAU <- sample.mean
  }
  LAMBDA <- MLIST$lambda
  ALPHA  <- matrix(0, ncol(LAMBDA), 1L)
  EETA   <- as.vector(ALPHA)
  EY     <- as.vector(TAU) + as.vector(LAMBDA %*% EETA)
  EY
}

# predict_eta_bartlett
# factor scores - normal case - Bartlett method
# NOTES: 1) classic 'Bartlett' method; for the linear/continuous
#           case, this is equivalent to 'ML'
#        2) usual formula is:
#               FSC = solve(lambda' theta.inv lambda) (lambda' theta.inv)
#           BUT to deal with singular THETA (with zeroes on the diagonal),
#           we use the 'GLS' version instead:
#               FSC = solve(lambda' sigma.inv lambda) (lambda' sigma.inv)
#           Reference: Bentler & Yuan (1997) 'Optimal Conditionally Unbiased
#                      Equivariant Factor Score Estimators'
#                      in Berkane (Ed) 'Latent variable modeling with
#                      applications to causality' (Springer-Verlag)
#        3) instead of solve(), use MASS::ginv, for special settings where
#           -by construction- (lambda' sigma.inv lambda) is singular
#           note: this will destroy the conditionally unbiased property
#                 of Bartlett scores
predict_eta_bartlett <- function(object = NULL, model = NULL, data = NULL,
                                 samplestats = NULL, implied = NULL,
                                 data.obs = NULL){

  # full object?
  if(inherits(object, "penfa")) {
    model       <- object@Model
    data        <- object@Data
    samplestats <- object@SampleStats
    implied     <- object@Implied
  }else{
    stopifnot(!is.null(model), !is.null(data), !is.null(samplestats), !is.null(implied))
  }

  if(is.null(data.obs)) {
    data.obs <- data@X
  } else {
    newdata.flag <- TRUE
  }

  LAMBDA    <- computeLAMBDA(model = model)
  Sigma.hat <- implied$cov
  Sigma.inv <- lapply(Sigma.hat, MASS::ginv)
  VETA      <- computePHI(model = model)
  EETA      <- computeEETA(model = model, samplestats = samplestats)
  EY        <- computeEY(model = model, samplestats = samplestats)
  FS        <- vector("list", length = data@ngroups)

  for(g in 1:data@ngroups){

      data.obs.g  <- data.obs[[g]]
      VETA.g      <- VETA[[g]]
      EETA.g      <- EETA[[g]]
      LAMBDA.g    <- LAMBDA[[g]]
      EY.g        <- EY[[g]]
      Sigma.inv.g <- Sigma.inv[[g]]
      nfac        <- length(EETA[[g]])
      if(nfac == 0L) {
        FS[[g]] <- matrix(0, data@nobs[[g]], nfac)
        next
      }

    # center data
    Yc <- t( t(data.obs.g) - EY.g)
    # factor score coefficient matrix 'C'
    FSC <- (MASS::ginv(t(LAMBDA.g) %*% Sigma.inv.g %*% LAMBDA.g) %*% t(LAMBDA.g) %*% Sigma.inv.g)
    # compute factor scores
    FS[[g]] <- t(FSC %*% t(Yc) + EETA.g)
  }
  FS
}
