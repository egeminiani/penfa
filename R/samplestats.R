#######################################################
# ---- Sample statistics & Unrestricted (h1) model ----
#######################################################
# Content:
# - samplestats_from_data : compute group-specific sample statistics and related quantities
# - samplestats_icov : compute the inverse of a covariance matrix
# - inv.chol : inverts a positive-definite symmetric matrix (e.g. cov matrix) using Choleski decomposition
#              and returns the log determinant as attribute
#
# - h1_implied_logl : compute the sample statistics for the unrestricted (h1) model and the logl
# - h1_logl : compute the group-specific logl for the unrestricted (h1) model
# - mvnorm_h1_loglik_samplestats : multivariate normal distribution, unrestricted (h1); everything
#                   evalued under the MLEs (Mu = ybar, Sigma = S); input are logdet, N and P
#
# Some of these functions are adaptations from the routines present in the lavaan package (https://CRAN.R-project.org/package=lavaan)
#
# --------------------------------------------------------------------------------------------------------------

# sample_stats_from_data
samplestats_from_data <- function(moddata = NULL, meanstructure = FALSE,
                                  debug = FALSE, verbose = FALSE) {
  # check moddata
  stopifnot(!is.null(moddata))

  X       <- moddata@X
  ngroups <- moddata@ngroups
  nobs    <- moddata@nobs

  cov         <- vector("list", length = ngroups)
  var         <- vector("list", length = ngroups)
  mean        <- vector("list", length = ngroups)
  icov        <- vector("list", length = ngroups)
  cov.log.det <- vector("list", length = ngroups)
  group.w     <- vector("list", length = ngroups)

  # compute some sample statistics per group
  for(g in 1:ngroups) {

    # check nobs
    if(nobs[[g]] < 2L) {
      if(nobs[[g]] == 0L) {
        stop("penfa ERROR: data set contains no observations",
             ifelse(ngroups > 1L, paste(" in group ", g, sep=""), ""))
      } else {
        stop("penfa ERROR: data set contains only a single observation",
             ifelse(ngroups > 1L, paste(" in group ", g, sep=""), ""))
      }
    }

    # group weight
    group.w[[g]] <- nobs[[g]] / sum(unlist(nobs))

    # 'transform' the sample cov (divided by n-1) to a sample cov divided by 'n'
    cov[[g]]  <-   (nobs[[g]]-1)/nobs[[g]] * stats::cov(X[[g]], use = "pairwise")
    var[[g]]  <-   diag(cov[[g]]) # also var is rescaled now
    mean[[g]] <-   colMeans(X[[g]], na.rm=TRUE)

    # icov and cov.log.det
    out              <- samplestats_icov(COV = cov[[g]], ngroups = ngroups, g = g, warn = TRUE)
    icov[[g]]        <- out$icov
    cov.log.det[[g]] <- out$cov.log.det

  } # ngroups

  # construct SampleStats object
  SampleStats <- new("penfaSampleStats",
                      mean         = mean,
                      cov          = cov,
                      var          = var,
                      group.w      = group.w,
                      nobs         = nobs,
                      ntotal       = sum(unlist(nobs)),
                      ngroups      = ngroups,
                      icov         = icov,
                      cov.log.det  = cov.log.det)
  SampleStats
}

# samplestats_icov
samplestats_icov <- function(COV = NULL, ngroups = 1L, g = 1L, warn = TRUE) {

  tmp <- try(inv.chol(COV, logdet = TRUE), silent = TRUE)
  if(inherits(tmp, "try-error")) {
    stop("penfa ERROR: sample covariance matrix is not positive-definite")
  } else {
    cov.log.det <- attr(tmp, "logdet")
    attr(tmp, "logdet") <- NULL
    icov <- tmp
  }
  list(icov = icov, cov.log.det = cov.log.det)
}

# inv.chol
inv.chol <- function(S, logdet=FALSE) {
  cS <- chol(S)
  S.inv <- chol2inv( cS )
  if(logdet) {
    diag.cS <- diag(cS)
    attr(S.inv, "logdet") <- sum(log(diag.cS*diag.cS))
  }
  S.inv
}



# h1_implied_logl
h1_implied_logl <- function(moddata = NULL, samplestats = NULL, options = NULL) {
  # complete data
  implied <- list(cov = samplestats@cov, mean  = samplestats@mean)
  logl    <- h1_logl(moddata = moddata, samplestats = samplestats, options = options)
  list(implied = implied, logl = logl)
}

# h1_logl
h1_logl <- function(moddata = NULL, samplestats = NULL, options = NULL) {

  # number of groups
  ngroups <- moddata@ngroups
  logl.group <- rep(as.numeric(NA), ngroups)

  logl.ok <- FALSE
  # check if everything is numeric, OR if we have factors with 2 levels only
  if(all(moddata@ov$type == "numeric")) {
    logl.ok <- TRUE
  }else {
    not.idx <- which(moddata@ov$type != "numeric")
    for(i in not.idx) {
      if(moddata@ov$type[i] == "factor") {
        logl.ok <- TRUE
      } else {
        logl.ok <- FALSE
        break
      }
    }
  }

  if(logl.ok) {
    for(g in seq_len(ngroups) ) {
      # we need logdet of covariance matrix, nobs and nvar
      logl.group[g] <- mvnorm_h1_loglik_samplestats(sample.cov.logdet = samplestats@cov.log.det[[g]],
                                                    sample.nvar       = NCOL(samplestats@cov[[g]]),
                                                    sample.nobs       = samplestats@nobs[[g]])
    } # g
  } # logl.ok is TRUE

  out <- list(loglik  = sum(logl.group), loglik.group = logl.group)
  out
}

# mvnorm_h1_loglik_samplestats
mvnorm_h1_loglik_samplestats <- function(sample.cov.logdet = NULL, sample.nvar = NULL, sample.nobs = NULL,
                                         sample.cov = NULL, Sinv.method = "eigen"){

  if(is.null(sample.nvar)) {
    P <- NCOL(sample.cov)
  } else {
    P <- sample.nvar # number of variables
  }

  N <- sample.nobs
  stopifnot(!is.null(P), !is.null(N))

  LOG.2PI <- log(2 * pi)
  # we need logdet
  if(is.null(sample.cov.logdet)) {
    sample.cov.inv <- matrix_symmetric_inverse(S = sample.cov, logdet = TRUE, Sinv.method = Sinv.method)
    logdet <- attr(sample.cov.inv, "logdet")
  } else {
    logdet <- sample.cov.logdet
  }

  loglik <- -N/2 * (P * LOG.2PI + logdet + P)
  loglik
}

