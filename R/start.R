###########################
# ---- Starting values ----
###########################
# Content:
# - get_start : provides starting values for model parameters
# - cfa_1fac_fabin : starting values for 1-factor model according to FABIN (Hagglund, 1982)
# - cfa_1fac_3ind : starting values for 1-factor model wih p = 3, no iterations, solved analytically
# - start_check_cov : sanity check: (user-specified) variances smaller than covariances
#
# Some of these functions are adaptations from the routines present in the lavaan package (https://CRAN.R-project.org/package=lavaan)
#
# --------------------------------------------------------------------------------------------------------------

# get_start
get_start <- function(partable = NULL, samplestats  = NULL, debug = FALSE) {

  # check arguments
  stopifnot(is.list(partable))

  start.user    <- NULL

  # check model list elements, if provided
  if(!is.null(start.user)) {
    if(is.null(start.user$lhs) || is.null(start.user$op) || is.null(start.user$rhs)) {
      stop("penfa ERROR: problem with start argument: model list does not contain all elements: lhs/op/rhs")
    }
    if(!is.null(start.user$est)) {
      # an est column; nothing to do
    } else if(!is.null(start.user$start)) {
      # no est column, but we use the start column
      start.user$est <- start.user$start
    } else if(!is.null(start.user$ustart)) {
      # not ideal
      start.user$est <- start.user$ustart
    } else {
      stop("penfa ERROR: problem with start argument: could not find est/start column in model list")
    }
  }

  # global settings
  # 0. everyting is zero
  start <- numeric( length(partable$ustart) )

  # 1. =~ factor loadings: factor loadings get 1
  start[ which(partable$op == "=~") ] <- 1.0

  # 2. ~~ variances and covariances of latent variables get 0.05
  lv.names    <- partable_vnames(partable, "lv") # all groups
  lv.var.idx <- which(partable$op == "~~" & partable$lhs %in% lv.names & partable$lhs == partable$rhs)
  start[lv.var.idx] <- 0.05

  # group-specific settings
  ngroups <- partable_ngroups(partable)
  if (is.null(partable$group) && ngroups == 1L) {
    partable$group <- rep(1L, length(partable$lhs))
  }

  for(g in 1:ngroups) {

    # group values
    group.values <- partable_group_values(partable)

    # info from user model for this group
    ov.names     <- partable_vnames(partable, "ov", group = group.values[g])
    ov.names.num <- ov.names
    lv.names    <- partable_vnames(partable, "lv",   group = group.values[g])

    # (residual ov) unique variances
    ov.var.idx <- which(partable$group == group.values[g] &
                        partable$op    == "~~"            &
                        partable$lhs %in% ov.names.num    &
                        partable$lhs == partable$rhs)
    sample.var.idx <- match(partable$lhs[ov.var.idx], ov.names)
    start[ov.var.idx] <- (1.0 - 0.50)* samplestats@var[[g]][sample.var.idx]
    # var is rescaled by N

    # 1-fac measurement model: loadings, phi, psi (MPLUS, CFA)
    # fabin3 estimator (2sls) of Hagglund (1982) per factor
    for(f in lv.names) {
      lambda.idx <- which(partable$lhs == f & partable$op == "=~" & partable$group == group.values[g])
      std.lv <- FALSE # standardized?
      var.f.idx <- which(partable$lhs == f & partable$op == "~~" &
                         partable$group == group.values[g] &
                         partable$rhs == f)
      if(length(var.f.idx) > 0L && partable$free[var.f.idx] == 0 && partable$ustart[var.f.idx] == 1) {
        std.lv <- TRUE
      }

      # no second order
      if(any(partable$rhs[lambda.idx] %in% lv.names)) next

      # get observed indicators for this latent variable
      ov.idx <- match(partable$rhs[lambda.idx], ov.names)
      if(length(ov.idx) > 0L && !any(is.na(ov.idx))) {
        COV <- samplestats@cov[[g]][ov.idx, ov.idx, drop = FALSE]

        # fabin for 1-factor
        fabin <- cfa_1fac_fabin(COV, std.lv = std.lv, lambda.only = TRUE)

        # factor loadings
        tmp <- fabin$lambda
        tmp[ !is.finite(tmp) ] <- 1.0 # just in case (e.g. 0/0)
        start[lambda.idx] <- tmp
      }

    } # end factors


    # 3g) intercepts/means
    ov.int.idx <- which(partable$group == group.values[g] & partable$op == "~1" &
                        partable$lhs %in% ov.names)
    sample.int.idx <- match(partable$lhs[ov.int.idx], ov.names)
    start[ov.int.idx] <- samplestats@mean[[g]][sample.int.idx]

    # 6b. exogenous lv variances if single indicator
    lv.x <- partable_vnames(partable, "lv", group = group.values[g]) # "lv.x"
    lv.x <- unique(unlist(lv.x))
    if(length(lv.x) > 0L) {
      for(ll in lv.x) {
        ind.idx <- which(partable$op == "=~" & partable$lhs == ll & partable$group == group.values[g])
        if(length(ind.idx) == 1L) {
          single.ind <- partable$rhs[ind.idx]
          single.fvar.idx <- which(partable$op == "~~" & partable$lhs == ll &
                                   partable$rhs == ll & partable$group == group.values[g])
          single.var.idx <- which(partable$op == "~~" & partable$lhs == single.ind &
                                  partable$rhs == single.ind & partable$group == group.values[g])
          single.var <- partable$ustart[single.var.idx[1]]
          if(is.na(single.var)) {
            single.var <- 1
          }
          ov.idx <- match(single.ind, ov.names)
          ov.var <- diag(samplestats@cov[[g]])[ov.idx]
          tmp <- (1 - (single.var/ov.var)) * ov.var
          # just in case
          if(is.na(tmp) || tmp < 0.05) {
            tmp <- 0.05
          }
          start[single.fvar.idx] <- tmp
        }
      }
    }


  } # groups

  # override if a user list with starting values is provided, only look at the 'est' column for now
  if(!is.null(start.user)) {
    if(is.null(partable$group)) {
      partable$group <- rep(1L, length(partable$lhs))
    }
    if(is.null(start.user$group)) {
      start.user$group <- rep(1L, length(start.user$lhs))
    }
    for(i in 1:length(partable$lhs)) {
      # find corresponding parameters
      lhs <- partable$lhs[i]
      op  <- partable$op[i]
      rhs <- partable$rhs[i]
      grp <- partable$group[i]

      start.user.idx <- which(start.user$lhs == lhs & start.user$op  ==  op &
                              start.user$rhs == rhs & start.user$group == grp)
      if(length(start.user.idx) == 1L && is.finite(start.user$est[start.user.idx])) {
        start[i] <- start.user$est[start.user.idx]
      }
    }
  }

  # override if the model syntax contains explicit starting values
  user.idx <- which(!is.na(partable$ustart))
  start[user.idx] <- partable$ustart[user.idx]

  # final check: no NaN or other non-finite values
  bad.idx <- which(!is.finite(start))
  if(length(bad.idx) > 0L) {
    cat("starting values:\n")
    print( start )
    warning("penfa WARNING: some starting values are non-finite; replacing them with 0.5; please provide better starting values.\n")
    start[ bad.idx ] <- 0.5
  }

  if(debug) {
    cat(" DEBUG: Start\n")
    print( start )
  }

  start
}

# cfa_1fac_fabin
cfa_1fac_fabin <- function(S, lambda.only = FALSE, std.lv = FALSE, extra = NULL) {

  # check arguments
  if(std.lv) {
    lambda.only = FALSE # we need phi
  }
  nvar <- NCOL(S)

  # catch nvar < 4
  if(nvar < 4L) {
    out <- cfa_1fac_3ind(sample.cov = S, std.lv = std.lv, warn.neg.triad = FALSE)
    return(out)
  }

  # 1. lambda
  lambda <- numeric( nvar ); lambda[1L] <- 1.0
  for(i in 2:nvar) {
    idx3 <- (1:nvar)[-c(i, 1L)]
    s23 <- S[i, idx3]
    S31 <- S13 <- S[idx3, 1L]
    S33 <- S[idx3,idx3]
    tmp <- try(solve(S33, S31), silent = TRUE) # GaussJordanPivot slighty more efficient
    if(inherits(tmp, "try-error")) {
      lambda[i] <- sum(s23 * S31) / sum(S13^2)
    } else {
      lambda[i] <- sum(s23 * tmp) / sum(S13 * tmp)
    }
  }

  if(lambda.only) {
    return(list(lambda = lambda, phi = as.numeric(NA), psi = rep(as.numeric(NA), nvar)))
  }

  # 2. psi
  D <- tcrossprod(lambda) / sum(lambda^2)
  psi <- solve(diag(nvar) - D*D, diag(S - (D %*% S %*% D)))

  # 3. phi
  S1 <- S - diag(psi)
  l2 <- sum(lambda^2)
  phi <- sum(colSums(as.numeric(lambda) * S1) * lambda) / (l2 * l2)

  # std.lv?
  if(std.lv) {
    # we allow for negative phi
    lambda <- lambda * sign(phi) * sqrt(abs(phi))
    phi <- 1
  }

  list(lambda = lambda, psi = psi, phi = phi)
}

# cfa_1fac_3ind
cfa_1fac_3ind <- function(sample.cov, std.lv = FALSE, warn.neg.triad = TRUE) {
  # check sample cov
  stopifnot(is.matrix(sample.cov))
  nRow <- NROW(sample.cov); nCol <- NCOL(sample.cov)
  stopifnot(nRow == nCol, nRow < 4L, nCol < 4L)
  nvar <- nRow

  # we expect a 3x3 sample covariance matrix
  # however, if we get a 2x2 (or 1x1 covariance matrix), do something anyway
  if(nvar == 1L) {
    # lambda = 1, psi = 0, phi = sample.cov[1,1]
    # lambda = 1, psi = 0, phi = 1
    sample.cov <- matrix(1, 3L, 3L) * 1.0
  } else if(nvar == 2L) {
    mean.2var <- mean(diag(sample.cov))
    max.var <- max(diag(sample.cov))
    extra <- c(mean.2var, sample.cov[2,1])
    sample.cov <- rbind( cbind(sample.cov, extra, deparse.level = 0),
                         c(extra, max.var), deparse.level = 0)
  }
  s11 <- sample.cov[1,1]; s22 <- sample.cov[2,2]; s33 <- sample.cov[3,3]
  stopifnot(s11 > 0, s22 > 0, s33 > 0)

  s21 <- sample.cov[2,1]; s31 <- sample.cov[3,1]; s32 <- sample.cov[3,2]
  # note: s21*s31*s32 should be positive!
  if(s21 * s31 * s32 < 0 && warn.neg.triad) {
    warning("penfa WARNING: product of the three covariances is negative!")
  }

  # first, we assume l1 = 1
  phi    <- (s21*s31)/s32
  l1     <- 1
  l2     <- s32/s31 # l2 <- s21/psi
  l3     <- s32/s21 # l3 <- s31/psi
  psi1   <- s11 - phi
  psi2   <- s22 - l2*l2*phi
  psi3   <- s33 - l3*l3*phi
  lambda <- c(l1, l2, l3)
  psi    <- c(psi1, psi2, psi3)

  # std.lv?
  if(std.lv) {
    # allow for negative phi
    lambda <- lambda * sign(phi) * sqrt(abs(phi))
    phi <- 1
  }

  # special cases
  if(nvar == 1L) {
    lambda <- lambda[1]
    psi <- psi[1]
  } else if(nvar == 2L) {
    lambda <- lambda[1:2]
    psi <- psi[1:2]
    phi <- phi / 2 # smaller works better?
  }

  list(lambda = lambda, psi = psi, phi = phi)
}

# start_check_cov
start_check_cov <- function(partable = NULL, start = partable$start) {

  nblocks <-  partable_ngroups(partable)

  for(g in 1:nblocks) {
    # block values
    block.values <- partable_group_values(partable)

    # collect all non-zero covariances
    cov.idx <- which(partable$op == "~~" & partable$block == block.values[g] &
                     partable$lhs != partable$rhs & start != 0)

    # for each covariance, use corresponding variances to standardize;
    # the end result should not exceed abs(1)
    for(cc in seq_along(cov.idx)) {
      this.cov.idx <- cov.idx[cc]

      # find corresponding variances
      var.lhs <- partable$lhs[this.cov.idx]
      var.rhs <- partable$rhs[this.cov.idx]

      var.lhs.idx <- which(partable$op == "~~" &
                           partable$block == block.values[g] &
                           partable$lhs == var.lhs &
                           partable$lhs == partable$rhs)

      var.rhs.idx <- which(partable$op == "~~" &
                           partable$block == block.values[g] &
                           partable$lhs == var.rhs &
                           partable$lhs == partable$rhs)

      var.lhs.value <- start[var.lhs.idx]
      var.rhs.value <- start[var.rhs.idx]

      block.txt <- ""
      if(nblocks > 1L) {
        block.txt <- paste(" [in group ", g, "]", sep = "")
      }

      # check for zero variances
      if(var.lhs.value == 0 || var.rhs.value == 0) {
        # this can only happen if it is user-specified
        # cov.idx free? set it to zero
        if(start[this.cov.idx] == 0) {
          # nothing to do
        } else if(partable$free[this.cov.idx] > 0L) {
          warning(
            "penfa WARNING: non-zero covariance element set to zero, due to fixed-to-zero variances\n",
            " variables involved are: ", var.lhs, " ", var.rhs, block.txt)
          start[this.cov.idx] <- 0
        } else {
          stop("penfa ERROR: please provide better fixed values for (co)variances;\n",
               " variables involved are: ", var.lhs, " ", var.rhs, block.txt)
        }
        next
      }

      # which is the smallest? abs() in case of negative variances
      if(abs(var.lhs.value) < abs(var.rhs.value)) {
        var.min.idx <- var.lhs.idx
        var.max.idx <- var.rhs.idx
      } else {
        var.min.idx <- var.rhs.idx
        var.max.idx <- var.lhs.idx
      }

      # check
      COR <- start[this.cov.idx] / sqrt(var.lhs.value * var.rhs.value)
      # NOTE: treat this as an unconditional COR

      if(!is.finite(COR)) {
        # force simple values
        warning("penfa WARNING: starting values imply NaN for a correlation value;\n",
                " variables involved are: ", var.lhs, " ", var.rhs, block.txt)
        start[var.lhs.idx] <- 1
        start[var.rhs.idx] <- 1
        start[this.cov.idx] <- 0
      } else if(abs(COR) > 1) {
        txt <- paste("penfa WARNING: starting values imply a correlation larger than 1;\n",
                     " variables involved are: ", var.lhs, " ", var.rhs, block.txt)

        # three ways to fix it: rescale cov12, var1 or var2

        # we prefer a free parameter, and not user-specified
        if(partable$free[this.cov.idx] > 0L &&
           is.na(partable$ustart[this.cov.idx])) {
          warning(txt)
          start[this.cov.idx] <- start[this.cov.idx] / (COR * 1.1)
        } else if(partable$free[var.min.idx] > 0L &&
                  is.na(partable$ustart[var.min.idx])) {
          warning(txt)
          start[var.min.idx] <- start[var.min.idx] * (COR * 1.1)^2
        } else if(partable$free[var.max.idx] > 0L &&
                  is.na(partable$ustart[var.max.idx])) {
          warning(txt)
          start[var.max.idx] <- start[var.max.idx] * (COR * 1.1)^2

          # not found? try just a free parameter
        } else if (partable$free[this.cov.idx] > 0L) {
          warning(txt)
          start[this.cov.idx] <- start[this.cov.idx] / (COR * 1.1)
        } else if( partable$free[var.min.idx] > 0L) {
          warning(txt)
          start[var.min.idx] <- start[var.min.idx] * (COR * 1.1)^2
        } else if( partable$free[var.max.idx] > 0L) {
          warning(txt)
          start[var.max.idx] <- start[var.max.idx] * (COR * 1.1)^2

          # nothing? warn (and fail later...)
        } else {
          warning(txt)
        }
      } # COR > 1
    } # cov.idx
  } # blocks (groups)

  start
}
