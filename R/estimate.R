################################
# ---- Penalized Estimation ----
################################
# Content:
# - model_estimate : penalized-likelihood estimation of factor models via trust-region algorithm
#
# objective function computation
# - model_objective : returns the value of the objective function (negative logl) for the unpenalized model
# - estimator.ML : computes the group-specific value of the objective function (without multiplication by N/2)
#
# gradient computation
# - model_gradient : returns the gradient vector of the objective function
# - computeOmega : computes the Omega matrix; the expression differs with or without meanstructure
# - derivative.F.LISREL : computes and vectorizes the first derivatives for each parameter matrix
#
# Fisher information computation
# - model_information : computes expected information matrix (UNIT, not total, information)
# - computeDelta : computes the Delta matrix
# - derivative.sigma.LISREL : computes dSigma/dx per model matrix
# - derivative.mu.LISREL : computes dMu/dx per model matrix
# - model_information_invert : takes the inverse of the information matrix
# - model_h1_information_expected : expected fisher information matrix of the unrestricted (h1) model
# - mvnorm_information_expected: unit expected information h0 Mu and vech(Sigma)
#
# Hessian matrix computation
# - model_hessian : computes the Hessian matrix of the model
# - hessian.analytical : computes the analytical Hessian (no meanstructure allowed)
# - dx_Lambda, dx_Phi, dx_Psi, dx_Lambda_Phi, dx_Lambda_Psi, dx_Phi_Psi: compute analytic second-order
#                 derivatives with respect to all parameter matrices
# - set_label_lambda_noconstr, set_label_lambda_constr, get_lambda_constr_row_idx, set_label_phivcov_constr,
#   set_label_phiv_constr, set_label_phiv_constr, set_label_phicov_constr, set_label_phi_noconstr,
#   set_label_phi_constr get_phi_constr_row_idx order_idx_matrix : set the labels for Lambda and Phi
# - get_label_lambda, get_label_phi, get_label_psi : get labels sets for Lambda, Phi and Psi
#
# Implied
# - model_implied : computes model-implied moments per group
# - model_loglik : computes the log-likelihood of the data, given the model
# - mvnorm_loglik_samplestats : logl with sample statistics (mean, cov, N) only as input
#
# Some of these functions are adaptations from the routines present in the lavaan package (https://CRAN.R-project.org/package=lavaan)
#----------------------------------------------------------------------------------------------------------------

model_estimate <- function(model = NULL, samplestats = NULL, moddata = NULL,
                           modpenalty = NULL,  options = NULL) {

  verbose       <- options$verbose
  debug         <- options$debug
  information   <- options$information
  N             <- samplestats@ntotal
  counter.npd   <- 0

  # starting values for optimizer
  start.x <- model_get_parameters(model)

  if(debug) {
   cat("start.x = ", start.x, "\n")
  }

  # Objective function
  objfun <- function(x) {

    # Penalty matrix
    info <- penalty_mat(x = x, model = model, modpenalty = modpenalty, options = options, N = N)
    S.h  <- info$S.h
    S.h.shrink <- info$Ss.shrink
    S.h.diff   <- info$Ss.diff

    # update GLIST (change `state') and make a copy
    GLIST <- model_x2GLIST(model, x = x)

    ### Value
    # objective function value for unpenalized model
    fx <- model_objective(model = model, GLIST = GLIST, samplestats = samplestats,
                          moddata = moddata, verbose = verbose, debug = debug)
    S.res <- fx
    S.h1  <- 0.5 * crossprod(x, S.h) %*% x
    # penalized objective function value
    fx.pen <- fx + S.h1

    if(debug) {
      cat("Objective function  = ", sprintf("%18.16f", fx.pen), "\n", sep="")
      cat("Current parameter values =\n"); print(x); cat("\n")
    }

    # Record how frequently the objective function becomes non-finite
    if(!is.finite(fx)){
      counter.npd <<- counter.npd + 1
      return(list(value = fx))
    }

    ### Gradient
    # gradient for unpenalized model
    dx <- model_gradient(model = model, GLIST = GLIST, samplestats = samplestats,
                         moddata = moddata, verbose = verbose)
    S.h2   <- S.h %*% x
    # penalized gradient
    dx.pen <- dx + S.h2

    if(debug) {
      cat("Gradient function =\n"); print(dx.pen); cat("\n")
    }

    ### Information matrix
    # Hessian matrix or Expected Fisher information for unpenalized model
    if(information == "fisher"){
      hessian <- model_information(model = model, GLIST = GLIST, samplestats = samplestats)
    }else if (information == "hessian"){
      hessian <- model_hessian(model = model, GLIST = GLIST, samplestats = samplestats)
    }
    # penalized Hessian/Fisher
    hessian.pen <- hessian + S.h

    if(debug) {
      cat(information, "matrix =\n"); print(hessian.pen); cat("\n")
    }

    # Collect output for trust-region algorithm
    output <- list(value       = fx.pen,
                   gradient    = dx.pen,
                   hessian     = hessian.pen,
                   counter.npd = counter.npd,
                   S.h1        = S.h1,
                   S.h2        = S.h2,
                   S.h         = S.h,
                   l           = S.res,
                   SS.shrink   = S.h.shrink,
                   SS.diff     = S.h.diff)
    return(output)
  }

  # Parameter scaling: if start.x > 1.0, rescale by using 1/start.x
  SCALE <- rep(1.0, length(start.x))
  idx   <- which(abs(start.x) > 1.0)
  if(length(idx) > 0L) {
    SCALE[idx] <- abs(1.0/start.x[idx])
  }
  if(debug) {
    cat("SCALE = ", SCALE, "\n")
  }

  control.trust <- list(rinit=1L, rmax=100L, iterlim=1000L,
                        fterm = sqrt(.Machine$double.eps),
                        mterm = sqrt(.Machine$double.eps))
  control.trust <- utils::modifyList(control.trust, options$control)

  if(debug) cat("Trust-region algorithm steps:\n")

  optim.out <- trust::trust(objfun   = objfun,
                            parinit  = start.x,
                            rinit    = control.trust$rinit,
                            rmax     = control.trust$rmax,
                            parscale = SCALE,
                            iterlim  = control.trust$iterlim,
                            fterm    = control.trust$fterm,
                            mterm    = control.trust$mterm,
                            minimize = TRUE,
                            blather  = verbose)

  iterations <- optim.out$iterations
  x          <- optim.out$argument

  attr(x, "converged")      <- optim.out$converged
  attr(x, "iterations")     <- iterations
  attr(x, "fx.pen")         <- optim.out$value
  attr(x, "dx.pen")         <- optim.out$gradient
  attr(x, "hessian.pen")    <- optim.out$hessian
  attr(x, "l")              <- optim.out$l
  attr(x, "S.h1")           <- optim.out$S.h1
  attr(x, "S.h2")           <- optim.out$S.h2
  attr(x, "S.h")            <- optim.out$S.h
  attr(x, "SS.shrink")      <- optim.out$SS.shrink
  attr(x, "SS.diff")        <- optim.out$SS.diff
  attr(x, "control")        <- control.trust
  attr(x, "npd")            <- optim.out$counter.npd; rm(counter.npd)

  if(verbose){
    attr(x, "argpath")   <- optim.out$argpath
    attr(x, "argtry")    <- optim.out$argtry
    attr(x, "type")      <- optim.out$steptype
    attr(x, "accept")    <- optim.out$accept
    attr(x, "radii")     <- optim.out$r
    attr(x, "rho")       <- optim.out$rho
    attr(x, "fx.val")    <- optim.out$valpath
    attr(x, "fx.valtry") <- optim.out$valtry
    attr(x, "change")    <- optim.out$preddiff
    attr(x, "stepnorm")  <- optim.out$stepnorm
  }

  x
}

# model_objective
model_objective <- function(model = NULL, GLIST = NULL, samplestats = NULL,
                            moddata = NULL, verbose  = FALSE, debug = FALSE) {

  # state or final?
  if(is.null(GLIST)) GLIST <- model@GLIST

  meanstructure <- model@meanstructure
  num.idx       <- model@num.idx

  # Compute Sigma.hat
  Sigma.hat <- computeSigmaHat(model = model, GLIST = GLIST)

  # Check if the matrix is positive definite
  is.pdef <- is.Pdef(mat = Sigma.hat, ngroups = samplestats@ngroups)

  if(meanstructure) {
    Mu.hat <- computeMuHat(model = model, GLIST = GLIST)
  }

  fx <- 0.0
  fx.group   <- numeric(samplestats@ngroups)
  logl.group <- rep(as.numeric(NA), samplestats@ngroups)

  # Compute Sigma.hat.inv if possible
  for(g in 1:samplestats@ngroups){
    if(is.pdef[g]){
      Sigma.hat.inv     <- compute.inv(Sigma.hat[[g]], pdef = is.pdef[g])
      Sigma.hat.log.det <- attr(Sigma.hat.inv, "logdet")
      attr(Sigma.hat[[g]], "po")      <- TRUE
      attr(Sigma.hat[[g]], "inv")     <- Sigma.hat.inv
      attr(Sigma.hat[[g]], "log.det") <- Sigma.hat.log.det
    }else{
      attr(Sigma.hat[[g]], "po") <- FALSE
      fx <- Inf
      return(fx)
    }
  }# group

  for(g in 1:samplestats@ngroups){
    # complete data; MLE
    group.fx <- estimator.ML(Sigma.hat        = Sigma.hat[[g]],
                             Mu.hat           = Mu.hat[[g]],
                             data.cov         = samplestats@cov[[g]],
                             data.mean        = samplestats@mean[[g]],
                             data.cov.log.det = samplestats@cov.log.det[[g]],
                             meanstructure    = meanstructure)

    group.fx <- samplestats@nobs[[g]] * 0.5 * group.fx # Ng/2 * f
    fx.group[g] <- group.fx
  } # g

  if(samplestats@ngroups > 1) {
    fx   <- sum(fx.group) # already multiplied by Ng
  } else { # single group
    fx <- fx.group[1]
  }

  fx.value <- as.numeric(fx)
  attr(fx, "fx.group") <- fx.group
  fx
}

# estimator.ML
estimator.ML <- function(Sigma.hat=NULL, Mu.hat=NULL, data.cov=NULL, data.mean=NULL,
                         data.cov.log.det=NULL, meanstructure=FALSE){

  if(!attr(Sigma.hat, "po")) return(Inf) # just to be sure

  Sigma.hat.inv     <- attr(Sigma.hat, "inv")
  Sigma.hat.log.det <- attr(Sigma.hat, "log.det")
  nvar <- ncol(Sigma.hat)

  if(!meanstructure) {
    fx <- (Sigma.hat.log.det + sum(data.cov * Sigma.hat.inv) + nvar*log(2*pi))
  } else {
    W.tilde <- data.cov + tcrossprod(data.mean - Mu.hat)
    fx <- (Sigma.hat.log.det + sum(W.tilde * Sigma.hat.inv) + nvar*log(2*pi))
  }

  # no negative values
  if(is.finite(fx) && fx < 0.0) fx <- 0.0
  fx
}

# model gradient
model_gradient <- function(model = NULL, GLIST  = NULL, samplestats = NULL, moddata  = NULL, verbose = FALSE,
                           m.el.idx = NULL, x.el.idx = NULL){

  nmat           <- model@nmat
  meanstructure  <- model@meanstructure
  num.idx        <- model@num.idx
  nx.free        <- model@nx.free

  # state or final?
  if(is.null(GLIST)) GLIST <- model@GLIST

  group.w <- unlist(samplestats@nobs) # Ng

  Sigma.hat <- computeSigmaHat(model = model, GLIST = GLIST)
  for(g in 1:samplestats@ngroups){
    is.pdef           <- is.Pdef(Sigma.hat[[g]])
    Sigma.hat.inv     <- compute.inv(Sigma.hat[[g]], pdef = is.pdef)
    Sigma.hat.log.det <- attr(Sigma.hat.inv, "logdet")
    attr(Sigma.hat[[g]], "po")      <- is.pdef
    attr(Sigma.hat[[g]], "inv")     <- Sigma.hat.inv
    attr(Sigma.hat[[g]], "log.det") <- Sigma.hat.log.det
  }

  if(meanstructure) {
    Mu.hat <- computeMuHat(model = model, GLIST = GLIST)
  }

  # MLE approach: using Omega and Omega.mu
  #  Omega = Sigma.inv %*% (S - Sigma) %*% t(Sigma.inv)
  if(meanstructure) {
    Omega <- computeOmega(Sigma.hat = Sigma.hat, Mu.hat = Mu.hat,
                          samplestats = samplestats, meanstructure=TRUE)
    Omega.mu <- attr(Omega, "mu")
  }else{
    Omega <- computeOmega(Sigma.hat = Sigma.hat, Mu.hat = NULL,
                          samplestats = samplestats, meanstructure = FALSE)
    Omega.mu <- vector("list", length = model@ngroups)
  }

  # compute DX (for all elements in every model matrix)
  DX <- vector("list", length=length(GLIST))

  for(g in 1:model@ngroups) {
    # which mm belong to group g?
    mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
    mm.names <- names( GLIST[mm.in.group] )
    DX.group <- derivative.F.LISREL(MLIST = GLIST[mm.in.group],
                                    Omega = Omega[[g]], Omega.mu = Omega.mu[[g]])

    # only save what we need
    DX[mm.in.group] <- DX.group[ mm.names ]

    # weight by group
    for(mm in mm.in.group) {
      DX[[mm]] <- group.w[g] * DX[[mm]] # Ng * DX
    }

  } # group

  # extract free parameters
  dx <- numeric( nx.free )
  for(g in 1:model@ngroups) {
    mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
    for(mm in mm.in.group) {
      m.free.idx  <- model@m.free.idx[[mm]]
      x.free.idx  <- model@x.free.idx[[mm]]
      dx[x.free.idx] <- DX[[mm]][m.free.idx]
    }
  }
  dx
}

#computeOmega
computeOmega <- function(Sigma.hat=NULL, Mu.hat=NULL,
                         samplestats=NULL, meanstructure=FALSE) {

  ngroups  <- length(Sigma.hat)
  Omega    <- vector("list", length = ngroups)
  Omega.mu <- vector("list", length = ngroups)

  for(g in 1:ngroups) {
    if(attr(Sigma.hat[[g]], "po") == FALSE) {
      # stop
      warning("penfa WARNING: model_gradient: Sigma.hat is not positive definite\n")
      Sigma.hat.inv <- MASS::ginv(Sigma.hat[[g]])
      Sigma.hat.log.det <- log(.Machine$double.eps)
    } else {
      Sigma.hat.inv <-  attr(Sigma.hat[[g]], "inv")
      Sigma.hat.log.det <- attr(Sigma.hat[[g]], "log.det")
    }

    if(meanstructure){
      diff <- samplestats@mean[[g]] - Mu.hat[[g]]
      W.tilde <- samplestats@cov[[g]] + tcrossprod(diff)
      # Browne 1995 eq 4.55
      Omega.mu[[g]] <- t(t(diff) %*% Sigma.hat.inv)
      Omega[[g]] <- (Sigma.hat.inv %*% (W.tilde - Sigma.hat[[g]]) %*% Sigma.hat.inv )
    } else {
      W.tilde <- samplestats@cov[[g]]
      Omega[[g]] <- (Sigma.hat.inv %*% (W.tilde - Sigma.hat[[g]]) %*% Sigma.hat.inv)
    }

  } # g

  if(meanstructure) attr(Omega, "mu") <- Omega.mu
  Omega
}

# derivative.F.LISREL
derivative.F.LISREL <- function(MLIST=NULL, Omega=NULL, Omega.mu=NULL) {

  LAMBDA <- MLIST$lambda
  PHI    <- MLIST$phi
  KAPPA  <- MLIST$kappa

  # meanstructure?
  meanstructure <- FALSE; if(!is.null(Omega.mu)) meanstructure <- TRUE

  # pre-compute some values
  tLAMBDA <- t(LAMBDA)
  Omega..LAMBDA <- Omega %*% LAMBDA

  # 1. LAMBDA
  if(meanstructure) {
    LAMBDA.deriv <- -1.0 * ( Omega.mu %*% t(KAPPA) + Omega..LAMBDA %*% PHI )
  }else{
    LAMBDA.deriv <- -1.0 * (Omega..LAMBDA %*% PHI)
  }

  # 3. PHI
  PHI.deriv <- -1.0 * ( tLAMBDA %*% Omega %*% LAMBDA )
  diag(PHI.deriv) <- 0.5 * diag(PHI.deriv)

  # 4. PSI
  PSI.deriv <- -1.0 * Omega
  diag(PSI.deriv) <- 0.5 * diag(PSI.deriv)

  if(meanstructure) {
    # 5. TAU
    TAU.deriv <- -1.0 * Omega.mu

    # 6. KAPPA
    KAPPA.deriv <- -1.0 * t( t(Omega.mu) %*% LAMBDA )
  } else {
    TAU.deriv <- NULL
    KAPPA.deriv <- NULL
  }

  list(lambda = LAMBDA.deriv, psi = PSI.deriv, phi = PHI.deriv, tau = TAU.deriv, kappa  = KAPPA.deriv)
}


# model_information
model_information <- function(model = NULL, GLIST = NULL, samplestats = NULL){

  # state or final?
  if(is.null(GLIST)) GLIST <- model@GLIST

  implied <- model_implied(model, GLIST = GLIST) # Sigma_hat and Muhat with the new estimates in GLIST

  # 1. Delta
  Delta <- computeDelta(model = model, GLIST. = GLIST) # pass GLIST!

  # 2. H1 information
  # information = Delta' I1 Delta, where I1 is the information of
  # the saturated model evaluated at the structured estimates
  A1 <- model_h1_information_expected(model = model, samplestats = samplestats, implied = implied)

  # 3. compute Information per group
  Info.group  <- vector("list", length=samplestats@ngroups)
  for(g in 1:samplestats@ngroups) {
    fg <- samplestats@nobs[[g]] # Ng

    # compute information for this group
    Info.group[[g]] <- fg * ( crossprod(Delta[[g]], A1[[g]]) %*% Delta[[g]] )
  }

  # 4. assemble over groups
  Information <- Info.group[[1]]
  if(samplestats@ngroups > 1) {
    for(g in 2:samplestats@ngroups) {
      Information <- Information + Info.group[[g]]
    }
  }

  Information
}

# computeDelta
computeDelta <- function (model = NULL, GLIST. = NULL, m.el.idx. = NULL, x.el.idx. = NULL) {

  nmat    <- model@nmat
  ngroups <- model@ngroups
  nvar    <- model@nvar
  num.idx <- model@num.idx

  # state or final?
  # pass GLIST in computeDelta otherwise estimates are not updated during optimization
  if (is.null(GLIST.))
    GLIST <- model@GLIST
  else GLIST <- GLIST.

  type <- "nonfree"
  m.el.idx <- m.el.idx.
  x.el.idx <- x.el.idx.
  if (is.null(m.el.idx) && is.null(x.el.idx))
    type <- "free"
  pstar <- integer(ngroups)
  for (g in 1:ngroups) {
    pstar[g] <- as.integer(nvar[g] * (nvar[g] + 1)/2)
    if (model@meanstructure) {
      pstar[g] <- nvar[g] + pstar[g]
    }
  }
  if (type == "free") {
    NCOL <- model@nx.free
    m.el.idx <- x.el.idx <- vector("list", length = length(GLIST))
    for (mm in 1:length(GLIST)) {
      m.el.idx[[mm]] <- model@m.free.idx[[mm]]
      x.el.idx[[mm]] <- model@x.free.idx[[mm]]
      if (model@isSymmetric[mm]) {
        dix <- duplicated(x.el.idx[[mm]])
        if (any(dix)) {
          m.el.idx[[mm]] <- m.el.idx[[mm]][!dix]
          x.el.idx[[mm]] <- x.el.idx[[mm]][!dix]
        }
      }
    }
  }
  else {
    NCOL <- sum(unlist(lapply(x.el.idx, function(x) length(unique(x)))))
  }

  Delta <- vector("list", length = ngroups)

  for (g in 1:ngroups) {
    Delta.group <- matrix(0, nrow = pstar[g], ncol = NCOL)
    mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]
    for (mm in mm.in.group) {
      mname <- names(model@GLIST)[mm]
      if (!length(m.el.idx[[mm]]))
        next
      DELTA <- dxSigma <- derivative.sigma.LISREL(m = mname, idx = m.el.idx[[mm]],
                                                  MLIST = GLIST[mm.in.group])
      if (model@meanstructure) {
        DELTA.mu <- derivative.mu.LISREL(m = mname, idx = m.el.idx[[mm]],
                                         MLIST = GLIST[mm.in.group])
        DELTA <- rbind(DELTA.mu, DELTA)
      }
      Delta.group[, x.el.idx[[mm]]] <- DELTA
    }
    Delta[[g]] <- Delta.group
  }
  Delta
}

# derivative.sigma.LISREL
derivative.sigma.LISREL <- function(m = "lambda",
                                    idx = seq_len(length(MLIST[[m]])),
                                    MLIST = NULL, vech  = TRUE){

  LAMBDA <- MLIST$lambda
  nvar   <- nrow(LAMBDA)
  nfac   <- ncol(LAMBDA)
  PHI    <- MLIST$phi

  # only lower.tri part of sigma
  v.idx <- matrix_vech_idx( nvar )
  pstar <- nvar*(nvar+1)/2

  # shortcut for tau and kappa : empty matrix
  if(m == "tau" || m == "kappa") {
    return( matrix(0.0, nrow=pstar, ncol=length(idx)) )
  }

  IB.inv <- diag(nrow(PHI))

  if(m == "lambda") {
    L1 <- LAMBDA %*% IB.inv %*% PHI %*% t(IB.inv)
  }

  if(m == "phi") {
    LAMBDA..IB.inv <- LAMBDA %*% IB.inv
  }

  if(m == "lambda") {
    KOL.idx <- matrix(1:(nvar*nfac), nvar, nfac, byrow = TRUE)[idx]
    DX <- (L1 %x% diag(nvar))[,idx, drop = FALSE] + (diag(nvar) %x% L1)[,KOL.idx, drop = FALSE]
  } else if(m == "phi") {

    DX <- (LAMBDA..IB.inv %x% LAMBDA..IB.inv) # Kronecker product of LAMBDA
    # symmetry correction, but keeping all duplicated elements since we depend on idx=m.el.idx
    lower.idx <- matrix_vech_idx(nfac, diagonal = FALSE)
    upper.idx <- matrix_vechru_idx(nfac, diagonal = FALSE)
    offdiagSum <- DX[,lower.idx] + DX[,upper.idx]
    DX[,c(lower.idx, upper.idx)] <- cbind(offdiagSum, offdiagSum)
    DX <- DX[,idx, drop = FALSE]
  } else if(m == "psi") {

    DX <- matrix(0, nvar*nvar, length(idx))
    DX[cbind(idx,seq_along(idx))] <- 1
    # symmetry correction not needed
  } else {
    stop("wrong model matrix names: ", m, "\n")
  }

  # vech?
  if(vech) {
    DX <- DX[v.idx,, drop=FALSE]	# eliminate rows of redundant elements
  }

  DX
}

# derivative.mu.LISREL
derivative.mu.LISREL <- function(m="kappa", idx=seq_len(length(MLIST[[m]])), MLIST=NULL) {

  LAMBDA <- MLIST$lambda
  nvar <- nrow(LAMBDA)
  nfac <- ncol(LAMBDA)

  # shortcut for empty matrices
  if(m == "phi" || m == "psi") {
    return( matrix(0.0, nrow=nvar, ncol=length(idx) ) )
  }

  # missing kappa
  if(is.null(MLIST$kappa))
    KAPPA <- matrix(0, nfac, 1L)
  else
    KAPPA  <- MLIST$kappa

  IB.inv <- diag(nrow(MLIST$phi))

  if(m == "tau") {
    DX <- diag(nvar)
  } else if(m == "lambda") {
    DX <- t(IB.inv %*% KAPPA) %x% diag(nvar)
  } else if(m == "kappa") {
    DX <- LAMBDA %*% IB.inv
  } else {
    stop("wrong model matrix names: ", m, "\n")
  }

  DX <- DX[, idx, drop=FALSE]
  DX
}

# model_information_invert
model_information_invert <- function(model = NULL,information = NULL,
                                     inverted = FALSE, check.pd = FALSE){

  if(check.pd) {
    eigvals <- eigen(information, symmetric = TRUE, only.values = TRUE)$values
    if(any(eigvals < -1 * .Machine$double.eps^(3/4))) {
      warning("penfa WARNING: information matrix is not positive definite; the model may not be identified")
    }
  }

  if(inverted) {
    information <- try( solve(information), silent = TRUE )
    # Now that is inverted, it is the variance covariance matrix
  }
  # inverted information
  information
}

# model_h1_information_expected
model_h1_information_expected <- function(object = NULL, model = NULL,
                                          samplestats = NULL, implied = NULL){
  # ML single level
  A1 <- vector("list", length=samplestats@ngroups)

  # structured --> compute model-implied statistics
  if(length(implied) == 0L) {
    implied <- model_implied(model)
  }

  for(g in 1:samplestats@ngroups) {
    # complete data
    if(model@meanstructure) {
      MEAN <- implied$mean[[g]]
    } else {
      MEAN <- samplestats@mean[[g]]
    }
    A1[[g]] <- mvnorm_information_expected(Sigma = implied$cov[[g]], meanstructure = model@meanstructure)
  } # g
  A1
}

# mvnorm_information_expected
mvnorm_information_expected <- function(Sigma = NULL, Sinv.method = "eigen",
                                        Sigma.inv = NULL, meanstructure = TRUE){

  if(is.null(Sigma.inv)) {
    # invert Sigma
    Sigma.inv <- matrix_symmetric_inverse(S = Sigma, logdet = FALSE, Sinv.method = Sinv.method)
  }

  I22 <- 0.5 * matrix_duplication_pre_post(Sigma.inv %x% Sigma.inv)

  if(meanstructure) {
    I11 <- Sigma.inv
    out <- matrix_bdiag(I11, I22)
  } else {
    out <- I22
  }

  out
}


# model_hessian
model_hessian <- function(model = NULL, GLIST = NULL, samplestats = NULL){

  nmat <- model@nmat
  ngroups <- samplestats@ngroups

  # compute Hessian per group
  Hess.group  <- vector("list", length=ngroups)
  lab.group   <- vector("list", length=ngroups)
  phiFirst    <- vector("logical",length=ngroups)

  for (g in 1:ngroups) {
    mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]

    # Get labels
    for(mm in mm.in.group){

      if(mm == mm.in.group[1]){ # lambda: 1 (for group 1), 4 (for group 1), etc
        lam.mat <- matrix(0, nrow = length(model@dimNames[[mm]][[1]]),
                          ncol = length(model@dimNames[[mm]][[2]]))
        lam.mat[model@m.free.idx[[mm]]] <- model@m.free.idx[[mm]]
        fix_lambda_item <- which(lam.mat == 0, arr.ind = T)
        fix_lambda <- ifelse(length(fix_lambda_item) > 0L, TRUE, FALSE)

        info.lab.lambda  <- get_label_lambda(restrict_lambda = fix_lambda, restrict_pos = fix_lambda_item,
                                             dim.row = length(model@dimNames[[mm]][[1]]),
                                             dim.col = length(model@dimNames[[mm]][[2]]))

      }else if(mm == mm.in.group[2]){ # psi: 2 (for group 1), 5 (for group 2), etc
        psi.mat <- matrix(0, nrow = length(model@dimNames[[mm]][[1]]), ncol = length(model@dimNames[[mm]][[2]]))
        psi.mat[model@m.free.idx[[mm]]] <- model@m.free.idx[[mm]]
        fix_psi_item <- which(diag(psi.mat) == 0, arr.ind = T)
        fix_psi      <- ifelse(length(fix_psi_item) > 0, TRUE, FALSE)

        info.lab.psi <- get_label_psi(restrict_psi = fix_psi, restrict_pos = fix_psi_item,
                                      dim = length(model@dimNames[[mm]][[1]]))

      }else if(mm == mm.in.group[3]){ # phi: 3 (for group 1), 6 (for group 2)
        phi.mat <- matrix(0, nrow = length(model@dimNames[[mm]][[1]]),
                          ncol = length(model@dimNames[[mm]][[2]]))
        phi.mat[model@m.free.idx[[mm]]] <- model@m.free.idx[[mm]]
        fix_phivar_item <- which(diag(phi.mat) == 0, arr.ind = T)
        fix_phicov_item <- which(phi.mat[lower.tri(phi.mat)] == 0, arr.ind = T)
        fix_phivar <- ifelse(length(fix_phivar_item) > 0, TRUE, FALSE)
        fix_phicov <- ifelse(length(fix_phicov_item) > 0, TRUE, FALSE)

        info.lab.phi <- get_label_phi(restrict_phivar = fix_phivar, restrict_phicov = fix_phicov,
                                      restrict_var_item = fix_phivar_item, restrict_cov_item = fix_phicov_item,
                                      dim = length(model@dimNames[[mm]][[1]]))
      }

    }# end matrices

    lab.group[[g]] <- list("label.load" = info.lab.lambda$label.load,
                           "label.phi" = info.lab.phi$label_phi,
                           "label.psi" = info.lab.psi$label_psi)

    # Detect which comes first, Phi or Psi
    copy        <- model@x.free.idx[mm.in.group]
    names(copy) <- names(GLIST)[mm.in.group] # the order of the slots is always lambda, psi, phi
    tmp         <- unlist(lapply(copy, function(x) x[1]))
    index       <- names(tmp)[order(tmp)]
    phi.idx     <- match(c("phi"), index)
    psi.idx     <- match(c("psi"), index)
    phiFirst[g] <- ifelse(phi.idx < psi.idx, TRUE, FALSE)

    # Compute Hessian
    # only pass parameter matrices of the correct group
    # pass an index which tells whether Phi or Psi comes first
    # for convenience only pass the arguments that we need, already filtered by group
    Hess.group[[g]] <- hessian.analytical(GLIST = GLIST[mm.in.group],
                                          N = samplestats@nobs[[g]], # Ng
                                          data.cov = samplestats@cov[[g]], # S_g
                                          Sigma = model_implied(model, GLIST = GLIST)$cov[[g]], # Sigma_g
                                          labels = lab.group[[g]],
                                          phiFirst = phiFirst[g])
  } # end group

  # assemble over groups
  Hessian <- Hess.group[[1]]
  if(ngroups > 1) {
    Hessian <- matrix_bdiag(Hess.group) # block diagonal matrix from a list of matrices
  }
  return(Hessian)
}

# hessian.analytical
hessian.analytical <- function(GLIST = NULL, N = NULL, data.cov = NULL,
                               Sigma = NULL, labels = NULL,
                               phiFirst = NULL){

  # # state or final?
  # if(is.null(GLIST)) GLIST <- model@GLIST

  Lambda <- GLIST$lambda
  Phi    <- GLIST$phi
  Psi    <- GLIST$psi
  Sigma.inv <- compute.inv(Sigma, pdef = is.Pdef(Sigma))

  lab.lambda <- labels$label.load
  lab.phi    <- labels$label.phi
  lab.psi    <- labels$label.psi
  pstar <- nrow(lab.lambda)
  rstar <- nrow(lab.phi)
  tstar <- nrow(lab.psi)

  M <- Sigma.inv %*% data.cov %*% Sigma.inv
  Q <-  Sigma.inv - M
  I <- diag(ncol(Lambda))
  tLambda <- t(Lambda)

  csi   <- Sigma.inv %*% Lambda
  eta   <- csi %*% Phi
  alpha <- tLambda %*% csi
  beta  <- Phi %*% alpha
  mu    <- beta %*% Phi

  gamma <- M %*% Lambda
  zeta  <- tLambda %*% gamma
  iota  <- Phi %*% zeta
  nu    <- gamma %*% Phi

  omicron <- Q %*% Lambda
  chi     <- tLambda %*% omicron
  omega   <- omicron %*% Phi

  tbeta <- t(beta)
  tiota <- t(iota)

  if(dim(lab.phi)[1]==0){
    delta_Phi <-  delta_Lambda_Phi <- delta_Phi_Psi <- NULL
    t_delta_Lambda_Phi <- t_delta_Phi_Psi <- NULL
  }else{
    delta_Phi <- dx_Phi(I = I, zeta = zeta, alpha = alpha, chi = chi, rstar = rstar, lab.phi = lab.phi)
    delta_Phi <- pmax(delta_Phi, t(delta_Phi), na.rm=TRUE)
  }

  if(dim(lab.psi)[1]==0){
    delta_Psi <-  delta_Lambda_Psi <- delta_Phi_Psi <- NULL
    t_delta_Lambda_Psi <- t_delta_Phi_Psi <- NULL
  }else{
    delta_Psi <- dx_Psi(Sigma.inv = Sigma.inv, M = M,tstar = tstar,lab.psi = lab.psi)
    delta_Psi <- pmax(delta_Psi, t(delta_Psi), na.rm=TRUE)
  }

  if(dim(lab.lambda)[1]==0){
    delta_Lambda <- delta_Lambda_Phi <- delta_Lambda_Psi <- NULL
    t_delta_Lambda_Phi <- t_delta_Lambda_Psi <- NULL
  }else{
    delta_Lambda <- dx_Lambda(Phi = Phi, Sigma.inv  = Sigma.inv, Q = Q, mu = mu, iota = iota, omega = omega,
                              eta = eta, nu = nu, pstar = pstar, lab.lambda = lab.lambda)
    delta_Lambda <- pmax(delta_Lambda, t(delta_Lambda), na.rm=TRUE)
  }

  if(dim(lab.lambda)[1]!=0 & dim(lab.phi)[1]!=0){
    delta_Lambda_Phi <- dx_Lambda_Phi(I = I, omicron = omicron, tbeta = tbeta, csi = csi, tiota = tiota,
                                      rstar = rstar, pstar = pstar, lab.lambda = lab.lambda, lab.phi = lab.phi)
    t_delta_Lambda_Phi <- t(delta_Lambda_Phi)
  }

  if(dim(lab.lambda)[1]!=0 & dim(lab.psi)[1]!=0){
    delta_Lambda_Psi <- dx_Lambda_Psi(Sigma.inv  = Sigma.inv, M = M, eta = eta, omega = omega, pstar = pstar,
                                      tstar = tstar,  lab.lambda = lab.lambda, lab.psi = lab.psi)

    t_delta_Lambda_Psi <- t(delta_Lambda_Psi)
  }

  if(dim(lab.phi)[1]!=0 & dim(lab.psi)[1]!=0){
    delta_Phi_Psi <- dx_Phi_Psi(I = I, csi = csi,gamma = gamma, omicron = omicron, tstar = tstar,
                                rstar = rstar, lab.psi = lab.psi, lab.phi = lab.phi)
    t_delta_Phi_Psi <- t(delta_Phi_Psi)
  }

  if(phiFirst){
    A <- rbind(delta_Lambda, delta_Lambda_Phi, delta_Lambda_Psi)
    B <- rbind(t_delta_Lambda_Phi, delta_Phi, delta_Phi_Psi)
    C <- rbind(t_delta_Lambda_Psi, t_delta_Phi_Psi, delta_Psi)
  }else{
    A <- rbind(delta_Lambda, delta_Lambda_Psi, delta_Lambda_Phi)
    B <- rbind(t_delta_Lambda_Psi, delta_Psi, t_delta_Phi_Psi)
    C <- rbind(t_delta_Lambda_Phi, delta_Phi_Psi, delta_Phi)
  }

  hessian <- N * cbind(A, B, C)

  return(hessian)
}

# Lambda
dx_Lambda <- function(Phi, Sigma.inv, Q, mu, iota, omega, eta, nu, pstar, lab.lambda){

  e <- matrix(NA, pstar, pstar)

  for(row in 1:pstar){
    for(column in 1:row){
      j <- lab.lambda[row, 1]
      i <- lab.lambda[row, 2]
      s <- lab.lambda[column, 1]
      t <- lab.lambda[column, 2]

      one   <- Q[t,i] * (Phi - mu)[s,j]
      two   <- Sigma.inv[t,i] * (iota %*% Phi)[s,j]
      three <- omega[i,s] * eta[t,j]
      four  <- eta[i,s] * nu[t,j]

      e[row, column] <- one + two - three + four
    }
  }
  rownames(e) <- colnames(e) <- paste0("F", lab.lambda[,1], "=~x", lab.lambda[,2])
  return(e)
}

# Phi
dx_Phi <- function(I, zeta, alpha, chi, rstar, lab.phi){

  e <- matrix(NA, rstar, rstar)

  for(row in 1:rstar){
    for(column in 1:row){
      g <- lab.phi[row, 1]
      h <- lab.phi[row, 2]
      l <- lab.phi[column, 1]
      q <- lab.phi[column, 2]

      one   <- 2 - I[l,q] - I[g,h] + I[l,q] * I[g,h]
      two   <- zeta[h,l] * alpha[q,g] - alpha[h,l] * chi[q,g]
      three <- 2 - I[l,q] - I[g,h]
      four  <- zeta[g,l] * alpha[q,h] - alpha[g,l] * chi[q,h]

      e[row, column] <- 0.5 * (one * two + three * four)

    }
  }
  rownames(e) <- colnames(e) <- paste0("F",lab.phi[,1], "~~F", lab.phi[,2])
  return(e)
}

# Psi
dx_Psi <- function(Sigma.inv, M, tstar, lab.psi){

  e <- matrix(NA, tstar, tstar)

  for(column in 1:tstar){
    for(row in column:tstar){
      i <- lab.psi[row, 1]
      t <- lab.psi[column, 1]

      e[row, column] <- 0.5 * Sigma.inv[i,t] %*% (2 * M - Sigma.inv)[i,t]
    }
  }
  rownames(e) <- colnames(e) <- paste0("x", lab.psi[,1], "~~x", lab.psi[,2])
  return(e)
}

# Lambda and Phi
dx_Lambda_Phi <- function(I, omicron, tbeta, csi, tiota, rstar, pstar, lab.lambda, lab.phi){

  e <- matrix(NA, rstar, pstar)

  for(row in 1:rstar){
    for(column in 1:pstar){
      j <- lab.lambda[column, 1]
      i <- lab.lambda[column, 2]
      h <- lab.phi[row, 1]
      g <- lab.phi[row, 2]

      one   <- 2 - I[g,h]
      two   <- omicron[i,g] * (I - tbeta)[h,j]
      three <- omicron[i,h] * (I - tbeta)[g,j]
      four  <- csi[i,g] * tiota[h,j]
      five  <- csi[i,h] * tiota[g,j]

      e[row, column] <- 0.5 * one * (two + three + four + five)
    }
  }
  rownames(e) <- paste0("F", lab.phi[,1], "~~F", lab.phi[,2])
  colnames(e) <- paste0("F", lab.lambda[,1], "=~x", lab.lambda[,2])
  return(e)
}

# Lambda and Psi
dx_Lambda_Psi <- function(Sigma.inv, M, eta, omega, pstar, tstar, lab.lambda, lab.psi){

  e <- matrix(NA, tstar, pstar)

  for(row in 1:tstar){
    for(column in 1:pstar){
      j <- lab.lambda[column, 1]
      i <- lab.lambda[column, 2]
      t <- lab.psi[row, 1]

      one <- M[i,t] * eta[t, j]
      two <- Sigma.inv[i,t] * omega[t, j]

      e[row, column] <- one - two
    }
  }
  rownames(e) <- paste0("x", lab.psi[,1], "~~x", lab.psi[,2])
  colnames(e) <- paste0("F", lab.lambda[,1], "=~x", lab.lambda[,2])
  return(e)
}

# Phi and Psi
dx_Phi_Psi <- function(I, csi, gamma, omicron, tstar, rstar, lab.psi, lab.phi){

  e <- matrix(NA, tstar, rstar)
  for(row in 1:tstar){
    for(column in 1:rstar){
      t <- lab.psi[row, 1]
      h <- lab.phi[column, 2]
      g <- lab.phi[column, 1]

      one   <- 2 - I[g, h]
      two   <- csi[t, h] * gamma[t, g]
      three <- csi[t, g] * omicron[t, h]

      e[row, column] <- 0.5 * one * (two - three)
    }
  }
  rownames(e) <- paste0("x", lab.psi[,1], "~~x", lab.psi[,2])
  colnames(e) <- paste0("F", lab.phi[,1], "~~F", lab.phi[,2])
  return(e)
}

# Lambda
set_label_lambda_noconstr <- function(r, p){
  return(cbind(rep(1:r, each = p), rep(1:p, r)))
}

set_label_lambda_constr <- function(r, p, lambda_constr){
  constr_row_idx<-numeric()

  for(row in 1:nrow(lambda_constr))
    constr_row_idx[row]<-get_lambda_constr_row_idx(r, p, lambda_constr[row,])

  label <- set_label_lambda_noconstr(r, p)
  label <- label[-c(constr_row_idx), , drop = FALSE]

  return(label)
}

get_lambda_constr_row_idx <- function(r, p, constraint){
  block <- p * (constraint[2] - 1)

  return( block + constraint[1] )
}

# Phi
set_label_phivcov_constr <- function(r, varconstr, covconstr){
  constr_row_idx<-numeric()

  for(row in 1:nrow(covconstr))
    constr_row_idx[row]<-get_phi_constr_row_idx(r, covconstr[row,])

  label <- set_label_phi_noconstr(r)
  label <- label[-c(varconstr, constr_row_idx), , drop = FALSE]

  return(label)
}

set_label_phiv_constr <- function(r, varconstr){
  label <- set_label_phi_noconstr(r)
  label <- label[-c(varconstr), , drop = FALSE]
  return(label)
}

set_label_phicov_constr <- function(r, covconstr){
  constr_row_idx<-numeric()

  for(row in 1:nrow(covconstr))
    constr_row_idx[row]<-get_phi_constr_row_idx(r, covconstr[row,])

  label <- set_label_phi_noconstr(r)
  label <- label[-c(constr_row_idx), , drop = FALSE]

  return(label)
}

set_label_phi_noconstr <- function(r){
  if(r > 0){
    label <- cbind(1:r, 1:r)
    if(r > 1){ # if r > 1, otherwise we do not need indices for covariances
      for(idx in 1:(r-1)){
        label <- rbind(label, cbind(idx, (idx+1):r))
      }
    }
  }else{
    warning("penfa WARNING: The number of factors was set as zero")
    label <- NULL
  }
  return(label)
}


set_label_phi_constr <- function(r, varconstr, covconstr){
  constr_row_idx<-numeric()

  for(row in 1:nrow(covconstr))
    constr_row_idx[row]<-get_phi_constr_row_idx(r, covconstr[row,])

  label <- set_label_phi_noconstr(r)
  label <- label[-c(varconstr, constr_row_idx), , drop = FALSE]

  return(label)
}

get_phi_constr_row_idx <- function(r,constraint){
  const_baseline <- r
  i=1
  while (i <= (constraint[1]-1)) {
    const_baseline <- const_baseline + r - i
    i<- i+1
  }
  return(const_baseline + constraint[2] - constraint[1])
}

order_idx_matrix <- function(idx_matrix){
  for(i in 1:nrow(idx_matrix)){
    indicator <- idx_matrix[i,1] < idx_matrix[i,2]
    if(indicator == TRUE){
      # do nothing
    } else{
      idx_matrix[i,] <- c(idx_matrix[i, 2], idx_matrix[i, 1])
    }
  }
  return(idx_matrix)
}

# Lambda
get_label_lambda <- function(restrict_lambda, restrict_pos, dim.row, dim.col){

  if(restrict_lambda){
    nresL <- nrow(restrict_pos)
    pstar <- dim.row * dim.col - nresL
    label <- set_label_lambda_constr(r = dim.col,
                                     p = dim.row,
                                     lambda_constr = restrict_pos)
  }else{
    pstar <- dim.row * dim.col
    label <- set_label_lambda_noconstr(r = dim.col, p = dim.row)
  }

  label.reversed <- cbind(label[,2], label[,1])

  res <- list("label.load" = label,
              "label.load.reversed" = label.reversed,
              "pstar" = pstar)
  return(res)

}

# Phi
get_label_phi <- function(restrict_phivar, restrict_phicov, restrict_var_item, restrict_cov_item, dim){

  if(restrict_phivar){
    nresV <- length(restrict_var_item)

    if(restrict_phicov){
      nresC <- nrow(restrict_cov_item)
      nres <- sum(nresV, nresC)
      rstar <- dim*(dim+1)/2 - nres
      label <- set_label_phivcov_constr(r = dim,
                                        varconstr = restrict_var_item,
                                        covconstr = restrict_cov_item)
    }else{
      nres <- nresV
      rstar <- dim*(dim+1)/2 - nres
      label <- set_label_phiv_constr(r = dim, varconstr = restrict_var_item)
    }
  }else{
    if(restrict_phicov){
      nresC <- nrow(restrict_cov_item)
      nres <- nresC
      rstar <- dim*(dim+1)/2 - nres
      label <- set_label_phicov_constr(r = dim, covconstr = restrict_cov_item)
    }else{
      rstar <- dim*(dim+1)/2
      label <- set_label_phi_noconstr(r = dim)
    }
  }

  res <- list("label_phi" = label, "rstar" = rstar)
  return(res)
}

# Psi
get_label_psi <- function(restrict_psi, restrict_pos, dim){

  if(restrict_psi){
    label <- cbind(1:dim, 1:dim)[- restrict_pos, ]
    tstar <- dim - length(restrict_pos)
  }else{
    label <- cbind(1:dim, 1:dim)
    tstar <- dim
  }

  res <- list("label_psi" = label, "tstar" = tstar)
  return(res)
}


# model_implied
model_implied <- function(model = NULL, GLIST = NULL) {

  stopifnot(inherits(model, "penfaModel"))

  # state or final?
  if(is.null(GLIST)) GLIST <- model@GLIST

  # model-implied variance/covariance matrix ('sigma hat')
  Sigma.hat <- computeSigmaHat(model = model, GLIST = GLIST)

  # model-implied mean structure ('mu hat')
  if(model@meanstructure) {
    Mu.hat <- computeMuHat(model = model,  GLIST = GLIST)
  } else {
    Mu.hat <- vector("list", length = model@ngroups)
  }

  implied <- list(cov = Sigma.hat, mean = Mu.hat)

  implied
}

# model_loglik
model_loglik <- function(moddata = NULL,
                         samplestats = NULL,
                         implied = NULL,
                         model   = NULL,
                         options = NULL) {

  ngroups <- moddata@ngroups
  logl.group <- rep(as.numeric(NA), ngroups)
  logl.ok <- TRUE
  # samplestats filled in
  if(length(samplestats@ntotal) == 0L) {
    logl.ok <- FALSE
  }

  if(logl.ok) {
    for(g in seq_len(ngroups) ) {
      # single-level, complete data
      if(options$meanstructure) {
        Mu <- implied$mean[[g]]
      } else {
        Mu <- samplestats@mean[[g]]
      }
      logl.group[g] <- mvnorm_loglik_samplestats(sample.mean = samplestats@mean[[g]],
                                                 sample.cov  = samplestats@cov[[g]],
                                                 sample.nobs = samplestats@nobs[[g]],
                                                 Mu          = Mu,
                                                 Sigma       = implied$cov[[g]],
                                                 Sinv.method = "eigen",
                                                 Sigma.inv   = NULL)
    } # g
  } # logl.ok is TRUE

  logl <- sum(logl.group)
  npar <- model@nx.free

  out <- list(loglik = logl,
              loglik.group = logl.group,
              npar = npar,
              ntotal = samplestats@ntotal)
  out
}

# mvnorm_loglik_samplestats
mvnorm_loglik_samplestats <- function(sample.mean = NULL, sample.cov = NULL,
                                      sample.nobs  = NULL, Mu = NULL, Sigma = NULL,
                                      Sinv.method = "eigen", Sigma.inv = NULL){
  P <- length(sample.mean)
  N <- sample.nobs

  Mu <- as.numeric(Mu)
  sample.mean <- as.numeric(sample.mean)

  LOG.2PI <- log(2 * pi)
  if(is.null(Sigma.inv)) {
    Sigma.inv <- matrix_symmetric_inverse(S = Sigma, logdet = TRUE,
                                          Sinv.method = Sinv.method)
    logdet <- attr(Sigma.inv, "logdet")
  } else {
    logdet <- attr(Sigma.inv, "logdet")
    if(is.null(logdet)) {
      # compute - ln|Sigma.inv|
      ev <- eigen(Sigma.inv, symmetric = TRUE, only.values = TRUE)
      logdet <- -1 * sum(log(ev$values))
    }
  }

  # tr(Sigma^{-1} %*% S)
  DIST1 <- sum(Sigma.inv * sample.cov)
  # (ybar - mu)^T %*% Sigma.inv %*% (ybar - mu)
  Diff <- as.numeric(sample.mean - Mu)
  DIST2 <- sum(as.numeric(crossprod(Diff, Sigma.inv)) * Diff)

  loglik <- -N/2 * (P * LOG.2PI + logdet + DIST1 + DIST2)
  loglik
}
