########################
# ---- Post-checks ----
########################
# Content:
# - convcheck : checks convergence of the fitted model
# - posdef : inspection and possibly correction for positive definiteness of parameter covariance matrix
# - is_Admissible : checks the admissibility of the resulting solution
#
# - model_vcov_se : computes standard errors for model parameters
# - compute_CI : computes confidence intervals for the model parameters. They are based on the Bayesian variance covariance matrix
#                (i.e., inverse of penalized Hessian/expected Fisher information matrix) and they are of the
#                form param +/- z_alpha/2 * sqrt(diag(vcov)).
# - compute_edf : computes effective degrees of freedom (edf) as trace of F = (H + S)^{-1} H
# - compute_IC : computes information criteria based on edf
#
# Inspect*
# - object_inspect_phi, computePHI : inspect and extract group-specific PHI matrices from a fitted model
# - object_inspect_psi, computePSI : inspect and extract group-specific PSI matrices from a fitted model
# - computeLAMBDA : extracts group-specific LAMBDA matrices from a fitted model
# - computeTAU : extracts group-specific TAU vectors (intercepts of ov) from a fitted model
# - computeKAPPA : extracts group-specific KAPPA vectors (factor means) from a fitted model
# - object_inspect_gradient : inspects the gradient of a fitted model
# - inspect_coef, inspect_est : inspect model parameters
# - object_inspect_implied : displays the model-implied moments (covariance matrix and mean) of a 'penfa' object
#
# Some of these functions are adaptations from the routines present in the lavaan package (https://CRAN.R-project.org/package=lavaan)
#
# -----------------------------------------------------------------------------------------------------------------------

# convcheck
convcheck <- function(model = NULL, options = NULL, modoptim = NULL){

  verbose     <- options$verbose
  strategy    <- options$strategy
  information <- options$information
  threshold   <- options$optim.dx.tol

  if(information == "hessian"){
    type <- "Hessian matrix"
  }else if(information == "fisher"){
    type <- "Fisher information matrix"
  }

  max.grad <- max(abs(modoptim$dx.pen))
  # Only print model information when strategy = fixed
  # Not in presence of the automatic procedure
  if(verbose & strategy == "fixed")
    cat("\nLargest absolute gradient value:", format(round(max.grad,8), nsmall = 8))

  if(max.grad > threshold){
    condition1 <- FALSE
      cat("\n")
      warning("\npenfa WARNING: Not all elements of the gradient are zero.\n",
              "The optimizer may not have found a local solution\n\n")
  }else{
    condition1 <- TRUE
  }

  e.v <- eigen(modoptim$hessian.pen, symmetric = TRUE, only.values = TRUE)$values

  if(min(e.v) >= 0){
    condition2 <- TRUE
    if(verbose & strategy == "fixed")
      cat("\n", type, " is positive definite\n", sep="")
  }else{
    condition2 <- FALSE
      warning(paste0("\npenfa WARNING: ", type, " is not positive definite\n"))
  }

  if(verbose & strategy == "fixed"){
    cat("Eigenvalue range: [",min(e.v),", ",max(e.v),"]\n", sep = "")
    cat("Trust region iterations:", modoptim$iterations,"\n")
  }

  conditions <- c("first.order" = condition1, "second.order" = condition2)
  return(conditions)
}

# posdef
posdef <- function(omega, opt = 1){

  eS <- eigen(omega, symmetric = TRUE)
  e.val <- eS$values
  e.vec <- eS$vectors

  if(opt == 1){ # Pdef
    check.eigen <- any(e.val <= 0)

    if (check.eigen == TRUE){
      # Put negative or null eigenvalues to 0 and rescale
      n.e.val <- e.val[e.val <= 0]
      s <- sum(e.val[n.e.val]) * 2
      t <- s^2 * 100 + 1
      p <- min(e.val[(e.val <= 0) == FALSE])
      e.val[e.val <= 0] <- p * (s - n.e.val)^2/t
    }

  }else{  # Option 2: postVb

    epsilon     <- 0.0000001
    check.eigen <- FALSE

    # Take absolute value of negative eigen-values
    if(min(e.val) < sqrt(.Machine$double.eps) && sign( min( sign(e.val) ) ) == -1){
      check.eigen   <- TRUE
      e.val <- abs(e.val)
    }

    # Replace those close to zero with epsilon
    if(min(e.val) < sqrt(.Machine$double.eps) ) {
      check.eigen <- TRUE
      pep <- which(e.val < sqrt(.Machine$double.eps))
      e.val[pep] <- epsilon
    }
  }

  if(check.eigen){
    D <- diag(e.val, nrow = length(e.val), ncol =  length(e.val))
    D.inv <- diag(1/e.val, nrow = length(e.val), ncol =  length(e.val))
    res <- e.vec %*% D %*% t(e.vec)
    res.inv <- e.vec %*% D.inv %*% t(e.vec)
  }else{
    res <- omega
    res.inv <- e.vec %*% diag(1/e.val, nrow = length(e.val), ncol =  length(e.val)) %*% t(e.vec)
  }

  output <- list(check.eigen = check.eigen, res = res, res.inv = res.inv)

  return(output)
}

# is_Admissible
is_Admissible <- function(model = NULL, verbose = FALSE){

  stopifnot(inherits(model, "penfaModel"))
  ngroups     <- model@ngroups
  Lambda.all  <- computeLAMBDA(model)
  Phi.all     <- computePHI(model)
  Psi.all     <- computePSI(model)
  admis       <- logical(ngroups)

  for(g in seq_len(ngroups)){
    Lambda <- Lambda.all[[g]]
    Phi    <- Phi.all[[g]]
    Psi    <- Psi.all[[g]]
    r      <- ncol(Lambda)

    result.ok <- TRUE

    if(ngroups > 1){
      mex <- paste("in Group", g)
    }else{
      mex <- character(0)
    }

    # 1. No negative variances in Psi
    if(min(diag(Psi)) < 0.0 ){
      result.ok <- FALSE
      warning(paste0("penfa WARNING: Heywood case: negative unique variances produced ", mex, "\n"))
    }

    # 2. No negative variances in Phi
    if(result.ok & min(diag(Phi)) < 0.0){
      result.ok <- FALSE
      warning(paste("penfa WARNING: Negative factor variances produced", mex, "\n"))
    }

    # 3. Phi positive definite
    if(result.ok){
      eigvals <- eigen(Phi, symmetric=TRUE, only.values=TRUE)$values
      if(any(eigvals < -1 * .Machine$double.eps^(3/4))) {
        result.ok <- FALSE
        warning(paste0("penfa WARNING: Factor covariance matrix not positive definite ", mex, " \n"))
      }
    }

    # 4. Psi positive definite
    if(result.ok){
      eigvals <- eigen(Psi, symmetric = TRUE, only.values = TRUE)$values
      if(any(eigvals < -1 * .Machine$double.eps^(3/4))) {
        result.ok <- FALSE
        warning(paste0("penfa WARNING: Covariance matrix of unique factors not positive definite ", mex, "\n"))
      }
    }

    # 5. Lambda full column rank
    if(result.ok){
      if(qr(Lambda)$rank < r){
        result.ok <- FALSE
        warning(paste0("penfa WARNING: Factor loading matrix not of full column rank ", mex, " \n"))
      }
    }

    # 6. No zero rows in Lambda
    if(result.ok){
      if(any(abs(rowSums(Lambda)) < 0.15 & ncol(Lambda) > 1)){
        # if Lambda has 1 factor, the zero loading can be a fixed loading or
        # a small cross-loading
        result.ok <- FALSE
        warning(paste0("penfa WARNING: Factor loading matrix with null rows ", mex, " \n"))
      }
    }

    admis[g] <- result.ok

  }# end group

  # Is the whole (multigroup) solution admissible?
  admissible <- all(admis) # all(admis == TRUE)

  if(verbose){
    admissible.txt <- ifelse(admissible, 'admissible', 'not admissible')
    if(admissible){
      cat("Factor solution:", admissible.txt,"\n")
    }else{
      cat("Factor solution:", admissible.txt,"\n")
    }
  }

  return(admissible)
}

# model_vcov_se
model_vcov_se <- function(model, partable, VCOV = NULL){

  # 0. special case: zeros for fixed elements
  if(is.null(VCOV)) {
    se <- rep(as.numeric(NA), model@nx.user)
    se[ partable$free == 0L ] <- 0.0
    return(se)
  }

  # 1. free parameters only
  x.var <- diag(VCOV)
  # check for negative values (what to do: NA or 0.0?)
  x.var[x.var < 0] <- as.numeric(NA)

  # Now we can take square root
  x.se <- sqrt( x.var )
  # put it in the entries of GLIST to associate each se with the corresponding parameter
  GLIST <- model_x2GLIST(model = model, x = x.se, type = "free")

  # se for full parameter table
  se <- model_get_parameters(model = model, GLIST = GLIST, type = "user", extra = FALSE)

  # 2. fixed parameters -> se = 0.0
  se[ which(partable$free == 0L) ] <- 0.0

  res <- list(x.se = x.se, se = se)

  return(res)
}

# compute_CI
compute_CI <- function(param, VCOV, std.err,
                       level = 0.95, modpenalty, options){
  npar         <- ncol(VCOV)

  if(is.null(std.err) || is.null(VCOV) ){
    return(NULL)
  }

  alpha   <- (1 - level)/2; alpha <- c(alpha, 1 - alpha)
  ci      <- matrix(NA, nrow = npar, ncol = 2)
  z_alpha <- qnorm(alpha)
  bounds  <- std.err %o% z_alpha
  ci      <- param + bounds

  colnames(ci) <- sprintf("%.*g%%", 3, 100 * alpha)
  rownames(ci) <- names(param)
  return(ci)
}

# compute_edf
compute_edf <- function(VCOV = NULL, modoptim = NULL,
                        modpenalty = NULL, options = options){

  verbose  <- options$verbose
  strategy <- options$strategy

  if(is.null(VCOV)){
    return(NULL)
  }

  H.pen     <- modoptim$hessian.pen
  S.h       <- modpenalty@Sh.info$S.h
  H.unpen   <- H.pen - S.h
  H.pen.inv <- VCOV
  F_mat     <- H.pen.inv %*% H.unpen

  edf.single <- diag(F_mat)
  edf        <- sum(edf.single)

  if(verbose & strategy == "fixed"){
    # cat("\nSingle edfs:", round(diag(F_mat), 3)," \n")
    cat("Effective degrees of freedom:", edf, " \n")
  }

  edf.list <- list("edf.single" = edf.single, "edf" = edf, "influence.mat" = F_mat)

  return(edf.list)
}

# compute_IC
compute_IC <- function(modoptim = NULL, samplestats = NULL, dgf = NULL){

  compute_AIC <- function(logl, df){
    -2*(logl) + 2*df
  }

  # compute_AIC3 <- function(logl, df){
  #   -2*(logl) + 3*df
  # }
  #
  # compute_CAIC<- function(logl, df, N){
  #   -2*(logl) + log(N+1)*df
  # }

  compute_BIC<- function(logl, df, N){
    -2*(logl) + log(N)*df
  }

  # compute_ABIC<- function(logl, df, N){
  #   -2*(logl) + log((N+2)/24)*df
  # }
  #
  # compute_HBIC<- function(logl, df, N){
  #   -2*(logl) + log(N/(2*pi))*df
  # }

  loglik  <- modoptim$logl.unpen
  nobs    <- samplestats@ntotal # overall

  aic  <- compute_AIC (logl = loglik, df = dgf)
  # aic3 <- compute_AIC3(logl = loglik, df = dgf)
  # caic <- compute_CAIC(logl = loglik, df = dgf, N = nobs)
  bic  <- compute_BIC (logl = loglik, df = dgf, N = nobs)
  # abic <- compute_ABIC(logl = loglik, df = dgf, N = nobs)
  # hbic <- compute_HBIC(logl = loglik, df = dgf, N = nobs)

  IC <- list(  "AIC"  = aic
             # , "AIC3" = aic3
             # , "CAIC" = caic
             , "BIC"  = bic
             # , "ABIC" = abic
             # , "HBIC" = hbic
             )
  return(IC)
}

# object_inspect_phi
object_inspect_phi <- function (object){
  OUT     <- computePHI(model = object@Model)
  if (length(object@Data@group.label) > 0L) {
    names(OUT) <- unlist(object@Data@group.label)
  }
  OUT
}

# computePHI
computePHI <- function (model = NULL, GLIST = NULL) {
  if (is.null(GLIST))
    GLIST <- model@GLIST

  ngroups <- model@ngroups
  nmat    <- model@nmat
  PHI <- vector("list", length = ngroups)
  for (g in 1:ngroups) {
    mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]
    MLIST <- GLIST[mm.in.group]
    PHI.g <- MLIST$phi
    PHI[[g]] <- PHI.g
  }
  PHI
}

# object_inspect_psi
object_inspect_psi <- function(object){
  OUT <- computePSI(model = object@Model)
  if (length(object@Data@group.label) > 0L) {
    names(OUT) <- unlist(object@Data@group.label)
  }
  OUT
}

# computePSI
computePSI <- function (model = NULL, GLIST = NULL){
  if(is.null(GLIST))
    GLIST <- model@GLIST
  ngroups <- model@ngroups
  nmat    <- model@nmat
  PSI <- vector("list", length = ngroups)
  for (g in 1:ngroups) {
    mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]
    MLIST <- GLIST[mm.in.group]
    PSI.g <- MLIST$psi
    PSI[[g]] <- PSI.g
  }
  PSI
}

# computeLAMBDA
computeLAMBDA <- function(model = NULL, GLIST = NULL){

  if(is.null(GLIST))
    GLIST <- model@GLIST
  ngroups        <- model@ngroups
  nmat           <- model@nmat

  # return a list
  LAMBDA <- vector("list", length=ngroups)

  # compute LAMBDA for each group
  for(g in 1:ngroups) {
    # which mm belong to group g?
    mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
    MLIST <- GLIST[ mm.in.group ]

    LAMBDA.g <- MLIST$lambda
    LAMBDA[[g]] <- LAMBDA.g
  }
  LAMBDA
}

# computeTAU
computeTAU <- function(model = NULL, GLIST = NULL){

  if(is.null(GLIST))
    GLIST <- model@GLIST
  ngroups        <- model@ngroups
  nmat           <- model@nmat

  # return a list
  TAU <- vector("list", length=ngroups)

  # compute TAU for each group
  for(g in 1:ngroups) {
    # which mm belong to group g?
    mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
    MLIST <- GLIST[ mm.in.group ]

    TAU.g <- MLIST$tau
    TAU[[g]] <- TAU.g
  }
  TAU
}

# computeKAPPA
computeKAPPA <- function(model = NULL, GLIST = NULL){

  if(is.null(GLIST))
    GLIST <- model@GLIST
  ngroups        <- model@ngroups
  nmat           <- model@nmat

  # return a list
  KAPPA <- vector("list", length=ngroups)

  # compute KAPPA for each group
  for(g in 1:ngroups) {
    # which mm belong to group g?
    mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
    MLIST <- GLIST[ mm.in.group ]

    KAPPA.g <- MLIST$kappa
    KAPPA[[g]] <- KAPPA.g
  }
  KAPPA
}

# object_inspect_gradient
object_inspect_gradient <- function(object){
  model       <- object@Model
  moddata     <- object@Data
  samplestats <- object@SampleStats
  dx <- model_gradient(model = model, GLIST = NULL, samplestats = samplestats,
                       moddata = object@Data, verbose = FALSE)
  dx
}

# inspect_coef
inspect_coef <- function(object, type = "free", add.labels = FALSE, add.class = FALSE) {

  if(type == "user" || type == "all") {
    type <- "user"
    idx <- 1:length( object@ParTable$lhs )
  } else if(type == "free") {
    idx <- which(object@ParTable$free > 0L & !duplicated(object@ParTable$free))
  } else {
    stop("penfa ERROR: argument `type' must be one of free or user")
  }
  EST <- inspect_est(object)
  cof <- EST[idx]

  # labels?
  if(add.labels) {
    names(cof) <- partable_labels(object@ParTable, type = type)
  }

  # class
  if(add.class) {
    class(cof) <- c("penfa.vector", "numeric")
  }
  cof
}

# inspect_est
inspect_est <- function(object){
  if(inherits(object, "penfa") & !is.null(object@ParTable$est)) {
    OUT <- object@ParTable$est
  }else{
    # try generic coef()
    OUT <- coef(object, type = "user")
    if(is.matrix(OUT)) {
      OUT <- rowMeans(OUT)
    }
  }
  OUT
}

# object_inspect_implied
object_inspect_implied <- function(object, add.labels = FALSE,
                                   add.class = FALSE,
                                   drop.list.single.group = FALSE) {

  ov.names <- object@pta$vnames$ov
  implied  <- object@Implied
  model    <- object@Model
  ngroups  <- model@ngroups
  b        <- 1 # block

  OUT <- vector("list", length = ngroups)
  for(g in seq_len(ngroups)) {
    # covariance matrix
    OUT[[g]]$cov  <- implied$cov[[g]]
    if(add.labels && !is.null(OUT[[g]]$cov)) {
      rownames(OUT[[g]]$cov) <- colnames(OUT[[g]]$cov) <- ov.names[[b]]
    }
    if(add.class) {
      class(OUT[[g]]$cov) <- c("penfa.matrix.symmetric", "matrix")
    }

    # mean vector
    if(model@meanstructure){
      OUT[[g]]$mean <- as.numeric(implied$mean[[g]])
      if(add.labels) {
        names(OUT[[g]]$mean) <- ov.names[[b]]
      }
      if(add.class) {
        class(OUT[[g]]$mean) <- c("penfa.vector", "numeric")
      }
    }
  } # groups

  if(ngroups == 1L && drop.list.single.group) {
    OUT <- OUT[[1]]
  }else if(length(object@Data@group.label) > 0L){
    names(OUT) <- unlist(object@Data@group.label)
  }
  OUT
}
