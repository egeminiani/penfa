#################
# ---- Utils ----
#################
# Content:
#
#
# COMPUTE*
# - computeSigmaHat : returns a list of group-specific model-implied covariance matrices
# - computeSigmaHat.LISREL : returns a model-implied covariance matrix (i.e., Lambda Phi Lambda^T + Psi)
# - is.Pdef : checks if the input matrix is positive-definite via eigen-decomposition
# - compute.inv ; returns the inverse of an input matrix (either or not pdef)
# - computeMuHat : returns a list of group-specific model-implied mean vectors
# - computeMuHat.LISREL : returns a model-implied mean vector (i.e., tau + Lambda kappa)
#
#
# MATRIX*
# - matrix_symmetric_inverse : returns the inverse of a non-singular (not necessarily positive-definite) symmetric matrix
# - matrix_vech_idx : returns the *vector* indices of the lower triangular elements of a symmetric matrix of size 'n'
# - matrix_vechru_idx: returns the *vector* indices of the upper triangular elements of a symmetric matrix of size 'n'
#                      (ROW-WISE)
# - matrix_bdiag : constructs block diagonal matrix from a list of matrices
# - matrix_duplication_pre_post : computes D^T A D (without explicitly computing D); A square matrix and sqrt(ncol) integer
#
#
# MODEL*
# - model_get_parameters : extracts parameters from model
# - model_x2GLIST : creates a standalone GLIST, filled with (new) x values
# - model_set_parameters : sets model parameters and makes a copy of the model
#
# Some of these functions are adaptations from the routines present in the lavaan package (https://CRAN.R-project.org/package=lavaan)
#
# ---------------------------------------------------------------------------------------------------------------

# computeSigmaHat
computeSigmaHat <- function(model = NULL, GLIST = NULL, debug = FALSE) {

  # state or final?
  if(is.null(GLIST)) GLIST <- model@GLIST
  nmat           <- model@nmat
  nvar           <- model@nvar
  ngroups        <- model@ngroups

  # return a list
  Sigma.hat <- vector("list", length = ngroups)

  for(g in 1:ngroups) {
    # which mm belong to group g?
    mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
    MLIST <- GLIST[mm.in.group]
    Sigma.hat[[g]] <- computeSigmaHat.LISREL(MLIST = MLIST)
    if(debug){
      cat("[DEBUG] : Model-implied Sigma \n")
      print(Sigma.hat[[g]])
    }
  } # ngroups

  Sigma.hat
}

# computeSigmaHat.LISREL
computeSigmaHat.LISREL <- function (MLIST = NULL){
  LAMBDA <- MLIST$lambda
  nvar <- nrow(LAMBDA)
  PHI  <- MLIST$phi
  PSI  <- MLIST$psi
  VYx <- tcrossprod(LAMBDA %*% PHI, LAMBDA) + PSI
  VYx
}

# is.Pdef
is.Pdef <- function(mat = NULL, ngroups = 1L){
  is.pdef <- numeric(length(ngroups))
  for(g in 1:ngroups) {
    ev <- eigen(mat[[g]], symmetric=TRUE, only.values=TRUE)$values
    is.pdef[g] <- any(ev < sqrt(.Machine$double.eps)) || sum(ev) == 0
  }
  return(!is.pdef)
}

# compute.inv
compute.inv <- function(mat, pdef){
  if(!pdef) {
    mat.inv <-  MASS::ginv(mat)
  } else {
    mat.inv <- inv.chol(mat, logdet = TRUE)
  }
  return(mat.inv)
}

# computeMuHat
computeMuHat <- function(model = NULL, GLIST = NULL) {

  # state or final?
  if(is.null(GLIST)) GLIST <- model@GLIST
  nmat           <- model@nmat
  ngroups        <- model@ngroups
  meanstructure  <- model@meanstructure

  # return a list
  Mu.hat <- vector("list", length=ngroups)

  for(g in 1:ngroups) {
    # which mm belong to group g?
    mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
    if(!meanstructure) {
      Mu.hat[[g]] <- numeric(model@nvar[g])
    }else{
      Mu.hat[[g]] <- computeMuHat.LISREL(MLIST = GLIST[ mm.in.group ])
    }
  } # ngroups
  Mu.hat
}

# computeMuHat.LISREL
computeMuHat.LISREL <- function (MLIST = NULL){
  TAU    <- MLIST$tau
  KAPPA  <- MLIST$kappa
  LAMBDA <- MLIST$lambda
  if(is.null(KAPPA) || is.null(TAU))
    return(matrix(0, nrow(LAMBDA), 1L))
  Mu.hat <- TAU + LAMBDA %*% KAPPA
  Mu.hat
}


# matrix_symmetric_inverse
matrix_symmetric_inverse <- function(S, logdet = FALSE, Sinv.method = "eigen",
                                     zero.warn= FALSE) {

  # catch zero cols/rows
  zero.idx <- which(colSums(S) == 0 & diag(S) == 0 & rowSums(S) == 0)
  S.orig <- S
  if(length(zero.idx) > 0L) {
    if(zero.warn) {
      warning("penfa WARNING: matrix to be inverted contains zero cols/rows")
    }
    S <- S[-zero.idx, -zero.idx]
  }

  P <- NCOL(S)

  if(P == 0L) {
    S.inv <- matrix(0,0,0)
    if(logdet) {
      attr(S.inv, "logdet") <- 0
    }
    return(S.inv)
  } else if(P == 1L) {
    tmp <- S[1,1]
    S.inv <- matrix(1/tmp, 1, 1)
    if(logdet) {
      attr(S.inv, "logdet") <- log(tmp)
    }
  } else if(P == 2L) {
    a11 <- S[1,1]; a12 <- S[1,2]; a21 <- S[2,1]; a22 <- S[2,2]
    tmp <- a11*a22 - a12*a21
    if(tmp == 0) {
    } else {
      S.inv <- matrix(c(a22/tmp, -a21/tmp, -a12/tmp, a11/tmp), 2, 2)
      if(logdet) {
        attr(S.inv, "logdet") <- log(tmp)
      }
    }
  } else if(Sinv.method == "eigen") {

    EV <- eigen(S, symmetric = TRUE)
    # V %*% diag(1/d) %*% V^{-1}, where V^{-1} = V^T
    S.inv <- tcrossprod(EV$vector / rep(EV$values, each = length(EV$values)), EV$vector)
    if(logdet) {
      if(all(EV$values >= 0)) {
        attr(S.inv, "logdet") <- sum(log(EV$values))
      } else {
        attr(S.inv, "logdet") <- as.numeric(NA)
      }
    }
  } else if(Sinv.method == "solve") {
    S.inv <- solve(S)
    if(logdet) {
      ev <- eigen(S, symmetric = TRUE, only.values = TRUE)
      if(all(ev$values >= 0)) {
        attr(S.inv, "logdet") <- sum(log(ev$values))
      } else {
        attr(S.inv, "logdet") <- as.numeric(NA)
      }
    }
  } else if(Sinv.method == "chol") {
    # this will break if S is not positive definite
    cS <- chol(S)
    S.inv <- chol2inv(cS)
    if(logdet) {
      diag.cS <- diag(cS)
      attr(S.inv, "logdet") <- sum(log(diag.cS * diag.cS))
    }
  } else {
    stop("penfa ERROR: Inversion method must be either `eigen', `solve' or `chol'")
  }

  if(length(zero.idx) > 0L) {
    logdet <- attr(S.inv, "logdet")
    tmp <- S.orig
    tmp[-zero.idx, -zero.idx] <- S.inv
    S.inv <- tmp
    attr(S.inv, "logdet") <- logdet
    attr(S.inv, "zero.idx") <- zero.idx
  }
  S.inv
}

# matrix_vech_idx
matrix_vech_idx <- function(n = 1L, diagonal = TRUE) {
  n <- as.integer(n)
  ROW <- matrix(seq_len(n), n, n)
  COL <- matrix(seq_len(n), n, n, byrow = TRUE)
  if(diagonal) which(ROW >= COL) else which(ROW > COL)
}

# matrix_vechru_idx
matrix_vechru_idx <- function(n = 1L, diagonal = TRUE) {
  n <- as.integer(n)
  ROW <- matrix(seq_len(n), n, n)
  COL <- matrix(seq_len(n),   n, n, byrow = TRUE)
  tmp <- matrix(seq_len(n*n), n, n, byrow = TRUE)
  if(diagonal) tmp[ROW >= COL] else tmp[ROW > COL]
}

# matrix_bdiag: input multiple arguments (coerced to list) or a list with multiple arguments
matrix_bdiag <- function(...) {

  if(nargs() == 0L) return(matrix(0,0,0))
  dots <- list(...)

  # create list of matrices
  if(is.list(dots[[1]])) {
    mlist <- dots[[1]]
  } else {
    mlist <- dots
  }
  if(length(mlist) == 1L) return(mlist[[1]])

  # more than 1 matrix
  nmat  <- length(mlist)
  nrows <- sapply(mlist, NROW); crows <- cumsum(nrows)
  ncols <- sapply(mlist, NCOL); ccols <- cumsum(ncols)
  trows <- sum(nrows)
  tcols <- sum(ncols)
  x <- numeric(trows * tcols)

  for(m in seq_len(nmat)) {
    if(m > 1L) {
      rcoffset <- trows*ccols[m-1] + crows[m-1]
    } else {
      rcoffset <- 0L
    }
    m.idx <- ( rep((0:(ncols[m] - 1L))*trows, each=nrows[m]) +
                 rep(1:nrows[m], ncols[m]) + rcoffset )
    x[m.idx] <- mlist[[m]]
  }

  attr(x, "dim") <- c(trows, tcols)
  x
}

# matrix_duplication_pre_post
matrix_duplication_pre_post <- function(A = matrix(0,0,0)) {
  # number of columns
  n2 <- NCOL(A)

  # square A only, n2 = n^2
  stopifnot(NROW(A) == n2, sqrt(n2) == round(sqrt(n2)))

  # dimension
  n <- sqrt(n2)

  # dup idx
  idx1 <- matrix_vech_idx(n); idx2 <- matrix_vechru_idx(n)

  OUT <- A[idx1, , drop = FALSE] + A[idx2, , drop = FALSE]
  u <- which(idx1 %in% idx2); OUT[u,] <- OUT[u,] / 2.0
  OUT <- OUT[, idx1, drop = FALSE] + OUT[, idx2, drop = FALSE]
  OUT[,u] <- OUT[,u] / 2.0

  OUT
}

# model_get_parameters
model_get_parameters <- function(model = NULL, GLIST = NULL, type = "free",extra = TRUE) {
  # type == "free": only non-redundant free parameters (x)
  # type == "user": all parameters listed in User model

  # state or final?
  if(is.null(GLIST)) GLIST <- model@GLIST

  if(type == "free") {
    N <- model@nx.free
  } else if(type == "user") {
    N <- model@nx.user
  }

  x <- numeric(N)

  for(mm in 1:length(model@GLIST)) {
    if(type == "free") {
      m.idx <- model@m.free.idx[[mm]]
      x.idx <- model@x.free.idx[[mm]]
    } else if(type == "user") {
      m.idx <- model@m.user.idx[[mm]]
      x.idx <- model@x.user.idx[[mm]]
    }
    x[x.idx] <- GLIST[[mm]][m.idx]
  }
  x
}

# model_x2GLIST
model_x2GLIST <- function(model = NULL, x = NULL, type = "free", m.el.idx = NULL, x.el.idx = NULL) {
  GLIST <- model@GLIST
  for(mm in 1:length(GLIST)) {
    # skip empty matrix
    if(nrow(GLIST[[mm]]) == 0L)
      next

    M.EL.IDX <- model@m.free.idx[[mm]]
    X.EL.IDX <- model@x.free.idx[[mm]]

    # assign
    GLIST[[mm]][M.EL.IDX] <- x[X.EL.IDX]
  }
  GLIST
}

# model_set_parameters
model_set_parameters <- function(model = NULL, x = NULL) {
  tmp <- model@GLIST
  for(mm in 1:length(model@GLIST)) {
    m.free.idx <- model@m.free.idx[[mm]]
    x.free.idx <- model@x.free.idx[[mm]]
    tmp[[mm]][m.free.idx] <- x[x.free.idx]
  }

  model@GLIST <- tmp
  model
}
