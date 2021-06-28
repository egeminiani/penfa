################################
# ---- Model representation ----
################################
# Content:
# - get_model: constructor of the model representation
# - representation.LISREL : specifies the required model matrices according to the LISREL representation
#
# Some of these functions are adaptations from the routines present in the lavaan package (https://CRAN.R-project.org/package=lavaan)
#
# --------------------------------------------------------------------------------------------------------------

# get_model
get_model <- function(partable = NULL, options = NULL) {

  # global info from user model
  ngroups <- partable_ngroups(partable)
  meanstructure <- any(partable$op == "~1")

  # select model matrices
  REP <- representation.LISREL(partable, target = NULL, extra = TRUE)
  if(options$debug){
    cat(" [DEBUG] : LISREL representation\n")
    print(REP)
  }

  # check for non-existing parameters
  bad.idx <- which(REP$mat == "" )
  if(length(bad.idx) > 0L) {
    label <- paste(partable$lhs[bad.idx[1]], partable$op[bad.idx[1]],
                   partable$rhs[bad.idx[1]], sep = " ")
    stop("penfa ERROR: parameter is not defined: ", label)
  }

  # prepare nG-sized slots
  nG          <- sum(unlist(attr(REP, "mmNumber")))
  GLIST       <- vector(mode="list", nG)
  names(GLIST)<- unlist(attr(REP, "mmNames"))
  dimNames    <- vector(mode="list", length=nG)
  isSymmetric <- logical(nG)
  mmSize      <- integer(nG)
  m.free.idx  <- m.user.idx <- vector(mode="list", length=nG)
  x.free.idx  <- x.user.idx <- vector(mode="list", length=nG)
  nvar        <- integer(ngroups)
  nmat        <- unlist(attr(REP, "mmNumber"))
  num.idx     <- vector("list", length = ngroups)

  offset <- 0L
  for(g in 1:ngroups) {

    # observed and latent variables for this group
    ov.names     <- partable_vnames(partable, "ov", group = g)
    ov.num       <- partable_vnames(partable, "ov", group = g)
    nvar[g]      <- length(ov.names)
    num.idx[[g]] <- which(ov.names %in% ov.num)

    # model matrices for this group
    mmNumber    <- attr(REP, "mmNumber")[[g]]
    mmNames     <- attr(REP, "mmNames")[[g]]
    mmSymmetric <- attr(REP, "mmSymmetric")[[g]]
    mmDimNames  <- attr(REP, "mmDimNames")[[g]]
    mmRows      <- attr(REP, "mmRows")[[g]]
    mmCols      <- attr(REP, "mmCols")[[g]]

    for(mm in 1:mmNumber) {
      # offset in GLIST
      offset <- offset + 1L

      # matrix size, symmetric, dimNames
      if(mmSymmetric[mm]) {
        N <- mmRows[mm]
        mm.size <- as.integer(N*(N+1)/2)
      } else {
        mm.size <- as.integer(mmRows[mm] * mmCols[mm])
      }
      mmSize[offset]      <- mm.size
      isSymmetric[offset] <- mmSymmetric[mm]
      dimNames[[offset]]  <- mmDimNames[[mm]]

      # select elements for this matrix
      idx <- which(partable$group == g & REP$mat == mmNames[mm])

      # create empty `pattern' matrix
      tmp <- matrix(0L, nrow=mmRows[mm], ncol=mmCols[mm])

      # 1. first assign free values only, to get vector index -> used in model_objective
      tmp[ cbind(REP$row[idx], REP$col[idx]) ] <- partable$free[idx]

      if(mmSymmetric[mm]) {
        # NOTE: assume everything is in the UPPER tri
        T <- t(tmp); tmp[lower.tri(tmp)] <- T[lower.tri(T)]
      }
      m.free.idx[[offset]] <-     which(tmp > 0)
      x.free.idx[[offset]] <- tmp[which(tmp > 0)]

      # 3. general mapping between user and GLIST
      tmp[ cbind(REP$row[idx], REP$col[idx]) ] <- partable$id[idx]
      if(mmSymmetric[mm]) {
        T <- t(tmp); tmp[lower.tri(tmp)] <- T[lower.tri(T)]
      }
      m.user.idx[[offset]] <-     which(tmp > 0)
      x.user.idx[[offset]] <- tmp[which(tmp > 0)]

      # 4. now assign starting/fixed values
      # create empty matrix
      tmp <- matrix(0.0, nrow=mmRows[mm], ncol=mmCols[mm])
      tmp[ cbind(REP$row[idx], REP$col[idx]) ] <- partable$start[idx]
      if(mmSymmetric[mm]) {
        T <- t(tmp); tmp[lower.tri(tmp)] <- T[lower.tri(T)]
      }

      # assign matrix to GLIST
      GLIST[[offset]] <- tmp
    } # mm

  } # g

  # which free parameters are observed variances?
  ov.names <- partable_vnames(partable, "ov")
  x.free.var.idx <- partable$free[partable$free & partable$lhs %in% ov.names &
                                  partable$op == "~~" & partable$lhs == partable$rhs]

  Model <- new("penfaModel",
                GLIST          = GLIST,
                dimNames       = dimNames,
                isSymmetric    = isSymmetric,
                mmSize         = mmSize,
                meanstructure  = meanstructure,
                ngroups        = ngroups,
                nmat           = nmat,
                nvar           = nvar,
                num.idx        = num.idx,
                nx.free        = max(partable$free),
                nx.user        = max(partable$id),
                m.free.idx     = m.free.idx,
                x.free.idx     = x.free.idx,
                x.free.var.idx = x.free.var.idx,
                m.user.idx     = m.user.idx,
                x.user.idx     = x.user.idx)

  if(options$debug){
    cat(" [DEBUG] : Model\n")
    print( str(Model) )
    cat(" [DEBUG] : Parameter list \n")
    print( Model@GLIST)
  }
  Model
}

# representation.LISREL
representation.LISREL <- function(partable = NULL, target = NULL,
                                  extra = FALSE, remove.nonexisting = TRUE){

  # prepare target list
  if(is.null(target)) target <- partable
  stopifnot(!is.null(target$group))

  # prepare output
  N <- length(target$lhs)
  tmp.mat <- character(N); tmp.row <- integer(N); tmp.col <- integer(N)
  # global settings
  meanstructure <- any(partable$op == "~1")

  # number of groups
  ngroups <- partable_ngroups(partable)

  if(extra) {
    REP.mmNames     <- vector("list", ngroups)
    REP.mmNumber    <- vector("list", ngroups)
    REP.mmRows      <- vector("list", ngroups)
    REP.mmCols      <- vector("list", ngroups)
    REP.mmDimNames  <- vector("list", ngroups)
    REP.mmSymmetric <- vector("list", ngroups)
  }

  for(g in 1:ngroups) {

    # info from user model per group
    ov.names   <- partable_vnames(partable, "ov",  group=g); nvar <- length(ov.names)
    lv.names   <- partable_vnames(partable, "lv",  group=g); nfac <- length(lv.names)

    # 1a. "=~" regular indicators
    idx <- which(target$group == g & target$op == "=~" & !(target$rhs %in% lv.names))
    tmp.mat[idx] <- "lambda"
    tmp.row[idx] <- match(target$rhs[idx], ov.names)
    tmp.col[idx] <- match(target$lhs[idx], lv.names)

    # 3a. "~~" ov
    idx <- which(target$group == g & target$op == "~~" & !(target$lhs %in% lv.names))
    tmp.mat[idx] <- "psi"
    tmp.row[idx] <- match(target$lhs[idx], ov.names)
    tmp.col[idx] <- match(target$rhs[idx], ov.names)

    # 3b. "~~" lv
    idx <- which(target$group == g & target$op == "~~" & target$rhs %in% lv.names)
    tmp.mat[idx] <- "phi"
    tmp.row[idx] <- match(target$lhs[idx], lv.names)
    tmp.col[idx] <- match(target$rhs[idx], lv.names)

    # 4a. "~1" ov
    idx <- which(target$group == g & target$op == "~1" & !(target$lhs %in% lv.names))
    tmp.mat[idx] <- "tau"
    tmp.row[idx] <- match(target$lhs[idx], ov.names)
    tmp.col[idx] <- 1L

    # 4b. "~1" lv
    idx <- which(target$group == g & target$op == "~1" & target$lhs %in% lv.names)
    tmp.mat[idx] <- "kappa"
    tmp.row[idx] <- match(target$lhs[idx], lv.names)
    tmp.col[idx] <- 1L

    # catch lower-elements in psi/phi
    idx.lower <- which(tmp.mat %in% c("psi","phi") & tmp.row > tmp.col)
    if(length(idx.lower) > 0L) {
      tmp <- tmp.row[idx.lower]
      tmp.row[idx.lower] <- tmp.col[idx.lower]
      tmp.col[idx.lower] <- tmp
    }

    if(extra) {
      # mRows
      mmRows <- list(tau = nvar, lambda = nvar,  psi = nvar, kappa = nfac, phi = nfac)

      # mCols
      mmCols <- list(tau = 1L, lambda = nfac, psi = nvar, kappa = 1L, phi = nfac)

      # dimNames for LISREL model matrices
      mmDimNames <- list(tau    = list( ov.names,    "intercept"),
                         lambda = list( ov.names,       lv.names),
                         psi    = list( ov.names,       ov.names),
                         kappa  = list( lv.names,    "intercept"),
                         phi    = list( lv.names,       lv.names))

      # isSymmetric
      mmSymmetric <- list(tau = FALSE, lambda = FALSE, psi = TRUE, kappa  = FALSE, phi = TRUE)

      # which mm's do we need?
      IDX <- which(target$group == g)
      mmNames <- c("lambda", "psi", "phi") # (always include lambda, psi and phi)

      if(meanstructure) {
        mmNames <- c(mmNames, "tau", "kappa")
      }

      REP.mmNames[[g]]     <- mmNames
      REP.mmNumber[[g]]    <- length(mmNames)
      REP.mmRows[[g]]      <- unlist(mmRows[ mmNames ])
      REP.mmCols[[g]]      <- unlist(mmCols[ mmNames ])
      REP.mmDimNames[[g]]  <- mmDimNames[ mmNames ]
      REP.mmSymmetric[[g]] <- unlist(mmSymmetric[ mmNames ])
    } # extra

  } # ngroups

  REP <- list(mat = tmp.mat, row = tmp.row, col = tmp.col)

  if(extra) {
    attr(REP, "mmNames")     <- REP.mmNames
    attr(REP, "mmNumber")    <- REP.mmNumber
    attr(REP, "mmRows")      <- REP.mmRows
    attr(REP, "mmCols")      <- REP.mmCols
    attr(REP, "mmDimNames")  <- REP.mmDimNames
    attr(REP, "mmSymmetric") <- REP.mmSymmetric
  }

  REP
}

