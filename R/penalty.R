##################################
# ---- Penalty specification -----
##################################
# Content:
# - get_penfaPenalty : constructor for the penfaPenalty class
# - get_qstar : returns the indices of the parameters to penalize (for both shrinkage & differences)
# - penalty_mat : combines the penalty matrices deriving from the shrinkage & differences
# - Spens : returns the penalty matrix resulting from the local approximation of a certain penalty function
# - diff_mat : computes the matrix D of the pairwise differences
# - automatic : function for automatic multiple tuning parameter estimation (input: model with fixed tuning
#               parameters; output: optimal model with optimal vector of tuning parameters)
#
# ------------------------------------------------------------------------------------------------------------

# get_penfaPenalty
get_penfaPenalty <- function(model = NULL, options = NULL){

  penalty    <- list("shrink" = rep(options$pen.shrink, length(options$eta$shrink)),
                     "diff"   = rep(options$pen.diff,   length(options$eta$diff)))
  tuning     <- options$eta
  extra      <- list()

  if(any(unlist(penalty) == "scad")){
    extra[["a.scad"]] <- options$a.scad
  }
  if(any(unlist(penalty) == "mcp")){
    extra[["a.mcp"]] <- options$a.mcp
  }
  if(any(unlist(penalty) == "alasso")){
    extra[["a.alasso"]] <- options$a.alasso
  }

  pmat <- list("shrink" = names(options$eta$shrink), "diff" = names(options$eta$diff))

  if(any(unlist(penalty) == "alasso") & !is.null(options$weights))
    extra[["weights"]]  <- options$weights

  pen.idx    <- get_qstar(model = model, options = options)
  modpenalty <- new("penfaPenalty",
                    strategy  = options$strategy,
                    penalty   = penalty,
                    tuning    = tuning,
                    pmat      = pmat,
                    pen.idx   = pen.idx,
                    Sh.info   = list(),  # fill in the matrices after optimization
                    automatic = list(),  # fill in results of automatic procedure later on
                    extra     = extra)

  return(modpenalty)
}

# get_qstar
get_qstar <- function(model = NULL, options = NULL){

  GLIST <- model@GLIST
  # The indices of the parameters to penalize

  # Shrinkage
  mat.shrink.name   <- names(options$eta[["shrink"]])
  mat.shrink.name   <- mat.shrink.name[mat.shrink.name!="none"] # only penalized matrices
  qshrink           <- vector("list", length = length(mat.shrink.name))
  if(length(mat.shrink.name) > 0){
    names(qshrink)    <- mat.shrink.name
    for(i in 1:length(mat.shrink.name)){
      idx.name        <- which(names(GLIST) == mat.shrink.name[i])
      qshrink[[i]]    <- unlist(model@x.free.idx[idx.name])
    }
  }

  # Difference
  mat.diff.name   <- names(options$eta[["diff"]])
  mat.diff.name   <- mat.diff.name[mat.diff.name!="none"]
  qdiff           <- vector("list", length = length(mat.diff.name))
  if(length(mat.diff.name) > 0){
    names(qdiff)    <- mat.diff.name
    for(i in 1:length(mat.diff.name)){
      idx.name      <- which(names(GLIST) == mat.diff.name[i])
      qdiff[[i]]    <- unlist(model@x.free.idx[idx.name])
    }
  }

  qstar <- list("shrink" = qshrink, "diff" = qdiff)
  return(qstar)
}

# penalty_mat
penalty_mat <- function(x = NULL, model = NULL, modpenalty = NULL,
                        options = NULL, N = NULL){

  ngroups       <- model@ngroups
  pen.shrink    <- modpenalty@penalty$shrink
  pen.diff      <- modpenalty@penalty$diff
  eta.shrink    <- modpenalty@tuning$shrink
  eta.diff      <- modpenalty@tuning$diff
  m             <- model@nx.free
  qstar         <- modpenalty@pen.idx
  eps           <- options$cbar

  # Shrinkage
  Ss.shrink        <- vector("list", length = length(eta.shrink))
  names(Ss.shrink) <- names(eta.shrink)

  # For each penalty matrix to shrink
  for(i in 1:length(eta.shrink)){

    if(pen.shrink[i] == "none" | eta.shrink[i] == 0 ){
      Ss.shrink[[i]] <- diag(0, nrow = m, ncol = m)
    }else{
      # Matrix that penalizes all parameters
      tmp.mat <- N * Spens(params = x, model = model, modpenalty = modpenalty,
                           penalty = pen.shrink[i], eta = eta.shrink[i],
                           method = "shrink", eps = eps)
      # Only keep the elements of the matrix we want to penalize
      diag(tmp.mat)[-qstar$shrink[[i]] ] <- 0
      Ss.shrink[[i]] <- tmp.mat
    }
  }

  # Difference
  Ss.diff        <- vector("list", length = length(eta.diff))
  names(Ss.diff) <- names(eta.diff)

  # For each penalty matrix to compute difference penalty
  for(i in 1:length(eta.diff)){

    if(pen.diff[i] == "none" | eta.diff[i] == 0){
      Ss.diff[[i]] <- diag(0, nrow = m, ncol = m)
    }else{
      # Matrix that has differences between all parameters
      tmp.mat <- N * Spens(params = x, model = model, modpenalty = modpenalty,
                           penalty = pen.diff[i], eta = eta.diff[i],
                           method = "diff", eps = eps)
      # Only keep the differences of the elements that we want to penalize
      tmp.mat[-qstar$diff[[i]], ] <- tmp.mat[,-qstar$diff[[i]] ] <- 0
      Ss.diff[[i]] <- tmp.mat
    }
  }

  # Sum all shrinkage penalty matrices
  # Sum all difference matrices
  # and sum together these two elements
  S.h <- Reduce('+', Ss.shrink) + Reduce('+', Ss.diff)

  mat <- list("S.h" = S.h, "Ss.shrink" = Ss.shrink, "Ss.diff" = Ss.diff)
  return(mat)
}

# Spens
Spens <- function(params = NULL, model = NULL, modpenalty = NULL, penalty = "lasso",
                  eta = 0.001, method = "shrink", eps = 1e-08){

  # Get the R matrix
  if(method == "shrink"){
    R         <- diag(1, nrow = length(params))
  }else if(method == "diff"){
    nmat    <- model@nmat
    ngroups <- model@ngroups
    m       <- numeric(ngroups)
    for(g in seq_len(ngroups)){
      mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
      # dimension of group-specific parameter vector
      m[g] <- length(unique(unlist(model@x.free.idx[mm.in.group])))
    }
    # assume same number of parameters across groups
    R <- diff_mat(m1  = m[1], G = ngroups)
  } # end method shrink/diff

  Rtheta    <- as.numeric(R %*% params)
  Rtheta.sq <- Rtheta^2
  RtR       <- t(R) %*% R

  # Ridge
  if(penalty == "ridge"){
    factor <- rep(1, length(Rtheta))
    S      <- eta * (factor * RtR)
  }

  # Lasso
  if(penalty == "lasso"){
    factor <- 1/sqrt(Rtheta.sq + eps)
    S <- eta * (factor * RtR) * as.integer(Rtheta != 0)
  }

  # Adaptive lasso
  if(penalty == "alasso") {
    a.alasso <- modpenalty@extra$a.alasso
    w.alasso <- modpenalty@extra$weights
    if(is.null(w.alasso))
      w.alasso <- 1
    w.alasso <- as.numeric(R %*% w.alasso) # if diff, they are differences of MLE
    w.al     <- 1/abs(w.alasso)^a.alasso
    factor   <- w.al/sqrt(Rtheta.sq + eps)
    S        <- eta * (factor * RtR) * as.integer(Rtheta != 0)
  }

  # Scad
  if (penalty == "scad") {
    a.scad <- modpenalty@extra$a.scad
    theta  <- abs(Rtheta)
    p      <- length(Rtheta)
    f1     <- sapply(theta, function(theta) {
      max(a.scad * eta - theta, 0)/((a.scad - 1) * eta)
    })
    f.d    <- eta * ((theta <= eta) + f1 * (theta > eta))
    factor <- f.d/(sqrt(Rtheta.sq + eps) + 1e-06)
    S      <- factor * RtR
  }

  # Mcp
  if(penalty == "mcp"){
    a.mcp  <- modpenalty@extra$a.mcp
    theta  <- abs(Rtheta)
    f1     <- sapply(theta, function(theta) { eta - theta/a.mcp  } )
    f.d    <- f1 * (theta <= a.mcp * eta) + 0 * (theta > a.mcp * eta)
    factor <- f.d/(sqrt(Rtheta.sq + eps) + 1e-06)
    S      <- factor * RtR
  }

  return(S)
}

# diff_mat
diff_mat <- function(m1 = NULL, G = NULL){

  binom <- choose(G, 2)
  h1    <- diag(m1)
  h2    <- -diag(m1)

  ncol <- m1* G
  nrow <- m1 * binom

  Dmat <- matrix(0, nrow = nrow, ncol = ncol)

  seqs.col <- seq(1, ncol, by = m1)
  first.idx.col <- seqs.col
  last.idx.col  <- c(seqs.col[-1]-1,ncol)
  coords.col <- cbind(first.idx.col, last.idx.col)

  seqs.row <- seq(1, nrow, by = m1)
  first.idx.row <- seqs.row
  last.idx.row <- c(seqs.row[-1]-1, nrow)
  coords.row <- cbind(first.idx.row, last.idx.row)
  combs <- t(combn(G, 2))

  for(r in 1:binom){
    Dmat[coords.row[r,1] : coords.row[r, 2], coords.col[combs[r,1],1] : coords.col[combs[r,1], 2] ] <- h1
    Dmat[coords.row[r,1] : coords.row[r, 2], coords.col[combs[r,2],1] : coords.col[combs[r,2], 2] ] <- h2
  }

  return(Dmat)
}

# automatic
automatic <- function(model  = NULL,
                      data   = NULL,
                      syntax = NULL,
                      group  = NULL,
                      iterlimsp = 50,
                      tolsp     = 1e-07){


  gamma    <- model@Options$gamma
  fit      <- model               # the model with the tuning given in input
  options  <- fit@Options

  # Extract information on which penalization is there: shrinkage, difference or both?
  # Remove also tuning = 0, otherwise magic gives error
  # Separate the checks for shrinkage from those for differences, because they can be vector arguments
  flag <- names(unlist(fit@Penalize@penalty)[unlist(fit@Penalize@penalty)!="none" & unlist(fit@Penalize@tuning)!=0])
  # Remove numbers and duplicated elements
  flag <- unique(gsub('[[:digit:]]+', '', flag))

  # Supported penalties for automatic procedure: ridge, lasso, and alasso
  condition <- all(unlist(fit@Penalize@penalty)[unlist(fit@Penalize@penalty)!= "none"] %in% c("ridge", "lasso", "alasso"))
  if(condition == FALSE){
    stop("\n penfa ERROR: The automatic procedure only supports lasso, alasso, or ridge penalties. \nRespecify the model with a different penalty. \n")
  }


  iter.if     <- model@Optim$iterations
  bs.mgfit    <- wor.c <- magpp <- NULL
  stoprule.SP <- 1
  conv.sp     <- TRUE
  iter.inner  <- iter.sp <- 0

  # The input tuning must be != 0, otherwise magic gives error
  if(length(flag) != 0){ # Penalized model, go on with procedure

    if(options$verbose){
      cat("\nAutomatic procedure: \n")
    }



    while(stoprule.SP > tolsp){  # start tuning loop

      fit2   <- list()
      fito   <- fit@Optim$fx.unpen
      o.ests <- fit2$argument <- c(fit@Optim$x)

      fit2$gradient <- fit@Optim$dx.pen
      fit2$S.h2     <- fit@Penalize@Sh.info$S.h2
      fit2$hessian  <- fit@Optim$hessian.pen
      fit2$S.h      <- fit@Penalize@Sh.info$S.h

      # In fit@Options there is the updated tuning

      # Step 0: Get the penalty matrices
      if(length(flag)==1){
        if(flag == "shrink"){
          tuning <- fit@Penalize@tuning$shrink
          # Get a list of separate penalty matrices for each to-be-penalized parameter matrix
          Ss <- penalty_mat(x = o.ests, model = fit@Model, modpenalty = fit@Penalize,
                            options = fit@Options, N = fit@SampleStats@ntotal)$Ss.shrink
          Ss <- mapply("/", Ss, fit@Penalize@tuning$shrink, SIMPLIFY = FALSE)
        }else if(flag == "diff"){
          tuning <- fit@Penalize@tuning$diff
          # Get a list of separate penalty matrices for each to-be-penalized parameter matrix
          Ss   <- penalty_mat(x = o.ests, model = fit@Model, modpenalty = fit@Penalize,
                              options = fit@Options, N = fit@SampleStats@ntotal)$Ss.diff
          Ss   <- mapply("/", Ss, fit@Penalize@tuning$diff, SIMPLIFY = FALSE)
        }

      }else if(length(flag)==2){

        tuning <- unlist(fit@Penalize@tuning)

        # Get all the penalty matrices
        all.pmat  <- penalty_mat(x = o.ests, model = fit@Model, modpenalty = fit@Penalize,
                                 options = fit@Options, N = fit@SampleStats@ntotal)
        # Extract the shrinkage matrices and divide them by their corresponding tunings
        Ss.shrink <- mapply("/", all.pmat$Ss.shrink, fit@Penalize@tuning$shrink, SIMPLIFY = FALSE)

        # Extract the difference matrices and divide them by their corresponding tunings
        Ss.diff   <- mapply("/", all.pmat$Ss.diff, fit@Penalize@tuning$diff, SIMPLIFY = FALSE)

        # or (jointly)
        # Ss <- c(all.pmat$Ss.shrink, all.pmat$Ss.diff)
        # Ss <- mapply("/", Ss, tuning, SIMPLIFY = FALSE)

        # Combine
        Ss <- c(Ss.shrink, Ss.diff)
      }

      wor.c     <- GJRM::working.comp(fit = fit2)

      # Step 1: Update the tuning
      bs.mgfit  <- try(mgcv::magic(y     = wor.c$Z,
                                   X     = wor.c$X,
                                   sp    = tuning,
                                   S     = Ss,
                                   off   = rep(1, length(Ss)),
                                   rank  = NULL,
                                   gcv   = FALSE,
                                   gamma = gamma),
                       silent = TRUE)

      if(class(bs.mgfit)=="try-error") {
        conv.sp <- FALSE
        break
      }

      if(any(is.na(bs.mgfit$sp))){ # magic produced Nas for the tunings?
        # restart the process from 0.01 (hoping it works from there)
        tuning <- rep(0.01, times = length(bs.mgfit$sp))
      }else{ # tuning ok, keep the updated values
        tuning <- bs.mgfit$sp
      }

      iter.sp <- iter.sp + 1
      if(options$verbose){
        cat("Iteration ", iter.sp, ":",
            sapply(tuning, function(x){ format(round(x, 8), nsmall = 8) }),
            "\n")
      }

      tmp         <- list()
      tmp$verbose <- FALSE; tmp$strategy <- "fixed"

      ## updated tuning!
      if(length(flag) == 1){
        if(flag == "shrink"){                   # either shrinkage or difference
          tmp$eta <- list("shrink" = tuning)
          names(tmp$eta[["shrink"]]) <- fit@Penalize@pmat$shrink
        }else if(flag == "diff"){
          tmp$eta <- list("diff" = tuning)
          names(tmp$eta[["diff"]])   <- fit@Penalize@pmat$diff
        }
      }else if(length(flag) == 2){            # both shrinkage and difference
        tmp$eta <- list("shrink" = tuning[1:length(options$eta$shrink)],
                        "diff"   = tuning[-c(1:length(options$eta$shrink))])
        names(tmp$eta[["shrink"]]) <- fit@Penalize@pmat$shrink
        names(tmp$eta[["diff"]])   <- fit@Penalize@pmat$diff
      }

      tmp$user.start <- TRUE; tmp$start.val <- o.ests # but old parameters
      tmp.arg <- utils::modifyList(fit@Options, tmp)

      # Step 2: Refit the model to update the parameters
      fit <- try(do.call(penfa, args = c(tmp.arg, list(model = syntax, data = data,
                                                       group = group))),
                 silent = TRUE)

      if(class(fit) == "try-error" | is.null(fit@Optim$fx.unpen) | !fit@Vcov$admissibility) {
        # errors? Try with default starting values
        conv.sp <- FALSE
        tmp     <- list(); tmp$strategy <- "fixed"; tmp$user.start <- FALSE
        tmp.arg <- utils::modifyList(fit@Options, tmp)
        fit     <- try(do.call(penfa, args = c(tmp.arg, list(model = syntax,
                                                             data = data, group = group))),
                       silent = TRUE)
        if(class(fit) == "try-error" | is.null(fit@Optim$fx.unpen))
          stop("penfa ERROR: Error in the automatic tuning procedure. Try with a different starting value for the tuning.")
      }

      iter.inner <- iter.inner + fit@Optim$iterations

      if(iter.sp >= iterlimsp){
        conv.sp <- FALSE
        break
      }

      # Stop rule
      stoprule.SP <- abs(fit@Optim$fx.unpen - fito)/(0.1 + abs(fit@Optim$fx.unpen))

    } # end tuning fitting loop


    magpp <- try(mgcv::magic.post.proc(wor.c$X, bs.mgfit),
                 silent = TRUE)
    if(inherits(magpp, "try-error")){
      stop("penfa ERROR: Error in the automatic tuning procedure.")
    }

  }else{
    # Unpenalized model, return eta = 0
    warning("\n penfa WARNING: Cannot execute the automatic procedure. \n Choose a penalty different from \"none\" or a non-null starting value for the tuning.")
    tuning <- rep(0, length(unlist(fit@Penalize@tuning)))
    fit@Options$strategy <- "fixed" # necessary for Verbose output in penfa Step 17
  }

  res <- list(model      = fit,
              tuning     = tuning,
              iterlim    = iterlimsp,
              iter.inner = iter.inner,
              iter       = iter.sp,
              tol        = tolsp,
              gamma      = gamma,
              conv       = conv.sp,
              R          = bs.mgfit$R)

  return(res)
}
