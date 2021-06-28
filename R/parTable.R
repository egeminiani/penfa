##############################
# ---- Creating Partable -----
##############################
# Content:
# - ParTable : parameter table
# - partable_labels : get the labels for each parameter
# - partable_check  : checks the completeness and consistency of the partable
# - partable_attributes : returns the 'attributes' of a partable, and generates a new set if necessary
# - partable_vnames : returns variable names of a partable
# - partable_block_values : the block values  (not necessarily integers) of a parameter table
# - partable_nblocks : guess the number of blocks from a partable
# - partable_covariance_reorder : check the order of covariances (only fill the upper.tri, thus 'switch' lhs & rhs if appearing in wrong order)
# - partable_group_values : the group values of a parameter table
# - partable_ngroups : guess number of groups from a partable
# - partable_flat : flatten the partable
#
# Some of these functions are adaptations from the routines present in the lavaan package (https://CRAN.R-project.org/package=lavaan)
#
# ---------------------------------------------------------------------------------------------------------------

# ParTable
ParTable <- function(model           = NULL,
                     meanstructure   = FALSE,
                     int.ov.free     = FALSE,
                     int.lv.free     = FALSE,
                     orthogonal      = FALSE,
                     std.lv          = FALSE,
                     ngroups         = 1L,
                     auto.fix.first  = FALSE,
                     auto.fix.single = FALSE,
                     debug           = FALSE,
                     warn            = TRUE,
                     as.data.frame.  = TRUE){

  FLAT <- model
  # user-specified *modifiers* are returned as attributes
  MOD  <- attr(FLAT, "modifiers"); attr(FLAT, "modifiers") <- NULL

  if(any(FLAT$op == "~1")) meanstructure <- TRUE		# check for meanstructure

  # Parameter list without modifiers
  LIST <- partable_flat(FLAT,
                        meanstructure = meanstructure,
                        int.ov.free = int.ov.free,
                        int.lv.free = int.lv.free,
                        orthogonal = orthogonal,
                        std.lv = std.lv,
                        auto.fix.first = auto.fix.first,
                        auto.fix.single = auto.fix.single,
                        ngroups = ngroups)

  # Apply user-specified modifiers
  if(length(MOD)) {
    for(el in 1:length(MOD)) {
      idx <- which(LIST$mod.idx == el) # for each group

      # perhaps the corresponding element was duplicated, and removed
      if(length(idx) == 0L) {
        next
      }
      MOD.fixed <- MOD[[el]]$fixed
      MOD.start <- MOD[[el]]$start

      if (ngroups > 1L && length(idx) > 1L) {
        if (length(MOD.fixed) == 1L)
          MOD.fixed <- rep(MOD.fixed, ngroups)
        if (length(MOD.start) == 1L)
          MOD.start <- rep(MOD.start, ngroups)
      }

      # apply modifiers
      nidx <- length(idx)

      if((!is.null(MOD.fixed) && nidx != length(MOD.fixed)) ||
         (!is.null(MOD.start) && nidx != length(MOD.start))){
            el.idx <- which(LIST$mod.idx == el)[1L]
            stop("penfa ERROR: wrong number of arguments in modifier of element ",
                  LIST$lhs[el.idx], LIST$op[el.idx], LIST$rhs[el.idx])
      }

      if (!is.null(MOD.fixed)) {
        # two options: constant or NA
        na.idx <- which(is.na(MOD.fixed))
        not.na.idx <- which(!is.na(MOD.fixed))

        # constant
        LIST$ustart[idx][not.na.idx] <- MOD.fixed[not.na.idx]
        LIST$free[  idx][not.na.idx] <- 0L

        # NA* modifier
        LIST$free[  idx][na.idx] <- 1L # e.g. a factor loading
        LIST$ustart[idx][na.idx] <- as.numeric(NA)
      }
      if(!is.null(MOD.start)) {
        LIST$ustart[idx] <- MOD.start
      }
    }
  }
  # remove mod.idx column
  LIST$mod.idx <- NULL

  # count free parameters
  idx.free <- which(LIST$free > 0)
  LIST$free[idx.free] <- seq_along(idx.free)

  # data.frame?
  if(as.data.frame.) { LIST <- as.data.frame(LIST, stringsAsFactors = FALSE)  }
  LIST
}

# partable_labels
partable_labels <- function(partable, type = "user") {

  if (length(partable$lhs) == 0L)
    return(character(0L))

  # default labels
  label <- paste(partable$lhs, partable$op, partable$rhs, sep="")

  # groups
  if (is.character(partable$group)) {
    group.label <- unique(partable$group)
    group.label <- group.label[nchar(group.label) > 0L]
    ngroups <- length(group.label)
  }
  else {
    ngroups <- partable_ngroups(partable)
    group.label <- 1:ngroups
  }
  if (ngroups > 1L) {
    for (g in 2:ngroups) {
      label[partable$group == group.label[g]] <- paste(label[partable$group == group.label[g]], ".g", g, sep = "")
    }
  }

  user.idx <- which(nchar(partable$label) > 0L)
  # user-specified labels -- override everything
  label[user.idx] <- partable$label[user.idx]

  if (type == "user") {
    idx <- 1:length(label)
  }
  else if (type == "free") {
    idx <- which(partable$free > 0L & !duplicated(partable$free))
  }
  else {
    stop("penfa ERROR: argument `type' must be one of free or user")
  }
  label[idx]
}

# partable_check
partable_check <- function(partable, warn = TRUE) {

  check <- TRUE

  # check for empy table
  if(length(partable$lhs) == 0) return(check)

  # get observed/latent variables
  ov.names <- partable_vnames(partable, "ov")
  lv.names <- partable_vnames(partable, "lv")
  all.names <- c(ov.names, lv.names)

  # we should have a (residual) variance for *each* ov/lv
  var.idx <- which(partable$op == "~~" & partable$lhs == partable$rhs)
  missing.idx <- which(is.na(match(all.names, partable$lhs[var.idx])))
  if(length(missing.idx) > 0L) {
    check <- FALSE
    if(warn) {
      warning("penfa WARNING: parameter table does not contain (residual) variances for one or
              more variables: [", paste(all.names[missing.idx], collapse = " "), "]")
    }
  }

  # meanstructure?
  meanstructure <- any(partable$op == "~1")
  # if meanstructure, check for missing intercepts
  if(meanstructure) {
    # we should have a intercept for *each* ov/lv
    int.idx <- which(partable$op == "~1")
    missing.idx <- which(is.na(match(all.names, partable$lhs[int.idx])))
    if(length(missing.idx) > 0L) {
      check <- FALSE
      if(warn) {
        warning("penfa WARNING: parameter table does not contain intercepts for one or
                more variables: [", paste(all.names[missing.idx], collapse = " "), "]")
      }
    }
  }

  # do we have added intercepts (user = 0) that are fixed to zero?
  ov.ind <- unique(partable$rhs[partable$op == "=~"])
  lv.names <- unique(partable$lhs[partable$op == "=~"])
  int.fixed <- which(partable$op == "~1" & partable$user == 0L & partable$free == 0L &
                     partable$ustart == 0L & !partable$block == 1L & !partable$lhs %in% lv.names &
                    !partable$lhs %in% ov.ind)

  if(length(int.fixed) > 0L) {
    check <- FALSE
    if(warn) {
      warning("penfa WARNING: missing intercepts are set to zero: [", paste(partable$lhs[int.fixed],  collapse = " "), "]")
    }
  }

  # return check code
  check
}

# partable_attributes
partable_attributes <- function(partable, pta = NULL) {

  if(is.null(pta)) {
    # attached to partable?
    pta <- attributes(partable)
    if(!is.null(pta$vnames) && !is.null(pta$nvar)) {
      return(pta)
    } else {
      pta <- list()
    }
  }

  # vnames
  pta$vnames <- partable_vnames(partable, type="all")

  # vidx
  OV <- pta$vnames$ov
  LV <- pta$vnames$lv
  nblocks <- length(pta$vnames$ov)
  pta$vidx <- lapply(names(pta$vnames), function(v) {
    lapply(seq_len(nblocks), function(g) {
      if(grepl("lv", v)) {
        match(pta$vnames[[v]][[g]], LV[[g]])
      } else {
        match(pta$vnames[[v]][[g]], OV[[g]])
      }
    })
  })
  names(pta$vidx) <- names(pta$vnames)
  # meanstructure
  pta$meanstructure <- any(partable$op == "~1")

  # nblocks
  pta$nblocks <- nblocks

  # ngroups
  pta$ngroups <- partable_ngroups(partable)

  # nvar
  pta$nvar <- lapply(pta$vnames$ov, length)

  # nfac
  pta$nfac <- lapply(pta$vnames$lv, length)

  pta
}

# partable_vnames
# - the 'type' argument determines the status of the variable: "ov" (observed) ,"lv" (latent) or "all"
# - the 'group' argument either selects a single group (if group is an integer) or returns a list per group
partable_vnames <- function(partable, type = NULL, ..., warn = FALSE) {

  # check for empty table
  if(length(partable$lhs) == 0) return(character(0L))
  # dotdotdot
  dotdotdot <- list(...)

  type.list <- c("ov", "lv")  # observed (ov) or latent (lv) variables

  # sanity check
  stopifnot(is.list(partable), !missing(type), type %in% c(type.list, "all"))

  if(length(type) == 1L && type == "all") {
    type <- type.list
  }

  # always need `block' column -- create one if missing
  if(is.null(partable$block)) {
    partable$block <- rep(1L, length(partable$lhs))
  }

  nblocks <- 1       # nblocks -- block column is integer only
  block.select <- 1  # per default, use full partable

  # check for ... selection argument(s)
  ndotdotdot <- length(dotdotdot)
  if(ndotdotdot > 0L) {
    dot.names <- names(dotdotdot)
    block.select <- rep(TRUE, length(partable$lhs))
    for(dot in seq_len(ndotdotdot)) {
      # selection variable?
      block.var <- dot.names[dot]
      block.val <- dotdotdot[[block.var]]
      # do we have this 'block.var' in partable?
      if(is.null(partable[[block.var]])) {

        # treat "group = 1" special
        if(block.var == "group" && block.val == 1L) {
          partable$group <- rep(1L, length(partable$lhs))
          # remove block == 0
          idx <- which(partable$block == 0L)
          if(length(idx) > 0L) {
            partable$group[idx] <- 0L
          }
          block.select <- ( block.select & partable[[block.var]] %in% block.val )
        } else {
          stop("penfa ERROR: selection variable `", block.var, " not found in the parameter table.")
        }

      } else {
        if(!all(block.val %in% partable[[block.var]])) {
          stop("penfa ERROR: ", block.var , " column does not contain value `", block.val, "'")
        }
        block.select <- ( block.select & partable[[block.var]] %in% block.val )
      }
    } # dot
    block.select <- unique(partable$block[block.select])

    if(length(block.select) == 0L) {
      warnings("penfa WARNING: no blocks selected.")
    }
  }

  # output: list per block
  OUT     <- vector("list", length = nblocks)
  OUT$ov  <- vector("list", length = nblocks)
  OUT$lv  <- vector("list", length = nblocks)

  for(b in block.select) {
    # compute lv.names (regular latent variables only, defined by =~)
    lv.names <- unique( partable$lhs[ partable$block == b  & partable$op == "=~" ] )
    if("lv" %in% type) {
      OUT$lv[[b]] <- lv.names # store lv
    }

    # v.ind -- indicators of latent variables
    if(!(length(type) == 1L && type %in% "lv")) {
      v.ind <- unique( partable$rhs[ partable$block == b & partable$op == "=~" ] )
    }

    # ov.*
    if(!(length(type) == 1L && type %in% "lv")) {
      # indicators, which are not latent variables themselves
      ov.ind <- v.ind[ !v.ind %in% lv.names ]
    }

    # observed variables
    if(!(length(type) == 1L && type %in% "lv")) {

      # orphaned covariances
      ov.cov <- c(partable$lhs[ partable$block == b & partable$op == "~~" & !partable$lhs %in% lv.names ],
                  partable$rhs[ partable$block == b & partable$op == "~~" & !partable$rhs %in% lv.names ])
      # orphaned intercepts
      ov.int   <- partable$lhs[ partable$block == b & partable$op == "~1" & !partable$lhs %in% lv.names ]
      ov.extra <- unique(c(ov.cov, ov.int)) # must be in this order!
      ov.names <- c(ov.ind, ov.extra[ !ov.extra %in% ov.ind ])
    }

    if("ov" %in% type) {
      OUT$ov[[b]] <- ov.names # store ov
    }

  }

  if(length(type) == 1L) {
    OUT <- OUT[[type]]
    if(ndotdotdot == 0L) {
      OUT <- unique(unlist(OUT))
    } else if(length(block.select) == 1L) {
      OUT <- OUT[[block.select]]
    } else {
      OUT <- OUT[block.select]
    }
  } else {
    OUT <- OUT[type]
  }

  OUT
}

# partable_block_values
partable_block_values <- function(partable) {

  if(is.null(partable$block)) {
    block.values <- 1L
  } else {
    # always integers
    tmp <- partable$block[ partable$block > 0L ] # non-zero only
    block.values <- unique(na.omit(tmp))
  }
  block.values
}

# partable_nblocks
partable_nblocks <- function(partable) {
  length( partable_block_values(partable) )
}

# parable_covariance_reorder
partable_covariance_reorder <- function(partable, ov.names = NULL, lv.names = NULL) {

  # shortcut
  cov.idx <- which(partable$op == "~~" & partable$lhs != partable$rhs)
  if(length(cov.idx) == 0L) {
    # nothing to do
    return(partable)
  }

  # get names
  if(is.null(ov.names)) {
    ov.names <- partable_vnames(partable, "ov")
  } else {
    ov.names <- unlist(ov.names)
  }
  if(is.null(lv.names)) {
    lv.names <- partable_vnames(partable, "lv")
  } else {
    lv.names <- unlist(lv.names)
  }
  lv.ov.names <- c(lv.names, ov.names)

  # identify wrong ordering
  lhs.idx <- match(partable$lhs[ cov.idx ], lv.ov.names)
  rhs.idx <- match(partable$rhs[ cov.idx ], lv.ov.names)
  swap.idx <- cov.idx[ lhs.idx > rhs.idx ]

  if(length(swap.idx) == 0L) {
    # nothing to do
    return(partable)
  }

  # swap
  tmp <- partable$lhs[ swap.idx ]
  partable$lhs[ swap.idx ] <- partable$rhs[ swap.idx ]
  partable$rhs[ swap.idx ] <- tmp

  partable
}

# partable_group_values
partable_group_values <- function(partable) {

  if(is.null(partable$group)) {
    group.values <- 1L
  } else if(is.numeric(partable$group)) {
    tmp <- partable$group[ partable$group > 0L ]
    group.values <- unique(na.omit(tmp))
  } else { # character
    tmp <- partable$group[nchar(partable$group) > 0L]
    group.values <- unique(na.omit(tmp))
  }
  group.values
}

# partable_ngroups
partable_ngroups <- function(partable) {
  length( partable_group_values(partable) )
}

# partable_flat
partable_flat <- function(FLAT = NULL, meanstructure = FALSE, int.ov.free = FALSE, int.lv.free = FALSE,
                          orthogonal = FALSE, std.lv = FALSE, auto.fix.first = FALSE,
                          auto.fix.single  = FALSE, ngroups = 1L) {

  # extract `names' of various types of variables:
  lv.names     <- partable_vnames(FLAT, type="lv")
  ov.names     <- partable_vnames(FLAT, type="ov")
  lhs          <- rhs <- character(0)

  # 2. default (residual) variances and covariances
  # a) (residual) VARIANCES (all ov's except exo, and all lv's)
  ov.var <- ov.names # partable_vnames(FLAT, type="ov.nox")
  lhs <- c(lhs, ov.var, lv.names)
  rhs <- c(rhs, ov.var, lv.names)

  # b) `independent` latent variable COVARIANCES (lv.names.x)
  if(length(lv.names) > 1L) {
    tmp <- utils::combn(lv.names, 2)
    lhs <- c(lhs, tmp[1,]) # to fill upper.tri
    rhs <- c(rhs, tmp[2,])
  }

  # create 'op'
  op <- rep("~~", length(lhs))
  if (meanstructure) {
    ov.int  <- ov.names
    int.lhs <- c(ov.int, lv.names)
    lhs     <- c(lhs, int.lhs)
    rhs     <- c(rhs, rep("",  length(int.lhs)))
    op      <- c(op,  rep("~1",length(int.lhs)))
  }

  DEFAULT <- data.frame(lhs=lhs, op=op, rhs=rhs, mod.idx=rep(0L, length(lhs)), stringsAsFactors=FALSE)

  # 4. USER: user-specified elements
  lhs     <- FLAT$lhs
  op      <- FLAT$op
  rhs     <- FLAT$rhs
  mod.idx <- FLAT$mod.idx

  lv.names <- partable_vnames(FLAT, type="lv")     # latent variables
  ov.names <- partable_vnames(FLAT, type="ov")     # observed variables
  USER 		 <- data.frame(lhs=lhs, op=op, rhs=rhs, mod.idx=mod.idx, stringsAsFactors=FALSE)

  # check for duplicated elements in USER
  TMP <- USER[,1:3]
  idx <- which(duplicated(TMP))
  if(length(idx) > 0L) {
    txt <- sapply(1:length(idx), function(i) {
      paste("    ", TMP[idx[i],"lhs"], TMP[idx[i], "op"], TMP[idx[i],"rhs"]) })
    warning("penfa WARNING: duplicated elements in model syntax have been ignored:\n",
            paste(txt, collapse = "\n"))
    USER <- USER[-idx,]
  }

  # check for duplicated elements in DEFAULT
  TMP <- rbind(DEFAULT[,1:3], USER[,1:3])

  idx <- which(duplicated(TMP, fromLast=TRUE)) # idx should be in DEFAULT
  if(length(idx)) {
    for(i in idx) {
      flat.idx <- which(USER$lhs   == DEFAULT$lhs[i] & USER$op    == DEFAULT$op[i]  &  USER$rhs   == DEFAULT$rhs[i])
      if(length(flat.idx) != 1L) {
        cat("[DEBUG] idx in TMP: i = ", i, "\n"); print(TMP[i,])
        cat("[DEBUG] idx in DEFAULT: i = ", i, "\n"); print(DEFAULT[i,])
        cat("[DEBUG] flat.idx:"); print(flat.idx)
      }
    }
    DEFAULT <- DEFAULT[-idx,]
  }

  # After the removal of all duplicated elements, we can construct the LIST for a single group
  lhs     <- c(USER$lhs, DEFAULT$lhs)
  op      <- c(USER$op,  DEFAULT$op)
  rhs     <- c(USER$rhs, DEFAULT$rhs)
  user    <- c(rep(1L, length(USER$lhs)), rep(0L, length(DEFAULT$lhs)))
  mod.idx <- c(USER$mod.idx, DEFAULT$mod.idx)
  free    <- rep(1L,  length(lhs))
  ustart <- rep(as.numeric(NA), length(lhs))

  # 1. fix metric of regular latent variables
  if(std.lv) {
    # fix metric by fixing the variance of the latent variable
    lv.var.idx <- which(op == "~~" & lhs %in% lv.names & lhs == rhs)
    ustart[lv.var.idx] <- 1.0
    free[lv.var.idx] <- 0L
  }
  if (auto.fix.first) {
    mm.idx    <- which(op == "=~")
    first.idx <- mm.idx[which(!duplicated(lhs[mm.idx]))]
    ustart[first.idx] <- 1
    free[first.idx]   <- 0L
  }

  # 2. fix residual variance of single indicators to zero
  if(auto.fix.single) {
    mm.idx <- which(op == "=~")
    T <- table(lhs[mm.idx])
    if(any(T == 1L)) {
      # ok, we have a LV with only a single indicator
      lv.names.single <- names(T)[T == 1L]
      # get corresponding indicator if unique
      lhs.mm <- lhs[mm.idx]; rhs.mm <- rhs[mm.idx]
      single.ind <- rhs.mm[which(lhs.mm %in% lv.names.single & !(duplicated(rhs.mm) | duplicated(rhs.mm, fromLast=TRUE)))]
      # is the indicator unique?
      if(length(single.ind)) {
        var.idx <- which(op == "~~" & lhs %in% single.ind & rhs %in% single.ind & lhs == rhs & user == 0L)
        ustart[var.idx] <- 0.0
        free[var.idx] <- 0L
      }
    }
  }

  # 3. orthogonal=TRUE?
  if(orthogonal) {
    lv.cov.idx <- which(op == "~~" & lhs %in% lv.names & lhs != rhs & user == 0L)
    ustart[lv.cov.idx] <- 0.0
    free[lv.cov.idx] <- 0L
  }

  if (meanstructure) {
    if (int.ov.free == FALSE) {
      ov.int.idx <- which(op == "~1" & lhs %in% ov.names & user == 0L)
      ustart[ov.int.idx] <- 0
      free[ov.int.idx] <- 0L
    }
    if (int.lv.free == FALSE) {
      lv.int.idx <- which(op == "~1" & lhs %in% lv.names & user == 0L)
      ustart[lv.int.idx] <- 0
      free[lv.int.idx] <- 0L
    }
  }

  group <- rep(1L, length(lhs))
  if (ngroups > 1) {
    group   <- rep(1:ngroups, each = length(lhs))
    user    <- rep(user, times = ngroups)
    lhs     <- rep(lhs, times = ngroups)
    op      <- rep(op, times = ngroups)
    rhs     <- rep(rhs, times = ngroups)
    free    <- rep(free, times = ngroups)
    ustart  <- rep(ustart, times = ngroups)
    mod.idx <- rep(mod.idx, times = ngroups)
  }

  # construct LIST
  LIST <- list(id = seq_along(lhs), lhs = lhs, op = op, rhs = rhs, user = user)
  LIST[["group"]] <- group

  # other columns
  LIST2 <- list(mod.idx = mod.idx, free = free, ustart = ustart)
  LIST <- c(LIST, LIST2)
}
