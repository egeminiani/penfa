#####################
# ---- penfaData ----
#####################
# Content:
# - penfaData : constructor for the 'penfaData' class describing how the data look like
# - penfa_data_full: handles the full data set
# - dataframe_vartable : returns a data frame with the data set statistics
#
# Some of these functions are adaptations from the routines present in the lavaan package (https://CRAN.R-project.org/package=lavaan)
#
# -------------------------------------------------------------------------------------------------------------

# penfaData
penfaData <- function(data              = NULL,             # data.frame
                      group             = NULL,             # multiple-group analysis
                      ov.names          = NULL,             # variables in model
                      options           = penfaOptions()    # options
){

  # get info from options
  group.label <- options$group.label
  if(is.null(group.label)) {
    group.label <- character(0L)
  }
  std.ov <- options$std.ov
  if(is.null(std.ov)) {
    std.ov <- FALSE
  }
  warn <- options$warn
  if(is.null(warn)) {
    warn <- TRUE
  }

  # Only possible scenario (for now): the data is a full data.frame (or a matrix)
  if(!is.null(data)){
    # catch matrix
    if(!is.data.frame(data)) {
      if(is.matrix(data)) {     # matrix?
        if(nrow(data) == ncol(data)) { # covariance matrix?
          if(data[2,1] == data[1,2] && warn) {
            warning("penfa WARNING: the data look like a covariance matrix; please provide the full data set")
          }
        }
        data <- as.data.frame(data, stringsAsFactors = FALSE)
      } else {
        stop("penfa ERROR: data object of class ", class(data))
      }
    }

    if(is.null(ov.names)) {     # no ov.names
      ov.names <- names(data)
      # remove group variable, if provided
      if(length(group) > 0L) {
        group.idx <- which(ov.names == group)
        ov.names <- ov.names[-group.idx]
      }
    }

    modData <- penfa_data_full(data = data, group = group, group.label = group.label,
                               ov.names = ov.names, std.ov = std.ov, warn = warn)
  }

  # construct penfaData object
  modData2 <- new("penfaData",
                  ngroups     = modData$ngroups,
                  group       = modData$group,
                  group.label = modData$group.label,
                  std.ov      = modData$std.ov,
                  nobs        = modData$nobs,
                  norig       = modData$norig,
                  ov.names    = modData$ov.names,
                  ov          = modData$ov,
                  case.idx    = modData$case.idx,
                  X           = modData$X)
  modData2
}

# penfa_data_full
penfa_data_full <- function(data          = NULL,          # data.frame
                            group         = NULL,          # multiple-group analysis
                            group.label   = NULL,          # custom group label
                            ov.names      = NULL,          # variables in model
                            std.ov        = FALSE,         # standardize ov's?
                            warn          = TRUE           # produce warnings?
){
  # number of groups and group labels
  if(!is.null(group) && length(group) > 0L) {
    if(!(group %in% names(data))) {
      stop("penfa ERROR: grouping variable ", sQuote(group), " not found;\n  ",
           "variable names found in data frame are:\n  ", paste(names(data), collapse=" "))
    }
    # note: by default, use the order as in the data, not as in levels(data[,group])
    if(length(group.label) == 0L) {
      group.label <- unique(as.character(data[[group]]))
      if(warn && any(is.na(group.label))) {
        warning("penfa WARNING: group variable ", sQuote(group), " contains missing values\n", sep="")
      }
      group.label <- group.label[!is.na(group.label)]
    } else {
      group.label <- unique(as.character(group.label))
      # check if user-provided group labels exist
      LABEL <- unique(as.character(data[[group]]))
      idx <- match(group.label, LABEL)
      if(warn && any(is.na(idx))) {
        warning("penfa WARNING: some group.labels do not appear ",
                "in the grouping variable: ", paste(group.label[which(is.na(idx))], collapse=" "))
      }
      group.label <- group.label[!is.na(idx)]
      # any groups left?
      if(length(group.label) == 0L)
        stop("penfa ERROR: no group levels left; check the group.label argument")
    }
    ngroups <- length(group.label)
  } else {
    if(warn && length(group.label) > 0L)
      warning("penfa WARNING: `group.label' argument will be ignored if `group' argument is missing")
    ngroups <- 1L
    group.label <- character(0L)
    group <- character(0L)
  }

  # check ov.names vs ngroups
  if(ngroups > 1L) {
    if(is.list(ov.names)) {
      if(length(ov.names) != ngroups)
        stop("penfa ERROR: ov.names assumes ", length(ov.names), " groups; data contains ", ngroups, " groups")
    } else {
      tmp <- ov.names
      ov.names <- vector("list", length = ngroups)
      ov.names[1:ngroups] <- list(tmp)
    }
  } else {
    if(is.list(ov.names)) {
      if(length(ov.names) > 1L)
        stop("penfa ERROR: model syntax defines multiple groups; data suggests a single group")
    } else {
      ov.names <- list(ov.names)
    }
  }

  # check if all ov.names can be found in the data.frame
  for(g in 1:ngroups) {
    # does the data contain all the observed variables needed in the user-specified model for this group
    ov.all <- unique(ov.names[[g]])

    # handle interactions
    ov.int.names <- ov.all[ grepl(":", ov.all) ]
    n.int <- length(ov.int.names)
    if(n.int > 0L) {
      ov.names.noint <- ov.all[!ov.all %in% ov.int.names]
      for(iv in seq_len(n.int)) {
        NAMES <- strsplit(ov.int.names[iv], ":", fixed = TRUE)[[1L]]
        if(all(NAMES %in% ov.names.noint)) {
          # add this interaction term to the data.frame, unless it already exists
          if(is.null(data[[ ov.int.names[iv] ]])) {
            data[[ ov.int.names[iv] ]] <- data[[NAMES[1L]]] * data[[NAMES[2L]]]
          }
        }
      }
    }

    # check for missing observed variables
    idx.missing <- which(!(ov.all %in% names(data)))

    if(length(idx.missing)) {
      stop("penfa ERROR: missing observed variables in dataset: ",
           paste(ov.all[idx.missing], collapse=" "))
    }
  }

  # here, we know for sure all ov.names exist in the data.frame
  # create varTable
  ov <- dataframe_vartable(frame = data, ov.names = ov.names, as.data.frame. = FALSE)

  # check for mix small/large variances
  if(!std.ov && warn && any(ov$type == "numeric")) {
    num.idx <- which(ov$type == "numeric")
    if(length(num.idx) > 0L) {
      min.var <- min(ov$var[num.idx])
      max.var <- max(ov$var[num.idx])
      rel.var <- max.var/min.var
      if(warn && rel.var > 1000) {
        warning("penfa WARNING: some observed variances are (at least) a factor 1000 times larger than others")
      }
    }
  }

  # check for really large variances (perhaps -999999 for missing?)
  if(!std.ov && warn && any(ov$type == "numeric")) {
    num.idx <- which(ov$type == "numeric")
    if(length(num.idx) > 0L) {
      max.var <- max(ov$var[num.idx])
      if(warn && max.var > 1000000) {
        warning("penfa WARNING: some observed variances are larger than 1000000\n")
      }
    }
  }

  # prepare empty lists, group-based
  case.idx <- vector("list", length = ngroups)
  norig    <- vector("list", length = ngroups)
  nobs     <- vector("list", length = ngroups)
  X        <- vector("list", length = ngroups)

  # collect information per group
  for(g in 1:ngroups) {

    # extract variables in correct order
    ov.idx  <- ov$idx[match(ov.names[[g]],   ov$name)]
    all.idx <- unique(c(ov.idx))

    # extract cases per group
    if(ngroups > 1L || length(group.label) > 0L) {
        # missing == listwise
        case.idx[[g]] <- which(data[[group]] == group.label[g] & complete.cases(data[all.idx]))
        nobs[[g]] <- length(case.idx[[g]])
        norig[[g]] <- length(which(data[[group]] == group.label[g]))
    } else {
      # missing == "listwise"
        case.idx[[g]] <- which(complete.cases(data[all.idx]))
        nobs[[g]] <- length(case.idx[[g]])
        norig[[g]] <- nrow(data)
    }

    # extract data
    X[[g]] <- data.matrix( data[case.idx[[g]], ov.idx, drop = FALSE] )
    dimnames(X[[g]]) <- NULL

    # standardize observed variables? numeric only!
    if(std.ov) {
      num.idx <- which(ov$name %in% ov.names[[g]] & ov$type == "numeric" )
      if(length(num.idx) > 0L) {
        X[[g]][,num.idx] <- scale(X[[g]][,num.idx,drop = FALSE])[,,drop = FALSE]
      }
    }

    # warn if we have a small number of observations (but NO error)
    if( warn && nobs[[g]] < (nvar <- length(ov.idx)) ) {
      txt <- ""
      if(ngroups > 1L) txt <- paste(" in group ", g, sep="")
      warning("penfa WARNING: small number of observations (nobs < nvar)", txt,
              "\n  nobs = ", nobs[[g]], " nvar = ", nvar)
    }

  } # groups

  out <- list(ngroups = ngroups, group = group, group.label = group.label,
              std.ov = std.ov, nobs = nobs, norig = norig, ov.names = ov.names,
              ov = ov, case.idx = case.idx, X = X)
  out
}

# dataframe_vartable
dataframe_vartable <- function(frame = NULL, ov.names = NULL,
                               as.data.frame. = FALSE) {

  if(missing(ov.names)){
    var.names <- names(frame)
  } else{
    ov.names <- unlist(ov.names, use.names=FALSE)
    var.names <- unique(c(ov.names))
  }
  nvar    <- length(var.names)
  var.idx <- match(var.names, names(frame))
  nobs    <- integer(nvar)
  type    <- character(nvar)
  mean    <- numeric(nvar)
  var     <- numeric(nvar)

  for(i in seq_len(nvar)) {
    x <- frame[[var.idx[i]]]
    type.x <- class(x)[1L]

    # correct for matrix with 1 column
    if(type.x == "matrix" && ncol(x) == 1L) {
      type.x <- "numeric"
    }

    # correct for integers
    if(type.x == "integer") {
      type.x <- "numeric"
    }

    type[i] <- type.x
    nobs[i] <- sum(!is.na(x))
    mean[i] <- ifelse(type.x == "numeric", mean(x, na.rm=TRUE), as.numeric(NA))
    var[i]  <- ifelse(type.x == "numeric", var(x, na.rm=TRUE),  as.numeric(NA))
  }

  VAR <- list(name=var.names, idx=var.idx, nobs=nobs,
              type=type, mean=mean, var=var)

  if(as.data.frame.) {
    VAR <- as.data.frame(VAR, stringsAsFactors=FALSE,
                         row.names=1:length(VAR$name))
  }
  VAR
}
