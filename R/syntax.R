############################################
# ---- Options and Parsing model syntax ----
############################################
# Content:
# - penfaOptions : set the default options
# - ParseModelString : parser of the model syntax
# - syntax_parse_rhs : parser of the right-hand side of the formula
# - syntax_get_modifier : get the modifiers (e.g. fixed elements, starting values, etc.) in the syntax
#
# Some of these functions are adaptations from the routines present in the lavaan package (https://CRAN.R-project.org/package=lavaan)
# -------------------------------------------------------------------------------------------------------------

#' \code{penfa} Options
#'
#' @description The default options internally used by the \code{\link{penfa}}
#'   function. These options can be changed by passing "name = value" arguments to
#'   the \code{penfa} function call, where they are being added to the "..."
#'   argument.
#'
#' @param opt List of default options. See below for details.
#'
#' @details The following section details the full list of options currently
#'   accepted by the \code{penfa} function.
#'
#'   Model features:
#'
#'\describe{
#'  \item{\code{meanstructure}:}{Logical. If \code{TRUE}, a meanstructure is
#'  requested. It should be used in conjunction with \code{int.ov.free} and
#'  \code{int.lv.free} or intercept-like formulas in the model syntax.
#'  Default to \code{FALSE}.}
#'  \item{\code{int.ov.free}:}{Logical. If \code{FALSE}, the intercepts of the
#'  observed variables are fixed to zero. Default to \code{FALSE}.}
#'  \item{\code{int.lv.free}:}{Logical. If \code{FALSE}, the intercepts of the
#'  common factors are fixed to zero. Default to \code{FALSE}.}
#'  \item{\code{orthogonal}:}{Logical. If \code{TRUE}, all covariances among the
#'  common factors are set to zero. Default to \code{FALSE}.}
#'  \item{\code{std.lv}:}{Logical. If \code{TRUE}, the factor variances are fixed
#'  to 1.0. Default to \code{FALSE}.}
#'  \item{\code{auto.fix.first}:}{Logical. If \code{TRUE}, the factor loading of
#'  the first indicator is set to 1.0 for every factor. Default to \code{FALSE}.}
#'  \item{\code{auto.fix.single}:}{Logical. If \code{TRUE}, the residual variance
#'  (if included) of an observed indicator is set to zero if it is the only
#'  indicator of a common factor. Default to \code{FALSE}.}
#'}
#'
#' Data options:
#' \describe{
#'   \item{\code{std.ov}:}{Logical. If \code{TRUE}, all observed variables are
#'   standardized before entering the analysis. Default to \code{FALSE}.}
#' }
#'
#' Estimation and optimization:
#' \describe{
#'   \item{\code{information}:}{Character. If \code{"fisher"}, the penalized
#'   expected Fisher information matrix is used as second-order derivatives
#'   in the trust-region algorithm and for computing the standard errors of the
#'   model parameters. If \code{"hessian"}, the penalized Hessian matrix
#'   is used. Default to \code{"fisher"}.}
#'   \item{\code{control}:}{A list containing control parameters passed to the
#'   trust-region optimizer. See the manual page of \code{trust} from the
#'   \code{trust} package for an overview of its control parameters. Default
#'   values for these parameters are \code{rinit=1L}, \code{rmax=100L},
#'   \code{iterlim=1000L}, \cr \code{fterm = sqrt(.Machine$double.eps)},
#'   \code{mterm = sqrt(.Machine$double.eps)}. }
#'   \item{\code{optim.dx.tol}}{Numeric. The tolerance value used when checking
#'   the size of the elements of the gradient of the objective function. Default
#'   equal to 100.}
#'  }
#'
#'
#' Penalization:
#' \describe{
#'  \item{\code{a.scad}}{Numeric. The shape parameter for the scad penalty.
#'  Default to 3.7, as recommended by Fan & Li (2001).}
#'  \item{\code{a.mcp}}{Numeric. The shape parameter of the mcp penalty. Default
#'  to 3.}
#'  \item{\code{a.alasso}}{Numeric. The exponent in the adaptive weights for the
#'  alasso penalty. Default to 1.}
#'  \item{\code{weights}}{Numeric. Only valid when either \code{pen.shrink} or
#'  \code{pen.diff} is equal to "alasso". An optional vector of values provided
#'  by the user representing a consistent estimate for each model parameter. The
#'  vector is then internally used for computing the adaptive weights. If
#'  unspecified, the maximum likelihood estimates (MLE) from the unpenalized
#'  model are used. }
#'  \item{\code{cbar}}{Numeric. Numerical constant used in the local approximation
#'  of the penalty functions. Default to 1e-08.}
#'
#'
#'  Automatic procedure:
#'  \item{\code{gamma}}{Numeric. The value of the influence factor used in the
#'  automatic tuning parameter procedure. Default to 4.}
#'  \item{\code{user.start}}{Logical whether the user has provided a vector of
#'  starting values for the model parameter estimates. }
#'  \item{\code{start.val}}{Numeric. An optional vector of parameter estimates to
#'  be used as starting values for the model parameters. This option is also
#'  internally used by the automatic procedure.}
#'  }
#'
#' Verbosity options:
#' \describe{
#'   \item{\code{verbose}:}{Logical. If \code{TRUE}, some information on the
#'   estimation process (e.g., convergence and admissibility checks, effective
#'   degrees of freedom) are printed out. Default to \code{TRUE}.}
#'  \item{\code{warn}:}{Logical. If \code{TRUE}, some warnings are printed out
#'  during the iterations. Default to \code{TRUE}.}
#'  \item{\code{debug}:}{Logical. If \code{TRUE}, debugging information is
#'  printed out. Default to \code{FALSE}.}
#'}
#'
#' @return A list of default options internally used by the \code{\link{penfa}}
#'   function.
#'
#' @export
#'
#'
penfaOptions <- function(opt = list(
  # Options for the factor model
  meanstructure      = FALSE,
  int.ov.free        = FALSE,
  int.lv.free        = FALSE,
  orthogonal         = FALSE,
  std.lv             = FALSE,
  auto.fix.first     = FALSE,
  auto.fix.single    = FALSE,

  # Option for the data
  std.ov             = FALSE,

  # Estimation & optimization options
  information        = "fisher",
  control            = list(),
  optim.dx.tol       = 100,

  # Options for penalization
  a.scad             = 3.7,
  a.mcp              = 3,
  a.alasso           = 1,
  weights            = NULL,
  cbar               = 1e-08,

  # Options for automatic procedure
  gamma              = 4,
  user.start         = FALSE,
  start.val          = c(),

  # Options for debugging
  verbose            = TRUE,
  warn               = TRUE,
  debug              = FALSE)) {

  opt
}


# ParseModelString
ParseModelString <- function(model.syntax = '', as.data.frame. = FALSE, warn = TRUE, debug = FALSE){

  # check for empty syntax
  if(length(model.syntax) == 0) {
    stop("penfa ERROR: empty model syntax")
  }

  # Remove comments prior to split
  # Match from comment character to newline, but don't eliminate newline
  model.syntax <- gsub("[#!].*(?=\n)", "", model.syntax, perl = TRUE)
  # Replace semicolons with newlines prior to split
  model.syntax <- gsub(";", "\n", model.syntax, fixed = TRUE)
  # Remove all whitespace prior to split
  model.syntax <- gsub("[ \t]+", "", model.syntax, perl = TRUE)
  # Remove any occurrence of >= 2 consecutive newlines to eliminate
  # blank statements; this retains a blank newline at the beginning,
  # if such exists, but parser will not choke because of start.idx
  model.syntax <- gsub("\n{2,}", "\n", model.syntax, perl = TRUE)

  # Break up in lines
  model <- unlist( strsplit(model.syntax, "\n") )
  # Check for multi-line formulas: they contain no operator symbol
  # but before we do that, we remove all strings between double quotes
  # to avoid confusion with for example equal("f1=~x1") statements
  model.simple <- gsub("\\\".[^\\\"]*\\\"", "LABEL", model)

  start.idx <- grep("[~=<>:|%]", model.simple)
  # Check for empty start.idx: no operator found
  if (length(start.idx) == 0L) {
    stop("penfa ERROR: model does not contain syntax (no operator found)")
  }
  # Check for non-empty string, without an operator in the first lines
  if (start.idx[1] > 1L) {
    # two possibilities:
    # - we have an empty line (ok)
    # - the element contains no operator (warn)
    for (el in 1:(start.idx[1] - 1L)) {
      if (nchar(model.simple[el]) > 0L) {
        warning("penfa WARNING: no operator found in this syntax line: ",
                model.simple[el], "\n", "  This syntax line will be ignored!")
      }
    }
  }

  end.idx <- c( start.idx[-1]-1, length(model) )
  model.orig    <- model
  model <- character( length(start.idx) )
  for(i in 1:length(start.idx)) {
    model[i] <- paste(model.orig[start.idx[i]:end.idx[i]], collapse="")
  }

  # ok, in all remaining lines, we should have an operator outside the ""
  model.simple <- gsub("\\\".[^\\\"]*\\\"", "LABEL", model)
  idx.wrong <- which(!grepl("[~=<>:|%]", model.simple))
  if(length(idx.wrong) > 0) {
    cat("Missing operator in formula(s):\n")
    print(model[idx.wrong])
    stop("penfa ERROR: syntax error in model syntax")
  }

  # but perhaps we have a '+' as the first character?
  idx.wrong <- which(grepl("^\\+", model))
  if(length(idx.wrong) > 0) {
    cat("Some formula(s) start with a plus (+) sign:\n")
    print(model[idx.wrong])
    stop("penfa ERROR: syntax error in model syntax")
  }

  # Main operation: flatten formulas into single bivariate pieces
  # with a left-hand-side (lhs), an operator (eg "=~"), and a right-hand-side (rhs)
  # both lhs and rhs can have a modifier
  FLAT.lhs         <- character(0)
  FLAT.op          <- character(0)
  FLAT.rhs         <- character(0)
  FLAT.rhs.mod.idx <- integer(0)
  FLAT.block       <- integer(0)    # keep track of groups using ":" operator
  FLAT.fixed       <- character(0)  # only for display purposes!
  FLAT.start       <- character(0)  # only for display purposes!
  FLAT.idx <- 0L
  MOD.idx  <- 0L
  MOD <- vector("list", length=0L)
  BLOCK <- 1L
  BLOCK_OP <- FALSE
  for(i in 1:length(model)) {
    x <- model[i]
    if(debug) {
      cat("formula to parse:\n"); print(x); cat("\n")
    }

    # 1. which operator is used?
    line.simple <- gsub("\\\".[^\\\"]*\\\"", "LABEL", x)

    # Admissible operators for factor analysis: =~, ~~, ~
    if(grepl("=~", line.simple, fixed=TRUE)) {        # "=~" operator?
      op <- "=~"
    } else if(grepl("~~", line.simple, fixed=TRUE)) { # "~~" operator?
      op <- "~~"
    } else if(grepl("~", line.simple, fixed=TRUE)) {  # "~" operator?
      op <- "~"
    } else {
      stop("unknown operator in ", model[i])
    }

    # 2. split by operator (only the *first* occurence!)
    # check first if equal/label modifier has been used on the LEFT
    if(substr(x,1,6) == "label(")
      stop("label modifier can not be used on the left-hand side of the operator")
    op.idx <- regexpr(op, x)

    lhs <- substr(x, 1L, op.idx-1L)
    rhs <- substr(x, op.idx+attr(op.idx, "match.length"), nchar(x))

    # check if first character of rhs is '+'; if so, remove silently
    if(substr(rhs, 1, 1) == "+") {
      rhs <- substr(rhs, 2, nchar(rhs))
    }

    # 3. parse left hand

    # first check if all lhs names are valid (in R); see ?make.names and ?reserved
    LHS <- strsplit(lhs, split = "+", fixed = TRUE)[[1]]
    # remove modifiers
    LHS <- gsub("^\\S*\\*", "", LHS)
    if( !all(make.names(LHS) == LHS) ) {
      stop("penfa ERROR: left hand side (lhs) of this formula:\n  ", lhs, " ", op, " ", rhs,
           "\n    contains either a reserved word (in R) or an illegal charachter: ",
           dQuote(LHS[!make.names(LHS) == LHS]), "\n    see ?reserved for a list of reserved words in R",
           "\n    please use a variable name that is not a reserved word in R",
           "\n    and use only characters, digits, or the dot symbol.")
    }

    lhs.formula <- as.formula(paste("~",lhs))
    out <- syntax_parse_rhs(rhs=lhs.formula[[2L]], op = op)
    lhs.names <- names(out)
    if (sum(sapply(out, length)) > 0L) {
      warning("penfa WARNING: left-hand side of formula below contains modifier:\n", x, "\n")
    }

    # 4. syntax_parse_rhs (as rhs of a single-sided formula)

    # before we do this, replace '0.2?' by 'start(0.2)*'
    rhs <- gsub('\\(?([-]?[0-9]*\\.?[0-9]*)\\)?\\?',"start(\\1)\\*", rhs)
    rhs.formula <- as.formula(paste("~",rhs))
    out <- syntax_parse_rhs(rhs = rhs.formula[[2L]], op = op)

    if(debug) print(out)

    # for each lhs element
    for(l in 1:length(lhs.names)) {
      # for each rhs element
      for(j in 1:length(out)) {
        # catch intercepts
        if(names(out)[j] == "intercept") {
          if(op == "~") {
            rhs.name <- ""
          } else {
            # either number (1), or reserved name?
            stop("penfa ERROR: right-hand side of formula contains an invalid variable name:\n    ", x)
          }
        } else if(names(out)[j] == "..zero.." && op == "~") {
          rhs.name <- ""
        } else if(names(out)[j] == "..constant.." && op == "~") {
          rhs.name <- ""
        } else {
          rhs.name <- names(out)[j]
        }

        # catch lhs = rhs and op = "=~"
        if(op == "=~" && lhs.names[l] == names(out)[j]) {
          stop("penfa ERROR: latent variable `", lhs.names[l], "' cannot be measured by itself")
        }

        # check if we not already have this combination (in this group)
        # 1. asymmetric (=~, ~, ~1) # e.g. factor loadings
        if(op != "~~") {
          idx <- which(FLAT.lhs == lhs.names[l] & FLAT.op == op & FLAT.block == BLOCK & FLAT.rhs == rhs.name)
          if(length(idx) > 0L) {
            stop("penfa ERROR: duplicate model element in: ", model[i])
          }
        } else {
          # 2. symmetric (~~) # e.g. factor variances
          idx <- which(FLAT.lhs == rhs.name & FLAT.op == "~~" & FLAT.block == BLOCK & FLAT.rhs == lhs.names[l])
          if(length(idx) > 0L) {
            stop("penfa ERROR: duplicate model element in: ", model[i])
          }
        }

        FLAT.idx <- FLAT.idx + 1L
        FLAT.lhs[FLAT.idx] <- lhs.names[l]
        FLAT.op[ FLAT.idx] <- op
        FLAT.rhs[FLAT.idx] <- rhs.name
        FLAT.block[FLAT.idx] <- BLOCK
        FLAT.fixed[FLAT.idx] <- ""
        FLAT.start[FLAT.idx] <- ""

                mod <- list()
        rhs.mod <- 0L
        if(length(out[[j]]$fixed) > 0L) {
          mod$fixed <- out[[j]]$fixed
          FLAT.fixed[FLAT.idx] <- paste(mod$fixed, collapse=";")
          rhs.mod <- 1L
        }
        if(length(out[[j]]$start) > 0L) {
          mod$start <- out[[j]]$start
          FLAT.start[FLAT.idx] <- paste(mod$start, collapse=";")
          rhs.mod <- 1L
        }
        if(op == "=~" && rhs == "0") {
          mod$fixed <- 0
          FLAT.rhs[FLAT.idx] <- FLAT.lhs[FLAT.idx]
          FLAT.fixed[FLAT.idx] <- paste(mod$fixed, collapse=";")
          rhs.mod <- 1L
        }
        FLAT.rhs.mod.idx[FLAT.idx] <- rhs.mod
        if(rhs.mod > 0L) {
          MOD.idx <- MOD.idx + 1L
          MOD[[MOD.idx]] <- mod
        }
      } # rhs elements
    } # lhs elements
  } # model elements

  # enumerate modifier indices
  mod.idx <- which(FLAT.rhs.mod.idx > 0L)
  FLAT.rhs.mod.idx[ mod.idx ] <- 1:length(mod.idx)

  FLAT <- list(lhs=FLAT.lhs, op=FLAT.op, rhs=FLAT.rhs, mod.idx=FLAT.rhs.mod.idx,
               block=FLAT.block, fixed=FLAT.fixed, start=FLAT.start #,
               #lower=FLAT.lower, upper=FLAT.upper, label=FLAT.label, prior=FLAT.prior
               )

  # change op for intercepts (for convenience only)
  int.idx <- which(FLAT$op == "~" & FLAT$rhs == "")
  if(length(int.idx) > 0L) {
    FLAT$op[int.idx] <- "~1"
  }

  # reorder covariances here
  FLAT <- partable_covariance_reorder(FLAT)
  if(as.data.frame.) {
    FLAT <- as.data.frame(FLAT, stringsAsFactors=FALSE)
  }
  attr(FLAT, "modifiers") <- MOD
  FLAT
}

# syntax_parse_rhs
syntax_parse_rhs <- function(rhs, op = "") {

  # fill in rhs list
  out <- list()
  repeat {
    if(length(rhs) == 1L) { # last one and only a single element
      out <- c(vector("list", 1L), out)
      NAME <- all.vars(rhs)
      if(length(NAME) > 0L) {
        names(out)[1L] <- NAME
      } else { # intercept or zero?
        if(as.character(rhs) == "1") {
          names(out)[1L] <- "intercept"
        } else if(as.character(rhs) == "0") {
          names(out)[1L] <- "..zero.."
          out[[1L]]$fixed <- 0
        } else {
          names(out)[1L] <- "..constant.."
          out[[1L]]$fixed <- 0
        }
      }
      break
    } else if(rhs[[1L]] == "*") { # last one, but with modifier
      out <- c(vector("list", 1L), out)
      NAME <- all.vars(rhs[[3L]])

      if(length(NAME) > 0L) { # not an intercept
        # catch interaction term
        rhs3.names <- all.names(rhs[[3L]])
        if(rhs3.names[1L] == ":") {
          NAME <- paste(NAME[1L], ":", NAME[2L], sep = "")
        }
        names(out)[1L] <- NAME
      } else { # intercept
        names(out)[1L] <- "intercept"
      }
      i.var <- all.vars(rhs[[2L]], unique=FALSE)
      if(length(i.var) > 0L) {
        # modifier are unquoted labels
        out[[1L]]$label <- i.var
      } else {
        # modifer is something else
        out[[1L]] <- syntax_get_modifier(rhs[[2L]])
      }
      break
    } else if(rhs[[1L]] == ":") { # last one, but interaction term
      out <- c(vector("list", 1L), out)
      NAME <- all.vars(rhs)
      NAME <- paste(NAME[1L], ":", NAME[2L], sep = "")
      names(out)[1L] <- NAME
      break
    } else if(rhs[[1L]] == "+") { # not last one!

      # three possibilities:
      # 1. length(rhs[[3]] == 3), and rhs[[3L]][[1]] == "*" -> modifier
      # 2. length(rhs[[3]] == 3), and rhs[[3L]][[1]] == ":" -> interaction
      # 3. length(rhs[[3]] == 1) -> single element
      out <- c(vector("list", 1L), out)

      # modifier or not?
      if(length(rhs[[3L]]) == 3L && rhs[[3L]][[1]] == "*") {
        # modifier!!
        NAME <- all.vars(rhs[[3L]][[3]])

        if(length(NAME) > 0L) { # not an intercept
          # catch interaction term
          rhs3.names <- all.names(rhs[[3L]][[3]])
          if(rhs3.names[1L] == ":") {
            NAME <- paste(NAME[1L], ":", NAME[2L], sep = "")
          }
          names(out)[1L] <- NAME
        } else { # intercept
          names(out)[1L] <- "intercept"
        }
        i.var <- all.vars(rhs[[3]][[2L]], unique = FALSE)
        if(length(i.var) > 0L) {
          # modifier are unquoted labels
          out[[1L]]$label <- i.var
        } else {
          # modifer is something else
          out[[1L]] <- syntax_get_modifier(rhs[[3]][[2L]])
        }

        # interaction term?
      } else if(length(rhs[[3L]]) == 3L && rhs[[3L]][[1]] == ":") {
        # interaction term, without modifier
        NAME <- all.vars(rhs[[3L]])
        NAME <- paste(NAME[1L], ":", NAME[2L], sep = "")
        names(out)[1L] <- NAME

      } else { # no modifier!!
        NAME <- all.vars(rhs[[3]])
        if(length(NAME) > 0L) {
          names(out)[1L] <- NAME
        } else { # intercept or zero?
          if(as.character(rhs[[3]]) == "1") {
            names(out)[1L] <- "intercept"
          } else if(as.character(rhs[[3]]) == "0") {
            names(out)[1L] <- "..zero.."
            out[[1L]]$fixed <- 0
          } else {
            names(out)[1L] <- "..constant.."
            out[[1L]]$fixed <- 0
          }
        }
      }

      # next element
      rhs <- rhs[[2L]]
    } else {
      stop("penfa ERROR: problems when parsing this line: ", rhs, "\n")
    }
  }

  # if multiple elements, check for duplicated elements and merge if found
  if(length(out) > 1L) {
    rhs.names <- names(out)
    while( !is.na(idx <- which(duplicated(rhs.names))[1L]) ) {
      dup.name <- rhs.names[ idx ]
      orig.idx <- match(dup.name, rhs.names)
      merged <- c( out[[orig.idx]], out[[idx]] )
      if(!is.null(merged)) # be careful, NULL will delete element
        out[[orig.idx]] <- merged
      out <- out[-idx]
      rhs.names <- names(out)
    }
  }
  out
}

# syntax_get_modifier
syntax_get_modifier <- function(mod) {
  if(length(mod) == 1L) {
    # three possibilites: 1) numeric, 2) NA, or 3) quoted character
    if( is.numeric(mod) )
      return( list(fixed=mod) )
    if( is.na(mod) )
      return( list(fixed=as.numeric(NA)) )
    if( is.character(mod) )
      return( list(label=mod) )
  } else if(mod[[1L]] == "start") {
    cof <- unlist(lapply(as.list(mod)[-1], eval, envir=NULL, enclos=NULL))
    return( list(start=cof) )
  } else if(mod[[1L]] == "c") {
    # vector: we allow numeric and character only
    cof <- unlist(lapply(as.list(mod)[-1], eval, envir=NULL, enclos=NULL))
    if(all(is.na(cof))) {
      return( list(fixed=rep(as.numeric(NA), length(cof))) )
    } else if(is.numeric(cof))
      return( list(fixed=cof) )
    else if(is.character(cof)) {
      cof[is.na(cof)] <- "" # catch 'NA' elements in a label
      return( list(label=cof) )
    } else {
      stop("penfa ERROR: cannot parse modifier:", mod, "\n")
    }
  } else {
    # unknown expression
    # as a final attempt, we will evaluate it and coerce it
    # to either a numeric or character (vector)
    cof <- try( eval(mod, envir=NULL, enclos=NULL), silent=TRUE)
    if(inherits(cof, "try-error")) {
      stop("penfa ERROR: evaluating modifier failed: ",
           paste(as.character(mod)[[1]], "()", sep = ""), "\n")
    } else if(is.numeric(cof)) {
      return( list(fixed=cof) )
    } else if(is.character(cof)) {
      return( list(label=cof) )
    } else {
      stop("penfa ERROR: cannot parse modifier: ", paste(as.character(mod)[[1]], "()", sep = ""), "\n")
    }
  }
}

