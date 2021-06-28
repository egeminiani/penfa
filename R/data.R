####################
# ---- Data set ----
####################

#' Data set for cross-cultural analysis
#'
#' @description A data set for cross-cultural analysis containing the
#'   standardized ratings to 12 items concerning organizational citizenship
#'   behavior. Employees from different countries were asked to rate their
#'   attitudes towards helping other employees and giving suggestions for
#'   improved work conditions. The items are thought to measure two latent
#'   factors: **helping** behavior (first seven items) and **voice** behavior
#'   (last five items). See below for details.
#'
#' @details The original data come from the
#'   [`ccpsyc`](https://github.com/Jo-Karl/ccpsyc/) package. For convenience,
#'   the following pre-processing has been applied:
#'   * The data were filtered to only include employees from Lebanon and Taiwan.
#'   * The answers, originally on a 7-point Likert scale, were standardized.
#'   * The items were renamed as described above.
#'
#' @format A data frame with 767 rows and 13 variables:
#'
#' \describe{
#'   \item{country}{Character. Country of origin of the employee: Lebanon
#'   ("LEB") or Taiwan ("TAIW").}
#'   \item{h1}{Numeric. Standardized ratings to the item "*I volunteer to do
#'   things for this organization*".}
#'   \item{h2}{Numeric.  Standardized ratings to the item "*I help orient new
#'   employees in this organization*".}
#'   \item{h3}{Numeric. Standardized ratings to the item "*I attend functions
#'   that help this organization*".}
#'   \item{h4}{Numeric. Standardized ratings to the item "*I help others in
#'   this work group with their work for the benefit of the group*".}
#'   \item{h5}{Numeric. Standardized ratings to the item "*I get involved to
#'   benefit this organization*".}
#'   \item{h6}{Numeric. Standardized ratings to the item "*I help others in this
#'   organization learn about the work*".}
#'   \item{h7}{Numeric. Standardized ratings to the item "*I help others in this
#'   organization with their work responsibilities*".}
#'   \item{v1}{Numeric. Standardized ratings to the item "*I develop and make
#'   recommendations concerning issues that affect this organization*".}
#'   \item{v2}{Numeric. Standardized ratings to the item "*I speak up
#'   and encourage other in this organization to get involved in issues that
#'   affect the group*".}
#'   \item{v3}{Numeric. Standardized ratings to the item "*I
#'   communicate my opinions about work issues to others in this organization
#'   even if my opinion is different and others in the organization disagree
#'   with me*".}
#'   \item{v4}{Numeric. Standardized ratings to the item "*I keep well informed
#'   about issues where my opinion might be useful to this organization*".}
#'   \item{v5}{Numeric. Standardized ratings to the item "*I speak up in
#'   this organization with ideas for new projects or changes in procedures*".}
#'   }
#'
#'
#'
#' @source
#'
#' The original data set is available from the
#' [`ccpsyc`](https://github.com/Jo-Karl/ccpsyc/tree/master/data/) package.
#' Please refer to [Fischer and Karl
#' (2019)](https://www.frontiersin.org/articles/10.3389/fpsyg.2019.01507/full)
#' and [Fischer et al.
#' (2019)](https://link.springer.com/article/10.1057/s41267-017-0132-6) for a
#' description and analysis of these data.
"ccdata"
