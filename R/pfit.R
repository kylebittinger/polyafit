#' S3 interface to polyafit
#'
#' @param x A matrix or data frame with one row per observation and column per
#'   feature
#' @param min_counts Minimum number of observed counts for a feature to be
#'   included in the fit
#' @param ... Additional arguments passed to other methods
#' @export
pfit <- function (x, min_counts = 5, ...) {
  UseMethod("pfit")
}

#' @rdname pfit
#' @examples
#' x <- matrix(c(56, 42, 122, 100, 8, 15, 21, 14, 82, 98), nrow = 2)
#' p <- pfit(x)
#' p$par
#' p$theta
#' @export
pfit.matrix <- function (x, min_counts = 5) {
  is_included <- colSums(x) >= min_counts
  fitobj <- optim_polya(x[,is_included])
  params <- fitobj$par
  theta <- 1 / sum(params)
  obj <- list(
    data = x,
    min_counts = min_counts,
    is_included = is_included,
    parameters = params,
    theta = theta)
  class(obj) <- "pfit"
  obj
}

#' @rdname pfit
#' @param rownames_in The column of the data frame that gives the rownames
#'   for the data. Can be specified as an integer index or with the column
#'   name in quotes. Set to \code{NULL} to use the existing rownames of the
#'   data frame.
#' @export
pfit.data.frame <- function (x, min_counts = 5, rownames_in = 1) {
  x <- make_rownames(x, rownames_in)
  x <- as.matrix(x)
  pfit.matrix(x, min_counts)
}

# Convert a column into rownames
make_rownames <- function (df, rownames_in = 1) {
  if (!is.null(rownames_in)) {
    if (length(rownames_in) != 1) {
      stop("rownames_in must have length 1")
    }
    if (is.character(rownames_in)) {
      rownames_in <- match(rownames_in, colnames(df))
    }
    df <- as.data.frame(df)
    df_rownames <- df[,rownames_in]
    df <- df[,-rownames_in, drop = FALSE]
    rownames(df) <- df_rownames
  }
  df
}

#' Test features for enrichment
#' @param p A \code{pfit} object
#' @return A data frame with three columns:
#'   \item{observation_idx}{The index of the observation, i.e. the row in the
#'     input data.}
#'   \item{feature_idx}{The index of the feature, i.e. the column in the input
#'     data.}
#'   \item{p_value}{The p-value for feature enrichment against a null
#'     hypothesis of all features arising from a single Dirichlet-Multinomial
#'     distribution.}
#' @export
pfit_enrichment <- function (p) {
  get_pval <- function (idx) {
    x <- p$data[idx,p$is_included]
    cdf_val <- ppolya_marginal(x, p$parameters, log.p = FALSE)
    1 - cdf_val
  }
  feature_idx = which(p$is_included)
  names(feature_idx) <- NULL
  tibble::tibble(observation_idx = seq_len(nrow(p$data))) %>%
    dplyr::group_by(observation_idx) %>%
    dplyr::summarise(
      feature_idx = feature_idx,
      p_value = get_pval(observation_idx),
      .groups = "drop")
}
