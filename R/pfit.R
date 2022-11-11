#' S3 interface to polyafit
#'
#' @param x A matrix or data frame with one row per observation and column per
#'   feature
#' @param min_counts Minimum number of observed counts for a feature to be
#'   included in the fit
#' @param ... Additional arguments passed to other methods
#' @return An object of class "pfit", with the following attributes:
#'   \item{data}{The data matrix. Row and column names will be created
#'     automatically if the original data matrix lacked them.}
#'   \item{is_included}{Boolean vector indicating if the feature was included
#'     in the model.}
#'   \item{params}{Best-fit model parameters.}
#'   \item{theta}{The overdispersion parameter, \code{1 / sum(params)}.}
#' @export
pfit <- function (x, min_counts = 5, ...) {
  UseMethod("pfit")
}

#' @rdname pfit
#' @examples
#' x <- matrix(c(56, 42, 122, 100, 8, 15, 21, 14, 82, 98), nrow = 2)
#' p <- pfit(x)
#' p$params
#' p$theta
#' @export
pfit.matrix <- function (x, min_counts = 5) {
  if (is.null(rownames(x))) {
    rownames(x) <- paste0("observation", seq_len(nrow(x)))
  }
  if (is.null(colnames(x))) {
    colnames(x) <- paste0("feature", seq_len(ncol(x)))
  }
  is_included <- colSums(x) >= min_counts
  fitobj <- optim_polya(x[,is_included])
  params <- fitobj$par
  theta <- 1 / sum(params)
  obj <- list(
    data = x,
    is_included = is_included,
    params = params,
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
#'
#' Tests each feature for enrichment. The null hypothesis is that all the
#' features across all observations arise from a single Dirichlet-Multinomial
#' distribution.
#'
#' @param p A \code{pfit} object
#' @return A data frame with three columns, giving the observation, the
#'   feature, and the p-value for enrichment.
#' @export
feature_enrichment <- function (p) {
  p$data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("observation") %>%
    tidyr::pivot_longer(
      cols = -observation, names_to = "feature", values_to = "counts") %>%
    dplyr::group_by(observation) %>%
    dplyr::filter(p$is_included) %>%
    dplyr::mutate(
      p.value = 1 - ppolya_marginal(counts, p$params, log.p = FALSE)) %>%
    dplyr::ungroup() %>%
    dplyr::select(observation, feature, p.value)
}

#' Plot a \code{pfit} object
#' @param p A \code{pfit} object
#' @return A ggplot object
#' @export
plot.pfit <- function (p) {
  ref_observation <- rownames(p$data)[1]
  props <- sweep(p$data, 1, rowSums(p$data), "/")
  props <- props[,p$is_included]
  rownames(props)[1] <- "ref_prop"
  lower_limit <- min(props) / 2
  props %>%
    t() %>% # each row is a feature
    as.data.frame() %>%
    tibble::rownames_to_column(var = "feature") %>%
    tidyr::pivot_longer(
      cols = -c(feature, ref_prop),
      names_to = "observation", values_to = "prop") %>%
    ggplot2::ggplot() +
    ggplot2::geom_abline(
      slope = 1, intercept = 0, linetype = "dashed", color = "#333333") +
    ggplot2::geom_point(ggplot2::aes(x = ref_prop, y = prop)) +
    ggplot2::facet_wrap(~ observation) +
    ggplot2::scale_x_log10(limits = c(lower_limit, 1)) +
    ggplot2::scale_y_log10(limits = c(lower_limit, 1)) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      x = paste("Proportion in", ref_observation), y = "Proportion in other")
}
