#' Arithmetic mean estimate with control variates
#'
#' Computes the Monte Carlo estimate of the arithmetic mean using
#' first order zero-variance control variates (Mira et al, 2013) in the
#' hopes of obtaining reduced-variance estimates compared to `mean()`.
#' The control variates in this setting are the `gradients`.
#' A sufficient condition for the control variates to have zero expectation
#' is that the distribution for the transformed variables has tails that decay
#' faster than polynomially.
#'
#' @param x (rvar) An [`rvar`].
#' @param gradients the gradients of the natural logarithm of the probability
#' distribution function for the transformed variables with respect to the
#' transformed variables, as used for gradient-based methods
#'
#' @return The estimate of `mean(x)` using control variates for variance reduction.
#'
#' @template ref-mira-cv-2013
#'
#' @examples
#' x <- as_draws_df(example_draws())
#' x$.gradient_mu <- rnorm(ndraws(x))
#' x$.gradient_sigma <- rnorm(ndraws(x))
#'
#' grads <- gradients(x)
#'
#' mean_control_variates(x$mu, gradients = grads)
#' summarise_draws(x, mean, mean_control_variates, .args = list(gradients = grads))
#'
#' @export
mean_control_variates <- function(x, gradients = NULL, ...) {
  # TODO: make generic
  # TODO: we can add an efficient 'matrix' method as per github
  # but we need this default method to reliably work with summarize_draws
  if (is.null(gradients)) {
    return(mean(x, ...))
  }
  x <- as.vector(x)
  # Getting optimal coefficients for control variates
  coefs <- as.matrix(lm(x ~ gradients)$coefficients[-1])
  # Obtaining the estimate with control variates
  mean(x - gradients %*% coefs)
}


#' Standard deviation estimate with control variates
#'
#' Computes the Monte Carlo estimate of the standard deviation using
#' first order zero-variance control variates (Mira et al, 2013) in the
#' hopes of obtaining reduced-variance estimates compared to `sd()`.
#' Estimates are based on estimating `E[X^2]` and `E[x]` separately
#' with first order zero-variance control variates.
#' The control variates in this setting are the `gradients`.
#' A sufficient condition for the control variates to have zero expectation
#' is that the distribution for the transformed variables has tails that decay
#' faster than polynomially.
#'
#' @param x (rvar) An [`rvar`].
#' @param gradients the gradients of the natural logarithm of the probability
#' distribution function for the transformed variables with respect to the
#' transformed variables, as used for gradient-based methods
#'
#' @return The estimate of `sd(x)` using control variates for variance reduction.
#'
#' @template ref-mira-cv-2013
#'
#' @examples
#' x <- as_draws_df(example_draws())
#' x$.gradient_mu <- rnorm(ndraws(x))
#' x$.gradient_sigma <- rnorm(ndraws(x))
#'
#' grads <- gradients(x)
#'
#' sd_control_variates(x$mu, gradients = grads)
#' summarise_draws(x, sd, sd_control_variates, .args = list(gradients = grads))
#'
#' @export
sd_control_variates <- function(x, gradients = NULL, ...) {
  # TODO: make generic
  # TODO: we can add an efficient 'matrix' method as per github
  # but we need this default method to reliably work with summarize_draws
  if (is.null(gradients)) {
    return(var(x, ...))
  }
  x <- as.vector(x)

  # Getting optimal coefficients for control variates (E[x])
  coefs1 <- as.matrix(lm(x ~ gradients)$coefficients[-1])
  # Obtaining the estimate with control variates (E[x])
  E_x <- mean(x - gradients %*% coefs1)

  # Getting optimal coefficients for control variates (E[x^2])
  coefs2 <- as.matrix(lm(x^2 ~ gradients)$coefficients[-1])
  # Obtaining the estimate with control variates (E[x^2])
  E_x_squared <- mean(x^2 - gradients %*% coefs2)

  sqrt(E_x_squared - E_x^2)
}
