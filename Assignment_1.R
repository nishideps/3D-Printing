#' Log-Exponential density
#'
#' Compute the density or log-density for a Log-Exponential (LogExp)
#' distribution
#'
#' @param x vector of quantiles
#' @param rate vector of rates
#' @param log logical; if TRUE, the log-density is returned

dlogexp <- function(x, rate = 1, log = FALSE) {
  result <- log(rate) + x - rate * exp(x)
  if (!log) {
    exp(result)
  }
  result
}

#' Log-Sum-Exp
#'
#' Convenience function for computing log(sum(exp(x))) in a
#' numerically stable manner
#'
#' @param x numerical vector

log_sum_exp <- function(x) {
  max_x <- max(x, na.rm = TRUE)
  max_x + log(sum(exp(x - max_x)))
}


#' wquantile 
#'
#' Calculates empirical sample quantiles with optional weights, for given probabilities. 
#' Like in quantile(), the smallest observation corresponds to a probability of 0 and the largest to a probability of 1. 
#' Interpolation between discrete values is done when type=7, as in quantile(). 
#' Use type=1 to only generate quantile values from the raw input samples.
#'
#' @param x numeric vector whose sample quantiles are wanted
#' NA and NaN values are not allowed in numeric vectors unless na.rm is TRUE
#' @param probs numeric vector of probabilities with values in [0,1]
#' @param na.rm logical; if true, any NA and NaN's are removed from x before the quantiles are computed
#' @param type numeric, 1 for no interpolation, or 7, for interpolated quantiles. Default is 7
#' @param weights	 numeric vector of non-negative weights, the same length as x, or NULL. The weights are normalised to sum to 1. If NULL, then wquantile(x) behaves the same as quantile(x), with equal weight for each sample value

wquantile <- function (x, probs = seq(0, 1, 0.25), na.rm = FALSE, type = 7, 
                       weights = NULL, ...) 
{
  if (is.null(weights) || (length(weights) == 1)) {
    weights <- rep(1, length(x))
  }
  stopifnot(all(weights >= 0))
  stopifnot(length(weights) == length(x))
  if (length(x) == 1) {
    return(rep(x, length(probs)))
  }
  n <- length(x)
  q <- numeric(length(probs))
  reorder <- order(x)
  weights <- weights[reorder]
  x <- x[reorder]
  wecdf <- pmin(1, cumsum(weights)/sum(weights))
  if (type == 1) {
  }
  else {
    weights2 <- (weights[-n] + weights[-1])/2
    wecdf2 <- pmin(1, cumsum(weights2)/sum(weights2))
  }
  for (pr_idx in seq_along(probs)) {
    pr <- probs[pr_idx]
    if (pr <= 0) {
      q[pr_idx] <- x[1]
    }
    else if (pr >= 1) {
      q[pr_idx] <- x[n]
    }
    else {
      if (type == 1) {
        j <- 1 + pmax(0, pmin(n - 1, sum(wecdf <= pr)))
        q[pr_idx] <- x[j]
      }
      else {
        j <- 1 + pmax(0, pmin(n - 2, sum(wecdf2 <= pr)))
        g <- (pr - c(0, wecdf2)[j])/(wecdf2[j] - c(0, 
                                                   wecdf2)[j])
        q[pr_idx] <- (1 - g) * x[j] + g * x[j + 1]
      }
    }
  }
  q
}

#' Compute empirical weighted cumulative distribution
#'
#' Version of `ggplot2::stat_ecdf` that adds a `weights` property for each
#' observation, to produce an empirical weighted cumulative distribution function.
#' The empirical cumulative distribution function (ECDF) provides an alternative
#' visualisation of distribution. Compared to other visualisations that rely on
#' density (like [geom_histogram()]), the ECDF doesn't require any
#' tuning parameters and handles both continuous and discrete variables.
#' The downside is that it requires more training to accurately interpret,
#' and the underlying visual tasks are somewhat more challenging.
#'
# @inheritParams layer
# @inheritParams geom_point
#' @param na.rm If `FALSE` (the default), removes missing values with
#'    a warning.  If `TRUE` silently removes missing values.
#' @param n if NULL, do not interpolate. If not NULL, this is the number
#'   of points to interpolate with.
#' @param pad If `TRUE`, pad the ecdf with additional points (-Inf, 0)
#'   and (Inf, 1)
#' @section Computed variables:
#' \describe{
#'   \item{x}{x in data}
#'   \item{y}{cumulative density corresponding x}
#' }
#' @seealso wquantile
#' @export
#' @examples
#' library(ggplot2)
#'
#' n <- 100
#' df <- data.frame(
#'   x = c(rnorm(n, 0, 10), rnorm(n, 0, 10)),
#'   g = gl(2, n),
#'   w = c(rep(1/n, n), sort(runif(n))^sqrt(n))
#' )
#' ggplot(df, aes(x, weights = w)) + stat_ewcdf(geom = "step")
#'
#' # Don't go to positive/negative infinity
#' ggplot(df, aes(x, weights = w)) + stat_ewcdf(geom = "step", pad = FALSE)
#'
#' # Multiple ECDFs
#' ggplot(df, aes(x, colour = g, weights = w)) + stat_ewcdf()
#' ggplot(df, aes(x, colour = g, weights = w)) +
#'   stat_ewcdf() +
#'   facet_wrap(vars(g), ncol = 1)

stat_ewcdf <- function(mapping = NULL, data = NULL,
                       geom = "step", position = "identity",
                       ...,
                       n = NULL,
                       pad = TRUE,
                       na.rm = FALSE,
                       show.legend = NA,
                       inherit.aes = TRUE) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = StatEwcdf,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      n = n,
      pad = pad,
      na.rm = na.rm,
      ...
    )
  )
}


#' @title StatEwcdf ggproto object
#' @name StatEwcdf
#' @rdname StatEwcdf
#' @aliases StatEwcdf
#' @format NULL
#' @usage NULL
#' @export
#' @importFrom ggplot2 aes after_stat has_flipped_aes Stat
NULL

StatEwcdf <- ggplot2::ggproto(
  "StatEwcdf", ggplot2::Stat,
  required_aes = c("x|y", "weights"),
  dropped_aes = c("weights"),     
  
  default_aes = ggplot2::aes(y = ggplot2::after_stat(y)),
  
  setup_params = function(data, params) {
    params$flipped_aes <-
      ggplot2::has_flipped_aes(data,
                               params,
                               main_is_orthogonal = FALSE,
                               main_is_continuous = TRUE)
    
    has_x <- !(is.null(data$x) && is.null(params$x))
    has_y <- !(is.null(data$y) && is.null(params$y))
    if (!has_x && !has_y) {
      rlang::abort("stat_ewcdf() requires an x or y aesthetic.")
    }
    has_weights <- !(is.null(data$weights) && is.null(params$weights))
    #    if (!has_weights) {
    #      rlang::abort("stat_ewcdf() requires a weights aesthetic.")
    #    }
    
    params
  },
  
  compute_group = function(data, scales, n = NULL, pad = TRUE, flipped_aes = FALSE) {
    data <- flip_data(data, flipped_aes)
    # If n is NULL, use raw values; otherwise interpolate
    if (is.null(n)) {
      x <- unique(data$x)
    } else {
      x <- seq(min(data$x), max(data$x), length.out = n)
    }
    
    if (pad) {
      x <- c(-Inf, x, Inf)
    }
    if (is.null(data$weights)) {
      data_ecdf <- ecdf(data$x)(x)
    } else {
      data_ecdf <-
        spatstat.geom::ewcdf(
          data$x,
          weights = data$weights / sum(abs(data$weights)) 
        )(x)
    }
    
    df_ecdf <- vctrs::new_data_frame(list(x = x, y = data_ecdf), n = length(x))
    df_ecdf$flipped_aes <- flipped_aes
    ggplot2::flip_data(df_ecdf, flipped_aes)
  }
)

### Classical Estimation

#' neg_log_like
#' 
#' Compute the negative log likelihood for a specified model, either A or B that
#' follows:
#' Model A: yi ∼ Normal(β1 + β2xi, exp(β3 + β4xi))
#' Model B: yi ∼ Normal(β1 + β2xi, exp(β3) + exp(β4)xi^2))
#' 
#' @param beta a numeric vector of model parameters
#' @param data a data.frame containing the required variables
#' @param model either A or B
#' @returns the negated log-likelihood for the specified model
 
neg_log_like <- function(beta, data, model) {
  
  xi = data$CAD_Weight
  yi = data$Actual_Weight
  
  if (model == "A"){
    mu = beta[1] + beta[2]*xi
    var <- exp(beta[3] + beta[4]*xi)
  } else if (model == "B") {
    mu <- beta[1] + beta[2]*xi
    var <- exp(beta[3]) + exp(beta[4])*xi^2
  }
  
  neg_llike <- -sum(dnorm(x = yi, mean = mu, sd = sqrt(var), log = TRUE))
  
  return(neg_llike)
}

#' filament_1_estimate
#' 
#' @param data a data frame with columns CAD_Weight and Actual_Weight
#' @param model the model choice (either A or B)
#' @returns a data frame the best set of parameters, their confidence intervals and and the Hessian

filament1_estimate <- function(data, model){
  
  if (model == "A"){
    betaA <- c(-0.1, 1.07, -2, 0.05) # initial values for (β1, β2, β3, β4)
    opt <- optim(par = betaA, fn = neg_log_like, data = data, model = model, hessian = TRUE)
  } else if (model == "B") {
    betaB <- c(-0.15, 1.07, -13.5, -6.5)  # initial values for (β1, β2, β3, β4)
    opt <- optim(par = betaB, fn = neg_log_like, data = data, model = model, hessian = TRUE)
  }
 
  parameters <- opt$par
  hessian <- opt$hessian
  df <- list("parameters" = parameters, "hessian" = hessian)
  
  z_alpha <- qnorm(1 - 0.1 / 2)
  se <- sqrt(diag(solve(hessian)))
  Parameter <- parameters
  Actual <- parameters
  Lower <- parameters - z_alpha * se
  Upper  <- parameters + z_alpha * se
  
  return(data.frame(Actual, Lower, Upper))
}


### Bayesian Estimation


#' log_prior_density
#' 
#' @param theta a numeric vector equal to θ parameters
#' @param params a numeric vector of γ parameters
#' 
#' @returns the logarithm of the joint prior density p(θ) for the four θi parameters
#' 
#' @examples 
#' theta <- c(0, 0, 0, 0)
#' params <- c(1, 1, 1, 1)
#' log_prior_density(theta, params)

log_prior_density <- function(theta, params) {
  
  prior_t1 <- dnorm(theta[1], mean = 0, sd = sqrt(params[1]), log = TRUE)
  prior_t2 <- dnorm(theta[2], mean = 1, sd = sqrt(params[2]), log = TRUE)
  prior_t3 <- dlogexp(theta[3], rate = params[3], log = TRUE)
  prior_t4 <- dlogexp(theta[4], rate = params[4], log = TRUE)
  
  return(prior_t1 + prior_t2 + prior_t3 + prior_t4)
}

#' log_like
#' 
#' @param theta a numeric vector equal to θ parameters
#' @param x a numeric vector of x values as defined in the model
#' @param y a numeric vector of y values as defined in the model
#' 
#' @returns the log likelihood for the model defined above
#' 
#' @example 
#' x <- seq(100)
#' y <- seq(100)
#' log_like(c(0,0,0,0), x, y)

log_like <- function(theta, x, y) {
  
  mu = theta[1] + theta[2]*x
  var = exp(theta[3]) + exp(theta[4])*x^2
  
  llike = sum(dnorm(y, mean = mu, sd = sqrt(var), log = TRUE))
  
  return(llike)
}

#' log_posterior_density
#' 
#' This function calculates the logarithm of the posterior density p(θ|y) for the Bayesian model, 
#' where θ represents the parameter vector, based on observed data x and y, and given γ parameters.
#' 
#' @param theta a numeric vector equal to θ parameters
#' @param x a numeric vector of predictor x values as defined in the model
#' @param y a numeric vector of observed y values as defined in the model
#' @param params a numeric vector of γ parameters
#' 
#' @returns the logarithm of the posterior density p(θ|y), apart from some 
#' unevaluated normalisation constant
#' 
#' @seealso log_prior_density, log_like
#' 
#' @examples
#' # Compute log posterior density for given parameters and data
#' log_posterior_density(c(0.1, 0.2, -1.5, 0.5), x_data, y_data, c(1, 1, 1, 1))

log_posterior_density <- function(theta, x, y, params) {
  
  prior = log_prior_density(theta, params)
  llike = log_like(theta, x, y)
  
  return(llike + prior)
}

#' posterior_mode
#' 
#' @param theta_start initial values of (θ1, θ2, θ3, θ4)
#' @param x a numeric vector of predictor x values as defined in the model
#' @param y a numeric vector of observed y values as defined in the model
#' @param params a numeric vector of γ parameters
#' 
#' @returns a list containing the posterior mode location, the Hessian of the
#' log density at the mode and the inverse of the negated Hessian at the mode
#' 
#' @seealso log_posterior_density

posterior_mode <- function(theta_start, x, y, params) {
  
  opt <- optim(theta_start, log_posterior_density, x = x, y = y, params = params, method = "BFGS", control = list(fnscale = -1), hessian = TRUE)
  
  parameters <- opt$par
  hessian <- opt$hessian
  
  return(list("mode" = parameters, "hessian" = hessian, "S" = solve(-hessian)))
}

#' do_importance
#' 
#' This function generates samples from a multivariate normal distribution and 
#' computes the normalised log-importance weights based on the parameters and 
#' data. The function outputs a data frame containing the samples and their 
#' corresponding normalized log-importance-weights.
#' 
#' @param N the number of samples to generate
#' @param mu a numeric vector of the mean for the importance distribution.
#' @param S the covariance matrix for the importance distribution.
#' @param x a numeric vector of predictor x values as defined in the model
#' @param y a numeric vector of observed y values as defined in the model
#' @param params a numeric vector of γ parameters
#'
#' @returns a data frame with the beta values and the normalised log importance 
#' weights
#' 
#' @seealso log_posterior_density, log_sum_exp
#' 
#' @example 
#' # Generate importance samples with N = 10000
#' do_importance(10000, c(0, 0, 0, 0), diag(4), x_data, y_data, c(1, 1, 1, 1))

do_importance <- function(N, mu, S, x, y, params) {
  
  samples <- rmvnorm(n = N, mean = mu, sigma = S)
  
  log_weights <- log_posterior_density(samples, x, y, params) - dmvnorm(samples, mean = mu, sigma = S, log = TRUE)

  normalise_log_weights <- log_weights - log_sum_exp(log_weights)
  
  Beta1 = samples[, 1]
  Beta2 = samples[, 2]
  Beta3 = exp(samples[, 3])
  Beta4 = exp(samples[, 4])
  
  B_i <- cbind(Beta1, Beta2, Beta3, Beta4, log_weights = normalise_log_weights)
  
  return(data.frame(B_i))
}

#' make_CI
#' 
#' This function creates confidencene intervals given predictor x values and
#' a confidence interval 
#' 
#' @param x a numeric vector of predictor x values as defined in the model
#' @param weights a numeric vector of weights from the samples
#' @param prob number that represents the credible interval we are looking for
#' 
#' @returns a data frame containing the upper and lower bounds of the credible interval
#' 
#' @seealso wquantile
#' 
#' @examples 
#' make_CI(Value, exp(log_weights), 0.90)

make_CI <- function(x, weights, prob) {
  
  quantiles <- wquantile(x, probs = c((1 - prob) / 2, 1 - (1 - prob) / 2), weights = weights)
  CI <- data.frame(Lower = quantiles[[1]], Upper = quantiles[[2]])
  
  return(CI)
}
