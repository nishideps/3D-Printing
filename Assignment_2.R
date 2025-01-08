#' Nisha Depala, s2149899
#' Add your own function definitions on this file.

#' neg_log_lik
#
#' @description Evaluate the negated log-likelihood for model A and B
#' @param beta A vector with the beta parameters
#' @param data A `data.frame` with the same variables as the `filament1` data set.
#' Must have columns `CAD_Weight` and `Actual_Weight`
#' @param model Either "A" for a log-linear variance model, or "B" for a proportional
#' scaling error model

neg_log_lik <- function(beta, data, model){
  
  mu <- beta[1] + beta[2]*data[["CAD_Weight"]]
  
  # distinguish between the two models to find the particular standard deviation for the betas
  if(model == "A") {
    sigma <- sqrt(exp(beta[3] + beta[4]*data[["CAD_Weight"]]))
  }else{
    sigma <- sqrt(exp(beta[3])+exp(beta[4]) * (data[["CAD_Weight"]]^2))
  }
  - sum(dnorm(data[["Actual_Weight"]],
              mean = mu,
              sd=sigma,
              log = TRUE))
  
}

#' filament_estimate
#
#' @description Estimate filament models with different variance structure
#' @param data A `data.frame` with the same variables as the `filament1` data set.
#' Must have columns `CAD_Weight` and `Actual_Weight`
#' @param model Either "A" for a log-linear variance model, or "B" for a proportional
#' scaling error model
#' @return An estimation object suitable for use with [filament1_predict()]

filament1_estimate <- function(data, model) {
  model <- match.arg(model, c("A", "B"))
  if (model == "A") {
    beta_start <- c(-0.1, 1.07, -2, 0.05)
  } else {
    beta_start <- c(-0.15, 1.07, -13.5, -6.5)
  }
  opt <- optim(beta_start,
               neg_log_lik,
               data = data,
               model = model,
               hessian = TRUE,
               method = "Nelder-Mead",
               control = list(maxit = 5000)
  )
  fit <- list(
    model = model,
    par = opt$par,
    hessian = opt$hessian
  )
  class(fit) <- c("filament1_estimate", "list")
  fit
}

#' filament1_aux_EV
#' 
#' @description Evaluate the expectation and variance for model A and B
#' @param beta A vector with the beta parameters
#' @param data A `data.frame` containing the required predictors, including `CAD_Weight`
#' @param model Either "A" for a log-linear variance model, or "B" for a proportional
#' scaling error model
#' @param Sigma_beta : If not NULL, an estimate of the covariance matrix for
#                 the uncertainty of estimated betas
#' @return A list with four elements:
#     E : E(y|beta,x)
#     V : Var(y|beta,x)
#     VE : Var(E(y|beta,x)|x) or NULL
#     EV : E(Var(y|beta,x)|x) or NULL

filament1_aux_EV <- function(beta, data, model = c("A", "B"),
                             Sigma_beta = NULL) {
  
  model <- match.arg(model)
  if (model == "A") {
    
    ZE.0 <- model.matrix( ~ 1 + CAD_Weight, data = data)
    ZV.0 <- model.matrix( ~ 1 + CAD_Weight, data = data)
    ZE = cbind(ZE.0, ZV.0 * 0) 
    ZV = cbind(ZE.0 * 0, ZV.0)
    
    VE <- EV <- NULL
    if (!is.null(Sigma_beta)) {
      # E(Var(y|beta,x)|x)
      EV <- exp(ZV %*% beta + rowSums(ZV * (ZV %*% Sigma_beta)) / 2)
      # Var(E(y|beta,x)|x)
      VE <- rowSums(ZE * (ZE %*% Sigma_beta))
    }
    out <- list(
      E = ZE %*% beta,
      V = exp(ZV %*% beta),
      VE = VE,
      EV = EV
    )
  } else {
    
    ZE.0 <- model.matrix( ~ 1 + CAD_Weight, data = data)
    ZV.0 <- model.matrix( ~ 1 + I(CAD_Weight^2), data = data)
    ZE = cbind(ZE.0, ZV.0 * 0) 
    ZV = cbind(ZE.0 * 0, ZV.0)
    
    VE <- EV <- NULL
    if (!is.null(Sigma_beta)) {
      # E(Var(y|beta,x)|x)
      # (pmin: Ignore large Sigma_beta values)
      EV <- ZV %*% exp(beta + pmin(0.5^2, diag(Sigma_beta)) / 2)
      # Var(E(y|beta,x)|x)
      VE <- rowSums(ZE * (ZE %*% Sigma_beta))
    }
    out <- list(
      E = ZE %*% beta,
      V = ZV %*% exp(beta),
      VE = VE,
      EV = EV
    )
  }
  out
}

#' filament1_predict
#'
#' @description Computes predictive distributions and 95% confidence intervals for a dataset.
#' @param new_data The new dataset we want to make predictions for.
#' @param filament_fit Estimated parameters for a given model using filament1_estimate.
#'
#' @return a data frame with 4 variables summarising the predictive distribution for each row of the new data:
#'    mean: E(y|beta,x)
#'    sd: sqrt(Var(y|beta,x))
#'    lwr: lower predictive interval
#'    upr: upper predictive interval
#' 
#' @seealso filament1_aux_EV,  filament1_estimate
#' 
#' @examples 
#' # Assuming 'new_data' is a data frame with predictor variables,
#' # and 'filament_fit' is the model fit object obtained from filament1_estimate
#' A_fit <- filament1_estimate(filament1, "A")
#' pred_A <- filament1_predict(filament1, A_fit)

filament1_predict <- function(new_data, filament_fit){
  # extract model and parameters from filament_fit
  model <- filament_fit$model
  beta <- filament_fit$par
  sig_b <- solve(filament_fit$hessian)
  
  # Calculate expected value and variance
  aux_ev <- filament1_aux_EV(beta, new_data, model, sig_b)
  
  # Extract predictive mean and standard deviation
  mean <- aux_ev$E
  sd <- sqrt(aux_ev$EV + aux_ev$VE)
  
  # Calculate z for 95% confidence interval
  z <- qnorm(0.975)
  
  # Calculate lower and upper predictive intervals
  lwr <- mean - z * sd
  upr <- mean + z * sd
  
  result <- data.frame(mean, sd, lwr, upr)
  
  return(result)
}

#' se
#' 
#' @description Computed the squared error between predicted and actual values
#' 
#' @param pred A numeric vector of predicted values.
#' @param act A numeric vector of actual values.
#' 
#' @return A numeric vector of squared errors between predicted and actual values.
#' 
#' @examples 
#' # Assuming 'pred' and 'act' are vectors
#' squared_error <- se(pred, act)
#' 
#' @export

se <- function(pred, act){
  se <- (pred - act)^2
  return(se)
}

#' ds
#' 
#' @description Computes the Dawid-Sebastiani score
#' 
#' @param pred_mean A numeric vector of predicted mean.
#' @param pred_sd A numeric vector of predicted standard deviation.
#' @param act_mean A numeric vector of actual values.
#' 
#' @return A numeric vector containing the Dawid-Sebastiani score.
#' 
#' @examples 
#' # Assuming 'pred_mean', 'pred_sd' and 'act_mean' are vectors
#' dawid_sebastiani <- (pred_mean, pred_sd, act_mean)
#' 
#' @export

ds <- function(pred_mean, pred_sd, act){
  ds <- (pred_mean - act)^2 / (pred_sd^2) + 2 * log(pred_sd)
  return(ds)
}

#' leave1out
#' 
#' @description Performs leave-one-out cross validation for a given model and dataset.
#'
#' @param data The dataset for which leave-one-out cross validation is performed.
#' @param model The model used for estimation and prediction, "A" or "B".
#'
#' @return A modified dataset with columns added for predictive mean, predictive standard deviation, squared error, and deviance statistic for each observation.
#'
#' @examples
#' # Assuming 'data' is a data frame with predictor and response variables, and 'model' 
#' result <- leave1out(data, model)
#' head(result)
#' 
#' @seealso filament1_estimate, filament1_predict.
#'
#' @export

leave1out <- function(data, model) {
  # Number of observations
  n <- nrow(data)
  
  # Add columns for predictive mean and standard deviation
  data <- data %>%
    mutate(mean = NA_real_, sd = NA_real_)
  
  # Execute leave1out cross validation
  for (i in 1:n) {
    # Exclude i-th observation in fit
    fit <- filament1_estimate(data[-i, , drop = FALSE], model)
    
    # Predict i-th observation
    pred <- filament1_predict(data[i, , drop = FALSE], fit)
    
    data[i, "mean"] <- pred$mean
    data[i, "sd"] <- pred$sd
  }
  
  # Calcuate squared error and Dawid-Sebastiani statistic
  data <- data %>%
    mutate(se = se(mean, data$Actual_Weight),
           ds = ds(mean, sd, data$Actual_Weight))
  
  return(data)
}

#' arch_loglike
#' 
#' @description Evaluates the combined log-likelihood for a collection of y observations.
#' 
#' @param y Vector of the number of left and right femurs found respectively.
#' @param N Number of people buried in total.
#' @param phi The probability of finding a femur.
#' 
#' @return The log-likelihood for each row-pair, (N, phi).
#' 
#' # Assuming 'y' and 'N' are vectors of observed and total counts respectively, and 'phi' is the parameter of the ARCH model
#' log_likelihood <- arch_loglike(y, N, phi)
#' 
#' @export

arch_loglike <- function(y, N, phi){
  # Compute sum of log-gamma functions
  y_1 <- sum(lgamma(y + 1))
  N_y_1 <- sum(lgamma(N-y+1))
  N_1 <- lgamma(N+1)
  
  log_likelihood <- - y_1 - N_y_1 + 2*N_1 + sum(y)*log(phi) + (2*N - sum(y)) * log(1-phi)
  
  return(log_likelihood)
}

#' estimate
#' 
#' @param y Vector of the number of left and right femurs found respectively.
#' @param xi Parameter of the geometric distribution representing the probability of success in the first trial.
#' @param a Shape parameter of the beta distribution for the model.
#' @param b Shape parameter of the beta distribution for the model.
#' @param K Number of Monte Carlo samples to generate.
#' 
#' @return A list containing the estimated parameters:
#'   - p_y: Estimated probability of observing the data y.
#'   - E_N_y: Estimated expected value of the total count N.
#'   - E_phi_y: Estimated expected value of the probability of success phi.
#' 
#' @example 
#' estimate(y=c(237,256), xi=1/1001, a=0.5, b=0.5, K=10000)
#' 
#' @seealso arch_loglike
#' 
#' @export

estimate <- function(y, xi, a, b, K) {
  # Generate Monte Carlo samples for N (total counts) and phi (probability of success)
  N_samples <- rgeom(K, xi)
  phi_samples <- rbeta(K, a, b)
  
  # Calculate likelihood for each sample
  likelihood <- numeric(K)
  for (i in 1:K) {
    likelihood[i] <- exp(arch_loglike(y, N_samples[i], phi_samples[i]))
  }
  
  # Estimate parameters based on Monte Carlo samples
  p_y <- 1/K * sum(likelihood)
  E_N_y <- 1/(K*p_y) * sum(N_samples * likelihood)
  E_phi_y <- 1/(K*p_y) * sum(phi_samples * likelihood)
  
  return(list(p_y = p_y, E_N_y = E_N_y, E_phi_y = E_phi_y))
}
