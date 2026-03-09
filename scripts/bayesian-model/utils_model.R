# Function to convert mean and CV to Gamma parameters
get_gamma_param_from_mean_CV <- function(mu, cv){
  alpha <- 1 / (cv * cv)
  beta <- alpha / mu
  return(c('alpha' = alpha, 
           'beta' = beta))
}

get_proba_beta_prime_above_x <- function(x, alpha, beta){
  return(
    1. -  Rbeta(x, alpha, beta, lower=TRUE, log = F)
  )
}

# Function that creates the model and functions to extract information about posterior
create_model <- function(data, prior_param){
  
  # Prior distribution for dN
  get_prior_density_dN <- function(){
    f <- function(x){
      return(
        dgamma(x, shape = prior_param$alpha_N, rate = prior_param$beta_N)
      )
    }
    return(f)
  }
  # Prior distribution for dS
  get_prior_density_dS <- function(){
    f <- function(x){
      return(
        dgamma(x, shape = prior_param$alpha_S, rate = prior_param$beta_S)
      )
    }
    return(f)
  }
  # Prior distribution for omega
  get_log_prior_density_omega <- function(){
    f <- function(w){
      return(
        log(prior_param$beta_N) - log(prior_param$beta_S) - lbeta(prior_param$alpha_N, prior_param$alpha_S) + 
          (prior_param$alpha_N - 1.) * log(w * prior_param$beta_N / prior_param$beta_S) - (prior_param$alpha_N + prior_param$alpha_S) * log(1 + w * prior_param$beta_N / prior_param$beta_S)
      )
    }
    return(f)
  }
  
  
  # Posterior distribution for dN
  get_posterior_density_dN <- function(){
    f <- function(x){
      return(
        dgamma(x, 
               shape = prior_param$alpha_N + data$N, 
               rate = prior_param$beta_N + data$LN)
      )
    }
    return(f)
  }
  # Posterior distribution for dS
  get_posterior_density_dS <- function(){
    f <- function(x){
      return(
        dgamma(x, shape = prior_param$alpha_S + data$S, rate = prior_param$beta_S + data$LS)
      )
    }
    return(f)
  }
  
  # Quantile distribution for dN
  get_quantile_posterior_dN <- function(proba_values){
    return(
      qgamma(p = proba_values, 
             shape = prior_param$alpha_N + data$N, 
             rate = prior_param$beta_N + data$LN)
    )
  }
  # Quantile distribution for dS
  get_quantile_posterior_dS <- function(proba_values){
    return(
      qgamma(p = proba_values, shape = prior_param$alpha_S + data$S, rate = prior_param$beta_S + data$LS)
    )
  }
  # Mean and mode of the posterior
  mean_posterior_dN <- (prior_param$alpha_N + data$N) / (prior_param$beta_N + data$LN)
  mode_posterior_dN <- ifelse((prior_param$alpha_N + data$N) >= 1,
                              (prior_param$alpha_N + data$N - 1)/(prior_param$beta_N + data$LN),
                              0)
  mean_posterior_dS <- (prior_param$alpha_S + data$S) / (prior_param$beta_S + data$LS)
  mode_posterior_dS <- ifelse((prior_param$alpha_S + data$S) >= 1,
                              (prior_param$alpha_S + data$S - 1)/(prior_param$beta_S + data$LS),
                              0)
  
  
  # Posterior distribution for omega
  beta_1 <- prior_param$beta_N + data$LN
  beta_2 <- prior_param$beta_S + data$LS
  alpha_1 <- prior_param$alpha_N + data$N
  alpha_2 <- prior_param$alpha_S + data$S
  beta_ratio <- beta_1 / beta_2
  
  get_log_posterior_density_omega <- function(){
    f <- function(w){
      return(
        log(beta_ratio) - lbeta(alpha_1, alpha_2) + 
          (alpha_1 - 1.) * log(w * beta_ratio) - (alpha_1 + alpha_2) * log(1 + w * beta_ratio)
      )
    }
    return(f)
  }
  
  # Quantile distribution for omega
  get_quantile_posterior_omega <- function(proba_values){
    quantile_proba <- 
      (prior_param$beta_S + data$LS) / 
      (prior_param$beta_N + data$LN) *
      qbetapr(proba_values, prior_param$alpha_N + data$N, prior_param$alpha_S + data$S, scale = 1)
    
    return(quantile_proba)
  }
  
  #  Mean and mode of the posterior for omega
  mean_posterior_omega <- 
    (prior_param$beta_S + data$LS) / 
    (prior_param$beta_N + data$LN) * 
    (prior_param$alpha_N + data$N) / 
    (prior_param$alpha_S + data$S - 1)
  
  mode_posterior_omega <- 
    (prior_param$beta_S + data$LS) / 
    (prior_param$beta_N + data$LN) * 
    ifelse(prior_param$alpha_N + data$N >= 1,
           (prior_param$alpha_N + data$N - 1) / (prior_param$alpha_S + data$S + 1),
           0)
  
  posterior_proba_omega_greater_1 <- get_proba_beta_prime_above_x(x = beta_1 / (beta_2 + beta_1),
                                                                  alpha = alpha_1,
                                                                  beta = alpha_2)
  
  prior_proba_omega_greater_1 <- get_proba_beta_prime_above_x(x = prior_param$beta_N / (prior_param$beta_S + prior_param$beta_N),
                                                              alpha = prior_param$alpha_N,
                                                              beta = prior_param$alpha_S)
  
  # Marginal likelihood
  marginal_likelihood_N <- dnbinom(x = data$N, 
                                   size = prior_param$alpha_N, 
                                   prob = prior_param$beta_N / (prior_param$beta_N + data$LN))
  
  marginal_likelihood_S <- dnbinom(x = data$S, 
                                   size = prior_param$alpha_S, 
                                   prob = prior_param$beta_S / (prior_param$beta_S + data$LS))
  
  # Simulate from posterior (for posterior predictive checks)
  simulate_from_posterior_N <- function(n_rep){
    return(rnbinom(n = n_rep, 
                   size = prior_param$alpha_N + data$N, 
                   prob = (prior_param$beta_N + data$LN) / (prior_param$beta_N + 2 * data$LN)))
  }
  simulate_from_posterior_S <- function(n_rep){
    return(rnbinom(n = n_rep, 
                   size = prior_param$alpha_S + data$S, 
                   prob = (prior_param$beta_S + data$LS) / (prior_param$beta_S + 2 * data$LS)))
  }
  
  return(
    list(
      # Prior distribution
      'get_prior_density_dN' = get_prior_density_dN(),
      'get_prior_density_dS' = get_prior_density_dS(),
      'get_log_prior_density_omega' = get_log_prior_density_omega(),
      
      # Posterior distribution
      'get_posterior_density_dN' = get_posterior_density_dN(),
      'get_posterior_density_dS' = get_posterior_density_dS(),
      'get_log_posterior_density_omega' = get_log_posterior_density_omega(),
      
      'get_quantile_posterior_dN' = get_quantile_posterior_dN,
      'get_quantile_posterior_dS' = get_quantile_posterior_dS,
      'get_quantile_posterior_omega' = get_quantile_posterior_omega,
      
      'mean_posterior_dN' = mean_posterior_dN,
      'mean_posterior_dS' = mean_posterior_dS,
      'mean_posterior_omega' = mean_posterior_omega,
      'mode_posterior_dN' = mode_posterior_dN,
      'mode_posterior_dS' = mode_posterior_dS,
      'mode_posterior_omega' = mode_posterior_omega,
      
      # Probability that the posterior and prior for omega are > 1
      'posterior_proba_omega_greater_1' = posterior_proba_omega_greater_1,
      'prior_proba_omega_greater_1' = prior_proba_omega_greater_1,
      
      # Marginal likelihood
      'marginal_likelihood_N' = marginal_likelihood_N,
      'marginal_likelihood_S' = marginal_likelihood_S,
      
      # Simulations from the posterior
      'simulate_from_posterior_N' = simulate_from_posterior_N,
      'simulate_from_posterior_S' = simulate_from_posterior_S
    )
  )
}
