# Some useful functions for IGF-1 data analysis script

# Convert multiple columns to factor 
convert_to_factor <- function(data, cols) {
  for (i in cols) {
    data[[i]] <- factor(data[[i]])
  }
  return(data)
}

# Rescale variables
standardize_columns <- function(data, cols) {
  for (col in cols) {
    data[[paste0(col, "_sc")]] <- as.numeric(scale(data[[col]]))
  }
  return(data)
}

# # Extract variance components from rstanarm model 
# extract_random_variances <- function(posterior_draws) {
#   
#   # Check if the input is a dataframe
#   if (!is.data.frame(posterior_draws)) {
#     stop("Input must be a dataframe of posterior draws.")
#   }
#   
#   #Create a dataframe to store parameters and variances
#   random_vars <- data.frame(Parameter = character(),
#                             Variance = numeric(),
#                             stringsAsFactors = FALSE)
#   
#   # Extract random effect variance
#   for (col in names(posterior_draws)) {
#     if (startsWith(col, "Sigma")) {
#       var_name <- substr(col, 6, nchar(col))
#       variance <- median(posterior_draws[[col]])
#       random_vars <- rbind(random_vars, data.frame(Parameter = var_name, Variance = variance))
#     }
#   }
#     # Extract residual variance
#     residual_variance <- median(posterior_draws$sigma^2)
#     random_vars <- rbind(random_vars, data.frame(Parameter = "Residual", Variance = residual_variance))
#   
#     # Return data frame of variances
#     return(random_vars)
# }
# 
 # Add proportion of variance explained in the dataframe 
add_prop_var <- function(data) {
   var_sum = sum(data$vcov)
   prop_variance = data$vcov/var_sum
   data <- cbind(data, prop_variance)
   return(data)
 }


# Convert multiple columns to factor 
convert_to_num_fac <- function(data, cols) {
  for (i in cols) {
    data[[i]] <- factor(data[[i]])
    data[[i]] <- as.numeric(data[[i]])
  }
  return(data)
}

# Function to prepare data for stan modelling
prepare_data_list <- function(data) {
  data_list <- list(
    N = nrow(data),
    BirthYear = data$BirthYear,
    MumID = data$MumID,
    ELISARunDate = data$ELISARunDate,
    PlateNumber = data$PlateNumber,
    IGF1_sc = data$IGF1_sc,
    # StorageTime_sc = data$StorageTime_sc,
    SexF = data$SexF,
    Twin = data$Twin,
    PopSize_sc = data$PopSize_sc,
    num_years = length(unique(data$BirthYear)),
    num_mums = length(unique(data$MumID)),
    num_days = length(unique(data$ELISARunDate)),
    num_plates = length(unique(data$PlateNumber))
  )
  if ("Weight" %in% names(data)) {
    data_list$Weight = data$Weight
  }
  if ("Survival" %in% names(data)) {
    data_list$Survival = data$Survival
  }
  if ("Survival" %in% names(data) & "Weight_sc" %in% names(data)) {
    data_list$Survival = data$Survival
    data_list$Weight_sc = data$Weight_sc
  }
  if ("ForeLeg" %in% names(data)) {
    data_list$ForeLeg = data$ForeLeg
  }
  if ("BredAsAYearling" %in% names(data)) {
    data_list$BredAsAYearling = data$BredAsAYearling
    data_list$Survival = data$Survival
  }
  if ("DaysSinceBirth_sc" %in% names(data)) {
    data_list$Weight = data$Weight
    data_list$BirthWt_sc = data$BirthWt_sc
    data_list$DaysSinceBirth_sc = data$DaysSinceBirth_sc
  }
  return(data_list)
}

# Function to run stan model
run_stan_model <- function(model_file, data_list, seed = 123, chains = 4, parallel_chains = 4, refresh = 500   
                           # iter_sampling=4000,
                           # iter_warmup = 1000, 
                           # adapt_delta = 0.99999
                           ) {
  mod <- cmdstan_model(model_file)
  fit_mod <- mod$sample(
    data = data_list,
    seed = seed,
    chains = chains,
    parallel_chains = parallel_chains,
    refresh = refresh
    #iter_sampling = iter_sampling,
    #iter_warmup = iter_warmup,
    #adapt_delta = adapt_delta
  )
  return(fit_mod)
}

# Function to run stan model (bernoulli)
run_stan_bernoulli_model <- function(model_file, data_list, seed = 123, chains = 4, parallel_chains = 4, refresh = 500,   
                           iter_sampling=4000,
                           iter_warmup = 2000, 
                           adapt_delta = 0.99999
) {
  mod <- cmdstan_model(model_file)
  fit_mod <- mod$sample(
    data = data_list,
    seed = seed,
    chains = chains,
    parallel_chains = parallel_chains,
    refresh = refresh,
    iter_sampling = iter_sampling,
    iter_warmup = iter_warmup,
    adapt_delta = adapt_delta
  )
  return(fit_mod)
}

# Function to run stan diagnostics
run_cmdstan_diagnose <- function(fit_mod_list) {
  for (fit_mod in fit_mod_list) {
    fit_mod$cmdstan_diagnose()
  }
}

# Function to generate and plot ppc
generate_and_plot_ppc <- function(fit_mod, variable_name, data, observed_data, trait, n_samples = 200) {
  posterior_draws <- fit_mod$draws(variables = variable_name, format = "matrix")
  ppc_data <- ppc_dens_overlay(as.numeric(observed_data), posterior_draws[1:n_samples, ]) + labs(subtitle = trait)
  return(ppc_data)
}

# Function to process post dist and get IGF_true values
process_fit_mod <- function(fit_mod, temp_data) {
  post_IGFTrue <- fit_mod %>%
    spread_draws(IGF_true[i]) %>%
    summarise_draws()
  
  temp_data$IGF_true <- post_IGFTrue$mean
  return(temp_data)
}

# Violin plotting function
plot_violin <- function(temp_data, trait) {
  temp_data %>%
    select(IGF1_sc, IGF_true) %>%
    melt() %>%
    mutate(SampleID = rep(1:nrow(temp_data), times=2)) %>%
    ggplot(aes(x=variable,y=value)) + 
    geom_violinhalf(aes(group=variable, fill=variable, colour=variable), alpha=0.5, flip=c(1)) +
    geom_line(aes(group=SampleID), alpha=0.2) +
    scale_fill_manual(values=c("#B4A7D6","#C27BA0"), labels= c("Observed IGF-1", "True IGF-1")) +
    scale_colour_manual(values=c("#B4A7D6","#C27BA0"), labels= c("Observed IGF-1", "True IGF-1")) +
    ylab("Normalized IGF-1 Concentration") +
    xlab("") +
    scale_x_discrete(labels= c("Observed IGF-1", "True IGF-1")) +
    theme(legend.position="none") +
    labs(subtitle=trait)
}


# Model processing to plot y~x and 500 draws from posterior 
process_fit_linear_model <- function(fit_mod, temp_data, variables, sample_variable) {
  post_mod <- fit_mod$draws(variables = variables, format = "df")
  
  mu.link <- function(IGF1_true)
      post_mod$alpha +
      post_mod$beta_SexF * 1 +
      post_mod$beta_Twin * 0 +
      post_mod$beta_PopSize * mean(temp_data$PopSize_sc) +
      post_mod$beta_IGF * IGF1_true
  
  igf_seq <- modelr::seq_range(temp_data$IGF_true, n = 200)
  mu <- sapply(igf_seq, mu.link)
  
  mu_mean <- apply(mu, 2, mean)
  mu_CI <- apply(mu, 2, PI, prob = 0.95)
  summary_mu <- rbind(mu_mean, mu_CI)
  
  mu_df <- as.data.frame(mu)
  mu_500 <- sample_n(mu_df, 500)
  mu_500$draw <- seq(1, 500, 1)
  mu_500_melt <- mu_500 %>% melt(id.vars = c("draw")) %>% arrange(draw)
  summ_draw_mu <- cbind(igf_seq, mu_500_melt)
  summ_draw_mu <- summ_draw_mu %>%
    dplyr::rename(IGF_true = igf_seq, !!as.name(sample_variable) := variable)
  mu_median <- apply(mu, 2, median) %>%
    as.data.frame()
  mu_median <- cbind(igf_seq, mu_median)
  colnames(mu_median) <- c("IGF_true", "value")
  
  return(list(summary_mu = summary_mu, summary_500draws_mu = summ_draw_mu, mu_median = mu_median))
}


# Model processing to plot y~x and 500 draws from posterior (growth model)
process_fit_growth_model <- function(fit_mod, temp_data, variables, sample_variable) {
  post_mod <- fit_mod$draws(variables = variables, format = "df")
  
  mu.link <- function(IGF1_true)
    post_mod$alpha +
    post_mod$beta_SexF * 1 +
    post_mod$beta_Twin * 0 +
    post_mod$beta_PopSize * mean(temp_data$PopSize_sc) +
    post_mod$beta_IGF * IGF1_true +
    post_mod$beta_DSB * mean(temp_data$DaysSinceBirth_sc) +
    post_mod$beta_BirthWt * mean(temp_data$BirthWt_sc) 
  
  igf_seq <- modelr::seq_range(temp_data$IGF_true, n = 200)
  mu <- sapply(igf_seq, mu.link)
  
  mu_mean <- apply(mu, 2, mean)
  mu_CI <- apply(mu, 2, PI, prob = 0.95)
  summary_mu <- rbind(mu_mean, mu_CI)
  
  mu_df <- as.data.frame(mu)
  mu_500 <- sample_n(mu_df, 500)
  mu_500$draw <- seq(1, 500, 1)
  mu_500_melt <- mu_500 %>% melt(id.vars = c("draw")) %>% arrange(draw)
  summ_draw_mu <- cbind(igf_seq, mu_500_melt)
  summ_draw_mu <- summ_draw_mu %>%
    dplyr::rename(IGF_true = igf_seq, !!as.name(sample_variable) := variable)
  mu_median <- apply(mu, 2, median) %>%
    as.data.frame()
  mu_median <- cbind(igf_seq, mu_median)
  colnames(mu_median) <- c("IGF_true", "value")
  
  return(list(summary_mu = summary_mu, summary_500draws_mu = summ_draw_mu, mu_median = mu_median))
}


# Model processing to plot y~x and 500 draws from posterior (binomial model)
process_fit_bernoulli_model <- function(fit_mod, temp_data, variables, sample_variable) {
  post_mod <- fit_mod$draws(variables = variables, format = "df")
  
  mu.link <- function(IGF1_true)
    inv_logit(post_mod$alpha +
    post_mod$beta_SexF * 1 +
    post_mod$beta_Twin * 0 +
    post_mod$beta_PopSize * mean(temp_data$PopSize_sc) +
    post_mod$beta_IGF * IGF1_true)
  
  igf_seq <- modelr::seq_range(temp_data$IGF_true, n = 200)
  mu <- sapply(igf_seq, mu.link)
  
  mu_mean <- apply(mu, 2, mean)
  mu_CI <- apply(mu, 2, PI, prob = 0.95)
  summary_mu <- rbind(mu_mean, mu_CI)
  
  mu_df <- as.data.frame(mu)
  mu_500 <- sample_n(mu_df, 500)
  mu_500$draw <- seq(1, 500, 1)
  mu_500_melt <- mu_500 %>% melt(id.vars = c("draw")) %>% arrange(draw)
  summ_draw_mu <- cbind(igf_seq, mu_500_melt)
  summ_draw_mu <- summ_draw_mu %>%
    dplyr::rename(IGF_true = igf_seq, !!as.name(sample_variable) := variable)
  mu_median <- apply(mu, 2, median) %>%
    as.data.frame()
  mu_median <- cbind(igf_seq, mu_median)
  colnames(mu_median) <- c("IGF_true", "value")
  
  return(list(summary_mu = summary_mu, summary_500draws_mu = summ_draw_mu, mu_median = mu_median))
}

# Model processing to plot y~x and 500 draws from posterior (binomial model)
process_fit_bern_wt_model <- function(fit_mod, temp_data, variables, sample_variable) {
  post_mod <- fit_mod$draws(variables = variables, format = "df")
  
  mu.link <- function(IGF1_true)
    inv_logit(post_mod$alpha +
                post_mod$beta_SexF * 1 +
                post_mod$beta_Twin * 0 +
                post_mod$beta_PopSize * mean(temp_data$PopSize_sc) +
                post_mod$beta_Weight * mean(temp_data$Weight_sc) +
                post_mod$beta_IGF * IGF1_true)
  
  igf_seq <- modelr::seq_range(temp_data$IGF_true, n = 200)
  mu <- sapply(igf_seq, mu.link)
  
  mu_mean <- apply(mu, 2, mean)
  mu_CI <- apply(mu, 2, PI, prob = 0.95)
  summary_mu <- rbind(mu_mean, mu_CI)
  
  mu_df <- as.data.frame(mu)
  mu_500 <- sample_n(mu_df, 500)
  mu_500$draw <- seq(1, 500, 1)
  mu_500_melt <- mu_500 %>% melt(id.vars = c("draw")) %>% arrange(draw)
  summ_draw_mu <- cbind(igf_seq, mu_500_melt)
  summ_draw_mu <- summ_draw_mu %>%
    dplyr::rename(IGF_true = igf_seq, !!as.name(sample_variable) := variable)
  mu_median <- apply(mu, 2, median) %>%
    as.data.frame()
  mu_median <- cbind(igf_seq, mu_median)
  colnames(mu_median) <- c("IGF_true", "value")
  
  return(list(summary_mu = summary_mu, summary_500draws_mu = summ_draw_mu, mu_median = mu_median))
}

# Another function pooling all these 3 together
process_fit_model <- function(fit_mod, temp_data, variables, sample_variable, model_type) {
  post_mod <- fit_mod$draws(variables = variables, format = "df")
  
  mu.link <- switch(
    model_type,
    "linear" = function(IGF1_true)
      post_mod$alpha +
      post_mod$beta_SexF * 1 +
      post_mod$beta_Twin * 0 +
      post_mod$beta_PopSize * mean(temp_data$PopSize_sc) +
      post_mod$beta_IGF * IGF1_true,
    "growth" = function(IGF1_true)
      post_mod$alpha +
      post_mod$beta_SexF * 1 +
      post_mod$beta_Twin * 0 +
      post_mod$beta_PopSize * mean(temp_data$PopSize_sc) +
      post_mod$beta_IGF * IGF1_true +
      post_mod$beta_DSB * mean(temp_data$DaysSinceBirth_sc) +
      post_mod$beta_BirthWt * mean(temp_data$BirthWt_sc),
    "bernoulli" = function(IGF1_true)
      inv_logit(post_mod$alpha +
                  post_mod$beta_SexF * 1 +
                  post_mod$beta_Twin * 0 +
                  post_mod$beta_PopSize * mean(temp_data$PopSize_sc) +
                  post_mod$beta_IGF * IGF1_true)
  )
  
  igf_seq <- modelr::seq_range(temp_data$IGF_true, n = 200)
  mu <- sapply(igf_seq, mu.link)
  
  mu_mean <- apply(mu, 2, mean)
  mu_CI <- apply(mu, 2, PI, prob = 0.95)
  summary_mu <- rbind(mu_mean, mu_CI)
  
  mu_df <- as.data.frame(mu)
  mu_500 <- sample_n(mu_df, 500)
  mu_500$draw <- seq(1, 500, 1)
  mu_500_melt <- mu_500 %>% melt(id.vars = c("draw")) %>% arrange(draw)
  summ_draw_mu <- cbind(igf_seq, mu_500_melt)
  summ_draw_mu <- summ_draw_mu %>%
    dplyr::rename(IGF_true = igf_seq, !!as.name(sample_variable) := variable)
  mu_median <- apply(mu, 2, median) %>%
    as.data.frame()
  mu_median <- cbind(igf_seq, mu_median)
  colnames(mu_median) <- c("IGF_true", "value")
  
  return(list(summary_mu = summary_mu, summary_500draws_mu = summ_draw_mu, mu_median = mu_median))
}

summarize_model <- function(model, variables) {
  model$summary(
    variables = variables,
    posterior::default_summary_measures()[1:4],
    quantiles = ~ quantile2(., probs = c(0.025, 0.975)),
    posterior::default_convergence_measures()
  )
}

# bla_bla <- function(fit_mod, temp_data, variables, sample_variable) {
#   post_mod <- fit_mod$draws(variables = variables, format = "df")
#   
#   mu.link <- function(IGF1_true)
#     post_mod$alpha +
#     # post_mod$beta_SexF * 1 +
#     # post_mod$beta_Twin * 0 +
#     # post_mod$beta_PopSize * mean(temp_data$PopSize_sc) +
#     post_mod$beta_IGF * IGF1_true 
#     # post_mod$beta_DSB * mean(temp_data$DaysSinceBirth_sc) +
#     # post_mod$beta_BirthWt * mean(temp_data$BirthWt_sc)
#   
#   igf_seq <- modelr::seq_range(temp_data$IGF_true, n = 200)
#   mu <- sapply(igf_seq, mu.link)
#   
#   mu_mean <- apply(mu, 2, mean)
#   mu_CI <- apply(mu, 2, PI, prob = 0.95)
#   summary_mu <- rbind(mu_mean, mu_CI)
#   
#   mu_df <- as.data.frame(mu)
#   mu_500 <- sample_n(mu_df, 500)
#   mu_500$draw <- seq(1, 500, 1)
#   mu_500_melt <- mu_500 %>% melt(id.vars = c("draw")) %>% arrange(draw)
#   summ_draw_mu <- cbind(igf_seq, mu_500_melt)
#   summ_draw_mu <- summ_draw_mu %>%
#     dplyr::rename(IGF_true = igf_seq, !!as.name(sample_variable) := variable)
#   mu_median <- apply(mu, 2, median) %>%
#     as.data.frame()
#   mu_median <- cbind(igf_seq, mu_median)
#   colnames(mu_median) <- c("IGF_true", "value")
#   
#   return(list(summary_mu = summary_mu, summary_500draws_mu = summ_draw_mu, mu_median = mu_median))
# }
