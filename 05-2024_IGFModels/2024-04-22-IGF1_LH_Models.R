# Load libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(posterior, parameters, insight, flextable, see, tidybayes, rethinking, rstanarm, cmdstanr, tidyverse, lmerTest, ggeffects, plyr, reshape2, rptR, viridis, cowplot, bayesplot, patchwork, lubridate)

# Set wd
setwd("C:/Users/sanja/OneDrive/Desktop/PostDocStuff/IGF-1")

# Read Data
igf_lh_data <- read.csv('2024MARCH_IGF1DataCollection/2024May_FinalIGF-LHDataset.csv',  header = T, stringsAsFactors = F, fileEncoding="UTF-8-BOM")
source("03-2024-IGF1DataAnalysis/05-2024_IGFModels/useful_functions.R")

# Only include Ashworth Data and data that passed intra-assay CV < 10% (14 samples failed i.e. 2% of samples)
igf_lh_data <- igf_lh_data %>%
  filter(Lab =="Ashworth" & IntraCVPass == 1)

# # Temp dataset for stan models 
temp <- igf_lh_data

# Subset relevant columns
temp <- temp %>%
  select(ID, Weight, ForeLeg, Survival, BredAsAYearling, IGF1, SexF, Twin, PopSize, BirthYear, MumID, ELISARunDate, PlateNumber, DaysSinceBirth, BirthWt) 


# Split data into different subsets (weight, foreleg, survival, repro)
temp_s <- temp %>% select(-ForeLeg, -Weight, -BredAsAYearling, -BirthWt, -DaysSinceBirth) %>% drop_na()
temp_sw <- temp %>% select(-ForeLeg, -BredAsAYearling, -BirthWt, -DaysSinceBirth) %>% drop_na()
temp_r <- temp %>% filter(Survival == 1) %>% select(-ForeLeg, -Weight, -Survival, -BirthWt, -DaysSinceBirth) %>% drop_na()

# # Convert cols to numeric
cols_n <- c("ELISARunDate", "PlateNumber", "BirthYear", "MumID")
temp_s <- convert_to_num_fac(temp_s, cols_n)
temp_sw <- convert_to_num_fac(temp_sw, cols_n)
temp_r <- convert_to_num_fac(temp_r, cols_n)

# Rescale variables
temp_s <- standardize_columns(temp_s, c("IGF1", "PopSize"))
temp_sw <- standardize_columns(temp_sw, c("IGF1", "PopSize", "Weight"))
temp_r <- standardize_columns(temp_r, c("IGF1", "PopSize"))

# Prepare data list to pass to stan model
data_survival <- prepare_data_list(temp_s)
data_survwt <- prepare_data_list(temp_sw)
data_repro <- prepare_data_list(temp_r)

# Run the model using cmdstanr
file_s <- here::here("03-2024-IGF1DataAnalysis/05-2024_IGFModels/StanModel_IGF_Survival.stan")
fit_mod_surv <- run_stan_model(file_s, data_survival)

file_sw <- here::here("03-2024-IGF1DataAnalysis/05-2024_IGFModels/StanModel_IGF_SurvivalWeight.stan")
fit_mod_survwt <- run_stan_model(file_sw, data_survwt)

file_r <- here::here("03-2024-IGF1DataAnalysis/05-2024_IGFModels/StanModel_IGF_Reproduction.stan")
fit_mod_repro <- run_stan_model(file_r, data_repro)

# Check the model
fit_mod_list <- list(
                     fit_mod_surv,
                     fit_mod_survwt,
                     fit_mod_repro
                     )
run_cmdstan_diagnose(fit_mod_list)

fit_mod_surv$cmdstan_diagnose()
fit_mod_survwt$cmdstan_diagnose()
fit_mod_repro$cmdstan_diagnose()

# Check posterior dist
color_scheme_set("mix-teal-pink")
(p4 <- generate_and_plot_ppc(fit_mod_surv, "Survival_rep", temp_s, temp_s$Survival, "Survival"))
(p5 <- generate_and_plot_ppc(fit_mod_survwt, "Survival_rep", temp_sw, temp_sw$Survival, "Survival"))
(p6 <- generate_and_plot_ppc(fit_mod_repro, "BredAsAYearling_rep", temp_r, temp_r$BredAsAYearling, "BredAsAYearling"))


# Exrtact samples 
# Add "true" IGF values estimated in model to original dataset with observed values to compare against
temp_s <- process_fit_mod(fit_mod_surv, temp_s)
temp_sw <- process_fit_mod(fit_mod_survwt, temp_sw)
temp_r <- process_fit_mod(fit_mod_repro, temp_r)

# Plotting
# Scatter-plot of IGF_true and observed IGF values 
plot_violin(temp_s, "Survival")
plot_violin(temp_sw, "Survival (Model accounting for variation in August weight)")
plot_violin(temp_r, "BredAsAYearling")

# Compare means and variances of IGF_obs vs IGF_true
var(temp_s$IGF1_sc)
var(temp_s$IGF_true)
cor.test(temp_s$IGF1_sc, temp_s$IGF_true)

var(temp_sw$IGF1_sc)
var(temp_sw$IGF_true)
cor.test(temp_sw$IGF1_sc, temp_sw$IGF_true)

var(temp_r$IGF1_sc)
var(temp_r$IGF_true)
cor.test(temp_r$IGF1_sc, temp_r$IGF_true)

#
# Get posterior draws of weight and compare with observed
# bayesplot::mcmc_intervals(fit_mod_weight$draws(), pars = c("beta_IGF", "beta_SexF"))
post_survival <- process_fit_bernoulli_model(fit_mod_surv, temp_s, 
                                  c("alpha", "beta_IGF", "beta_SexF", "beta_Twin", "beta_PopSize"), 
                                  "Survival")

post_survwt <- process_fit_bern_wt_model(fit_mod_survwt, temp_sw, 
                                             c("alpha", "beta_IGF", "beta_SexF", "beta_Twin", "beta_PopSize", "beta_Weight"), 
                                             "Survival")

post_repro <- process_fit_bernoulli_model(fit_mod_repro, temp_r, 
                                             c("alpha", "beta_IGF", "beta_SexF", "beta_Twin", "beta_PopSize"), 
                                             "BredAsAYearling")

#Plotting
ggplot(data= temp_s, aes(x=IGF_true, y=Survival)) +
  geom_point(colour="#078282FF", alpha=0.6) +
  geom_line(data=post_survival$summary_500draws_mu,  aes(y=value, group=draw), colour="#339E66FF", alpha=1/15) +
  geom_line(data=post_survival$mu_median,  aes(y=value), size=1.2, colour="#078282FF", alpha=0.8) +
  theme_cowplot() +
  ylab("First-Year Survival Probability") + 
  xlab("Normalized True IGF-1 Concentration (ng/ml)") +
  panel_bg(fill = "gray95", color = NA) +
  grid_lines(color = "white") +
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=10))

ggplot(data= temp_sw, aes(x=IGF_true, y=Survival)) +
  geom_point(colour="#078282FF", alpha=0.6) +
  geom_line(data=post_survwt$summary_500draws_mu,  aes(y=value, group=draw), colour="#339E66FF", alpha=1/15) +
  geom_line(data=post_survwt$mu_median,  aes(y=value), size=1.2, colour="#078282FF", alpha=0.8) +
  theme_cowplot() +
  ylab("First-Year Survival Probability") + 
  xlab("Normalized True IGF-1 Concentration (ng/ml)") +
  panel_bg(fill = "gray95", color = NA) +
  grid_lines(color = "white") +
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=10))

ggplot(data= temp_r, aes(x=IGF_true, y=BredAsAYearling)) +
  geom_point(colour="#078282FF", alpha=0.6) +
  geom_line(data=post_repro$summary_500draws_mu,  aes(y=value, group=draw), colour="#339E66FF", alpha=1/15) +
  geom_line(data=post_repro$mu_median,  aes(y=value), size=1.2, colour="#078282FF", alpha=0.8) +
  theme_cowplot() +
  ylab("Annual Reproduction") + 
  xlab("Normalized True IGF-1 Concentration") +
  panel_bg(fill = "gray95", color = NA) +
  grid_lines(color = "white") +
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=10))


# Model summary
# List the predictors for each model
predictors_surv <- c("alpha", "beta_IGF", "beta_SexF", "beta_Twin", "beta_PopSize", "sigma_e1", "sigma_t", "sigma_u", "sigma_v", "sigma_w")
predictors_survwt <- c("alpha", "beta_IGF", "beta_SexF", "beta_Twin", "beta_PopSize", "beta_Weight", "sigma_e1", "sigma_t", "sigma_u", "sigma_v", "sigma_w")
predictors_repro <- c("alpha", "beta_IGF", "beta_SexF", "beta_Twin", "beta_PopSize", "sigma_e1", "sigma_t", "sigma_u", "sigma_v", "sigma_w")

# Call the function for each model
full_post_surv <- summarize_model(fit_mod_surv, predictors_surv)
full_post_survwt <- summarize_model(fit_mod_survwt, predictors_survwt)
full_post_repro <- summarize_model(fit_mod_repro, predictors_repro)

full_post_df <- bind_rows(full_post_weight, full_post_foreleg, full_post_growth)
full_post_df

# Combine model summary tables
ft <- flextable(full_post_df) 
ft <- colformat_double(x=ft, j=c(2:9), digits = 3)  
ft <- autofit(ft)
ft
ft %>% save_as_docx(path="./IGF1_Writeup/Tables/IGF1_MorphoModelsTable.docx")



# Check survival model
# mod_surv <- cmdstan_model(file_s)
# fit_mod_surv <- mod_surv$sample(
#   data = data_survival,
#   seed = 123,
#   chains = 4,
#   parallel_chains = 4,
#   refresh = 500,
#   iter_sampling=6000,
#   iter_warmup = 2000,
#   adapt_delta = 0.99999
# )



# Random bits of code: 
# Weight_rep_post <- fit_mod_weight$draws(variables="Weight_rep", format="matrix")
# p1 <- ppc_dens_overlay(as.numeric(temp_w$Weight_sc), Weight_rep_post[1:100, ])
# ForeLeg_rep_post <- fit_mod_foreleg$draws(variables="ForeLeg_rep", format="matrix")
# p2 <- ppc_dens_overlay(as.numeric(temp_f$ForeLeg_sc), ForeLeg_rep_post[1:100, ])
# Survival_rep_post <- fit_mod_survival$draws(variables="Survival_rep", format="matrix")
# p3 <- ppc_dens_overlay(as.numeric(temp_s$Survival), Survival_rep_post[1:100, ])
# Repro_rep_post <- fit_mod_repro$draws(variables="Repro_rep", format="matrix")
# p4 <- ppc_dens_overlay(as.numeric(temp_s$OffspringNumber), Repro_rep_post[1:100, ])
# 

# Extract samples
# stanfit_weight <- rstan::read_stan_csv(fit_mod_weight$output_files())
# post_IGFTrue <- as.data.frame(rstan::extract(stanfit_weight, pars = c('IGF_true'), permuted = TRUE, include = TRUE)$IGF_true)
# 
# chk <- as.data.frame(rstan::extract(stanfit_weight, pars = c('IGF_true'), permuted = TRUE))
# 
# chk <- rstan::extract(stanfit_weight, pars = c('IGF_true'), permuted = FALSE) %>% reshape2::melt()
# 
# post_WeightRep <- as.data.frame(rstan::extract(stanfit_weight, pars = c('Weight_rep'), permuted = FALSE, include = TRUE)$Weight_rep)

# Add "true" IGF values estimated in model to original dataset with observed values to compare against
# temp_w$IGF_true <- apply( post_IGFTrue , 2 , mean )
# temp_w$Weight_pred <- apply(post_WeightRep, 2, mean)

# post_w_IGFTrue <- fit_mod_weight %>%
#   spread_draws(IGF_true[i]) %>%
#   summarise_draws()
# 
# post_f_IGFTrue <- fit_mod_foreleg %>%
#   spread_draws(IGF_true[i]) %>%
#   summarise_draws()
# 
# post_s_IGFTrue <- fit_mod_survival %>%
#   spread_draws(IGF_true[i]) %>%
#   summarise_draws()
# 
# post_r_IGFTrue <- fit_mod_repro %>%
#   spread_draws(IGF_true[i]) %>%
#   summarise_draws()
# Add "true" IGF values estimated in model to original dataset with observed values to compare against
# temp_w$IGF_true <- post_w_IGFTrue$mean
# temp_f$IGF_true <- post_f_IGFTrue$mean
# temp_s$IGF_true <- post_s_IGFTrue$mean
# temp_r$IGF_true <- post_r_IGFTrue$mean
# 
# # Plotting
# # Scatter-plot of IGF_true and observed IGF values 
# temp_w %>%
#   select(IGF1_sc, IGF_true) %>%
#   melt() %>%
#   mutate(SampleID = rep(1:680, times=2)) %>%
#   ggplot(aes(x=variable,y=value)) + 
#   geom_violin(aes(group=variable, fill=variable, colour=variable), alpha=0.5) +
#   geom_line(aes(group=SampleID), alpha=0.2) +
#   scale_fill_manual(values=c("#B4A7D6","#C27BA0")) +
#   scale_colour_manual(values=c("#B4A7D6","#C27BA0")) +
#   coord_flip()
# 
# temp_f %>%
#   select(IGF1_sc, IGF_true) %>%
#   melt() %>%
#   mutate(SampleID = rep(1:680, times=2)) %>%
#   ggplot(aes(x=variable,y=value)) + 
#   geom_violin(aes(group=variable, fill=variable, colour=variable), alpha=0.5) +
#   geom_line(aes(group=SampleID), alpha=0.2) +
#   scale_fill_manual(values=c("#B4A7D6","#C27BA0")) +
#   scale_colour_manual(values=c("#B4A7D6","#C27BA0")) +
#   coord_flip()
# 
# temp_s %>%
#   select(IGF1_sc, IGF_true) %>%
#   melt() %>%
#   mutate(SampleID = rep(1:680, times=2)) %>%
#   ggplot(aes(x=variable,y=value)) + 
#   geom_violin(aes(group=variable, fill=variable, colour=variable), alpha=0.5) +
#   geom_line(aes(group=SampleID), alpha=0.2) +
#   scale_fill_manual(values=c("#B4A7D6","#C27BA0")) +
#   scale_colour_manual(values=c("#B4A7D6","#C27BA0")) +
#   coord_flip()
# 
# temp_r %>%
#   select(IGF1_sc, IGF_true) %>%
#   melt() %>%
#   mutate(SampleID = rep(1:680, times=2)) %>%
#   ggplot(aes(x=variable,y=value)) + 
#   geom_violin(aes(group=variable, fill=variable, colour=variable), alpha=0.5) +
#   geom_line(aes(group=SampleID), alpha=0.2) +
#   scale_fill_manual(values=c("#B4A7D6","#C27BA0")) +
#   scale_colour_manual(values=c("#B4A7D6","#C27BA0")) +
#   coord_flip()
# 
# post.mod_weight = fit_mod_weight$draws(variables = c("alpha", "beta_IGF", "beta_SexF", "beta_Twin", "beta_PopSize"), format = "df")
# 
# mu.wt.link <- function(IGF1_true) post.mod_weight$alpha + post.mod_weight$beta_SexF*1 + post.mod_weight$beta_Twin*0 + post.mod_weight$beta_PopSize*mean(temp_w$PopSize_sc) + post.mod_weight$beta_IGF*(IGF1_true)
# igf.seq <- modelr::seq_range(temp_w$IGF_true, n = 200)
# mu.wt <- sapply( igf.seq , mu.wt.link )
# 
# # Get means and 95% CIs
# mu.wt.mean <- apply( mu.wt , 2 , mean )
# mu.wt.CI <- apply( mu.wt , 2 , PI , prob=0.95)
# summ.mu <- rbind(mu.wt.mean, mu.wt.CI)
# 
# # Get 500 draws and put in a df
# mu.wt_df <- as.data.frame(mu.wt)
# mu.wt_100 <- sample_n(mu.wt_df, 500)
# mu.wt_100$draw <- seq(1,500,1)
# mu.wt_100.melt <- mu.wt_100 %>% melt(id.vars = c("draw")) %>% arrange(draw)
# summ.draw.mu <- cbind(igf.seq, mu.wt_100.melt)
# summ.draw.mu <- summ.draw.mu %>%
#   dplyr::rename(IGF_true = igf.seq, Weight_sc=variable)
# mu.wt.median <- apply( mu.wt , 2 , median ) %>%
#   as.data.frame() 
# mu.wt.median <- cbind(igf.seq, mu.wt.median)
# colnames(mu.wt.median) <- c("IGF_true", "value")

ggplot(temp, aes(x=DaysSinceBirth, y=Weight)) +
  geom_point()

temp_50dsb <- temp %>%
  filter(DaysSinceBirth > 50)

ggplot(temp_50dsb, aes(x=DaysSinceBirth, y=Weight)) +
  geom_point() +
  geom_smooth()


library(mgcv)
gam1 <- gam(Weight ~ s(DaysSinceBirth, k=6), data=temp)
summary(gam1)

lm1 <- lm(Weight ~ DaysSinceBirth, data=temp)
summary(lm1)
AIC(lm1)

temp_nom <- temp %>%
  drop_na(DaysSinceBirth)

lm2 <- lmer(Weight ~ factor(SexF) + factor(Twin) + (1|BirthYear) + (1|MumID) + PopSize + BirthWt + DaysSinceBirth, REML=F, data=temp_nom)
summary(lm2)
AIC(lm2)
