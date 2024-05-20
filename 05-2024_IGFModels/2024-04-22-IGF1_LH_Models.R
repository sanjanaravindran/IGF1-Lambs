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
temp_rw <- temp %>% filter(Survival == 1) %>% select(-ForeLeg, -Survival, -BirthWt, -DaysSinceBirth) %>% drop_na()

# # Convert cols to numeric
cols_n <- c("ELISARunDate", "PlateNumber", "BirthYear", "MumID")
temp_s <- convert_to_num_fac(temp_s, cols_n)
temp_sw <- convert_to_num_fac(temp_sw, cols_n)
temp_r <- convert_to_num_fac(temp_r, cols_n)
temp_rw <- convert_to_num_fac(temp_rw, cols_n)

# Rescale variables
temp_s <- standardize_columns(temp_s, c("IGF1", "PopSize"))
temp_sw <- standardize_columns(temp_sw, c("IGF1", "PopSize", "Weight"))
temp_r <- standardize_columns(temp_r, c("IGF1", "PopSize"))
temp_rw <- standardize_columns(temp_rw, c("IGF1", "PopSize", "Weight"))

# Prepare data list to pass to stan model
data_survival <- prepare_data_list(temp_s)
data_survwt <- prepare_data_list(temp_sw)
data_repro <- prepare_data_list(temp_r)
data_reprowt <- prepare_data_list(temp_rw)

# Run the model using cmdstanr
file_s <- here::here("03-2024-IGF1DataAnalysis/05-2024_IGFModels/StanModel_IGF_Survival.stan")
file_sw <- here::here("03-2024-IGF1DataAnalysis/05-2024_IGFModels/StanModel_IGF_SurvivalWeight.stan")
file_r <- here::here("03-2024-IGF1DataAnalysis/05-2024_IGFModels/StanModel_IGF_Reproduction.stan")
file_rw <- here::here("03-2024-IGF1DataAnalysis/05-2024_IGFModels/StanModel_IGF_ReproductionWeight.stan")

fit_mod_surv <- run_stan_model(file_s, data_survival)
fit_mod_survwt <- run_stan_model(file_sw, data_survwt)
fit_mod_repro <- run_stan_model(file_r, data_repro)
fit_mod_reprowt <- run_stan_model(file_rw, data_reprowt)

# Check the model
fit_mod_list <- list(
                     fit_mod_surv,
                     fit_mod_survwt,
                     fit_mod_repro,
                     fit_mod_reprowt
                     )
run_cmdstan_diagnose(fit_mod_list)

fit_mod_surv$cmdstan_diagnose()
fit_mod_survwt$cmdstan_diagnose()
fit_mod_repro$cmdstan_diagnose()
fit_mod_reprowt$cmdstan_diagnose()

# Check posterior dist
color_scheme_set("mix-teal-pink")
(p4 <- generate_and_plot_ppc(fit_mod_surv, "Survival_rep", temp_s, temp_s$Survival, "Survival"))
(p5 <- generate_and_plot_ppc(fit_mod_survwt, "Survival_rep", temp_sw, temp_sw$Survival, "Survival"))
(p6 <- generate_and_plot_ppc(fit_mod_repro, "BredAsAYearling_rep", temp_r, temp_r$BredAsAYearling, "BredAsAYearling"))
(p7 <- generate_and_plot_ppc(fit_mod_reprowt, "BredAsAYearling_rep", temp_rw, temp_rw$BredAsAYearling, "BredAsAYearling"))


# Exrtact samples 
# Add "true" IGF values estimated in model to original dataset with observed values to compare against
temp_s <- process_fit_mod(fit_mod_surv, temp_s)
temp_sw <- process_fit_mod(fit_mod_survwt, temp_sw)
temp_r <- process_fit_mod(fit_mod_repro, temp_r)
temp_rw <- process_fit_mod(fit_mod_reprowt, temp_rw)

# Plotting
# Scatter-plot of IGF_true and observed IGF values 
plot_violin(temp_s, "Survival")
plot_violin(temp_sw, "Survival (Model accounting for variation in August weight and Population Size)")
plot_violin(temp_r, "Reproduction")
plot_violin(temp_rw, "Reproduction (Model accounting for variation in August weight and Population Size)")

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

var(temp_rw$IGF1_sc)
var(temp_rw$IGF_true)
cor.test(temp_rw$IGF1_sc, temp_rw$IGF_true)

#
# Get posterior draws of weight and compare with observed
# bayesplot::mcmc_intervals(fit_mod_weight$draws(), pars = c("beta_IGF", "beta_SexF"))
post_survival <- process_fit_bernoulli_model(fit_mod_surv, temp_s, 
                                  c("alpha", "beta_IGF", "beta_SexF", "beta_Twin"), 
                                  "Survival")

post_survwt <- process_fit_bern_wt_model(fit_mod_survwt, temp_sw, 
                                             c("alpha", "beta_IGF", "beta_SexF", "beta_Twin", "beta_PopSize", "beta_Weight"), 
                                             "Survival")

post_repro <- process_fit_repro_model(fit_mod_repro, temp_r, 
                                             c("alpha", "beta_IGF", "beta_SexF"), 
                                             "BredAsAYearling")

post_reprowt <- process_fit_reprowt_model(fit_mod_reprowt, temp_rw, 
                                         c("alpha", "beta_IGF", "beta_SexF", "beta_PopSize", "beta_Weight"), 
                                         "BredAsAYearling")

#Plotting
(plot_surv <- ggplot(data= temp_s, aes(x=IGF_true, y=Survival)) +
  geom_point(colour="#5e8d83", alpha=0.6) +
  geom_line(data=post_survival$summary_500draws_mu,  aes(y=value, group=draw), colour="#339E66FF", alpha=1/15) +
  geom_line(data=post_survival$mu_median,  aes(y=value), size=1.2, colour="#078282FF", alpha=0.8) +
  theme_cowplot() +
  ylab("First-Year Overwinter Survival") + 
  xlab("") +
  panel_bg(fill = "gray95", color = NA) +
  grid_lines(color = "white") +
  theme(legend.text = element_text(size=14),
        legend.title = element_text(size=14)))

(plot_survwt <- ggplot(data= temp_sw, aes(x=IGF_true, y=Survival)) +
  geom_point(colour="#ccb0be", alpha=0.6) +
  geom_line(data=post_survwt$summary_500draws_mu,  aes(y=value, group=draw), colour="#72668a", alpha=1/15) +
  geom_line(data=post_survwt$mu_median,  aes(y=value), size=1.2, colour="#374971", alpha=0.8) +
  theme_cowplot() +
  ylab("First-Year Overwinter Survival") + 
  xlab("") +
  panel_bg(fill = "gray95", color = NA) +
  grid_lines(color = "white") +
  ggtitle("Model accounting for condition and density effects") +
  theme(legend.text = element_text(size=14),
        legend.title = element_text(size=14),
        plot.title = element_text(size = 16, face="plain")))

(plot_repro <- ggplot(data= temp_r, aes(x=IGF_true, y=BredAsAYearling)) +
  geom_point(colour="#a67d65", alpha=0.6) +
  geom_line(data=post_repro$summary_500draws_mu,  aes(y=value, group=draw), colour="#a65e58", alpha=1/15) +
  geom_line(data=post_repro$mu_median,  aes(y=value), size=1.2, colour="#400101", alpha=0.8) +
  theme_cowplot() +
  ylab("First-Year Reproduction") + 
  xlab("Normalized True IGF-1 Concentration") +
  panel_bg(fill = "gray95", color = NA) +
  grid_lines(color = "white") +
  theme(legend.text = element_text(size=14),
        legend.title = element_text(size=14)))

(plot_reprowt <- ggplot(data= temp_rw, aes(x=IGF_true, y=BredAsAYearling)) +
    geom_point(colour="#ff9b54", alpha=0.6) +
    geom_line(data=post_reprowt$summary_500draws_mu,  aes(y=value, group=draw), colour="#fa9500", alpha=1/15) +
    geom_line(data=post_reprowt$mu_median,  aes(y=value), size=1.2, colour="#ff7f51", alpha=0.8) +
    theme_cowplot() +
    ylab("First-Year Reproduction") + 
    xlab("Normalized True IGF-1 Concentration") +
    panel_bg(fill = "gray95", color = NA) +
    grid_lines(color = "white") +
    ggtitle("Model accounting for condition and density effects") +
    theme(legend.text = element_text(size=14),
          legend.title = element_text(size=14),
          plot.title = element_text(size = 16, face="plain")))

# Save plots
plot2 <- cowplot::plot_grid(plot_surv, 
                            plot_survwt , 
                            plot_repro ,
                            plot_reprowt ,
                            nrow = 2, ncol=2,
                            labels = c("A", "B", "C", "D"),
                            align = "vh")

plot2

ggsave("./IGF1_Writeup/Figures/IGF1_LHModels.tiff", plot2, dpi=600, width=14, height=10, bg="white" )


# Model summary
# List the predictors for each model
predictors_surv <- c("alpha", "beta_IGF", "beta_SexF", "beta_Twin", "sigma_e1", "sigma_t", "sigma_u", "sigma_v", "sigma_w")
predictors_survwt <- c("alpha", "beta_IGF", "beta_SexF", "beta_Twin", "beta_PopSize", "beta_Weight", "sigma_e1", "sigma_t", "sigma_u", "sigma_v", "sigma_w")
predictors_repro <- c("alpha", "beta_IGF", "beta_SexF",  "sigma_e1", "sigma_t", "sigma_u", "sigma_v", "sigma_w")
predictors_reprowt <- c("alpha", "beta_IGF", "beta_SexF", "beta_PopSize", "beta_Weight", "sigma_e1", "sigma_t", "sigma_u", "sigma_v", "sigma_w")

# Call the function for each model
full_post_surv <- summarize_model(fit_mod_surv, predictors_surv)
full_post_survwt <- summarize_model(fit_mod_survwt, predictors_survwt)
full_post_repro <- summarize_model(fit_mod_repro, predictors_repro)
full_post_reprowt <- summarize_model(fit_mod_reprowt, predictors_reprowt)

full_post_df <- bind_rows(full_post_surv, full_post_survwt, full_post_repro, full_post_reprowt)
full_post_df

# Combine model summary tables
ft <- flextable(full_post_df) 
ft <- colformat_double(x=ft, j=c(2:9), digits = 3)  
ft <- autofit(ft)
ft
ft %>% save_as_docx(path="./IGF1_Writeup/Tables/IGF1_LHModelsTable.docx")

