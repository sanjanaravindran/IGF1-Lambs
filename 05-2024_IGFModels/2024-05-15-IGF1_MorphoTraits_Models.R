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

# Split data into different subsets (weight, foreleg, growth)
temp_w <- temp %>% select(-ForeLeg, -Survival, -BredAsAYearling, -BirthWt, -DaysSinceBirth) %>% drop_na()
temp_f <- temp %>% select(-Weight, -Survival, -BredAsAYearling, -BirthWt, -DaysSinceBirth) %>% drop_na()
temp_g <- temp %>% filter(DaysSinceBirth > 20) %>% select(-ForeLeg, -Survival, -BredAsAYearling) %>% drop_na()

# # Convert cols to numeric
cols_n <- c("ELISARunDate", "PlateNumber", "BirthYear", "MumID")
temp_w <- convert_to_num_fac(temp_w, cols_n)
temp_f <- convert_to_num_fac(temp_f, cols_n)
temp_g <- convert_to_num_fac(temp_g, cols_n)

# Rescale variables
temp_w <- standardize_columns(temp_w, c("Weight", "IGF1", "PopSize"))
temp_f <- standardize_columns(temp_f, c("ForeLeg", "IGF1", "PopSize"))
temp_g <- standardize_columns(temp_g, c("Weight", "IGF1", "PopSize", "DaysSinceBirth", "BirthWt"))


# Prepare data list to pass to stan model
data_weight <- prepare_data_list(temp_w)
data_foreleg <- prepare_data_list(temp_f)
data_growth <- prepare_data_list(temp_g)

# Run the model using cmdstanr
file_w <- here::here("03-2024-IGF1DataAnalysis/05-2024_IGFModels/StanModel_IGF_Weight.stan")
fit_mod_weight <- run_stan_model(file_w, data_weight)

file_f <- here::here("03-2024-IGF1DataAnalysis/05-2024_IGFModels/StanModel_IGF_ForeLeg.stan")
fit_mod_foreleg <- run_stan_model(file_f, data_foreleg)

file_g <- here::here("03-2024-IGF1DataAnalysis/05-2024_IGFModels/StanModel_IGF_Growth.stan")
fit_mod_growth <- run_stan_model(file_g, data_growth)


# Check the model diagnostics
fit_mod_list <- list(
  fit_mod_weight, 
  fit_mod_foreleg, 
  fit_mod_growth
)
run_cmdstan_diagnose(fit_mod_list)
# fit_mod_weight$cmdstan_diagnose()
# fit_mod_foreleg$cmdstan_diagnose()
# fit_mod_growth$cmdstan_diagnose()

# Check posterior dist
color_scheme_set("mix-teal-pink")
(p1 <- generate_and_plot_ppc(fit_mod_weight, "Weight_rep", temp_w, temp_w$Weight, "Weight"))
(p2 <- generate_and_plot_ppc(fit_mod_foreleg, "ForeLeg_rep", temp_f, temp_f$ForeLeg, "Foreleg"))
(p3 <- generate_and_plot_ppc(fit_mod_growth, "Weight_rep", temp_g, temp_g$Weight, "Weight (Model accounting for Birth Weight)"))

# Exrtact samples 
# Add "true" IGF values estimated in model to original dataset with observed values to compare against
temp_w <- process_fit_mod(fit_mod_weight, temp_w)
temp_f <- process_fit_mod(fit_mod_foreleg, temp_f)
temp_g <- process_fit_mod(fit_mod_growth, temp_g)

# Plotting
# Scatter-plot of IGF_true and observed IGF values 
plot_violin(temp_w, "Weight")
plot_violin(temp_f, "Foreleg Length")
plot_violin(temp_g, "Weight (Model accounts for variation in birth weight)")

# Compare means and variances of IGF_obs vs IGF_true
var(temp_w$IGF1_sc)
var(temp_w$IGF_true)
cor.test(temp_w$IGF_true, temp_w$IGF1_sc)

var(temp_f$IGF1_sc)
var(temp_f$IGF_true)
cor.test(temp_f$IGF_true, temp_f$IGF1_sc)

var(temp_g$IGF1_sc)
var(temp_g$IGF_true)
cor.test(temp_g$IGF_true, temp_g$IGF1_sc)

# Get posterior draws of weight and compare with observed
# bayesplot::mcmc_intervals(fit_mod_weight$draws(), pars = c("beta_IGF", "beta_SexF"))
post_weight <- process_fit_linear_model(fit_mod_weight, temp_w, 
                                        c("alpha", "beta_IGF", "beta_SexF", "beta_Twin", "beta_PopSize"), 
                                        "Weight")

post_foreleg <- process_fit_linear_model(fit_mod_foreleg, temp_f, 
                                         c("alpha", "beta_IGF", "beta_SexF", "beta_Twin", "beta_PopSize"), 
                                         "ForeLeg")

post_growth <- process_fit_growth_model(fit_mod_growth, temp_g, 
                                        c("alpha", "beta_IGF", "beta_SexF", "beta_Twin", "beta_PopSize", "beta_DSB", "beta_BirthWt"), 
                                        "Weight")

#Plotting
(plot_weight <- ggplot(data= temp_w, aes(x=IGF_true, y=Weight)) +
  geom_point(colour="#FFBA52FF", alpha=0.6) +
  geom_line(data=post_weight$summary_500draws_mu,  aes(y=value, group=draw), colour="#E683A9FF", alpha=1/15) +
  geom_line(data=post_weight$mu_median,  aes(y=value), size=1.2, colour="#FFBA52FF", alpha=0.8) +
  theme_cowplot() +
  ylab("August Body Weight (in kg)") + 
  xlab("") +
  scale_x_continuous(n.breaks=6) +
  panel_bg(fill = "gray95", color = NA) +
  grid_lines(color = "white") +
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=10)))


#Plotting
(plot_foreleg <- ggplot(data= temp_f, aes(x=IGF_true, y=ForeLeg)) +
  geom_point(colour="#949398FF", alpha=0.6) +
  geom_line(data=post_foreleg$summary_500draws_mu,  aes(y=value, group=draw), colour="#F4DF4EFF", alpha=1/15) +
  geom_line(data=post_foreleg$mu_median,  aes(y=value), size=1.2, colour="#949398FF", alpha=0.8) +
  theme_cowplot() +
  ylab("August Foreleg Length (in cm)") + 
  xlab("Normalized True IGF-1 Concentration") +
  scale_x_continuous(n.breaks=6) +
  panel_bg(fill = "gray95", color = NA) +
  grid_lines(color = "white") +
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=10)))


#Plotting
(plot_growth <- ggplot(data= temp_g, aes(x=IGF_true, y=Weight)) +
  geom_point(colour="#2C5F2D", alpha=0.6) +
  geom_line(data=post_growth$summary_500draws_mu,  aes(y=value, group=draw), colour="#97BC62FF", alpha=1/15) +
  geom_line(data=post_growth$mu_median,  aes(y=value), size=1.2, colour="#2C5F2D", alpha=0.8) +
  theme_cowplot() +
  ylab("August Body Weight (in kg)") + 
  xlab("") +
  scale_x_continuous(n.breaks=6) +
  panel_bg(fill = "gray95", color = NA) +
  grid_lines(color = "white") +
  ggtitle("Model accounting for variation in birth weight") +
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=10),
        plot.title = element_text(size = 10, face="plain")))

# Save plots
plot2 <- cowplot::plot_grid(plot_weight, 
                            plot_foreleg , 
                            plot_growth ,
                            nrow = 1,
                            labels = c("A", "B", "C"),
                            align = "vh")

plot2

ggsave("./IGF1_Writeup/Figures/IGF1_MorphoModels.tiff", plot2, dpi=600, width=12, height=4, bg="white" )


# Model summary
# List the predictors for each model
predictors_weight <- c("alpha", "beta_IGF", "beta_SexF", "beta_Twin", "beta_PopSize", "sigma_e", "sigma_e1", "sigma_t", "sigma_u", "sigma_v", "sigma_w")
predictors_foreleg <- c("alpha", "beta_IGF", "beta_SexF", "beta_Twin", "beta_PopSize", "sigma_e", "sigma_e1", "sigma_t", "sigma_u", "sigma_v", "sigma_w")
predictors_growth <- c("alpha", "beta_IGF", "beta_SexF", "beta_Twin", "beta_PopSize", "beta_DSB", "beta_BirthWt", "sigma_e", "sigma_e1", "sigma_t", "sigma_u", "sigma_v", "sigma_w")

# Call the function for each model
full_post_weight <- summarize_model(fit_mod_weight, predictors_weight)
full_post_foreleg <- summarize_model(fit_mod_foreleg, predictors_foreleg)
full_post_growth <- summarize_model(fit_mod_growth, predictors_growth)

full_post_df <- bind_rows(full_post_weight, full_post_foreleg, full_post_growth)
full_post_df

# Combine model summary tables
ft <- flextable(full_post_df) 
ft <- colformat_double(x=ft, j=c(2:9), digits = 3)  
ft <- autofit(ft)
ft
ft %>% save_as_docx(path="./IGF1_Writeup/Tables/IGF1_MorphoModelsTable.docx")

