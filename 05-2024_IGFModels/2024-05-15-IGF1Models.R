# Load libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidybayes, rstanarm, rethinking, cmdstanr, parameters, tidyverse, lmerTest, ggeffects, plyr, reshape2, rptR, viridis, cowplot, bayesplot, patchwork, ggpubr)

# Set wd
setwd("C:/Users/sanja/OneDrive/Desktop/PostDocStuff/IGF-1")

# Read Data
igf_lh_data <- read.csv('2024MARCH_IGF1DataCollection/2024May_FinalIGF-LHDataset.csv',  header = T, stringsAsFactors = F, fileEncoding="UTF-8-BOM")
source("03-2024-IGF1DataAnalysis/useful_functions.R")

# Only include Ashworth Data and data where Intra-Assay CV < 10
igf_lh_data <- igf_lh_data %>%
  filter(Lab =="Ashworth" & IntraCVPass==1)

# 
table(igf_lh_data$CapYear)
table(igf_lh_data$Sex_R)
table(igf_lh_data$Twin)

# Temp dataset for rstanarm model 
temp <- igf_lh_data

# Convert relevant columns to factor 
cols_f <- c("SexF", "Twin")
temp <- convert_to_factor(temp, cols_f)

# Rescale variables 
cols_sc <- c("IGF1", "PopSize", "StorageTime")
temp <- standardize_columns(temp, cols_sc)

# IGF-1 Model
mod_igf <- rstanarm::stan_lmer(IGF1 ~ SexF + Twin + PopSize_sc + (1|ELISARunDate) + (1|PlateNumber) + (1|BirthYear) + (1|MumID),
                     cores=4, 
                     seed=12345,
                     # QR=TRUE,
                     data=temp)

IGF <- temp$IGF1
IGF_rep <- posterior_predict(mod_igf)
(p0 <- ppc_dens_overlay(IGF, IGF_rep[1:200,]) + labs(subtitle = "IGF-1")) 

# Check model in lme4
mod_igf_lmer <-lmer(IGF1 ~ SexF + Twin + PopSize_sc + (1|ELISARunDate) + (1|PlateNumber) + (1|BirthYear) + (1|MumID),
                     # cores=4, 
                     # seed=12345,
                     # QR=TRUE,
                     data=temp)
summary(mod_igf_lmer)

# Model diagnostics
color_scheme_set("mix-teal-purple")
pp_check(mod_igf) + labs(subtitle="IGF-1")

rhat(mod_igf) %>%
  as.data.frame() %>%
  filter(.[[1]] > 1.1 & .[[1]] < 0.99)

mcmc_neff(neff_ratio(mod_igf), size = 2)

# Extract posterior draws for later use
posterior_igf <- as.array(mod_igf)

# Plot predictors

predictors <- as.array(mod_igf, pars = c("SexF1", "Twin1", "PopSize_sc"))
p1 <- bayesplot::mcmc_areas(predictors, point_est = c("median"), prob = 0.5, prob_outer = 0.95) +
  theme_cowplot() +
  vline_0(color = "darkgray", linetype = 2) +
  scale_x_continuous(limits=c(-100,10)) +
  scale_y_discrete(labels=c("Population Size", "Twin:1", "Sex:Female"), limits = rev) 
p1  

# Extract variance components (median; note: estimates differ a fair bit when you take mean instead)
posterior_igf_df <- as.data.frame(mod_igf)
med_igf_df <- median_hdci(posterior_igf_df)
p_value(mod_igf, method="hdi")

# Extract variance components (mean)
variance_table <- as.data.frame(VarCorr(mod_igf))
variance_table <- add_prop_var(variance_table)

# Getting model predictions (conditional effects)
# Sex effect
sex_eff <- temp %>%
  modelr::data_grid(PopSize_sc = mean(PopSize_sc), 
                    # StorageTime_sc = mean(StorageTime_sc), 
                    SexF = unique(SexF), 
                    Twin=0) 

sex_eff <- convert_to_factor(sex_eff, cols_f)

sex_preds <- add_epred_draws(mod_igf, newdata = sex_eff, re_formula = NA) %>% 
  group_by(SexF)  %>% 
  median_hdci(.epred, .width = .95) 

p2 <- ggplot(sex_preds, aes(x = SexF)) +
  geom_jitter(data = temp, aes(x=SexF, y = IGF1, colour=SexF), size = 1.5,  alpha=0.3, height=0) +
  geom_pointrange(aes(ymin = .lower, ymax=.upper, y=.epred)) +
  theme_cowplot() +
  scale_colour_brewer(name = "Sex=Female", palette="Set2") +
  ylab("IGF-1 Concentration\n(ng/ml)") + 
  xlab("") +
  scale_x_discrete(labels=c("Male", "Female")) +
  panel_bg(fill = "gray95", color = NA) +
  grid_lines(color = "white") +
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=10)) +
  theme(legend.position="none")+
  theme(axis.text.x = element_text(size = 14))

p2


# Twin effect
twin_eff <- temp %>%
  modelr::data_grid(PopSize_sc = mean(PopSize_sc),  
                    Twin = unique(Twin), 
                    SexF=1) 

twin_eff <- convert_to_factor(twin_eff, cols_f)

twin_preds <- add_epred_draws(mod_igf, newdata = twin_eff, re_formula = NA) %>% 
  group_by(Twin)  %>% 
  median_hdci(.epred, .width = .95) 

p3 <- ggplot(twin_preds, aes(x = Twin)) +
  geom_jitter(data = temp, aes(x=Twin, y = IGF1, colour=Twin), size = 1.5,  alpha=0.5, height=0) +
  geom_pointrange(aes(ymin = .lower, ymax=.upper, y=.epred)) +
  theme_cowplot() +
  scale_colour_brewer(name = "Twin", palette="Pastel1") +
  ylab("") + 
  xlab("") +
  scale_x_discrete(labels=c("Singleton", "Twin")) +
  panel_bg(fill = "gray95", color = NA) +
  grid_lines(color = "white") +
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=10))+
  theme(legend.position="none") +
  theme(axis.text.x = element_text(size = 14))

p3


# Pop effect
pop_eff <- temp %>%
  modelr::data_grid(PopSize_sc = modelr::seq_range(PopSize_sc, n = 200), 
                    # StorageTime_sc = mean(StorageTime_sc), 
                    SexF = 1, 
                    Twin=0) 

pop_eff <- convert_to_factor(pop_eff, cols_f)

pop_preds <- add_epred_draws(mod_igf, newdata = pop_eff, re_formula = NA, ndraws=500) %>% 
  group_by(PopSize_sc) 

pop_preds_med <- add_epred_draws(mod_igf, newdata = pop_eff, re_formula = NA) %>% 
  group_by(PopSize_sc) %>%
  median_hdci(.epred, .width = .95) 

p4 <- ggplot(pop_preds, aes(x = PopSize_sc, y = IGF1)) +
  geom_jitter(data = temp, size = 1.5,  alpha=0.2, height=0, width=0.1, color = "#3E285C") +
  geom_line(aes(y = .epred, group = .draw), alpha=1/15, color = "#3E285C")  +
  geom_line(data=pop_preds_med,  aes(y=.epred), linewidth=1.2, colour="#9457EB", alpha=0.8) +
  theme_cowplot() +
  scale_colour_brewer(name = "Population Size") +
  ylab("") + 
  xlab("Population Size") +
  panel_bg(fill = "gray95", color = NA) +
  grid_lines(color = "white") +
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=10))

p4

# Arrange the plots
plot1 <- cowplot::plot_grid(p1 , 
                            nrow = 1,
                            labels = c("A"),
                            align = "h")

plot1
ggsave("./IGF1_Writeup/Figures/IGF1_ForestPlot.tiff", plot1, dpi=600, width=4, height=3, bg="white")

plot2 <- cowplot::plot_grid(p2, 
                   p3 , 
                   p4 ,
                   nrow = 1,
                   labels = c("B", "C", "D"),
                   align = "h")

plot2

ggsave("./IGF1_Writeup/Figures/IGF1_SexLitSizePopSize.tiff", plot2, dpi=600, width=12, height=4, bg="white" )


# Get the tables

modelsummary(mod_igf, statistic = "conf.int", output="./IGF1_Writeup/Tables/IGFModelTable.docx")


# Extra figs 
# IGF histogram 
(plot_igf <- ggplot(data=igf_lh_data, aes(x=IGF1)) +
  geom_histogram(bins=30)) +
  theme(legend.title=element_blank(), 
        text = element_text(size = 16))
  
# Histogram per sex
means_sex <- ddply(igf_lh_data, "Sex_R", summarise, grp.mean=mean(IGF1))
means_sex

# Histogram
plot_igf_sex <- ggplot(data=igf_lh_data, aes(x=IGF1, color=Sex_R)) +
  geom_histogram(fill="white", alpha=0.5, position="dodge")+
  geom_vline(data=means_sex, aes(xintercept=grp.mean, color=Sex_R),linetype="dashed") +
  scale_color_manual(labels=c("Female", "Male"), values = c("purple","orange")) +
  theme(legend.title=element_blank(),
        text = element_text(size = 16))

# Histogram per twin
means_twin <- ddply(igf_lh_data, "Twin", summarise, grp.mean=mean(IGF1))
means_twin

# Histogram
plot_igf_twins <- ggplot(data=igf_lh_data, aes(x=IGF1, color=factor(Twin))) +
  geom_histogram(fill="white", alpha=0.5, position="dodge")+
  geom_vline(data=means_twin, aes(xintercept=grp.mean, color=factor(Twin)),linetype="dashed") +
  scale_color_manual(labels=c("Singleton", "Twin"), values = c("darkgreen", "brown")) +
  theme(legend.title=element_blank(), 
        text = element_text(size = 16))

plot_igf + plot_igf_sex+plot_igf_twins

plot3 <- cowplot::plot_grid(plot_igf,
                            plot_igf_sex, 
                            plot_igf_twins , 
                            nrow = 1,
                            labels = c("A", "B", "C"),
                            align = "h")

plot3

ggsave("./IGF1_Writeup/Figures/IGF1_Histograms.tiff", plot3, dpi=600, width=18, height=4, bg="white" )

# Extra code bits below

# Plot variance components as stacked bar plot
# 
# prop_var
# prop_var_table <- variance_table %>%
#   select(grp, prop_variance) %>%
#   arrange(desc(prop_variance))
# prop_var_table$Trait = c("IGF1")
# 
# prop_var <- ggplot(prop_var_table, aes(fill=factor(grp, levels=c("BirthYear", "PlateNumber", "MumID", "ELISARunDate", "Residual")), y=prop_variance, x=Trait)) + 
#   geom_bar(position="fill", stat="identity") +
#   scale_fill_viridis(option="D", discrete = TRUE)  +
#   theme_cowplot() +
#   labs(x="IGF-1 Levels", y = "Proportion of Variance") +
#   theme(axis.text=element_text(size=12), axis.title=element_text(size=12), axis.text.x=element_blank())+
#   theme(legend.title = element_blank(), legend.text=element_text(size=11)) 

