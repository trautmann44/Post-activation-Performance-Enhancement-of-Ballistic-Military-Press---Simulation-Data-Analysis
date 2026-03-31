rm(list = ls())

# needed libraries
packages_needed <- c("ggplot2", "tidyquant", "ggdist", "ggthemes",
                     "modeest", "faux", "dplyr", "tidyr", "pastecs",
                     "rstatix", "readxl", "writexl", "irr", "car",
                     "emmeans", "corrplot", "MuMIn", "fitdistrplus", 
                     "lme4", "flexplot", "ggeffects", "lmtest", "lmerTest",
                     "influence.ME", "emmeans", "simr", "performance",
                     "DHARMa", "broom.mixed", "purrr", "vroom", "readr",
                     "sjPlot", "modeest", "gridGraphics", "gridExtra")

lapply(packages_needed, FUN = require, character.only = T)


# Collected ES's of studies from Krzysztofik et al. (2021).
previous_studies <- data.frame("BBPT_ES" = c(0.38, 0.34, 0.03, 0.34, 0.35, 0.15, 0.09, 0.28, 0.35, 0.67, 0.70),
                           "sample_sizes" = c(16, 20, 8, 26, 23, 10, 8, 32, 9, 11, 10),
                           "strength_level" = c(1.18, 1.48, 1.46, 1.35, 1.27, 1.14, 1.40, 1.54, 1.09, 1.26, 1.26))
View(previous_studies)

# plot of ES's frequencies in bar plot
ggplot(previous_studies, aes(x = BBPT_ES)) +
  geom_histogram(binwidth = 0.1, fill = 'darkgrey', color = 'black') +
  xlab('Effect sizes presented by Krzysztofik et al. (2021)') +
  ylab('Frequency') +
  ggtitle("Histogram of ES's") +
  theme_minimal() +  # Use a minimal theme for simplicity
  theme(axis.title = element_text(size = 14, face = 'bold'),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 16, face = 'bold'))

# rain cloud plot
# library(tidyquant)
# library(ggdist)
# library(ggthemes)

xlabel <- c(seq(0, 35, by = 5), 5, 8)
xlabel2 <- c(0,5,8,10,15,20,25,30,35)

previous_studies %>% 
  ggplot(aes(x = sample_sizes)) +
  # add half-violin from {ggdist} package
  stat_halfeye(
    # adjust bandwidth
    adjust = 0.5,
    # move to the right
    justification = + 0.05,
    # remove the slub interval
    .width = 0,
    point_colour = NA, color = "black", fill = "grey"
  ) +
  geom_boxplot(
    width = 0.12,
    # removing outliers
    #outlier.color = NA,
    alpha = 0.5, fill = "red"
  ) + 
  stat_boxplot(
    geom = "errorbar",
    width = 0.15
  ) +
  stat_dots(
    # ploting on left side
    side = "right",
    # adjusting position
    justification = 0.15,
    # adjust grouping (binning) of observations
    binwidth = 0.25, position = "dodgejust", color = "black", fill = "darkgrey"
  ) + scale_x_continuous(breaks = xlabel2, name = "Sample sizes") + scale_y_continuous(limits = c(-0.3, 1.1), expand = c(0,0)) + 
  theme_bw() + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), 
        axis.line.y = element_blank(), 
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
  geom_vline(xintercept = 30, linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = 8, linetype = "longdash", linewidth = 0.7) + 
  geom_vline(xintercept = 10, linetype = "dashed", linewidth = 0.7)
  
# generate data with the mean of ES's 
mean(previous_studies$BBPT_ES) # The mean of all ES's is the same as presented by Krzysztofik et al. (2021).
mean(previous_studies$sample_sizes) # The mean of sample sizes of the previous studies presented by Krzysztofik et al. (2021).
min(previous_studies$sample_sizes)
max(previous_studies$sample_sizes)
table(previous_studies$sample_sizes)

# library(modeest)
mlv(previous_studies$sample_sizes, method = "mfv")

# calculating the mean/median values of strength level
mean(previous_studies$strength_level)
median(previous_studies$strength_level)

########################################################################################################
#################################### data set based on DeBruine article ################################
########################################################################################################

########################## PPO #############################

# set fixed effect parameters
PPO_beta_0 <- 550 # intercept; i.e., the grand mean of PPO
# setting the slope parameter based on the Cohens d 0.3
# (M2 = 0.3 ? 50 + 550 = 15 + 550 = 565)...565-550 = 15
PPO_beta_1 <- 15 # slope; i.e, effect of category BASE/PAPE 

# set random effect parameters
PPO_tau_0 <- 50 # by-subject random intercept sd
PPO_omega_0 <- 0 # by-item random intercept sd (expect that BASE and PAPE have both the same sd)

# set more random effect and error parameters
PPO_tau_1 <- 20 # by-subject random slope sd
PPO_rho <- 0.2 # correlation between intercept and slope
PPO_sigma <- 20 # residual (error) sd (twice the size of the by-subject random intercept SD)

# simulating the sampling process

# set number of subjects and items
PPO_n_subj <- 30 # number of subjects
PPO_n_BASE <- 12 # number of (3*4 BASE) 
PPO_n_PAPE <- 12 # number of (3*4 PAPE) 

# simulate a sample of items
# total number of items = BASE + PAPE
PPO_items <- data.frame(
  PPO_item_id = seq_len(PPO_n_BASE + PPO_n_PAPE),
  PPO_state = rep(c("BASE", "PAPE"), c(PPO_n_BASE, PPO_n_PAPE)),
  PPO_O_0i = rnorm(PPO_n_BASE + PPO_n_PAPE, mean = 0, sd = PPO_omega_0)
)

# effect-code category
PPO_items$PPO_X_i <- dplyr::recode(PPO_items$PPO_state, "BASE" = -0.5, "PAPE" = +0.5)


# simulate a sample of subjects
# calculate random intercept / random slope covariance
PPO_covar <- PPO_rho * PPO_tau_0 * PPO_tau_1

# put values into variance-covariance matrix
PPO_cov_mx  <- matrix(
  c(PPO_tau_0^2, PPO_covar,
    PPO_covar,   PPO_tau_1^2),
  nrow = 2, byrow = TRUE)

# generate the by-subject random effects
PPO_subject_rfx <- MASS::mvrnorm(n = PPO_n_subj,
                             mu = c(PPO_T_0s = 0, PPO_T_1s = 0),
                             Sigma = PPO_cov_mx)

# combine with subject IDs
PPO_subjects <- data.frame(PPO_subj_id = seq_len(PPO_n_subj),
                           PPO_subject_rfx)

# simulate a sample of subjects
# sample from a multivariate random distribution 
PPO_subjects <- faux::rnorm_multi(
  n = PPO_n_subj, 
  mu = 0, # means for random effects are always 0
  sd = c(PPO_tau_0, PPO_tau_1), # set SDs
  r = PPO_rho, # set correlation, see ?faux::rnorm_multi
  varnames = c("PPO_T_0s", "PPO_T_1s")
)

# add subject IDs
PPO_subjects$PPO_subj_id <- seq_len(PPO_n_subj)


# cross subject and item IDs; add an error term
# nrow(.) is the number of rows in the table
PPO_trials <- crossing(PPO_subjects, PPO_items)  %>%
  mutate(PPO_e_si = rnorm(nrow(.), mean = 0, sd = PPO_sigma)) %>%
  dplyr::select(PPO_subj_id, PPO_item_id, PPO_state, PPO_X_i, everything())

PPO_trials <- PPO_trials[order(PPO_trials$PPO_subj_id, decreasing = F),]


# calculate the response variable
dat_sim_PPO <- PPO_trials %>%
  mutate(PPO = PPO_beta_0 + PPO_T_0s + PPO_O_0i + (PPO_beta_1 + PPO_T_1s) * PPO_X_i + PPO_e_si) %>%
  dplyr::select(PPO_subj_id, PPO_item_id, PPO_state, PPO_X_i, PPO)

View(dat_sim_PPO)
################### PVO ####################

# set fixed effect parameters
PVO_beta_0 <- 1.90 # intercept; i.e., the grand mean of PPO
# setting the slope parameter based on the Cohens d 0.3
# (M2 = 0.3 ? 0.1 + 1.90 = 0.03 + 1.9 = 1.93)...1.93-1.90 = 0.03
PVO_beta_1 <- 0.03 # slope; i.e, effect of category BASE/PAPE

# set random effect parameters
PVO_tau_0 <- 0.1 # by-subject random intercept sd
PVO_omega_0 <- 0 # by-item random intercept sd (expect that BASE and PAPE have both the same sd)

# set more random effect and error parameters
PVO_tau_1 <- 0.04 # by-subject random slope sd
PVO_rho <- 0.2 # correlation between intercept and slope
PVO_sigma <- 0.04 # residual (error) sd (twice the size of the by-subject random intercept SD)

# simulating the sampling process

# set number of subjects and items
PVO_n_subj <- 30 # number of subjects
PVO_n_BASE <- 12 # number of (3*4 BASE) 
PVO_n_PAPE <- 12 # number of (3*4 PAPE) 

# simulate a sample of items
# total number of items = BASE + PAPE
PVO_items <- data.frame(
  PVO_item_id = seq_len(PVO_n_BASE + PVO_n_PAPE),
  PVO_state = rep(c("BASE", "PAPE"), c(PVO_n_BASE, PVO_n_PAPE)),
  PVO_O_0i = rnorm(PVO_n_BASE + PVO_n_PAPE, mean = 0, sd = PVO_omega_0)
)

# effect-code category
PVO_items$PVO_X_i <- dplyr::recode(PVO_items$PVO_state, "BASE" = -0.5, "PAPE" = +0.5)


# simulate a sample of subjects
# calculate random intercept / random slope covariance
PVO_covar <- PVO_rho * PVO_tau_0 * PVO_tau_1

# put values into variance-covariance matrix
PVO_cov_mx  <- matrix(
  c(PVO_tau_0^2, PVO_covar,
    PVO_covar,   PVO_tau_1^2),
  nrow = 2, byrow = TRUE)

# generate the by-subject random effects
PVO_subject_rfx <- MASS::mvrnorm(n = PVO_n_subj,
                             mu = c(PVO_T_0s = 0, PVO_T_1s = 0),
                             Sigma = PVO_cov_mx)

# combine with subject IDs
PVO_subjects <- data.frame(PVO_subj_id = seq_len(PVO_n_subj),
                           PVO_subject_rfx)

# simulate a sample of subjects
# sample from a multivariate random distribution 
PVO_subjects <- faux::rnorm_multi(
  n = PVO_n_subj, 
  mu = 0, # means for random effects are always 0
  sd = c(PVO_tau_0, PVO_tau_1), # set SDs
  r = PVO_rho, # set correlation, see ?faux::rnorm_multi
  varnames = c("PVO_T_0s", "PVO_T_1s")
)

# add subject IDs
PVO_subjects$PVO_subj_id <- seq_len(PVO_n_subj)

# cross subject and item IDs; add an error term
# nrow(.) is the number of rows in the table
PVO_trials <- crossing(PVO_subjects, PVO_items)  %>%
  mutate(PVO_e_si = rnorm(nrow(.), mean = 0, sd = PVO_sigma)) %>%
  dplyr::select(PVO_subj_id, PVO_item_id, PVO_state, PVO_X_i, everything())

PVO_trials <- PVO_trials[order(PVO_trials$PVO_subj_id, decreasing = F),]


# calculate the response variable
dat_sim_PVO <- PVO_trials %>%
  mutate(PVO = PVO_beta_0 + PVO_T_0s + PVO_O_0i + (PVO_beta_1 + PVO_T_1s) * PVO_X_i + PVO_e_si) %>%
  dplyr::select(PVO_subj_id, PVO_item_id, PVO_state, PVO_X_i, PVO)

View(dat_sim_PVO)
############### connect that into one simulated dataset ##########################

simulated_dataset <- cbind(dat_sim_PPO[, c('PPO_subj_id', 'PPO_item_id','PPO_state','PPO')], dat_sim_PVO$PVO)

colnames(simulated_dataset)[colnames(simulated_dataset) == "dat_sim_PPO$PPO"] ="PPO"
colnames(simulated_dataset)[colnames(simulated_dataset) == "dat_sim_PVO$PVO"] ="PVO"
colnames(simulated_dataset)[colnames(simulated_dataset) == "PPO_state"] ="state"
colnames(simulated_dataset)[colnames(simulated_dataset) == "PPO_subj_id"] ="subj_ID"
colnames(simulated_dataset)[colnames(simulated_dataset) == "PPO_item_id"] ="item_ID"

simulated_dataset$attempt_in_time_point <- rep(1:3)
simulated_dataset$time_point <- rep(c(0, 5, 7, 10), each = 3)

unique_subj_ids <- unique(simulated_dataset$subj_ID)

age_values <- rnorm(length(unique_subj_ids), mean = 28.42, sd = 7.0)
simulated_dataset$age <- age_values[match(simulated_dataset$subj_ID, unique_subj_ids)]

weight_values <- rnorm(length(unique_subj_ids), mean = 75.00, sd = 4.0)
simulated_dataset$weight <- weight_values[match(simulated_dataset$subj_ID, unique_subj_ids)]

height_values <- rnorm(length(unique_subj_ids), mean = 181.00, sd = 6.0)
simulated_dataset$height <- height_values[match(simulated_dataset$subj_ID, unique_subj_ids)]

one_rm_values <- rnorm(length(unique_subj_ids), mean = 60.00, sd = 5.0)
simulated_dataset$one_rm <- one_rm_values[match(simulated_dataset$subj_ID, unique_subj_ids)]

sex_values <- sample(c("M", "W"), length(unique_subj_ids), replace = TRUE)
simulated_dataset$sex <- sex_values[match(simulated_dataset$subj_ID, unique_subj_ids)]


simulated_dataset <- simulated_dataset[, c(1, 2, 12, 8:11, 3, 6, 7, 4, 5)]

View(simulated_dataset)

################## wrangling simulated data set into wide format data sets ################

raw_data_PPO <- simulated_dataset %>%
  unite("state_attempt", state, attempt_in_time_point, sep = "_") %>%
  pivot_wider(
    id_cols = c(subj_ID, sex, age, weight, height, one_rm),
    names_from = c("state_attempt", "time_point"),
    names_glue = "PPO{state_attempt}_{time_point}",
    values_from = PPO 
  )

raw_data_PVO <- simulated_dataset %>%
  unite("state_attempt", state, attempt_in_time_point, sep = "_") %>%
  pivot_wider(
    id_cols = subj_ID,
    names_from = c("state_attempt", "time_point"),
    names_glue = "PVO{state_attempt}_{time_point}",
    values_from = PVO
  )

simulated_wide_dataset_1 <- cbind(raw_data_PPO, raw_data_PVO[, c(2:25)])

writexl::write_xlsx(simulated_wide_dataset_1, 'simulated_wide_dataset_1.xlsx')
writexl::write_xlsx(simulated_dataset, 'simulated_dataset.xlsx')

simulated_dataset_PAPE_ballistic <- readxl::read_excel("E:/data/Statistics/Data/simulated_dataset.xlsx")
simulated_dataset_PAPE_ballistic_wide <- readxl::read_excel("E:/data/Statistics/Data/simulated_wide_dataset_1.xlsx")
# OR
simulated_dataset_PAPE_ballistic <- readxl::read_excel("C:/Users/malir/OneDrive - Univerzita Karlova/Plocha/Data/simulated_dataset.xlsx")
simulated_dataset_PAPE_ballistic_wide <- readxl::read_excel("C:/Users/malir/OneDrive - Univerzita Karlova/Plocha/Data/simulated_wide_dataset_1.xlsx")


# ----------------------------------------------------------------------------------------------------------- #
# ----------------------------------------------- Descriptive statistics ------------------------------------ #
# ----------------------------------------------------------------------------------------------------------- #

descriptive <- round(pastecs::stat.desc(simulated_dataset_PAPE_ballistic_wide[, c(3:54)], p = 0.95, norm = T), 2)
descriptive <- as.data.frame(descriptive)

(descriptive_body <- descriptive[c(4:6, 8:10, 13), c(1:52)])
print(descriptive_body)

(descriptive_body <- as.data.frame(t(descriptive_body)))
descriptive_body$variable <- rownames(descriptive_body)
descriptive_body <- descriptive_body[, c(8, 1:7)]
print(descriptive_body)

writexl::write_xlsx(descriptive_body, "descriptives_PAPE_ballistic.xlsx")

# ------------------------------------------------------------------------------------------------------------ #
# ----------------------------------------- data visualization - rain cloud plots ---------------------------- #
# ------------------------------------------------------------------------------------------------------------ #

# ------------------------------------------ rain cloud plots based on test ---------------------------------- #


# simulated_dataset_PAPE_isometric$state <- as.factor(simulated_dataset_PAPE_isometric$state)

(simulated_dataset_PAPE_ballistic_plot_PPO <- simulated_dataset_PAPE_ballistic %>% 
    filter(state %in% c('BASE', 'PAPE')) %>% 
    ggplot(aes(x = factor(state), y = PPO, fill = factor(state))) +
    # add half-violin from {ggdist} package
    stat_halfeye(
      # adjust bandwidth
      adjust = 1.0,
      # move to the right
      justification = 0.0,
      # remove the slub interval
      .width = 0.0,
      point_colour = NA
    ) +
    geom_boxplot(
      width = 0.12,
      # removing outliers
      outlier.color = NA,
      alpha = 0.5, position = position_dodge(width = 2.5)
    ) +
    stat_dots(
      # ploting on left side
      side = "left",
      # adjusting position
      justification = 1.1,
      # adjust grouping (binning) of observations
      binwidth = NA
    )+ theme(axis.line = element_line(colour = "black")) + theme_bw() + xlab("") + ylab("Peak power output (W)") + 
    scale_fill_manual(values = c("BASE" = "#AF46B4", "PAPE" = "#33a02c")))


(simulated_dataset_PAPE_ballistic_plot_PVO <- simulated_dataset_PAPE_ballistic %>% 
    filter(state %in% c('BASE', 'PAPE')) %>% 
    ggplot(aes(x = factor(state), y = PVO, fill = factor(state))) +
    # add half-violin from {ggdist} package
    stat_halfeye(
      # adjust bandwidth
      adjust = 1.0,
      # move to the right
      justification = 0.0,
      # remove the slub interval
      .width = 0.0,
      point_colour = NA
    ) +
    geom_boxplot(
      width = 0.12,
      # removing outliers
      outlier.color = NA,
      alpha = 0.5, position = position_dodge(width = 2.5)
    ) +
    stat_dots(
      # ploting on left side
      side = "left",
      # adjusting position
      justification = 1.1,
      # adjust grouping (binning) of observations
      binwidth = NA
    )+ theme(axis.line = element_line(colour = "black")) + theme_bw() + xlab("") + ylab("Peak velocity output (m/s)") + 
    scale_fill_manual(values = c("BASE" = "#AF46B4", "PAPE" = "#33a02c")))


grid.arrange(simulated_dataset_PAPE_ballistic_plot_PPO,
             simulated_dataset_PAPE_ballistic_plot_PVO, ncol=2)

# ------------------------------------ rain cloud plots based on time points -------------------------------- #

(simulated_dataset_PAPE_ballistic_plot_TP_PPO_0 <- simulated_dataset_PAPE_ballistic %>% 
    filter(state %in% c('BASE', 'PAPE'), time_point == 0) %>% 
    ggplot(aes(x = factor(state), y = PPO, fill = factor(state))) +
    # add half-violin from {ggdist} package
    stat_halfeye(
      # adjust bandwidth
      adjust = 1.0,
      # move to the right
      justification = 0.0,
      # remove the slub interval
      .width = 0.0,
      point_colour = NA
    ) +
    geom_boxplot(
      width = 0.12,
      # removing outliers
      outlier.color = NA,
      alpha = 0.5, position = position_dodge(width = 2.5)
    ) +
    stat_dots(
      # ploting on left side
      side = "left",
      # adjusting position
      justification = 1.1,
      # adjust grouping (binning) of observations
      binwidth = NA
    )+ theme(axis.line = element_line(colour = "black")) + theme_bw() + xlab("PPO_0") + ylab("Peak power output (W)") + 
    scale_fill_manual(values = c("BASE" = "#AF46B4", "PAPE" = "#33a02c")))


(simulated_dataset_PAPE_ballistic_plot_TP_PPO_5 <- simulated_dataset_PAPE_ballistic %>% 
    filter(state %in% c('BASE', 'PAPE'), time_point == 5) %>% 
    ggplot(aes(x = factor(state), y = PPO, fill = factor(state))) +
    # add half-violin from {ggdist} package
    stat_halfeye(
      # adjust bandwidth
      adjust = 1.0,
      # move to the right
      justification = 0.0,
      # remove the slub interval
      .width = 0.0,
      point_colour = NA
    ) +
    geom_boxplot(
      width = 0.12,
      # removing outliers
      outlier.color = NA,
      alpha = 0.5, position = position_dodge(width = 2.5)
    ) +
    stat_dots(
      # ploting on left side
      side = "left",
      # adjusting position
      justification = 1.1,
      # adjust grouping (binning) of observations
      binwidth = NA
    )+ theme(axis.line = element_line(colour = "black")) + theme_bw() + xlab("PPO_5") + ylab("Peak power output (W)") + 
    scale_fill_manual(values = c("BASE" = "#AF46B4", "PAPE" = "#33a02c")))


(simulated_dataset_PAPE_ballistic_plot_TP_PPO_7 <- simulated_dataset_PAPE_ballistic %>% 
    filter(state %in% c('BASE', 'PAPE'), time_point == 7) %>% 
    ggplot(aes(x = factor(state), y = PPO, fill = factor(state))) +
    # add half-violin from {ggdist} package
    stat_halfeye(
      # adjust bandwidth
      adjust = 1.0,
      # move to the right
      justification = 0.0,
      # remove the slub interval
      .width = 0.0,
      point_colour = NA
    ) +
    geom_boxplot(
      width = 0.12,
      # removing outliers
      outlier.color = NA,
      alpha = 0.5, position = position_dodge(width = 2.5)
    ) +
    stat_dots(
      # ploting on left side
      side = "left",
      # adjusting position
      justification = 1.1,
      # adjust grouping (binning) of observations
      binwidth = NA
    )+ theme(axis.line = element_line(colour = "black")) + theme_bw() + xlab("PPO_7") + ylab("Peak power output (W)") + 
    scale_fill_manual(values = c("BASE" = "#AF46B4", "PAPE" = "#33a02c")))


(simulated_dataset_PAPE_ballistic_plot_TP_PPO_10 <- simulated_dataset_PAPE_ballistic %>% 
    filter(state %in% c('BASE', 'PAPE'), time_point == 10) %>% 
    ggplot(aes(x = factor(state), y = PPO, fill = factor(state))) +
    # add half-violin from {ggdist} package
    stat_halfeye(
      # adjust bandwidth
      adjust = 1.0,
      # move to the right
      justification = 0.0,
      # remove the slub interval
      .width = 0.0,
      point_colour = NA
    ) +
    geom_boxplot(
      width = 0.12,
      # removing outliers
      outlier.color = NA,
      alpha = 0.5, position = position_dodge(width = 2.5)
    ) +
    stat_dots(
      # ploting on left side
      side = "left",
      # adjusting position
      justification = 1.1,
      # adjust grouping (binning) of observations
      binwidth = NA
    )+ theme(axis.line = element_line(colour = "black")) + theme_bw() + xlab("PPO_10") + ylab("Peak power output (W)") + 
    scale_fill_manual(values = c("BASE" = "#AF46B4", "PAPE" = "#33a02c")))


(simulated_dataset_PAPE_ballistic_plot_TP_PVO_0 <- simulated_dataset_PAPE_ballistic %>% 
    filter(state %in% c('BASE', 'PAPE'), time_point == 0) %>% 
    ggplot(aes(x = factor(state), y = PVO, fill = factor(state))) +
    # add half-violin from {ggdist} package
    stat_halfeye(
      # adjust bandwidth
      adjust = 1.0,
      # move to the right
      justification = 0.0,
      # remove the slub interval
      .width = 0.0,
      point_colour = NA
    ) +
    geom_boxplot(
      width = 0.12,
      # removing outliers
      outlier.color = NA,
      alpha = 0.5, position = position_dodge(width = 2.5)
    ) +
    stat_dots(
      # ploting on left side
      side = "left",
      # adjusting position
      justification = 1.1,
      # adjust grouping (binning) of observations
      binwidth = NA
    )+ theme(axis.line = element_line(colour = "black")) + theme_bw() + xlab("PVO_0") + ylab("Peak velocity output (m/s)") + 
    scale_fill_manual(values = c("BASE" = "#AF46B4", "PAPE" = "#33a02c")))


(simulated_dataset_PAPE_ballistic_plot_TP_PVO_5 <- simulated_dataset_PAPE_ballistic %>% 
    filter(state %in% c('BASE', 'PAPE'), time_point == 5) %>% 
    ggplot(aes(x = factor(state), y = PVO, fill = factor(state))) +
    # add half-violin from {ggdist} package
    stat_halfeye(
      # adjust bandwidth
      adjust = 1.0,
      # move to the right
      justification = 0.0,
      # remove the slub interval
      .width = 0.0,
      point_colour = NA
    ) +
    geom_boxplot(
      width = 0.12,
      # removing outliers
      outlier.color = NA,
      alpha = 0.5, position = position_dodge(width = 2.5)
    ) +
    stat_dots(
      # ploting on left side
      side = "left",
      # adjusting position
      justification = 1.1,
      # adjust grouping (binning) of observations
      binwidth = NA
    )+ theme(axis.line = element_line(colour = "black")) + theme_bw() + xlab("PVO_5") + ylab("Peak velocity output (m/s)") + 
    scale_fill_manual(values = c("BASE" = "#AF46B4", "PAPE" = "#33a02c")))


(simulated_dataset_PAPE_ballistic_plot_TP_PVO_7 <- simulated_dataset_PAPE_ballistic %>% 
    filter(state %in% c('BASE', 'PAPE'), time_point == 7) %>% 
    ggplot(aes(x = factor(state), y = PVO, fill = factor(state))) +
    # add half-violin from {ggdist} package
    stat_halfeye(
      # adjust bandwidth
      adjust = 1.0,
      # move to the right
      justification = 0.0,
      # remove the slub interval
      .width = 0.0,
      point_colour = NA
    ) +
    geom_boxplot(
      width = 0.12,
      # removing outliers
      outlier.color = NA,
      alpha = 0.5, position = position_dodge(width = 2.5)
    ) +
    stat_dots(
      # ploting on left side
      side = "left",
      # adjusting position
      justification = 1.1,
      # adjust grouping (binning) of observations
      binwidth = NA
    )+ theme(axis.line = element_line(colour = "black")) + theme_bw() + xlab("PVO_7") + ylab("Peak velocity output (m/s)") + 
    scale_fill_manual(values = c("BASE" = "#AF46B4", "PAPE" = "#33a02c")))


(simulated_dataset_PAPE_ballistic_plot_TP_PVO_10 <- simulated_dataset_PAPE_ballistic %>% 
    filter(state %in% c('BASE', 'PAPE'), time_point == 10) %>% 
    ggplot(aes(x = factor(state), y = PVO, fill = factor(state))) +
    # add half-violin from {ggdist} package
    stat_halfeye(
      # adjust bandwidth
      adjust = 1.0,
      # move to the right
      justification = 0.0,
      # remove the slub interval
      .width = 0.0,
      point_colour = NA
    ) +
    geom_boxplot(
      width = 0.12,
      # removing outliers
      outlier.color = NA,
      alpha = 0.5, position = position_dodge(width = 2.5)
    ) +
    stat_dots(
      # ploting on left side
      side = "left",
      # adjusting position
      justification = 1.1,
      # adjust grouping (binning) of observations
      binwidth = NA
    )+ theme(axis.line = element_line(colour = "black")) + theme_bw() + xlab("PVO_10") + ylab("Peak velocity output (m/s)") + 
    scale_fill_manual(values = c("BASE" = "#AF46B4", "PAPE" = "#33a02c")))


grid.arrange(simulated_dataset_PAPE_ballistic_plot_TP_PPO_0,
             simulated_dataset_PAPE_ballistic_plot_TP_PPO_5,
             simulated_dataset_PAPE_ballistic_plot_TP_PPO_7,
             simulated_dataset_PAPE_ballistic_plot_TP_PPO_10,
             simulated_dataset_PAPE_ballistic_plot_TP_PVO_0,
             simulated_dataset_PAPE_ballistic_plot_TP_PVO_5,
             simulated_dataset_PAPE_ballistic_plot_TP_PVO_7,
             simulated_dataset_PAPE_ballistic_plot_TP_PVO_10, ncol=2)


# ------------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------ ICC ------------------------------------------------- #
# ------------------------------------------------------------------------------------------------------------ #

########################################### ICC - among attempts in time points ###################

# PPO - ICC among attempts in time points - BASE
ICC_PPO_BASE_0 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PPOBASE_1_0', 'PPOBASE_2_0', 'PPOBASE_3_0')],
                           model = "twoway", type = "consistency", unit = "average")
ICC_PPO_BASE_0_summary <- cbind("ICC value" = round(ICC_PPO_BASE_0$value, 3),
                                "Lower 95%CI" = round(ICC_PPO_BASE_0$lbound, 3),
                                "Upper 95%CI" = round(ICC_PPO_BASE_0$ubound, 3),
                                "p-value" = ICC_PPO_BASE_0$p.value)
rownames(ICC_PPO_BASE_0_summary) <- "ICC_PPO_BASE_0_summary"

ICC_PPO_BASE_5 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PPOBASE_1_5', 'PPOBASE_2_5', 'PPOBASE_3_5')],
                           model = "twoway", type = "consistency", unit = "average")
ICC_PPO_BASE_5_summary <- cbind("ICC value" = round(ICC_PPO_BASE_5$value, 3),
                                "Lower 95%CI" = round(ICC_PPO_BASE_5$lbound, 3),
                                "Upper 95%CI" = round(ICC_PPO_BASE_5$ubound, 3),
                                "p-value" = ICC_PPO_BASE_5$p.value)
rownames(ICC_PPO_BASE_5_summary) <- "ICC_PPO_BASE_5_summary"

ICC_PPO_BASE_7 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PPOBASE_1_7', 'PPOBASE_2_7', 'PPOBASE_3_7')],
                           model = "twoway", type = "consistency", unit = "average")
ICC_PPO_BASE_7_summary <- cbind("ICC value" = round(ICC_PPO_BASE_7$value, 3),
                                "Lower 95%CI" = round(ICC_PPO_BASE_7$lbound, 3),
                                "Upper 95%CI" = round(ICC_PPO_BASE_7$ubound, 3),
                                "p-value" = ICC_PPO_BASE_7$p.value)
rownames(ICC_PPO_BASE_7_summary) <- "ICC_PPO_BASE_7_summary"

ICC_PPO_BASE_10 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PPOBASE_1_10', 'PPOBASE_2_10', 'PPOBASE_3_10')],
                           model = "twoway", type = "consistency", unit = "average")
ICC_PPO_BASE_10_summary <- cbind("ICC value" = round(ICC_PPO_BASE_10$value, 3),
                                "Lower 95%CI" = round(ICC_PPO_BASE_10$lbound, 3),
                                "Upper 95%CI" = round(ICC_PPO_BASE_10$ubound, 3),
                                "p-value" = ICC_PPO_BASE_10$p.value)
rownames(ICC_PPO_BASE_0_summary) <- "ICC_PPO_BASE_10_summary"


# PPO - ICC among attempts in time points - PAPE
ICC_PPO_PAPE_0 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PPOPAPE_1_0', 'PPOPAPE_2_0', 'PPOPAPE_3_0')],
                           model = "twoway", type = "consistency", unit = "average")
ICC_PPO_PAPE_0_summary <- cbind("ICC value" = round(ICC_PPO_PAPE_0$value, 3),
                                "Lower 95%CI" = round(ICC_PPO_PAPE_0$lbound, 3),
                                "Upper 95%CI" = round(ICC_PPO_PAPE_0$ubound, 3),
                                "p-value" = ICC_PPO_PAPE_0$p.value)
rownames(ICC_PPO_PAPE_0_summary) <- "ICC_PPO_PAPE_0_summary"

ICC_PPO_PAPE_5 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PPOPAPE_1_5', 'PPOPAPE_2_5', 'PPOPAPE_3_5')],
                           model = "twoway", type = "consistency", unit = "average")
ICC_PPO_PAPE_5_summary <- cbind("ICC value" = round(ICC_PPO_PAPE_5$value, 3),
                                "Lower 95%CI" = round(ICC_PPO_PAPE_5$lbound, 3),
                                "Upper 95%CI" = round(ICC_PPO_PAPE_5$ubound, 3),
                                "p-value" = ICC_PPO_PAPE_5$p.value)
rownames(ICC_PPO_PAPE_5_summary) <- "ICC_PPO_PAPE_5_summary"

ICC_PPO_PAPE_7 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PPOPAPE_1_7', 'PPOPAPE_2_7', 'PPOPAPE_3_7')],
                           model = "twoway", type = "consistency", unit = "average")
ICC_PPO_PAPE_7_summary <- cbind("ICC value" = round(ICC_PPO_PAPE_7$value, 3),
                                "Lower 95%CI" = round(ICC_PPO_PAPE_7$lbound, 3),
                                "Upper 95%CI" = round(ICC_PPO_PAPE_7$ubound, 3),
                                "p-value" = ICC_PPO_PAPE_7$p.value)
rownames(ICC_PPO_PAPE_7_summary) <- "ICC_PPO_PAPE_7_summary"

ICC_PPO_PAPE_10 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PPOPAPE_1_10', 'PPOPAPE_2_10', 'PPOPAPE_3_10')],
                           model = "twoway", type = "consistency", unit = "average")
ICC_PPO_PAPE_10_summary <- cbind("ICC value" = round(ICC_PPO_PAPE_10$value, 3),
                                "Lower 95%CI" = round(ICC_PPO_PAPE_10$lbound, 3),
                                "Upper 95%CI" = round(ICC_PPO_PAPE_10$ubound, 3),
                                "p-value" = ICC_PPO_PAPE_10$p.value)
rownames(ICC_PPO_PAPE_10_summary) <- "ICC_PPO_PAPE_10_summary"


# PVO - ICC among attempts in time points - BASE
ICC_PVO_BASE_0 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PVOBASE_1_0', 'PVOBASE_2_0', 'PVOBASE_3_0')],
                           model = "twoway", type = "consistency", unit = "average")
ICC_PVO_BASE_0_summary <- cbind("ICC value" = round(ICC_PVO_BASE_0$value, 3),
                                "Lower 95%CI" = round(ICC_PVO_BASE_0$lbound, 3),
                                "Upper 95%CI" = round(ICC_PVO_BASE_0$ubound, 3),
                                "p-value" = ICC_PVO_BASE_0$p.value)
rownames(ICC_PVO_BASE_0_summary) <- "ICC_PVO_BASE_0_summary"

ICC_PVO_BASE_5 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PVOBASE_1_5', 'PVOBASE_2_5', 'PVOBASE_3_5')],
                           model = "twoway", type = "consistency", unit = "average")
ICC_PVO_BASE_5_summary <- cbind("ICC value" = round(ICC_PVO_BASE_5$value, 3),
                                "Lower 95%CI" = round(ICC_PVO_BASE_5$lbound, 3),
                                "Upper 95%CI" = round(ICC_PVO_BASE_5$ubound, 3),
                                "p-value" = ICC_PVO_BASE_5$p.value)
rownames(ICC_PVO_BASE_5_summary) <- "ICC_PVO_BASE_5_summary"

ICC_PVO_BASE_7 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PVOBASE_1_7', 'PVOBASE_2_7', 'PVOBASE_3_7')],
                           model = "twoway", type = "consistency", unit = "average")
ICC_PVO_BASE_7_summary <- cbind("ICC value" = round(ICC_PVO_BASE_7$value, 3),
                                "Lower 95%CI" = round(ICC_PVO_BASE_7$lbound, 3),
                                "Upper 95%CI" = round(ICC_PVO_BASE_7$ubound, 3),
                                "p-value" = ICC_PVO_BASE_7$p.value)
rownames(ICC_PVO_BASE_7_summary) <- "ICC_PVO_BASE_7_summary"

ICC_PVO_BASE_10 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PVOBASE_1_10', 'PVOBASE_2_10', 'PVOBASE_3_10')],
                            model = "twoway", type = "consistency", unit = "average")
ICC_PVO_BASE_10_summary <- cbind("ICC value" = round(ICC_PVO_BASE_10$value, 3),
                                 "Lower 95%CI" = round(ICC_PVO_BASE_10$lbound, 3),
                                 "Upper 95%CI" = round(ICC_PVO_BASE_10$ubound, 3),
                                 "p-value" = ICC_PVO_BASE_10$p.value)
rownames(ICC_PVO_BASE_0_summary) <- "ICC_PVO_BASE_10_summary"


# PVO - ICC among attempts in time points - PAPE
ICC_PVO_PAPE_0 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PVOPAPE_1_0', 'PVOPAPE_2_0', 'PVOPAPE_3_0')],
                           model = "twoway", type = "consistency", unit = "average")
ICC_PVO_PAPE_0_summary <- cbind("ICC value" = round(ICC_PVO_PAPE_0$value, 3),
                                "Lower 95%CI" = round(ICC_PVO_PAPE_0$lbound, 3),
                                "Upper 95%CI" = round(ICC_PVO_PAPE_0$ubound, 3),
                                "p-value" = ICC_PVO_PAPE_0$p.value)
rownames(ICC_PVO_PAPE_0_summary) <- "ICC_PVO_PAPE_0_summary"

ICC_PVO_PAPE_5 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PVOPAPE_1_5', 'PVOPAPE_2_5', 'PVOPAPE_3_5')],
                           model = "twoway", type = "consistency", unit = "average")
ICC_PVO_PAPE_5_summary <- cbind("ICC value" = round(ICC_PVO_PAPE_5$value, 3),
                                "Lower 95%CI" = round(ICC_PVO_PAPE_5$lbound, 3),
                                "Upper 95%CI" = round(ICC_PVO_PAPE_5$ubound, 3),
                                "p-value" = ICC_PVO_PAPE_5$p.value)
rownames(ICC_PVO_PAPE_5_summary) <- "ICC_PVO_PAPE_5_summary"

ICC_PVO_PAPE_7 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PVOPAPE_1_7', 'PVOPAPE_2_7', 'PVOPAPE_3_7')],
                           model = "twoway", type = "consistency", unit = "average")
ICC_PVO_PAPE_7_summary <- cbind("ICC value" = round(ICC_PVO_PAPE_7$value, 3),
                                "Lower 95%CI" = round(ICC_PVO_PAPE_7$lbound, 3),
                                "Upper 95%CI" = round(ICC_PVO_PAPE_7$ubound, 3),
                                "p-value" = ICC_PVO_PAPE_7$p.value)
rownames(ICC_PVO_PAPE_7_summary) <- "ICC_PVO_PAPE_7_summary"

ICC_PVO_PAPE_10 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PVOPAPE_1_10', 'PVOPAPE_2_10', 'PVOPAPE_3_10')],
                            model = "twoway", type = "consistency", unit = "average")
ICC_PVO_PAPE_10_summary <- cbind("ICC value" = round(ICC_PVO_PAPE_10$value, 3),
                                 "Lower 95%CI" = round(ICC_PVO_PAPE_10$lbound, 3),
                                 "Upper 95%CI" = round(ICC_PVO_PAPE_10$ubound, 3),
                                 "p-value" = ICC_PVO_PAPE_10$p.value)
rownames(ICC_PVO_PAPE_10_summary) <- "ICC_PVO_PAPE_10_summary"


########################################### ICC - between BASE and PAPE #######################

# PPO - ICC between BASE and PAPE
# PPO 2 BASE - PAPE
ICC_PPO_BASE_vs_PAPE_1_0 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PPOBASE_1_0', 'PPOPAPE_1_0')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_PPO_BASE_vs_PAPE_1_0_summary <- cbind("ICC value" = round(ICC_PPO_BASE_vs_PAPE_1_0$value, 3),
                                          "Lower 95%CI" = round(ICC_PPO_BASE_vs_PAPE_1_0$lbound, 3),
                                          "Upper 95%CI" = round(ICC_PPO_BASE_vs_PAPE_1_0$ubound, 3),
                                          "p-value" = ICC_PPO_BASE_vs_PAPE_1_0$p.value)
rownames(ICC_PPO_BASE_vs_PAPE_1_0_summary) <- "ICC_PPO_BASE_vs_PAPE_1_0_summary"

ICC_PPO_BASE_vs_PAPE_2_0 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PPOBASE_2_0', 'PPOPAPE_2_0')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_PPO_BASE_vs_PAPE_2_0_summary <- cbind("ICC value" = round(ICC_PPO_BASE_vs_PAPE_2_0$value, 3),
                                          "Lower 95%CI" = round(ICC_PPO_BASE_vs_PAPE_2_0$lbound, 3),
                                          "Upper 95%CI" = round(ICC_PPO_BASE_vs_PAPE_2_0$ubound, 3),
                                          "p-value" = ICC_PPO_BASE_vs_PAPE_2_0$p.value)
rownames(ICC_PPO_BASE_vs_PAPE_2_0_summary) <- "ICC_PPO_BASE_vs_PAPE_2_0_summary"

ICC_PPO_BASE_vs_PAPE_3_0 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PPOBASE_3_0', 'PPOPAPE_3_0')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_PPO_BASE_vs_PAPE_3_0_summary <- cbind("ICC value" = round(ICC_PPO_BASE_vs_PAPE_3_0$value, 3),
                                          "Lower 95%CI" = round(ICC_PPO_BASE_vs_PAPE_3_0$lbound, 3),
                                          "Upper 95%CI" = round(ICC_PPO_BASE_vs_PAPE_3_0$ubound, 3),
                                          "p-value" = ICC_PPO_BASE_vs_PAPE_3_0$p.value)
rownames(ICC_PPO_BASE_vs_PAPE_3_0_summary) <- "ICC_PPO_BASE_vs_PAPE_3_0_summary"


# PPO 5 BASE - PAPE
ICC_PPO_BASE_vs_PAPE_1_5 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PPOBASE_1_5', 'PPOPAPE_1_5')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_PPO_BASE_vs_PAPE_1_5_summary <- cbind("ICC value" = round(ICC_PPO_BASE_vs_PAPE_1_5$value, 3),
                                          "Lower 95%CI" = round(ICC_PPO_BASE_vs_PAPE_1_5$lbound, 3),
                                          "Upper 95%CI" = round(ICC_PPO_BASE_vs_PAPE_1_5$ubound, 3),
                                          "p-value" = ICC_PPO_BASE_vs_PAPE_1_5$p.value)
rownames(ICC_PPO_BASE_vs_PAPE_1_5_summary) <- "ICC_PPO_BASE_vs_PAPE_1_5_summary"

ICC_PPO_BASE_vs_PAPE_2_5 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PPOBASE_2_5', 'PPOPAPE_2_5')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_PPO_BASE_vs_PAPE_2_5_summary <- cbind("ICC value" = round(ICC_PPO_BASE_vs_PAPE_2_5$value, 3),
                                          "Lower 95%CI" = round(ICC_PPO_BASE_vs_PAPE_2_5$lbound, 3),
                                          "Upper 95%CI" = round(ICC_PPO_BASE_vs_PAPE_2_5$ubound, 3),
                                          "p-value" = ICC_PPO_BASE_vs_PAPE_2_5$p.value)
rownames(ICC_PPO_BASE_vs_PAPE_2_5_summary) <- "ICC_PPO_BASE_vs_PAPE_2_5_summary"

ICC_PPO_BASE_vs_PAPE_3_5 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PPOBASE_3_5', 'PPOPAPE_3_5')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_PPO_BASE_vs_PAPE_3_5_summary <- cbind("ICC value" = round(ICC_PPO_BASE_vs_PAPE_3_5$value, 3),
                                          "Lower 95%CI" = round(ICC_PPO_BASE_vs_PAPE_3_5$lbound, 3),
                                          "Upper 95%CI" = round(ICC_PPO_BASE_vs_PAPE_3_5$ubound, 3),
                                          "p-value" = ICC_PPO_BASE_vs_PAPE_3_5$p.value)
rownames(ICC_PPO_BASE_vs_PAPE_3_5_summary) <- "ICC_PPO_BASE_vs_PAPE_3_5_summary"

# PPO 7 BASE - PAPE
ICC_PPO_BASE_vs_PAPE_1_7 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PPOBASE_1_7', 'PPOPAPE_1_7')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_PPO_BASE_vs_PAPE_1_7_summary <- cbind("ICC value" = round(ICC_PPO_BASE_vs_PAPE_1_7$value, 3),
                                          "Lower 95%CI" = round(ICC_PPO_BASE_vs_PAPE_1_7$lbound, 3),
                                          "Upper 95%CI" = round(ICC_PPO_BASE_vs_PAPE_1_7$ubound, 3),
                                          "p-value" = ICC_PPO_BASE_vs_PAPE_1_7$p.value)
rownames(ICC_PPO_BASE_vs_PAPE_1_7_summary) <- "ICC_PPO_BASE_vs_PAPE_1_7_summary"

ICC_PPO_BASE_vs_PAPE_2_7 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PPOBASE_2_7', 'PPOPAPE_2_7')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_PPO_BASE_vs_PAPE_2_7_summary <- cbind("ICC value" = round(ICC_PPO_BASE_vs_PAPE_2_7$value, 3),
                                          "Lower 95%CI" = round(ICC_PPO_BASE_vs_PAPE_2_7$lbound, 3),
                                          "Upper 95%CI" = round(ICC_PPO_BASE_vs_PAPE_2_7$ubound, 3),
                                          "p-value" = ICC_PPO_BASE_vs_PAPE_2_7$p.value)
rownames(ICC_PPO_BASE_vs_PAPE_2_7_summary) <- "ICC_PPO_BASE_vs_PAPE_2_7_summary"

ICC_PPO_BASE_vs_PAPE_3_7 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PPOBASE_3_7', 'PPOPAPE_3_7')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_PPO_BASE_vs_PAPE_3_7_summary <- cbind("ICC value" = round(ICC_PPO_BASE_vs_PAPE_3_7$value, 3),
                                          "Lower 95%CI" = round(ICC_PPO_BASE_vs_PAPE_3_7$lbound, 3),
                                          "Upper 95%CI" = round(ICC_PPO_BASE_vs_PAPE_3_7$ubound, 3),
                                          "p-value" = ICC_PPO_BASE_vs_PAPE_3_7$p.value)
rownames(ICC_PPO_BASE_vs_PAPE_3_7_summary) <- "ICC_PPO_BASE_vs_PAPE_3_7_summary"

# PPO 10 BASE - PAPE
ICC_PPO_BASE_vs_PAPE_1_10 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PPOBASE_1_10', 'PPOPAPE_1_10')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_PPO_BASE_vs_PAPE_1_10_summary <- cbind("ICC value" = round(ICC_PPO_BASE_vs_PAPE_1_10$value, 3),
                                          "Lower 95%CI" = round(ICC_PPO_BASE_vs_PAPE_1_10$lbound, 3),
                                          "Upper 95%CI" = round(ICC_PPO_BASE_vs_PAPE_1_10$ubound, 3),
                                          "p-value" = ICC_PPO_BASE_vs_PAPE_1_10$p.value)
rownames(ICC_PPO_BASE_vs_PAPE_1_10_summary) <- "ICC_PPO_BASE_vs_PAPE_1_10_summary"

ICC_PPO_BASE_vs_PAPE_2_10 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PPOBASE_2_10', 'PPOPAPE_2_10')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_PPO_BASE_vs_PAPE_2_10_summary <- cbind("ICC value" = round(ICC_PPO_BASE_vs_PAPE_2_10$value, 3),
                                          "Lower 95%CI" = round(ICC_PPO_BASE_vs_PAPE_2_10$lbound, 3),
                                          "Upper 95%CI" = round(ICC_PPO_BASE_vs_PAPE_2_10$ubound, 3),
                                          "p-value" = ICC_PPO_BASE_vs_PAPE_2_0$p.value)
rownames(ICC_PPO_BASE_vs_PAPE_2_10_summary) <- "ICC_PPO_BASE_vs_PAPE_2_10_summary"

ICC_PPO_BASE_vs_PAPE_3_10 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PPOBASE_3_10', 'PPOPAPE_3_10')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_PPO_BASE_vs_PAPE_3_10_summary <- cbind("ICC value" = round(ICC_PPO_BASE_vs_PAPE_3_10$value, 3),
                                          "Lower 95%CI" = round(ICC_PPO_BASE_vs_PAPE_3_10$lbound, 3),
                                          "Upper 95%CI" = round(ICC_PPO_BASE_vs_PAPE_3_10$ubound, 3),
                                          "p-value" = ICC_PPO_BASE_vs_PAPE_3_10$p.value)
rownames(ICC_PPO_BASE_vs_PAPE_3_10_summary) <- "ICC_PPO_BASE_vs_PAPE_3_10_summary"

# PVO - ICC between BASE and PAPE
# PVO 2 BASE - PAPE
ICC_PVO_BASE_vs_PAPE_1_0 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PVOBASE_1_0', 'PVOPAPE_1_0')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_PVO_BASE_vs_PAPE_1_0_summary <- cbind("ICC value" = round(ICC_PVO_BASE_vs_PAPE_1_0$value, 3),
                                          "Lower 95%CI" = round(ICC_PVO_BASE_vs_PAPE_1_0$lbound, 3),
                                          "Upper 95%CI" = round(ICC_PVO_BASE_vs_PAPE_1_0$ubound, 3),
                                          "p-value" = ICC_PVO_BASE_vs_PAPE_1_0$p.value)
rownames(ICC_PVO_BASE_vs_PAPE_1_0_summary) <- "ICC_PVO_BASE_vs_PAPE_1_0_summary"

ICC_PVO_BASE_vs_PAPE_2_0 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PVOBASE_2_0', 'PVOPAPE_2_0')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_PVO_BASE_vs_PAPE_2_0_summary <- cbind("ICC value" = round(ICC_PVO_BASE_vs_PAPE_2_0$value, 3),
                                          "Lower 95%CI" = round(ICC_PVO_BASE_vs_PAPE_2_0$lbound, 3),
                                          "Upper 95%CI" = round(ICC_PVO_BASE_vs_PAPE_2_0$ubound, 3),
                                          "p-value" = ICC_PVO_BASE_vs_PAPE_2_0$p.value)
rownames(ICC_PVO_BASE_vs_PAPE_2_0_summary) <- "ICC_PVO_BASE_vs_PAPE_2_0_summary"

ICC_PVO_BASE_vs_PAPE_3_0 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PVOBASE_3_0', 'PVOPAPE_3_0')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_PVO_BASE_vs_PAPE_3_0_summary <- cbind("ICC value" = round(ICC_PVO_BASE_vs_PAPE_3_0$value, 3),
                                          "Lower 95%CI" = round(ICC_PVO_BASE_vs_PAPE_3_0$lbound, 3),
                                          "Upper 95%CI" = round(ICC_PVO_BASE_vs_PAPE_3_0$ubound, 3),
                                          "p-value" = ICC_PVO_BASE_vs_PAPE_3_0$p.value)
rownames(ICC_PVO_BASE_vs_PAPE_3_0_summary) <- "ICC_PVO_BASE_vs_PAPE_3_0_summary"

# PPO 5 BASE - PAPE
ICC_PVO_BASE_vs_PAPE_1_5 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PVOBASE_1_5', 'PVOPAPE_1_5')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_PVO_BASE_vs_PAPE_1_5_summary <- cbind("ICC value" = round(ICC_PVO_BASE_vs_PAPE_1_5$value, 3),
                                          "Lower 95%CI" = round(ICC_PVO_BASE_vs_PAPE_1_5$lbound, 3),
                                          "Upper 95%CI" = round(ICC_PVO_BASE_vs_PAPE_1_5$ubound, 3),
                                          "p-value" = ICC_PVO_BASE_vs_PAPE_1_5$p.value)
rownames(ICC_PVO_BASE_vs_PAPE_1_5_summary) <- "ICC_PVO_BASE_vs_PAPE_1_5_summary"

ICC_PVO_BASE_vs_PAPE_2_5 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PVOBASE_2_5', 'PVOPAPE_2_5')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_PVO_BASE_vs_PAPE_2_5_summary <- cbind("ICC value" = round(ICC_PVO_BASE_vs_PAPE_2_5$value, 3),
                                          "Lower 95%CI" = round(ICC_PVO_BASE_vs_PAPE_2_5$lbound, 3),
                                          "Upper 95%CI" = round(ICC_PVO_BASE_vs_PAPE_2_5$ubound, 3),
                                          "p-value" = ICC_PVO_BASE_vs_PAPE_2_5$p.value)
rownames(ICC_PVO_BASE_vs_PAPE_2_5_summary) <- "ICC_PVO_BASE_vs_PAPE_2_5_summary"

ICC_PVO_BASE_vs_PAPE_3_5 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PVOBASE_3_5', 'PVOPAPE_3_5')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_PVO_BASE_vs_PAPE_3_5_summary <- cbind("ICC value" = round(ICC_PVO_BASE_vs_PAPE_3_5$value, 3),
                                          "Lower 95%CI" = round(ICC_PVO_BASE_vs_PAPE_3_5$lbound, 3),
                                          "Upper 95%CI" = round(ICC_PVO_BASE_vs_PAPE_3_5$ubound, 3),
                                          "p-value" = ICC_PVO_BASE_vs_PAPE_3_5$p.value)
rownames(ICC_PVO_BASE_vs_PAPE_3_5_summary) <- "ICC_PVO_BASE_vs_PAPE_3_5_summary"

# PPO 7 BASE - PAPE
ICC_PVO_BASE_vs_PAPE_1_7 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PVOBASE_1_7', 'PVOPAPE_1_7')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_PVO_BASE_vs_PAPE_1_7_summary <- cbind("ICC value" = round(ICC_PVO_BASE_vs_PAPE_1_7$value, 3),
                                          "Lower 95%CI" = round(ICC_PVO_BASE_vs_PAPE_1_7$lbound, 3),
                                          "Upper 95%CI" = round(ICC_PVO_BASE_vs_PAPE_1_7$ubound, 3),
                                          "p-value" = ICC_PVO_BASE_vs_PAPE_1_7$p.value)
rownames(ICC_PVO_BASE_vs_PAPE_1_7_summary) <- "ICC_PVO_BASE_vs_PAPE_1_7_summary"

ICC_PVO_BASE_vs_PAPE_2_7 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PVOBASE_2_7', 'PVOPAPE_2_7')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_PVO_BASE_vs_PAPE_2_7_summary <- cbind("ICC value" = round(ICC_PVO_BASE_vs_PAPE_2_7$value, 3),
                                          "Lower 95%CI" = round(ICC_PVO_BASE_vs_PAPE_2_7$lbound, 3),
                                          "Upper 95%CI" = round(ICC_PVO_BASE_vs_PAPE_2_7$ubound, 3),
                                          "p-value" = ICC_PVO_BASE_vs_PAPE_2_7$p.value)
rownames(ICC_PVO_BASE_vs_PAPE_2_7_summary) <- "ICC_PVO_BASE_vs_PAPE_2_7_summary"

ICC_PVO_BASE_vs_PAPE_3_7 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PVOBASE_3_7', 'PVOPAPE_3_7')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_PVO_BASE_vs_PAPE_3_7_summary <- cbind("ICC value" = round(ICC_PVO_BASE_vs_PAPE_3_7$value, 3),
                                          "Lower 95%CI" = round(ICC_PVO_BASE_vs_PAPE_3_7$lbound, 3),
                                          "Upper 95%CI" = round(ICC_PVO_BASE_vs_PAPE_3_7$ubound, 3),
                                          "p-value" = ICC_PVO_BASE_vs_PAPE_3_7$p.value)
rownames(ICC_PVO_BASE_vs_PAPE_3_7_summary) <- "ICC_PVO_BASE_vs_PAPE_3_7_summary"

# PPO 10 BASE - PAPE
ICC_PVO_BASE_vs_PAPE_1_10 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PVOBASE_1_10', 'PVOPAPE_1_10')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_PVO_BASE_vs_PAPE_1_10_summary <- cbind("ICC value" = round(ICC_PVO_BASE_vs_PAPE_1_10$value, 3),
                                          "Lower 95%CI" = round(ICC_PVO_BASE_vs_PAPE_1_10$lbound, 3),
                                          "Upper 95%CI" = round(ICC_PVO_BASE_vs_PAPE_1_10$ubound, 3),
                                          "p-value" = ICC_PVO_BASE_vs_PAPE_1_10$p.value)
rownames(ICC_PVO_BASE_vs_PAPE_1_10_summary) <- "ICC_PVO_BASE_vs_PAPE_1_10_summary"

ICC_PVO_BASE_vs_PAPE_2_10 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PVOBASE_2_10', 'PVOPAPE_2_10')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_PVO_BASE_vs_PAPE_2_10_summary <- cbind("ICC value" = round(ICC_PVO_BASE_vs_PAPE_2_10$value, 3),
                                          "Lower 95%CI" = round(ICC_PVO_BASE_vs_PAPE_2_10$lbound, 3),
                                          "Upper 95%CI" = round(ICC_PVO_BASE_vs_PAPE_2_10$ubound, 3),
                                          "p-value" = ICC_PVO_BASE_vs_PAPE_2_10$p.value)
rownames(ICC_PVO_BASE_vs_PAPE_2_10_summary) <- "ICC_PVO_BASE_vs_PAPE_2_10_summary"

ICC_PVO_BASE_vs_PAPE_3_10 <- irr::icc(simulated_dataset_PAPE_ballistic_wide[, c('PVOBASE_3_10', 'PVOPAPE_3_10')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_PVO_BASE_vs_PAPE_3_10_summary <- cbind("ICC value" = round(ICC_PVO_BASE_vs_PAPE_3_10$value, 3),
                                          "Lower 95%CI" = round(ICC_PVO_BASE_vs_PAPE_3_10$lbound, 3),
                                          "Upper 95%CI" = round(ICC_PVO_BASE_vs_PAPE_3_10$ubound, 3),
                                          "p-value" = ICC_PVO_BASE_vs_PAPE_3_10$p.value)
rownames(ICC_PVO_BASE_vs_PAPE_3_10_summary) <- "ICC_PVO_BASE_vs_PAPE_3_10_summary"

ICC_all <- rbind(ICC_PPO_BASE_0_summary, ICC_PPO_BASE_5_summary, ICC_PPO_BASE_7_summary, ICC_PPO_BASE_10_summary,
                 ICC_PPO_PAPE_0_summary, ICC_PPO_PAPE_5_summary, ICC_PPO_PAPE_7_summary, ICC_PPO_PAPE_10_summary,
                 ICC_PVO_BASE_0_summary, ICC_PVO_BASE_5_summary, ICC_PVO_BASE_7_summary, ICC_PVO_BASE_10_summary,
                 ICC_PVO_PAPE_0_summary, ICC_PVO_PAPE_5_summary, ICC_PVO_PAPE_7_summary, ICC_PVO_PAPE_10_summary,
                 ICC_PPO_BASE_vs_PAPE_1_0_summary, ICC_PPO_BASE_vs_PAPE_2_0_summary, ICC_PPO_BASE_vs_PAPE_3_0_summary,
                 ICC_PPO_BASE_vs_PAPE_1_5_summary, ICC_PPO_BASE_vs_PAPE_2_5_summary, ICC_PPO_BASE_vs_PAPE_3_5_summary,
                 ICC_PPO_BASE_vs_PAPE_1_7_summary, ICC_PPO_BASE_vs_PAPE_2_7_summary, ICC_PPO_BASE_vs_PAPE_3_7_summary,
                 ICC_PPO_BASE_vs_PAPE_1_10_summary, ICC_PPO_BASE_vs_PAPE_2_10_summary, ICC_PPO_BASE_vs_PAPE_3_10_summary,
                 ICC_PVO_BASE_vs_PAPE_1_0_summary, ICC_PVO_BASE_vs_PAPE_2_0_summary, ICC_PVO_BASE_vs_PAPE_3_0_summary,
                 ICC_PVO_BASE_vs_PAPE_1_5_summary, ICC_PVO_BASE_vs_PAPE_2_5_summary, ICC_PVO_BASE_vs_PAPE_3_5_summary,
                 ICC_PVO_BASE_vs_PAPE_1_7_summary, ICC_PVO_BASE_vs_PAPE_2_7_summary, ICC_PVO_BASE_vs_PAPE_3_7_summary,
                 ICC_PVO_BASE_vs_PAPE_1_10_summary, ICC_PVO_BASE_vs_PAPE_2_10_summary, ICC_PVO_BASE_vs_PAPE_3_10_summary)

ICC_all <- as.data.frame(ICC_all)

ICC_all$`p-value` <- ifelse(ICC_all$`p-value` <0.001, "<0.001", "<0.05") 
print(ICC_all)
range(ICC_all$`ICC value`)

# --------------------------------------------------- CORRPLOTS ------------------------------------------------- # 

# Visualization of correlations among variables using corr plot 
# PPO BASE and PAPE among all time points
cor_PPO_all <- corrplot(corr = cor(simulated_wide_dataset_1[7:30], method = "spearman"), method="number",
                               bg = "grey10",
                               addgrid.col = "gray50", 
                               tl.cex=0.70,
                               order="hclust", 
                               number.cex=0.5,
                               addCoef.col = "black",
                               tl.col="black",
                               tl.srt=45,
                               sig.level = 0.01,
                               insig = "blank",
                               diag=TRUE,
                               col = colorRampPalette(c("yellow","green","brown1"))(100))

# PVO BASE and PAPE among all time points
cor_PVO_all <- corrplot(corr = cor(simulated_wide_dataset_1[31:54], method = "spearman"), method="number",
                               bg = "grey10",
                               addgrid.col = "gray50", 
                               tl.cex=0.70,
                               order="hclust", 
                               number.cex=0.5,
                               addCoef.col = "black",
                               tl.col="black",
                               tl.srt=45,
                               sig.level = 0.01,
                               insig = "blank",
                               diag=TRUE,
                               col = colorRampPalette(c("yellow","green","brown1"))(100))

# ----------------------------------------------------------------------------------------------------------- #
# ------------------------------------------------ mixed effect models -------------------------------------- #
# ----------------------------------------------------------------------------------------------------------- #


#########################################################################################################
############################################# PPO model #################################################
#########################################################################################################

################################### DISTRIBUTION OF RESPONSE VARIABLE ###################################

# Check the normality by the Shapiro-Wilk test
shapiro_test(data = simulated_dataset_PAPE_ballistic$PPO)
# The PPO variable seems to be not normally distributed

# Visualization of the variable - density plot
ggplot(data = simulated_dataset_PAPE_ballistic, aes(x = PPO)) + geom_density(color = "black", fill = "orange")

# Find another probability distribution to fit the PPO using fitdistrplus package

descdist(simulated_dataset_PAPE_ballistic$PPO, boot = 5000)

(norm_PPO <- fitdist(data = simulated_dataset_PAPE_ballistic$PPO, "norm"))
(lognorm_PPO <- fitdist(data = simulated_dataset_PAPE_ballistic$PPO, "lnorm"))
(gamma_PPO <- fitdist(data = simulated_dataset_PAPE_ballistic$PPO, "gamma"))

summary(norm_PPO)
gofstat(norm_PPO)
summary(lognorm_PPO)
gofstat(lognorm_PPO)
summary(gamma_PPO)
gofstat(gamma_PPO)

# plotting the density and CDF of the variable
plotdist(simulated_dataset_PAPE_ballistic$PPO, histo = T, demp = T)

cdfcomp(list(norm_PPO, lognorm_PPO, gamma_PPO),
        xlogscale = F, ylogscale = F, 
        legendtext = c("normal", "lognormal", "gamma"))

par(mfrow = c(2, 2))
plot.legend <- c('normal', 'lognormal', 'gamma')
denscomp(list(norm_PPO, lognorm_PPO, gamma_PPO), legendtext = plot.legend)
qqcomp(list(norm_PPO, lognorm_PPO, gamma_PPO), legendtext = plot.legend)
cdfcomp(list(norm_PPO, lognorm_PPO, gamma_PPO), legendtext = plot.legend)
ppcomp(list(norm_PPO, lognorm_PPO, gamma_PPO), legendtext = plot.legend)
dev.off()

# Based on the CDF - the normal distribution is considered as the better option for the response variable

simulated_dataset_PAPE_ballistic$subj_ID <- as.numeric(simulated_dataset_PAPE_ballistic$subj_ID)
simulated_dataset_PAPE_ballistic$time_point <- as.factor(simulated_dataset_PAPE_ballistic$time_point)
simulated_dataset_PAPE_ballistic$time_point <- as.character(simulated_dataset_PAPE_ballistic$time_point)
simulated_dataset_PAPE_ballistic$state <- as.character(simulated_dataset_PAPE_ballistic$state)

simulated_dataset_PAPE_ballistic$state <- factor(simulated_dataset_PAPE_ballistic$state, levels = c("BASE", "PAPE"), ordered = T)

model_PPO_1_lme4 <- lme4::lmer(PPO ~ 1 + state + time_point + state:time_point + (1|subj_ID) + (1|subj_ID:attempt_in_time_point),
                      data = simulated_dataset_PAPE_ballistic)
model_PPO_1 <- lmerTest::lmer(PPO ~ 1 + state + time_point + state:time_point + (1|subj_ID) + (1|subj_ID:attempt_in_time_point),
                          data = simulated_dataset_PAPE_ballistic)
broom.mixed::tidy(model_PPO_1)


# show equation
equatiomatic::extract_eq(model = model_PPO_1_lme4)

# library(lmerTest)
# model_1a <- lmerTest::lmer(PPO ~ state + time_point + time_point:state + (1|ID), data = data_all_1)

# summary of the model
print(model_PPO_1)
summary(model_PPO_1)

# Variance estimation
VarCorr(model_PPO_1)

# R2M and R2C
# library(MuMIn)
r.squaredGLMM(model_PPO_1)

# Parameter estimates
# library(flexplot)
estimates(model_PPO_1_lme4)

# Visualize the model (sjPlot, flexplot packages)
# library(sjPlot)
sjPlot::plot_model(model_PPO_1, show.values = T, show.p = T)


simulated_dataset_PAPE_ballistic$state <- as.character(simulated_dataset_PAPE_ballistic$state)
simulated_dataset_PAPE_ballistic$state <- factor(simulated_dataset_PAPE_ballistic$state, levels = c("BASE", "PAPE"), ordered = T)

(simulated_dataset_flexplot_PPO <- flexplot(PPO ~ state + time_point, 
                                     data = simulated_dataset_PAPE_ballistic,
                                     model = model_PPO_1, jitter = NULL, spread = "stdev") + 
    theme(axis.text.x = element_text(angle = 0, hjust = 1), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + 
    coord_cartesian(ylim = c(450, 620)) + 
    scale_colour_manual(values = c("#4682B4", "#AF46B4", "#B47846", "#4BB446")) +
    scale_linetype_manual(values = c("solid", "solid", "solid", "solid")) + 
    scale_shape_manual(values = c(19, 19, 19, 19)))

#################### and with 95%CIs based on the ggplot 2 ############################################################

simulated_dataset_PAPE_ballistic$time_point <- as.numeric(simulated_dataset_PAPE_ballistic$time_point)

summary_data_PPO <- simulated_dataset_PAPE_ballistic %>%
  group_by(state, time_point) %>%
  summarise(mean_PPO = mean(PPO, na.rm = TRUE),
            sterr = sd(PPO, na.rm = TRUE) / sqrt(n()),
            lower_ci = mean_PPO - qt(0.975, df = n() - 1) * sterr,
            upper_ci = mean_PPO + qt(0.975, df = n() - 1) * sterr)

# Create the plot
(model_PPO_plot <- ggplot() +
    # 1. Individual participant data: gray points and connecting lines
    geom_line(data = simulated_dataset_PAPE_ballistic, 
              aes(x = time_point, y = PPO, group = interaction(subj_ID, state), color = state), 
              alpha = 0.15) +  # Connect time points per subject, colored by state
    geom_point(data = simulated_dataset_PAPE_ballistic, 
               aes(x = time_point, y = PPO, color = state), 
               size = 0.5, alpha = 0.2) +  # Individual data points
    
    # 2. Mean values with 95% confidence intervals
    geom_point(data = summary_data_PPO, 
               aes(x = time_point, y = mean_PPO, color = state), 
               size = 3, shape = 16) +  # Mean points
    geom_errorbar(data = summary_data_PPO, 
                  aes(x = time_point, ymin = lower_ci, ymax = upper_ci, color = state), 
                  width = 0.2, linewidth = 1.0) +  # 95% CI error bars
    geom_line(data = summary_data_PPO, 
              aes(x = time_point, y = mean_PPO, group = state, color = state), 
              linewidth = 1.0) +  # Line connecting means per state
    
    # Labels and theme
    xlab("Time Point") +
    ylab("Seated Ballistic Military Press Throw (PPO)") +
    scale_x_continuous(breaks = c(0, 5, 7, 10)) +  # Ensure time points appear correctly
    scale_color_manual(values = c("BASE" = "blue", "PAPE" = "red")) +  # Custom colors
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 1), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank()))  # Remove legend title

# ?pretty data frame? of predicted values for the response variable and its confidence interval
# library(ggeffects)
pred_model_PPO_1 <- ggpredict(model_PPO_1, terms = c("state"))
print(pred_model_PPO_1)

################################### MULTICOLLINEARITY #################################################

# --------------------------------- car package -------------------------------- #

(VIF_model_PPO_1_car <- car::vif(model_PPO_1_lme4))

# ---------------------------- performance package ----------------------------- #

(VIF_model_PPO_1_performance <- performance::check_collinearity(model_PPO_1_lme4, component = "all"))
plot(VIF_model_PPO_1_performance)

VIF_model_PPO_1_performance <- as.data.frame(VIF_model_PPO_1_performance)
write_xlsx(VIF_model_PPO_1_performance, "E:/data/Statistics/Data/VIF_model_PPO_1_performance.xlsx")

# -------------------------------- VIF visualization --------------------------- #

par(mar = c(5, 10, 4, 8) + 0.2)

barplot(VIF_model_PPO_1_performance$VIF, main = 'VIF values for the PPO model',beside = TRUE, horiz = T,
        col = 'orange', xlim = c(0,7), border = "black", axes = T,
        cex.names = 1.0, las = 1, xlab = "Variance Inflation Factor (performance)")
abline(v = 5, col = "black", lty = 2)

barplot(VIF_model_PPO_1_car, main = 'VIF values for the PPO model',beside = TRUE, horiz = T,
        col = 'orange', xlim = c(0,7), border = "black", axes = T,
        cex.names = 1.0, las = 1, xlab = "Variance Inflation Factor (car)")
abline(v = 5, col = "black", lty = 2)

########################################## POST HOCS ##################################################


# ------------------------------ effect sizes among individual pairs ------------------ #
# ------------------------------ post hoc --------------------------------------------- #

# at first see if the interaction is significant
# (https://stats.stackexchange.com/questions/602206/posthoc-analysis-lmm-for-interaction-containing-a-factor-with-3-levels)

model_PPO_1_null <- lmerTest::lmer(PPO ~ state + time_point + (1|subj_ID) + (1 | subj_ID:attempt_in_time_point),
                                   data = simulated_dataset_PAPE_ballistic)

anova(model_PPO_1, model_PPO_1_null)

#library(emmeans)
lsmeans(model_PPO_1, pairwise ~ state:time_point)

simulated_dataset_PAPE_ballistic$time_point <- as.character(simulated_dataset_PAPE_ballistic$time_point) # for post hocks should be as character
model_PPO_1 <- lme4::lmer(PPO ~ 1 + state + time_point + state:time_point + (1|subj_ID) + (1|subj_ID:attempt_in_time_point),
                          data = simulated_dataset_PAPE_ballistic)
str(simulated_dataset_PAPE_ballistic)

(emm_PPO <- emmeans(model_PPO_1, specs = pairwise ~ state:time_point))

######################################################################################## Assumptions

# Assumptions of linear mixed effect model
# 1 Linear relationship between dependent and independent variables
# 2 Homoscedasticity assumption
# 3 Independence of errors
# 4 Normal distributions of errors
# 5 No multicollinearity among independent variables (just fixed effect)
# 6 No outliers with abnormal influence on dependent variable

# ----------------------- 2 Homoscedasticity assumption ------------------------------- #

#library(DHARMa)
res_PPO <- simulateResiduals(fittedModel = model_PPO_1, n = 10000)
residuals(res_PPO, quantileFunction = qnorm)

plot_res <- plotQQunif(res_PPO)
plot_res1 <- plotResiduals(res_PPO)

# plotting the residuals vs fitted values
plot(model_PPO_1, resid(., type="pearson") ~ fitted(.), abline=0,
     ylab = "Paerson's residuals", xlab = "Fitted values")
# there is no pattern within the residuals

# ----------------------- 4  Normal distribution of errors ---------------------------- #

# car package - qqPLot function
model_1_residuals <- resid(model_PPO_1)
qqPlot(model_1_residuals, pch = 20)

# ----------------------- 6 Outliers detection in the model --------------------------- #

DHARMa::outliers(res_PPO, lowerQuantile = 0.025, upperQuantile = 0.975)

cooks_d_model_1 <- cooks.distance(model_PPO_1)
(cutoff_model_PPO_1 <- 4/(30))
plot(cooks_d_model_1, type = "b", ylab = "Cook's Distance",
     xlab = "Observation Index", pch = 20, col = "red", ylim = c(0, 0.15))
abline(h = cutoff_model_PPO_1,lty = 2)
# "As a rule of thumb, cases are regarded as too influential if the associated value for Cook?s Distance exceeds the cut-off
# value of (Van der Meer et al., 2010):" 4/n  (or 4/ (N-k-1)) in which n refers to the number of groups in the grouping factor 
# under evaluation. 

# example with special cases for LMM
# library(LMERConvenienceFunctions)
# example with special cases for LMM/GLMM
# library(influence.ME)


estex_model_PPO_1 <- influence(model = model_PPO_1, group = "subj_ID")

dfbetas(estex_model_PPO_1)

plot(estex_model_PPO_1, which = "dfbetas", xlab="DFbetaS", ylab="Participant ID")

cooks.distance.estex(estex_model_PPO_1, sort = TRUE)

plot(estex_model_PPO_1, which = "cook",
     cutoff = cutoff_model_PPO_1, sort = TRUE,
     xlab = "Cook's Distance",
     ylab = "Participant ID", cex.axis = 0.5)

sigtest(estex = estex_model_PPO_1)

#########################################################################################################
############################################# PVO model #################################################
#########################################################################################################

################################### DISTRIBUTION OF RESPONSE VARIABLE ###################################

# Check the normality by the Shapiro-Wilk test
shapiro_test(data = simulated_dataset_PAPE_ballistic$PVO)
# The PPO variable seems to be not normally distributed

# Visualization of the variable - density plot
ggplot(data = simulated_dataset_PAPE_ballistic, aes(x = PVO)) + geom_density(color = "black", fill = "orange")

# Find another probability distribution to fit the PPO using fitdistrplus package

# plotting the density and CDF of the variable
plotdist(simulated_dataset_PAPE_ballistic$PVO, histo = T, demp = T)

descdist(simulated_dataset_PAPE_ballistic$PVO, boot = 5000)

(norm_PVO <- fitdist(data = simulated_dataset_PAPE_ballistic$PVO, "norm"))
(lognorm_PVO <- fitdist(data = simulated_dataset_PAPE_ballistic$PVO, "lnorm"))
(gamma_PVO <- fitdist(data = simulated_dataset_PAPE_ballistic$PVO, "gamma"))

summary(norm_PVO)
gofstat(norm_PVO)
summary(lognorm_PVO)
gofstat(lognorm_PVO)
summary(gamma_PVO)
gofstat(gamma_PVO)

par(mfrow = c(2, 2))
plot.legend <- c('normal', 'lognormal', 'gamma')
denscomp(list(norm_PVO, lognorm_PVO, gamma_PVO), legendtext = plot.legend)
qqcomp(list(norm_PVO, lognorm_PVO, gamma_PVO), legendtext = plot.legend)
cdfcomp(list(norm_PVO, lognorm_PVO, gamma_PVO), legendtext = plot.legend)
ppcomp(list(norm_PVO, lognorm_PVO, gamma_PVO), legendtext = plot.legend)
dev.off()

# Based on the CDF - the normal distribution is considered as the better option for the response variable

simulated_dataset_PAPE_ballistic$subj_ID <- as.numeric(simulated_dataset_PAPE_ballistic$subj_ID)
simulated_dataset_PAPE_ballistic$time_point <- as.factor(simulated_dataset_PAPE_ballistic$time_point)
simulated_dataset_PAPE_ballistic$time_point <- as.character(simulated_dataset_PAPE_ballistic$time_point)

str(simulated_dataset_PAPE_ballistic)

model_PVO_1_lme4 <- lme4::lmer(PVO ~ 1 + state + time_point + state:time_point + (1|subj_ID) + (1|subj_ID:attempt_in_time_point),
                          data = simulated_dataset_PAPE_ballistic)
model_PVO_1 <- lmerTest::lmer(PVO ~ 1 + state + time_point + state:time_point + (1|subj_ID) + (1|subj_ID:attempt_in_time_point),
                          data = simulated_dataset_PAPE_ballistic)

broom.mixed::tidy(model_PVO_1)

# show equation
equatiomatic::extract_eq(model = model_PVO_1_lme4)

# library(lmerTest)
# model_1a <- lmerTest::lmer(PPO ~ state + time_point + time_point:state + (1|ID), data = data_all_1)

# summary of the model
print(model_PVO_1)
summary(model_PVO_1)

# Variance estimation
VarCorr(model_PVO_1)

# R2M and R2C
# library(MuMIn)
r.squaredGLMM(model_PVO_1)

# Parameter estimates
# library(flexplot)
estimates(model_PVO_1_lme4)

# Visualize the model (sjPlot, flexplot packages)
# library(sjPlot)
sjPlot::plot_model(model_PVO_1, show.values = T, show.p = T)


simulated_dataset_PAPE_ballistic$state <- as.character(simulated_dataset_PAPE_ballistic$state)
simulated_dataset_PAPE_ballistic$state <- factor(simulated_dataset_PAPE_ballistic$state, levels = c("BASE", "PAPE"), ordered = T)

(simulated_dataset_flexplot_PVO <- flexplot(PVO ~ state + time_point, 
                                        data = simulated_dataset_PAPE_ballistic,
                                        model = model_PVO_1, jitter = NULL, spread = "stdev") + 
    theme(axis.text.x = element_text(angle = 0, hjust = 1), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + 
    coord_cartesian(ylim = c(1.7, 2.1)) + 
    scale_colour_manual(values = c("#4682B4", "#AF46B4", "#B47846", "#4BB446")) +
    scale_linetype_manual(values = c("solid", "solid", "solid", "solid")) + 
    scale_shape_manual(values = c(19, 19, 19, 19)))

#################### and with 95%CIs based on the ggplot 2 ############################################################

simulated_dataset_PAPE_ballistic$time_point <- as.numeric(simulated_dataset_PAPE_ballistic$time_point)

summary_data_PVO <- simulated_dataset_PAPE_ballistic %>%
  group_by(state, time_point) %>%
  summarise(mean_PVO = mean(PVO, na.rm = TRUE),
            sterr = sd(PVO, na.rm = TRUE) / sqrt(n()),
            lower_ci = mean_PVO - qt(0.975, df = n() - 1) * sterr,
            upper_ci = mean_PVO + qt(0.975, df = n() - 1) * sterr)

# Create the plot
(model_PVO_plot <- ggplot() +
    # 1. Individual participant data: gray points and connecting lines
    geom_line(data = simulated_dataset_PAPE_ballistic, 
              aes(x = time_point, y = PVO, group = interaction(subj_ID, state), color = state), 
              alpha = 0.15) +  # Connect time points per subject, colored by state
    geom_point(data = simulated_dataset_PAPE_ballistic, 
               aes(x = time_point, y = PVO, color = state), 
               size = 0.5, alpha = 0.2) +  # Individual data points
    
    # 2. Mean values with 95% confidence intervals
    geom_point(data = summary_data_PVO, 
               aes(x = time_point, y = mean_PVO, color = state), 
               size = 3, shape = 16) +  # Mean points
    geom_errorbar(data = summary_data_PVO, 
                  aes(x = time_point, ymin = lower_ci, ymax = upper_ci, color = state), 
                  width = 0.2, linewidth = 1.0) +  # 95% CI error bars
    geom_line(data = summary_data_PVO, 
              aes(x = time_point, y = mean_PVO, group = state, color = state), 
              linewidth = 1.0) +  # Line connecting means per state
    
    # Labels and theme
    xlab("Time Point") +
    ylab("Seated Ballistic Military Press Throw (PVO)") +
    scale_x_continuous(breaks = c(0, 5, 7, 10)) +  # Ensure time points appear correctly
    scale_color_manual(values = c("BASE" = "blue", "PAPE" = "red")) +  # Custom colors
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 1), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank()))  # Remove legend title

# ?pretty data frame? of predicted values for the response variable and its confidence interval
# library(ggeffects)
pred_model_PVO_1 <- ggpredict(model_PVO_1, terms = c("state"))
print(pred_model_PVO_1)

################################### MULTICOLLINEARITY #################################################

# --------------------------------- car package -------------------------------- #

(VIF_model_PVO_1_car <- car::vif(model_PVO_1_lme4))

# ---------------------------- performance package ----------------------------- #

(VIF_model_PVO_1_performance <- performance::check_collinearity(model_PVO_1_lme4, component = "all"))
plot(VIF_model_PVO_1_performance)

VIF_model_PVO_1_performance <- as.data.frame(VIF_model_PVO_1_performance)
write_xlsx(VIF_model_PVO_1_performance, "E:/data/Statistics/Data/VIF_model_PVO_1_performance.xlsx")

# -------------------------------- VIF visualization --------------------------- #



par(mar = c(5, 10, 4, 8) + 0.2)

barplot(VIF_model_PVO_1_performance$VIF, main = 'VIF values for the PVO model',beside = TRUE, horiz = T,
        col = 'orange', xlim = c(0,7), border = "black", axes = T,
        cex.names = 1.0, las = 1, xlab = "Variance Inflation Factor (performance)")
abline(v = 5, col = "black", lty = 2)

barplot(VIF_model_PVO_1_car, main = 'VIF values for the PVO model',beside = TRUE, horiz = T,
        col = 'orange', xlim = c(0,7), border = "black", axes = T,
        cex.names = 1.0, las = 1, xlab = "Variance Inflation Factor (car)")
abline(v = 5, col = "black", lty = 2)

########################################## POST HOCS ##################################################


# ------------------------------ effect sizes among individual pairs ------------------ #
# ------------------------------ post hoc --------------------------------------------- #

# at first see if the interaction is significant
# (https://stats.stackexchange.com/questions/602206/posthoc-analysis-lmm-for-interaction-containing-a-factor-with-3-levels)

model_PVO_1_null <- lmerTest::lmer(PVO ~ state + time_point + (1|subj_ID) + (1 | subj_ID:attempt_in_time_point),
                                   data = simulated_dataset_PAPE_ballistic)

anova(model_PVO_1, model_PVO_1_null)

#library(emmeans)
lsmeans(model_PVO_1, pairwise ~ state:time_point)

simulated_dataset_PAPE_ballistic$time_point <- as.character(simulated_dataset_PAPE_ballistic$time_point) # for post hocks should be as character
model_PVO_1 <- lme4::lmer(PVO ~ 1 + state + time_point + state:time_point + (1|subj_ID) + (1|subj_ID:attempt_in_time_point),
                          data = simulated_dataset_PAPE_ballistic)
str(simulated_dataset_PAPE_ballistic)

(emm_PVO <- emmeans(model_PVO_1, specs = pairwise ~ state:time_point))

######################################################################################## Assumptions

# Assumptions of linear mixed effect model
# 1 Linear relationship between dependent and independent variables
# 2 Homoscedasticity assumption
# 3 Independence of errors
# 4 Normal distributions of errors
# 5 No multicollinearity among independent variables (just fixed effect)
# 6 No outliers with abnormal influence on dependent variable

# ----------------------- 2 Homoscedasticity assumption ------------------------------- #

#library(DHARMa)
res_PVO <- simulateResiduals(fittedModel = model_PVO_1, n = 10000)
residuals(res_PVO, quantileFunction = qnorm)

plot_res <- plotQQunif(res_PVO)
plot_res1 <- plotResiduals(res_PVO)

outliers(res_PVO, lowerQuantile = 0.025, upperQuantile = 0.975)

# plotting the residuals vs fitted values
plot(model_PVO_1, resid(., type="pearson") ~ fitted(.), abline=0,
     ylab = "Paerson's residuals", xlab = "Fitted values")
# there is no pattern within the residuals


# ----------------------- 4  Normal distribution of errors ---------------------------- #

# car package - qqPLot function
model_1_residuals <- resid(model_PVO_1)
qqPlot(model_1_residuals, pch = 20)

# ----------------------- 6 Outliers detection in the model --------------------------- #

cooks_d_model_1 <- cooks.distance(model_PVO_1)
(cutoff_model_PVO_1 <- 4/(30))
plot(cooks_d_model_1, type = "b", ylab = "Cook's Distance",
     xlab = "Observation Index", pch = 20, col = "red", ylim = c(0, 0.15))
abline(h = cutoff_model_PVO_1,lty = 2)

# example with special cases for LMM
library(LMERConvenienceFunctions)
# example with special cases for LMM/GLMM
library(influence.ME)

estex_model_PVO_1 <- influence(model = model_PVO_1, "subj_ID")

dfbetas(estex_model_PVO_1)

plot(estex_model_PVO_1, which = "dfbetas", xlab="DFbetaS", ylab="Participant ID")

cooks.distance.estex(estex_model_PVO_1, sort = TRUE)

plot(estex_model_PVO_1, which = "cook",
     cutoff = cutoff_model_PVO_1, sort = TRUE,
     xlab = "Cook's Distance",
     ylab = "Participant ID", cex.axis = 0.5)

sigtest(estex = estex_model_PVO_1)

###################################################################################################
############################################# POWER ANALYSIS ######################################
###################################################################################################

###################################################################################################
###################################### data set simulation of PPO #################################
###################################################################################################

my_sim_data_PPO <- function(
    n_subj = 30,   # number of subjects
    n_BASE = 12,   # number of BASE
    n_PAPE = 12,   # number of PAPE
    beta_0 = 550,   # grand mean
    beta_1 = 15,   # effect of state
    omega_0 = 0,   # by-item random intercept sd
    tau_0 = 50,   # by-subject random intercept sd
    tau_1 = 20,   # by-subject random slope sd
    rho = 0.2,   # correlation between intercept and slope
    sigma = 20) { # residual (standard deviation)
  
  items <- data.frame(
    item_ID = seq_len(n_BASE + n_PAPE),
    state = rep(c("BASE", "PAPE"), c(n_BASE, n_PAPE)),
    X_i = rep(c(-0.5, 0.5), c(n_BASE, n_PAPE)),
    O_0i = rnorm(n = n_BASE + n_PAPE, mean = 0, sd = omega_0))
  
  # variance-covariance matrix
  cov_mx  <- matrix(
    c(tau_0^2, rho * tau_0 * tau_1, rho * tau_0 * tau_1, tau_1^2),
    nrow = 2, byrow = TRUE)
  
  subjects <- data.frame(subj_ID = seq_len(n_subj),
                         MASS::mvrnorm(n = n_subj,
                                       mu = c(T_0s = 0, T_1s = 0),
                                       Sigma = cov_mx))
  
  crossing(subjects, items) %>%
    mutate(e_si = rnorm(nrow(.), mean = 0, sd = sigma),
           PPO = beta_0 + T_0s + O_0i + (beta_1 + T_1s) * X_i + e_si,
           attempt_in_time_point = rep(1:3, length.out = nrow(.)),
           time_point = rep(c(0, 5, 7, 10), each = 3, length.out = nrow(.))) %>%
    dplyr::select(subj_ID, item_ID, state, X_i, PPO, attempt_in_time_point, time_point)
}

my_sim_data_PPO()

#################################################################################################

single_run_PPO <- function(...) {
  dat_sim <- my_sim_data_PPO(...)
  mod_sim <- lme4::lmer(PPO ~ 1 + state + time_point + state:time_point + (1|subj_ID) + (1 | subj_ID:attempt_in_time_point),
                        data = dat_sim)
  table <- broom.mixed::tidy(mod_sim, effects = "fixed", conf.int = T)
  table <- as.data.frame(table)
  table <- table %>%
    dplyr::mutate(CI_includes_zero = ifelse(conf.low <= 0 & conf.high >= 0, 0, 1))
  table <- table %>%
    dplyr::mutate(CI_includes_zero_B = ifelse(conf.low <= 0 & conf.high >= 0, "Yes", "Not"))
  print(table)
}

single_run_PPO()

###################################################################################################
###################################### data set simulation of PVO #################################
###################################################################################################

my_sim_data_PVO <- function(
    n_subj = 30,   # number of subjects
    n_BASE = 12,   # number of BASE
    n_PAPE = 12,   # number of PAPE
    beta_0 = 1.90,   # grand mean
    beta_1 = 0.03,   # effect of state
    omega_0 = 0,   # by-item random intercept sd
    tau_0 = 0.1,   # by-subject random intercept sd
    tau_1 = 0.04,   # by-subject random slope sd
    rho = 0.2,   # correlation between intercept and slope
    sigma = 0.04) { # residual (standard deviation)
  
  items <- data.frame(
    item_ID = seq_len(n_BASE + n_PAPE),
    state = rep(c("BASE", "PAPE"), c(n_BASE, n_PAPE)),
    X_i = rep(c(-0.5, 0.5), c(n_BASE, n_PAPE)),
    O_0i = rnorm(n = n_BASE + n_PAPE, mean = 0, sd = omega_0))
  
  # variance-covariance matrix
  cov_mx  <- matrix(
    c(tau_0^2, rho * tau_0 * tau_1, rho * tau_0 * tau_1, tau_1^2),
    nrow = 2, byrow = TRUE)
  
  subjects <- data.frame(subj_ID = seq_len(n_subj),
                         MASS::mvrnorm(n = n_subj,
                                       mu = c(T_0s = 0, T_1s = 0),
                                       Sigma = cov_mx))
  
  crossing(subjects, items) %>%
    mutate(e_si = rnorm(nrow(.), mean = 0, sd = sigma),
           PVO = beta_0 + T_0s + O_0i + (beta_1 + T_1s) * X_i + e_si,
           attempt_in_time_point = rep(1:3, length.out = nrow(.)),
           time_point = rep(c(0, 5, 7, 10), each = 3, length.out = nrow(.))) %>%
    dplyr::select(subj_ID, item_ID, state, X_i, PVO, attempt_in_time_point, time_point)
}

my_sim_data_PVO()

#################################################################################################

single_run_PVO <- function(...) {
  dat_sim <- my_sim_data_PVO(...)
  mod_sim <- lme4::lmer(PVO ~ 1 + state + time_point + state:time_point + (1|subj_ID) + (1 | subj_ID:attempt_in_time_point),
                        data = dat_sim)
  table <- broom.mixed::tidy(mod_sim, effects = "fixed", conf.int = T)
  table <- as.data.frame(table)
  table <- table %>%
    dplyr::mutate(CI_includes_zero = ifelse(conf.low <= 0 & conf.high >= 0, 0, 1))
  table <- table %>%
    dplyr::mutate(CI_includes_zero_B = ifelse(conf.low <= 0 & conf.high >= 0, "Yes", "Not"))
  print(table)
}

single_run_PVO()

###################################################################################################
###################################### POWER ANALYSIS SIMULATION ##################################
###################################################################################################

######################### number of simulations ##########################

n_sims <- 100

# ------------------------------------- power PPO --------------------------------------------- #

simulations_PPO <- purrr::map_df(1:n_sims, ~ single_run_PPO())

readr::write_csv(simulations_PPO, "E:/data/Statistics/Data/simulations_PPO.csv")

simulations_PPO <- readr::read_csv("E:/data/Statistics/Data/simulations_PPO.csv")

simulations_PPO %>% 
  filter(effect == "fixed") %>%
  group_by(term) %>%
  summarize(
    mean_estimate = mean(estimate),
    mean_se = mean(std.error),
    power = mean(CI_includes_zero > 0),
    power_B = mean(CI_includes_zero_B == "Not") * 100
  )

# ------------------------------------- power PVO --------------------------------------------- #

simulations_PVO <- purrr::map_df(1:n_sims, ~ single_run_PVO())

readr::write_csv(simulations_PVO, "E:/data/Statistics/Data/simulations_PVO.csv")

simulations_PVO <- readr::read_csv("E:/data/Statistics/Data/simulations_PVO.csv")

simulations_PVO %>% 
  filter(effect == "fixed") %>%
  group_by(term) %>%
  summarize(
    mean_estimate = mean(estimate),
    mean_se = mean(std.error),
    power = mean(CI_includes_zero > 0),
    power_B = mean(CI_includes_zero_B == "Not") * 100
  )

# ---------------------------------------- power simr PPO ------------------------------------- #

str(simulated_dataset_PAPE_ballistic)
simulated_dataset_PAPE_ballistic$time_point <- as.numeric(simulated_dataset_PAPE_ballistic$time_point)
model_PPO_1_simr <- lme4::lmer(PPO ~ 1 + state + time_point + state:time_point + (1|subj_ID) + (1|subj_ID:attempt_in_time_point),
                               data = simulated_dataset_PAPE_ballistic)

power_model_1_1 <- powerSim(fit = model_PPO_1_simr, fixed('state'), nsim = 10)


print(power_model_1_1)
summary(power_model_1_1)

extended_model_1_1 <- simr::extend(model_PPO_1_simr, along = "subj_ID", n = 100)
pc_extended_model_1_1 <- powerCurve(extended_model_1_1,
                                    fixed('state'),
                                    along = 'subj_ID',
                                    nsim = 10, 
                                    alpha = 0.05, 
                                    progress = T)

print(summary(pc_extended_model_1_1))
plot(pc_extended_model_1_1)


# --------------------------------------------------------------------------------------------------------- #
# ---------------------------------------- Randomization of testing sessions ------------------------------ #
# --------------------------------------------------------------------------------------------------------- #

participants_1 <- data.frame(ID = 1:30)

# Define a function to generate a random permutation of integers from 1 to 2
generate_permutation <- function() {
  sample(1:2)
}

# Define a vector of strings corresponding to the integers 1 to 2
replacement_strings_control_PAPE <- c("control", "PAPE")

# Apply the function to generate two new columns
participants_1[, 2:3] <- t(replicate(nrow(participants_1), generate_permutation()))

participants_1 <- participants_1 %>%
  rename(testing_session_2 = V2, testing_session_3 = V3)

# print(participants_1)

# Replace integers with strings in the new columns
for (i in 2:3) {
  participants_1[[i]] <- replacement_strings_control_PAPE[participants_1[[i]]]
}

print(participants_1)
