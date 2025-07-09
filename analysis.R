###############################################

# R script for Statistical analyses
# Associated to "Trait-specific behavioural plasticity in response to repeated stressor exposure in wild crickets"
# Gilford, Emily R; Li, Ruonan; Rodríguez-Muñoz, Rolando; Tregenza, Tom; Kuijper, Bram

###############################################

# This R script is organised as follows

# 1. Dependencies and setup

# 2. Data Handling
### a) Sign reversal of post-emergence distance
### b) Within subject centering by ball type

# 3. Statistical Analysis 
### a) Model creation: rsmodel
### b) Posterior prediction plots: rsmodel1 (found in supplementary materials)
#### b_i) Emergence
#### b_ii) Escape speed
#### b_iii) Post-emergence distance
#### b_iii) Combined posterior predictions graphs
### c) Individual plasticity 

# 4. Result visualisation and plotting
### a) Figure 1 - Behavioural responses across repeated trials
### b) Figure 2 - Forest plot
### c) Figure 3 - Random slope estimates across stimuli

###############################################

#### 1. Dependencies and setup ####

## Set working directory
setwd("~/01 Science/3. Field/2024/Acute/data_code")

## Install and load packages
packages_list <- c("broom", "brms", "cowplot", "dplyr", "ggplot2",
                   "MASS","patchwork", "posterior", "purrr", "rstan", "tidyr")

install.packages(packages_list)

lapply(packages_list, function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
    library(package, character.only = TRUE)
  }
})


#### 2. Data handling ####

data_all <- read.csv("data.csv", header=TRUE)

## Column selection
data_all<- data_all %>%
  dplyr::select(idx, date, burrow, sex, ball, trial, temp, flee_dis, fs1, emergence, path_dist, n_trials_treatment)

## Converting ball type to stimulus intensity
data_all <- data_all %>%
  rename(stimulus = ball) %>%
  mutate(stimulus = factor(stimulus,
                           levels = c("cork", "steel"),
                           labels = c("weak", "strong")))

#### a) Sign reversal of post-emergence distance ####

## Smaller distances indicate habituation across all response variables
data_all$flee_disflip <- -data_all$flee_dis

### Variable naming notes:
## fs1 = Escape speed (m/s), calculated over first 1.5 cm of movement
## flee_dis = Post-emergence distance (cm), renamed as flee_disflip after sign reversal
## flee_disflip = Negative post-emergence distance, so that lower values reflect habituation
## emergence = Emergence time (s), latency to emerge after disturbance
## stimulus = Stimulus intensity (factor: "weak" = cork ball, "strong" = steel ball)

#### b) Within subject centering by stimulus type ####
data_all_centered <- data_all %>%
  group_by(burrow, stimulus) %>%
  mutate(
    fs1_centered = fs1 - mean(fs1, na.rm = TRUE),
    flee_disflip_centered = flee_disflip - mean(flee_disflip, na.rm = TRUE),
    emergence_centered = emergence - mean(emergence, na.rm = TRUE)
  ) %>%
  ungroup()

#### 3. Statistical Analysis ####
#### a) Model creation: rsmodel ####

### Model rationale:
## Although response variables were right-skewed, Box-Cox analyses suggested
## λ ≈ 1, indicating that transformations were unnecessary.
## Posterior predictive checks showed that Gaussian models provided a better fit
## than Gamma or log-normal alternatives, so Gaussian was used for all three 
## responses in the multivariate model.

## Default priors from the `brms` package were used, as the data were 
## sufficiently informative and posterior predictive checks showed good fit.

rsmodel <- brm(
  mvbind(fs1_centered, flee_disflip_centered, emergence_centered) ~ 
    trial * stimulus + temp + sex + path_dist + (trial * stimulus | p | burrow), 
  data = data_all_centered,
  family = gaussian(),
  chains = 4, cores = 4,
  control = list(adapt_delta = 0.99))

summary(rsmodel)

#### b) Posterior prediction plots: rsmodel ####
custom_colours <- c("weak" = "#D55E00", "strong" = "#0072B2")
theme_set(theme_classic(base_family = "serif"))

#### b_i) Emergence ####
newdata_emerg <- expand.grid(
  trial = 1:10,
  stimulus = c("weak", "strong"),
  temp = 0,
  sex = "female",
  path_dist = 0,
  burrow = NA
)

## Generate posterior predictions
emerg_pred <- fitted(
  rsmodel,
  newdata = newdata_emerg,
  resp = "emergencecentered",
  re_formula = NA,
  summary = TRUE
)

## Combine predictions with data
emerg_pred_df <- cbind(newdata_emerg, emerg_pred)

## Prep raw data
emergence_raw_df <- data_all_centered %>%
  dplyr::select(trial, stimulus, emergence_centered)

## Overlay plot
em_overlay_pp <- ggplot(emerg_pred_df, aes(x = trial, y = Estimate, colour = stimulus, fill = stimulus)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = 0.2, colour = NA) +
  geom_jitter(data = emergence_raw_df, aes(x = trial, y = emergence_centered, colour = stimulus),
              width = 0.2, alpha = 0.4, size = 1.5, inherit.aes = FALSE) +
  scale_color_manual(values = custom_colours) +
  scale_fill_manual(values = custom_colours) +
  scale_x_continuous(breaks = 1:10) +
  labs(
    x = "Trial",
    y = "Emergence time relative to mean",
    colour = "Stimulus intensity",
    fill = "Stimulus intensity"
  ) +
  theme(
    panel.grid.major = element_blank(),  # remove major gridlines
    panel.grid.minor = element_blank(),  # remove minor gridlines
    axis.line = element_line(color = "black"),  # restore axis lines
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

em_overlay_pp

ggsave("em_overlay_pp.png", em_overlay_pp, width = 10.813, height = 4.698, dpi = 300)


#### b_ii) Escape speed ####

newdata_fs1 <- expand.grid(
  trial = 1:10,
  stimulus = c("weak", "strong"),
  temp = 0,
  sex = "female",
  path_dist = 0,
  burrow = NA
)

## Generate posterior predictions
fs1_pred <- fitted(
  rsmodel,
  newdata = newdata_fs1,
  resp = "fs1centered",
  re_formula = NA,
  summary = TRUE
)

## Combine predictions with data
fs1_pred_df <- cbind(newdata_fs1, fs1_pred)

## Prep raw data
fs1_raw_df <- data_all_centered %>%
  dplyr::select(trial, stimulus, fs1_centered)

## Overlay plot
fs1_overlay_pp <- ggplot(fs1_pred_df, aes(x = trial, y = Estimate, colour = stimulus, fill = stimulus)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = 0.2, colour = NA) +
  geom_jitter(data = fs1_raw_df, aes(x = trial, y = fs1_centered, colour = stimulus),
              width = 0.2, alpha = 0.4, size = 1.5, inherit.aes = FALSE) +
  scale_color_manual(values = custom_colours) +
  scale_fill_manual(values = custom_colours) +
  scale_x_continuous(breaks = 1:10) +
  labs(
    x = "Trial",
    y = "Escape speed relative to mean",
    colour = "Stimulus intensity",
    fill = "Stimulus intensity"
  ) +
  theme(
    panel.grid.major = element_blank(),  # remove major gridlines
    panel.grid.minor = element_blank(),  # remove minor gridlines
    axis.line = element_line(color = "black"),  # restore axis lines
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

fs1_overlay_pp

ggsave("fs1_overlay_pp.png", fs1_overlay_pp, width = 10.813, height = 4.698, dpi = 300)

##### b_iii) Combined posterior predictions graphs ####
# Update legends and update plot themes
fs1_overlay_pp <- fs1_overlay_pp + theme(legend.position = "none")
em_overlay_pp <- em_overlay_pp + theme(legend.position = "none")
dfb_overlay_pp <- dfb_overlay_pp + theme(legend.position = "none")

em_overlay_pp <- em_overlay_pp +
  labs(y = "Emergence time\nrelative to mean") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(),
    axis.ticks.x = element_line()
  )

fs1_overlay_pp <- fs1_overlay_pp +
  labs(y = "Escape speed\nrelative to mean") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(),
    axis.ticks.x = element_line()
  )

dfb_overlay_pp <- dfb_overlay_pp +
  labs(y = "Post-emergence distance\nrelative to mean", x = "Trial")

legend_plot <- get_legend(
  em_overlay_pp +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 12)
    ) +
    guides(fill = guide_legend(nrow = 1), colour = guide_legend(nrow = 1))
)

# Combine plots, add legend separately
combined_plot <- (em_overlay_pp / fs1_overlay_pp / dfb_overlay_pp) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 14, face = "bold")
  )

final_plot <- plot_grid(
  combined_plot,
  legend_plot,
  ncol = 1,
  rel_heights = c(1, 0.08)
)

# View and save 
final_plot
ggsave("combined_figure.png", final_plot, width = 8, height = 13, dpi = 300)



#### b_iii) Post-emergence distance ####
newdata_dfb <- expand.grid(
  trial = 1:10,
  stimulus = c("weak", "strong"),
  temp = 0,
  sex = "female",
  path_dist = 0,
  burrow = NA
)

## Generate posterior predictions
dfb_pred <- fitted(
  rsmodel,
  newdata = newdata_dfb,
  resp = "fleedisflipcentered",
  re_formula = NA,
  summary = TRUE
)

## Combine predictions with data
dfb_pred_df <- cbind(newdata_dfb, dfb_pred)

## Prep raw data
dfb_raw_df <- data_all_centered %>%
  dplyr::select(trial, stimulus, flee_disflip_centered)

## Overlay plot
dfb_overlay_pp <- ggplot(dfb_pred_df, aes(x = trial, y = Estimate, colour = stimulus, fill = stimulus)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = 0.2, colour = NA) +
  geom_jitter(data = dfb_raw_df, aes(x = trial, y = flee_disflip_centered, colour = stimulus),
              width = 0.2, alpha = 0.4, size = 1.5, inherit.aes = FALSE) +
  scale_color_manual(values = custom_colours) +
  scale_fill_manual(values = custom_colours) +
  scale_x_continuous(breaks = 1:10) +
  labs(
    x = "Trial",
    y = "Post-emergence distance relative to mean",
    colour = "Stimulus intensity",
    fill = "Stimulus intensity"
  ) +
  theme(
    panel.grid.major = element_blank(),  # remove major gridlines
    panel.grid.minor = element_blank(),  # remove minor gridlines
    axis.line = element_line(color = "black"),  # restore axis lines
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

dfb_overlay_pp

ggsave("dfb_overlay_pp.png", dfb_overlay_pp, width = 10.813, height = 4.698, dpi = 300)


#### c) Individual plasticity ####

## Calculate trial slopes per burrow and stimulus type
get_trial_slope <- function(df, response) {
  df %>%
    group_by(burrow, stimulus) %>%
    nest() %>%
    mutate(model = purrr::map(data, ~ lm(reformulate("trial", response), data = .x)),
           tidied = purrr::map(model, broom::tidy)) %>%
    unnest(tidied) %>%
    filter(term == "trial") %>%
    dplyr::select(burrow, stimulus, estimate) %>%
    tidyr::pivot_wider(names_from = stimulus, values_from = estimate, names_prefix = paste0(response, "_"))
}

## Calculate slopes 
emergence_slopes <- get_trial_slope(data_all_centered, "emergence_centered")
fs1_slopes       <- get_trial_slope(data_all_centered, "fs1_centered")
dfb_slopes       <- get_trial_slope(data_all_centered, "flee_disflip_centered")

## Combine to one data frame
## filtered to only include individuals (n=17) which experienced both stimuli

paired_slopes <- emergence_slopes %>%
  full_join(fs1_slopes, by = "burrow") %>%
  full_join(dfb_slopes, by = "burrow") %>%
    filter(
    !is.na(emergence_centered_weak) & !is.na(emergence_centered_strong),
    !is.na(fs1_centered_weak) & !is.na(fs1_centered_strong),
    !is.na(flee_disflip_centered_weak) & !is.na(flee_disflip_centered_strong)
  )


## Correlations
cor_emergence <- cor(paired_slopes$emergence_centered_weak, paired_slopes$emergence_centered_strong, use = "complete.obs")
cor_fs1       <- cor(paired_slopes$fs1_centered_weak, paired_slopes$fs1_centered_strong, use = "complete.obs")
cor_dfb       <- cor(paired_slopes$flee_disflip_centered_weak, paired_slopes$flee_disflip_centered_strong, use = "complete.obs")

cor_emergence
cor_fs1
cor_dfb

## Paired-t tests
t_emergence <- t.test(paired_slopes$emergence_centered_weak, paired_slopes$emergence_centered_strong, paired = TRUE)
t_fs1       <- t.test(paired_slopes$fs1_centered_weak, paired_slopes$fs1_centered_strong, paired = TRUE)
t_dfb       <- t.test(paired_slopes$flee_disflip_centered_weak, paired_slopes$flee_disflip_centered_strong, paired = TRUE)

t_emergence
t_fs1
t_dfb

#### 4. Result visualisation and plotting ####

custom_colours <- c("weak" = "#D55E00", "strong" = "#0072B2")

#### a) Figure 1 - Behavioural responses across repeated trials ####

## Escape speed plot
## Calc means
fs_summary <- data_all %>%
  group_by(trial, stimulus) %>%
  summarize(mean_flee_speed = mean(fs1, na.rm = TRUE),
            se = sd(fs1, na.rm = TRUE)/sqrt(n()),
            ci_lower = mean_flee_speed - 1.96 * se,
            ci_upper = mean_flee_speed + 1.96 * se,
            .groups = 'drop')

fsmeanplot <- ggplot() +
  geom_point(data = data_all, aes(x = as.integer(trial), y = fs1, color = stimulus, group = interaction(burrow, stimulus)), 
             alpha = 0.2, size = 1) +
  geom_line(data = data_all, aes(x = as.integer(trial), y = fs1, color = stimulus, group = interaction(burrow, stimulus)), 
            alpha = 0.2, size = 0.7) +
  geom_errorbar(data = fs_summary, aes(x = as.integer(trial), ymin = ci_lower, ymax = ci_upper, color = stimulus), 
                width = 0.2, size = 0.7) +
  geom_line(data = fs_summary, aes(x = as.integer(trial), y = mean_flee_speed, color = stimulus, group = stimulus), 
            size = 1.2) +
  geom_point(data = fs_summary, aes(x = as.integer(trial), y = mean_flee_speed, color = stimulus), 
             size = 3) +
  labs(
    x = "Trial number",
    y = "Escape speed\n(m/s)"
  ) +
  scale_x_continuous(breaks = 1:10) + 
  scale_y_continuous(limits = c(0.0, 0.7), breaks = seq(0.0, 0.7, by = 0.1)) +
  scale_color_manual(values = custom_colours) +
  theme(
    axis.ticks = element_line(color = "black"),
    legend.position = "none",
    legend.title = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12))

print(fsmeanplot)

## DFB plots

## Calculating means
dfb_summary <- data_all %>%
  group_by(trial, stimulus) %>%
  summarize(mean_flee_dist = mean(flee_dis, na.rm = TRUE),
            se = sd(flee_dis, na.rm = TRUE)/sqrt(n()),
            ci_lower = mean_flee_dist - 1.96 * se,
            ci_upper = mean_flee_dist + 1.96 * se,
            .groups = 'drop')

## Plot of means and individual's DFB
dfbmeanplot <- ggplot() +
  geom_point(data = data_all, aes(x = as.integer(trial), y = flee_dis, color = stimulus, group = interaction(burrow, stimulus)), 
             alpha = 0.2, size = 1) +
  geom_line(data = data_all, aes(x = as.integer(trial), y = flee_dis, color = stimulus, group = interaction(burrow, stimulus)), 
            alpha = 0.2, size = 0.7) +
  geom_errorbar(data = dfb_summary, aes(x = as.integer(trial), ymin = ci_lower, ymax = ci_upper, color = stimulus), 
                width = 0.2, size = 0.7) +
  geom_line(data = dfb_summary, aes(x = as.integer(trial), y = mean_flee_dist, color = stimulus, group = stimulus), 
            size = 1.2) +
  geom_point(data = dfb_summary, aes(x = as.integer(trial), y = mean_flee_dist, color = stimulus), 
             size = 3) +
  labs(
    x = "Trial number",
    y = "Post-emergence\ndistance (cm)"
  ) +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(limits = c(1, 7), breaks = seq(1, 7, by = 1)) +
  scale_color_manual(values = custom_colours) +
  theme(
    axis.ticks = element_line(color = "black"),
    legend.position = "none",
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12))

print(dfbmeanplot)

## Emergence

## Calculating means
em_summary <- data_all %>%
  group_by(trial, stimulus) %>%
  summarize(
    mean_emergence = mean(emergence, na.rm = TRUE),
    sd = sd(emergence, na.rm = TRUE),
    n = sum(!is.na(emergence)),
    se = sd / sqrt(n),
    ci_lower = mean_emergence - qt(0.975, df = n - 1) * se,
    ci_upper = mean_emergence + qt(0.975, df = n - 1) * se,
    .groups = 'drop'
  )
## Plot of means and individual's emergence time
emmeanplot <- ggplot() +
  geom_point(data = data_all, aes(x = as.integer(trial), y = emergence, color = stimulus, group = interaction(burrow, stimulus)), 
             alpha = 0.2, size = 1) +
  geom_line(data = data_all, aes(x = as.integer(trial), y = emergence, color = stimulus, group = interaction(burrow, stimulus)), 
            alpha = 0.2, size = 0.7) +
  geom_errorbar(data = em_summary, 
                aes(x = as.integer(trial), ymin = ci_lower, ymax = ci_upper, color = stimulus), 
                width = 0.2, size = 0.8) +
  geom_line(data = em_summary, aes(x = as.integer(trial), y = mean_emergence, color = stimulus, group = stimulus), 
            size = 1.2) +
  geom_point(data = em_summary, aes(x = as.integer(trial), y = mean_emergence, color = stimulus), 
             size = 3) +
  labs(
    x = "Trial number",
    y = "Emergence time (s)"
  ) +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(limits = c(0, 550), breaks = seq(0, 550, by = 100)) +
  scale_color_manual(values = custom_colours, name = "Stimulus\nintensity") +
  theme(
    axis.ticks = element_line(color = "black"),
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

## Update other plots to match
fsmeanplot <- fsmeanplot + theme(
  axis.title = element_text(size = 14),
  axis.text = element_text(size = 12),
  legend.position = "none"
)

dfbmeanplot <- dfbmeanplot + theme(
  axis.title = element_text(size = 14),
  axis.text = element_text(size = 12),
  legend.position = "none"
)

## Combine plots
top <- emmeanplot
bottom <- plot_grid(fsmeanplot, dfbmeanplot, ncol = 2, labels = c("B", "C"), label_y = 1.05)

mean_indiv_resp_stimulustype_plot <- plot_grid(
  plot_grid(top, bottom, ncol = 1, rel_heights = c(2.05, 0.95), labels = c("A", ""), label_y = 1.02),
  ncol = 1
)

print(mean_indiv_resp_stimulustype_plot)

ggsave("emfs1dfb_plot_stimintens.png", mean_indiv_resp_stimulustype_plot, width = 12, height = 7, dpi = 300)


#### b) Figure 2 - Forest plot ####

## Axes switched, where y is x, and x is y

## Extract posterior draws
posterior <- as_draws_df(rsmodel)

## Summarise fixed effects
summary_df2 <- posterior %>%
  dplyr::select(starts_with("b_")) %>%
  pivot_longer(everything(), names_to = "Parameter") %>%
  group_by(Parameter) %>%
  summarise(
    mean = mean(value),
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975),
    .groups = "drop"
  ) %>%
  mutate(
    Response = case_when(
      grepl("fs1centered", Parameter) ~ "Flee speed",
      grepl("fleedisflipcentered", Parameter) ~ "Distance from burrow",
      grepl("emergencecentered", Parameter) ~ "Emergence time",
      TRUE ~ "Other"
    ),
    Parameter_clean = case_when(
      grepl("Intercept", Parameter) ~ "Intercept",
      grepl("trial:stimulusstrong", Parameter) ~ "Trial × Stimulus",
      grepl("trial", Parameter) ~ "Trial",
      grepl("stimulusstrong", Parameter) ~ "Stimulus: strong",
      grepl("temp", Parameter) ~ "Temperature",
      grepl("sexmale", Parameter) ~ "Sex: Male",
      grepl("sexunknown", Parameter) ~ "Sex: Unknown",
      grepl("path_dist", Parameter) ~ "Path Distance",
      TRUE ~ NA_character_
    )
  )

## Set factor level order
summary_df2$Parameter_clean <- factor(
  summary_df2$Parameter_clean,
  levels = rev(c("Intercept", "Trial", "Stimulus: strong",
                 "Temperature", "Sex: Male", "Sex: Unknown",
                 "Path Distance", "Trial × Stimulus"))
)

## Filter by response
summary_df2_flee_speed <- summary_df2 %>% filter(Response == "Flee speed")
summary_df2_distance <- summary_df2 %>% filter(Response == "Distance from burrow")
summary_df2_emergence <- summary_df2 %>% filter(Response == "Emergence time")

## Plot functions
make_forestplot <- function(df, title, y_label = NULL, show_y_labels = TRUE, y_limits = c(-1, NA)) {
  ggplot(df, aes(x = Parameter_clean, y = mean, ymin = lower, ymax = upper)) +
    geom_pointrange() +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    theme(
      axis.ticks = element_line(color = "black"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title.position = "plot",
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 12)
    ) +
    labs(
      title = title,
      x = if (show_y_labels) NULL else NULL,
      y = y_label
    ) +
    scale_y_continuous(limits = y_limits)
}

## Plot creation
forestplot_emergence2 <- make_forestplot(
  summary_df2_emergence,
  "Emergence time",
  y_label = NULL,
  show_y_labels = TRUE,
  y_limits = c(-50, NA)
)

forestplot_escape_speed2 <- make_forestplot(
  summary_df2_flee_speed,
  "Escape speed",
  y_label = "Estimate",
  show_y_labels = FALSE,
  y_limits = c(-0.1, NA)
) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())

forestplot_distance2 <- make_forestplot(
  summary_df2_distance,
  "Post-emergence distance",
  show_y_labels = FALSE,
  y_limits = c(-1, NA)
) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())

## Combine plots
forestplot_centered <- plot_grid(
  forestplot_emergence2,
  forestplot_escape_speed2,
  forestplot_distance2,
  ncol = 3,
  labels = c("A", "B", "C"),
  label_size = 14,
  align = "h"
)

print(forestplot_centered)

ggsave("forestplot_centered_stimulus.png", forestplot_centered, width = 10.813, height = 4.698, dpi = 300)


### c) Figure 3 - Random slope estimates across stimuli ####

## Filter for burrows that experienced both stimuli (n=17)
paired_slopes_filtered <- paired_slopes %>%
  filter(!is.na(emergence_centered_weak) & !is.na(emergence_centered_strong),
         !is.na(fs1_centered_weak) & !is.na(fs1_centered_strong),
         !is.na(flee_disflip_centered_weak) & !is.na(flee_disflip_centered_strong))

## Order x-axis
emergence_long <- paired_slopes_filtered %>%
  dplyr::select(burrow, emergence_centered_weak, emergence_centered_strong) %>%
  pivot_longer(cols = starts_with("emergence"),
               names_to = "stimulus",
               values_to = "slope") %>%
  mutate(stimulus = recode(stimulus,
                           emergence_centered_weak = "weak",
                           emergence_centered_strong = "strong"),
         stimulus = factor(stimulus, levels = c("weak", "strong")))

fs1_long <- paired_slopes_filtered %>%
  dplyr::select(burrow, fs1_centered_weak, fs1_centered_strong) %>%
  pivot_longer(cols = starts_with("fs1"),
               names_to = "stimulus",
               values_to = "slope") %>%
  mutate(stimulus = recode(stimulus,
                           fs1_centered_weak = "weak",
                           fs1_centered_strong = "strong"),
         stimulus = factor(stimulus, levels = c("weak", "strong")))

dfb_long <- paired_slopes_filtered %>%
  dplyr::select(burrow, flee_disflip_centered_weak, flee_disflip_centered_strong) %>%
  pivot_longer(cols = starts_with("flee"),
               names_to = "stimulus",
               values_to = "slope") %>%
  mutate(stimulus = recode(stimulus,
                           flee_disflip_centered_weak = "weak",
                           flee_disflip_centered_strong = "strong"),
         stimulus = factor(stimulus, levels = c("weak", "strong")))

## Plot creation per response variable

## Emergence
emps <- ggplot(emergence_long, aes(x = stimulus, y = slope, group = burrow)) +
  geom_point(aes(color = stimulus), size = 3) +
  geom_line(alpha = 0.4) +
  scale_color_manual(values = custom_colours) +
  labs(title = "Emergence time", x = NULL, y = "Slope of Trial Effect") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")

## Escape speed
fs1ps <- ggplot(fs1_long, aes(x = stimulus, y = slope, group = burrow)) +
  geom_point(aes(color = stimulus), size = 3) +
  geom_line(alpha = 0.4) +
  labs(title = "Escape speed", x = "Stimulus intensity", y = NULL) +
  scale_color_manual(values = custom_colours) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")

## Post-emergence distance
dfbps <- ggplot(dfb_long, aes(x = stimulus, y = slope, group = burrow)) +
  geom_point(aes(color = stimulus), size = 3) +
  geom_line(alpha = 0.4) +
  labs(title = "Post-emergence distance", x = NULL, y = NULL, color = "Stimulus") +
  scale_color_manual(values = custom_colours) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right")

## Combine plots
indiv_slopes_balltype_filtered <- plot_grid(
  emps, fs1ps, dfbps,
  ncol = 3,
  labels = c("A", "B", "C"),
  label_size = 14,
  align = "h"
)

indiv_slopes_balltype_filtered
ggsave("indiv_slopes_balltype_filtered.png", indiv_slopes_balltype_filtered, width = 12.188, height = 3.25, dpi = 300)


