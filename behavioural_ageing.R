
################################################################################
## Evaluate the effect of age on different behaviours: 
# a) Roost fidelity
# b) Probability of having a routine
# c) Strength of the movement routine
# d) Probability of occupying a popular roost
# e) Normalized degree
# f) Average strength
################################################################################


# ---- Set up ----
# Load the libraries
library(glmmTMB)
library(lme4)
library(lmerTest)
library(dplyr)
library(performance)
library(ggeffects)
library(ggplot2)


# For the models
functional_relationships <- list("age_scale", "poly(age_scale, 2, raw = TRUE)", 
                                 "expm1(age_scale)", "poly(age_scale, 3, raw = TRUE)")
func_rel_names <- list("Linear", "Quadratic", "Exponential", "Third-degree polynomial")

# For the figures
season_colours <- c("Breeding" = "#24586e", "Transient" = "#CE8B00", "Post-breeding" = "#A3BB98")


# ---- Upload datasets ----
# If necessary, download the datasets from Zenodo
roost_df <- read.csv("df_roost_movement.csv")
roost_df$age_scale <- scale(roost_df$age, center = TRUE, scale = TRUE)
roost.age.scale <- attr(roost_df$age_scale,"scaled:scale")
roost.age.center <- attr(roost_df$age_scale,"scaled:center")

routine_df <- read.csv("df_routine.csv")
routine_df$age_scale <- scale(routine_df$age, center = TRUE, scale = TRUE)
routine_df$seq_length_scale <- scale(routine_df$seq_length, center = FALSE, scale = TRUE)

routine.age.scale <- attr(routine_df$age_scale,"scaled:scale")
routine.age.center <- attr(routine_df$age_scale,"scaled:center")
routine.seq.scale <- attr(routine_df$seq_length_scale,"scaled:scale")

routine_index_df <- subset(routine_df, with_routine == 1)
routine_index_df$age_scale <- scale(routine_index_df$age, center = TRUE, scale = TRUE)
index.age.scale <- attr(routine_index_df$age_scale,"scaled:scale")
index.age.center <- attr(routine_index_df$age_scale,"scaled:center")

social_df <- read.csv("df_social.csv")
social_df$age_scale <- scale(social_df$age, center = TRUE, scale = TRUE)
social.age.scale <- attr(social_df$age_scale,"scaled:scale")
social.age.center <- attr(social_df$age_scale,"scaled:center")


# ---- a) Roost fidelity ----
roost_fidelity_models <- list()
roost_fidelity_mod_summary <- list()

for(i in 1:length(functional_relationships)){
  
  model_formula <- paste("roost_fid ~ season * ", functional_relationships[[i]], 
                         " + (1|ID) + (1|year)", sep ="")
  
  mod <- glmmTMB(formula = as.formula(model_formula), 
                 REML = F,                              # for model comparison
                 family = binomial(link="logit"),
                 data = roost_df)
  
  roost_fidelity_models[[i]] <- mod
  roost_fidelity_mod_summary[[i]] <- summary(mod)
  
  rm(mod)
  rm(model_formula)
  
}

roost_fidelity_comparison <- compare_performance(roost_fidelity_models[[1]], roost_fidelity_models[[2]], 
                                                 roost_fidelity_models[[3]], roost_fidelity_models[[4]])
roost_fidelity_comparison$Model <- func_rel_names

## Figure 3A ##
roost_fidel_mod_df <- ggpredict(roost_fidelity_models[[4]], terms = c("age_scale[all]", "season"))
roost_fidel_mod_df$age <- round(roost_fidel_mod_df$x * roost.age.scale + roost.age.center)
roost_fidel_mod_df <- rename(roost_fidel_mod_df, season = group)

roost_fidelity_figure <- 
  ggplot() +
  geom_ribbon(data = roost_fidel_mod_df,
              aes(x = age, y = predicted, 
                  ymin = conf.low, ymax = conf.high, 
                  fill = season), 
              alpha = 0.2, linewidth = 0.6) +
  geom_line(data = roost_fidel_mod_df, 
            aes(x = age, y = predicted, colour = season), 
            linewidth = 0.7) +
  scale_y_continuous(name = "Probability of roost fidelity",
                     limits = c(0,1)) +
  xlab("Age (years)") +
  scale_fill_manual(name = "Season", 
                    values =  season_colours) +
  scale_colour_manual(name = "Season", 
                      values = season_colours) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


# ---- b) Probability of having a routine ----
routine_models <- list()
routine_mod_summary <- list()

for(i in 1:length(functional_relationships)){
  
  model_formula <- paste("with_routine ~ log(seq_length_scale) + ", functional_relationships[[i]],
                         " + (1|ID) + (1|year)", sep ="")
  
  mod <- glmmTMB(formula = as.formula(model_formula),
                 REML = F,                              # for model comparison
                 family = binomial(link="logit"),
                 data = routine_df)
  
  routine_models[[i]] <- mod
  routine_mod_summary[[i]] <- summary(mod)
  
  rm(mod)
  rm(model_formula)
  
}

routine_comparison <- compare_performance(routine_models[[1]], routine_models[[2]], 
                                          routine_models[[3]], routine_models[[4]])
routine_comparison$Model <- func_rel_names


## Figure 3B ##
routine_mod_df <- ggpredict(routine_models[[1]], terms = c("age_scale[all]"))
routine_mod_df$age <- round(routine_mod_df$x * routine.age.scale + routine.age.center)

routine_figure <- 
  ggplot() +
  geom_ribbon(data = routine_mod_df,
              aes(x = age, y = predicted, 
                  ymin = conf.low, ymax = conf.high), 
              alpha = 0.2, linewidth = 0.6, fill = "#B46060") +
  geom_line(data = routine_mod_df, 
            aes(x = age, y = predicted), 
            linewidth = 0.7, colour = "#B46060") +
  scale_y_continuous(name = "Probability of having a routine",
                     limits = c(0,1)) +
  xlab("Age (years)") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


# ---- c) Routine index ----
index_models <- list()
index_mod_summary <- list()

for(i in 1:length(functional_relationships)){
  
  model_formula <- paste("R ~ ", functional_relationships[[i]],
                         " + (1|ID) + (1|year)", sep ="")
  
  mod <- glmmTMB(formula = as.formula(model_formula),
                 REML = F,                              # for model comparison
                 family = beta_family(link="logit"),
                 data = routine_index_df)
  
  index_models[[i]] <- mod
  index_mod_summary[[i]] <- summary(mod)
  
  rm(mod)
  rm(model_formula)
  
}

index_comparison <- compare_performance(index_models[[1]], index_models[[2]], 
                                        index_models[[3]], index_models[[4]])
index_comparison$Model <- func_rel_names


## Figure 3C ##
index_mod_df <- ggpredict(index_models[[3]], terms = c("age_scale[all]"))
index_mod_df$age <- round(index_mod_df$x * index.age.scale + index.age.center)

index_figure <- 
  ggplot() +
  geom_ribbon(data = index_mod_df,
              aes(x = age, y = predicted, 
                  ymin = conf.low, ymax = conf.high), 
              alpha = 0.2, linewidth = 0.6, fill = "#B46060") +
  geom_line(data = index_mod_df, 
            aes(x = age, y = predicted), 
            linewidth = 0.7, colour = "#B46060") +
  scale_y_continuous(name = "Index of routine (0-1)",
                     limits = c(0.3,0.7)) +
  xlab("Age (years)") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


# ---- d) Roost popularity ----
popularity_models <- list()
popularity_mod_summary <- list()

for(i in 1:length(functional_relationships)){
  
  model_formula <- paste("is_popular_roost ~ season * ", functional_relationships[[i]], 
                         " + (1|ID) + (1|year)", sep ="")
  
  mod <- glmmTMB(formula = as.formula(model_formula), 
                 REML = F,                              # for model comparison
                 family = binomial(link="logit"),
                 data = roost_df)
  
  popularity_models[[i]] <- mod
  popularity_mod_summary[[i]] <- summary(mod)
  
  rm(mod)
  rm(model_formula)
  
}

popularity_comparison <- compare_performance(popularity_models[[1]], popularity_models[[2]], 
                                             popularity_models[[3]], popularity_models[[4]])
popularity_comparison$Model <- func_rel_names

## Figure 3D ##
popularity_mod_df <- ggpredict(popularity_models[[4]], terms = c("age_scale[all]", "season"))
popularity_mod_df$age <- round(popularity_mod_df$x * roost.age.scale + roost.age.center)
popularity_mod_df <- rename(popularity_mod_df, season = group)

popularity_figure <- 
  ggplot() +
  geom_ribbon(data = popularity_mod_df,
              aes(x = age, y = predicted, 
                  ymin = conf.low, ymax = conf.high, 
                  fill = season), 
              alpha = 0.2, linewidth = 0.6) +
  geom_line(data = popularity_mod_df, 
            aes(x = age, y = predicted, colour = season), 
            linewidth = 0.7) +
  scale_y_continuous(name = "Probability of occupying a popular roost",
                     limits = c(0,1)) +
  xlab("Age (years)") +
  scale_fill_manual(name = "Season", 
                    values =  season_colours) +
  scale_colour_manual(name = "Season", 
                      values = season_colours) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


# ---- e) Normalized degree ----
degree_models <- list()
degree_mod_summary <- list()

for(i in 1:length(functional_relationships)){
  
  model_formula <- paste("norm_degree ~ season + ", functional_relationships[[i]], 
                         " + (1|ID) + (1|year)", sep ="")
  
  mod <- glmmTMB(formula = as.formula(model_formula), 
                 REML = F,                           
                 family = beta_family(link = "logit"),
                 data = social_df)
  
  degree_models[[i]] <- mod
  degree_mod_summary[[i]] <- summary(mod)
  
  rm(mod)
  rm(model_formula)
  
}

degree_comparison <- compare_performance(degree_models[[1]], degree_models[[2]], 
                                         degree_models[[3]], degree_models[[4]])
degree_comparison$Model <- func_rel_names


## Figure 3E ##
degree_mod_df <- ggpredict(degree_models[[1]], terms = c("age_scale[all]", "season"))
degree_mod_df$age <- round(degree_mod_df$x * social.age.scale + social.age.center)
degree_mod_df <- rename(degree_mod_df, season = group)

degree_figure <- 
  ggplot() +
  geom_ribbon(data = degree_mod_df,
              aes(x = age, y = predicted, 
                  ymin = conf.low, ymax = conf.high, 
                  fill = season), 
              alpha = 0.2, linewidth = 0.6) +
  geom_line(data = degree_mod_df, 
            aes(x = age, y = predicted, colour = season), 
            linewidth = 0.7) +
  scale_y_continuous(name = "Normalized degree",
                     limits = c(min(social_df$norm_degree), 
                                max(social_df$norm_degree))) +
  xlab("Age (years)") +
  scale_fill_manual(name = "Season", 
                    values =  season_colours) +
  scale_colour_manual(name = "Season", 
                      values = season_colours) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


# ---- f) Average strength ----
avg_strength_models <- list()
avg_strength_mod_summary <- list()

for(i in 1:length(functional_relationships)){
  
  model_formula <- paste("avg_strength ~ season * ", functional_relationships[[i]], 
                         " + (1|ID) + (1|year)", sep ="")
  
  mod <- glmmTMB(formula = as.formula(model_formula), 
                 REML = F,                               
                 family = gaussian(link = "identity"),
                 data = social_df)
  
  avg_strength_models[[i]] <- mod
  avg_strength_mod_summary[[i]] <- summary(mod)
  
  rm(mod)
  rm(model_formula)
  
}

avg_strength_comparison <- compare_performance(avg_strength_models[[1]], avg_strength_models[[2]], 
                                               avg_strength_models[[3]], avg_strength_models[[4]])
avg_strength_comparison$Model <- func_rel_names


## Figure 3F ##
avg_strength_mod_df <- ggpredict(avg_strength_models[[2]], terms = c("age_scale[all]", "season"))
avg_strength_mod_df$age <- round(avg_strength_mod_df$x * social.age.scale + social.age.center)
avg_strength_mod_df <- rename(avg_strength_mod_df, season = group)

avg_strength_figure <- 
  ggplot() +
  geom_ribbon(data = avg_strength_mod_df,
              aes(x = age, y = predicted, 
                  ymin = conf.low, ymax = conf.high, 
                  fill = season), 
              alpha = 0.2, linewidth = 0.6) +
  geom_line(data = avg_strength_mod_df, 
            aes(x = age, y = predicted, colour = season), 
            linewidth = 0.7) +
  scale_y_continuous(name = "Average strength") +
  xlab("Age (years)") +
  scale_fill_manual(name = "Season", 
                    values =  season_colours) +
  scale_colour_manual(name = "Season", 
                      values = season_colours) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
