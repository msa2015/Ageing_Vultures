
####################################################################################################
## Evaluate the effect of individual plasticity and selective disappearance in shaping the observed
## age-related changes in roost fidelity, roost popularity and average strength 
###################################################################################################

library(glmmTMB)
library(lme4)
library(lmerTest)
library(dplyr)
library(performance)
library(ggeffects)
library(ggplot2)

models <- c("Age", "Age + ID", "Longevity + ID", "Age + Longevity + ID")
model_colours <- c("Age" = "#f3df6c", "Age + ID" = "#ceab07", 
                   "Longevity + ID" = "#c9c9c7", "Age + Longevity + ID" = "black")

# ---- Upload datasets ----
# If necessary, download the datasets from Zenodo
roost_df <- read.csv("Data/df_roost_movement.csv")
roost_long_df <- droplevels(subset(roost_df, !is.na(longevity_5x)))
roost_long_df$age_scale <- scale(roost_long_df$age, center = TRUE, scale = TRUE)
roost_long_df$longevity_5x_scale <- scale(roost_long_df$longevity_5x, center = TRUE, scale = TRUE)

social_df <- read.csv("Data/df_social.csv")
social_long_df <- droplevels(subset(social_df, !is.na(longevity_5x)))
social_long_df$age_scale <- scale(social_long_df$age, center = TRUE, scale = TRUE)
social_long_df$longevity_5x_scale <- scale(social_long_df$longevity_5x, center = TRUE, scale = TRUE)


# ---- a) Roost fidelity and popularity ----
columns_estimate = c('analysis', 'model', 'variable','polynomial', 'estimate', 'conf.low', 'conf.high')
movement_estimates_df = data.frame(matrix(nrow = 22, 
                                          ncol = length(columns_estimate))) 
colnames(movement_estimates_df) = columns_estimate
movement_estimates_df$analysis[c(1:11)] <- "Roost fidelity"
movement_estimates_df$analysis[c(12:22)] <- "Roost popularity"

variables <- list("roost_fid", "is_popular_roost")

movement_mechanisms <- list("* poly(age_scale, 3, raw = TRUE)",
                            "* poly(age_scale, 3, raw = TRUE) + (1|ID)", 
                            "+ longevity_5x_scale + (1|ID)", 
                            "* poly(age_scale, 3, raw = TRUE) + longevity_5x_scale + (1|ID)")

movement_mech_models <- list()
movement_mech_mod_summary <- list()

j <- 1

for(i in 1:length(variables)){
  for(k in 1:length(models)){
    model_formula <- paste(variables[i], " ~ season ", movement_mechanisms[[k]], 
                           " + (1|year)", sep = "")
    
    mod <- glmmTMB(formula = as.formula(model_formula), 
                   REML = F,                              # for model comparison
                   family = binomial(link="logit"),
                   data = roost_long_df)
    
    mod_df <- broom.mixed::tidy(mod, conf.int = TRUE)
    
    movement_mech_models[[k]] <- mod
    movement_mech_mod_summary[[k]] <- summary(mod)
    
    if(k < 3){
      
      movement_estimates_df$variable[j] <- mod_df$term[4]
      movement_estimates_df$estimate[j] <- mod_df$estimate[4]
      movement_estimates_df$conf.low[j] <- mod_df$conf.low[4]
      movement_estimates_df$conf.high[j] <- mod_df$conf.high[4]
      movement_estimates_df$polynomial[j] <- 1
      
      movement_estimates_df$variable[j+1] <- mod_df$term[5]
      movement_estimates_df$estimate[j+1] <- mod_df$estimate[5]
      movement_estimates_df$conf.low[j+1] <- mod_df$conf.low[5]
      movement_estimates_df$conf.high[j+1] <- mod_df$conf.high[5]
      movement_estimates_df$polynomial[j+1] <- 2
      
      movement_estimates_df$variable[j+2] <- mod_df$term[6]
      movement_estimates_df$estimate[j+2] <- mod_df$estimate[6]
      movement_estimates_df$conf.low[j+2] <- mod_df$conf.low[6]
      movement_estimates_df$conf.high[j+2] <- mod_df$conf.high[6]
      movement_estimates_df$polynomial[j+2] <- 3
      
      movement_estimates_df$model[c(j:c(j+2))] <- models[k]
      
      j <- j + 3
    }
    
    if(k == 3){
      movement_estimates_df$variable[j] <- mod_df$term[4]
      movement_estimates_df$estimate[j] <- mod_df$estimate [4]
      movement_estimates_df$conf.low[j] <- mod_df$conf.low[4]
      movement_estimates_df$conf.high[j] <- mod_df$conf.high[4]
      movement_estimates_df$polynomial[j] <- 1
      
      movement_estimates_df$model[j] <- models[k]
      
      j <- j + 1
      
    }
    
    if(k == 4){
      movement_estimates_df$variable[j] <- mod_df$term[4]
      movement_estimates_df$estimate[j] <- mod_df$estimate [4]
      movement_estimates_df$conf.low[j] <- mod_df$conf.low[4]
      movement_estimates_df$conf.high[j] <- mod_df$conf.high[4]
      movement_estimates_df$polynomial[j] <- 1
      
      movement_estimates_df$variable[j+1] <- mod_df$term[5]
      movement_estimates_df$estimate[j+1] <- mod_df$estimate[5]
      movement_estimates_df$conf.low[j+1] <- mod_df$conf.low[5]
      movement_estimates_df$conf.high[j+1] <- mod_df$conf.high[5]
      movement_estimates_df$polynomial[j+1] <- 2
      
      movement_estimates_df$variable[j+2] <- mod_df$term[6]
      movement_estimates_df$estimate[j+2] <- mod_df$estimate[6]
      movement_estimates_df$conf.low[j+2] <- mod_df$conf.low[6]
      movement_estimates_df$conf.high[j+2] <- mod_df$conf.high[6]
      movement_estimates_df$polynomial[j+2] <- 3
      
      movement_estimates_df$variable[j+3] <- mod_df$term[7]
      movement_estimates_df$estimate[j+3] <- mod_df$estimate[7]
      movement_estimates_df$conf.low[j+3] <- mod_df$conf.low[7]
      movement_estimates_df$conf.high[j+3] <- mod_df$conf.high[7]
      movement_estimates_df$polynomial[j+3] <- 1
      
      movement_estimates_df$model[c(j:c(j+3))] <- models[k]
      
      j <- j + 4
      
    }
  }
  rm(mod)
  rm(model_formula)
}

# Figure 4
movement_estimates_df$variable <- ifelse(grepl("poly", movement_estimates_df$variable, fixed = TRUE), 
                                         "Age", "Longevity")
movement_estimates_df$variable_2 <- paste(movement_estimates_df$variable, 
                                          movement_estimates_df$polynomial, sep = "_")
movement_estimates_df$variable_2 <- factor(movement_estimates_df$variable_2,
                                           levels = c("Age_3", "Age_2", "Age_1", "Longevity_1"))

roost_mechanism_plot <- ggplot(subset(movement_estimates_df, analysis == "Roost fidelity"),
                               aes(x = variable_2, y = estimate, 
                                   colour = model, fill = model)) +
  geom_errorbar(data = subset(movement_estimates_df, analysis == "Roost fidelity"), 
                aes(x = variable_2, 
                    ymin = conf.low, ymax = conf.high, colour = model),
                width = 0, position = position_dodge(0.5), linewidth = 1.5, alpha = 0.8) +
  geom_point(position = position_dodge(0.5), size = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  ylab("Effect size") +
  coord_flip() +
  scale_x_discrete(name = "Model predictors",
                   labels = c("Age_1" = expression(Age ~ 1^st ~ polynomial),
                              "Age_2" = expression(Age ~ 2^nd ~ polynomial),
                              "Age_3" = expression(Age ~ 3^rd ~ polynomial),
                              "Longevity_1" = "Longevity")) +
  scale_colour_manual(name = "Model",
                      values = model_colours,
                      breaks = c("Age + Longevity + ID", "Age + ID", "Age"))+
  scale_fill_manual(values = model_colours,
                    breaks = c("Age + Longevity + ID", "Age + ID", "Age"))+
  theme_bw(base_size = 18)+
  theme(legend.position = "none",
        panel.grid = element_blank())


# Extended Data 2A
popularity_mechanism_plot <- ggplot(subset(movement_estimates_df, analysis == "Roost popularity"),
                                    aes(x = variable_2, y = estimate, 
                                        colour = model, fill = model)) +
  geom_errorbar(data = subset(movement_estimates_df, analysis == "Roost popularity"), 
                aes(x = variable_2, 
                    ymin = conf.low, ymax = conf.high, colour = model),
                width = 0, position = position_dodge(0.5), linewidth = 1.5, alpha = 0.8) +
  geom_point(position = position_dodge(0.5), size = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  ylab("Effect size") +
  coord_flip() +
  scale_x_discrete(name = "Model predictors",
                   labels = c("Age_1" = expression(Age ~ 1^st ~ polynomial),
                              "Age_2" = expression(Age ~ 2^nd ~ polynomial),
                              "Age_3" = expression(Age ~ 3^rd ~ polynomial),
                              "Longevity_1" = "Longevity")) +
  scale_colour_manual(name = "Model",
                      values = model_colours,
                      breaks = c("Age + Longevity + ID", "Age + ID", "Age"))+
  scale_fill_manual(values = model_colours,
                    breaks = c("Age + Longevity + ID", "Age + ID", "Age"))+
  theme_bw(base_size = 18)+
  theme(legend.position = "none",
        panel.grid = element_blank())


# ---- b) Average strength ----
columns_estimate = c('analysis', 'model', 'variable','polynomial', 'estimate', 'conf.low', 'conf.high', 'max.vif')
social_estimates_df = data.frame(matrix(nrow = 8,
                                        ncol = length(columns_estimate))) 
colnames(social_estimates_df) = columns_estimate
social_estimates_df$analysis <- "Average strength"

social_mechanisms <- list("+ poly(age_scale, 2, raw = TRUE)",
                          "+ poly(age_scale, 2, raw = TRUE) + (1|ID)", 
                          "+ longevity_5x_scale + (1|ID)", 
                          "+ poly(age_scale, 2, raw = TRUE) + longevity_5x_scale + (1|ID)")

social_mech_models <- list()
social_mech_mod_summary <- list()

j <- 1

for(i in 1:length(models)){
  model_formula <- paste("avg_strength ~ season ", social_mechanisms[[i]], 
                         " + (1|year)", sep = "")
  
  mod <- glmmTMB(formula = as.formula(model_formula), 
                 REML = F,                              # for model comparison
                 family = Gamma(link = "inverse"),
                 data = social_long_df)
  
  mod_df <- broom.mixed::tidy(mod, conf.int = TRUE)
  
  social_mech_models[[i]] <- mod
  social_mech_mod_summary[[i]] <- summary(mod)
  
  if(i < 3){
    
    social_estimates_df$variable[j] <- mod_df$term[4]
    social_estimates_df$estimate[j] <- mod_df$estimate[4]
    social_estimates_df$conf.low[j] <- mod_df$conf.low[4]
    social_estimates_df$conf.high[j] <- mod_df$conf.high[4]
    social_estimates_df$polynomial[j] <- 1
    
    social_estimates_df$variable[j+1] <- mod_df$term[5]
    social_estimates_df$estimate[j+1] <- mod_df$estimate[5]
    social_estimates_df$conf.low[j+1] <- mod_df$conf.low[5]
    social_estimates_df$conf.high[j+1] <- mod_df$conf.high[5]
    social_estimates_df$polynomial[j+1] <- 2
    
    social_estimates_df$model[c(j:c(j+1))] <- models[i]
    social_estimates_df$max.vif[c(j:c(j+1))] <- max(check_collinearity(mod)$VIF)
    
    j <- j + 2
    
  }
  
  if(i == 3){
    social_estimates_df$variable[j] <- mod_df$term[4]
    social_estimates_df$estimate[j] <- mod_df$estimate [4]
    social_estimates_df$conf.low[j] <- mod_df$conf.low[4]
    social_estimates_df$conf.high[j] <- mod_df$conf.high[4]
    social_estimates_df$polynomial[j] <- 1
    
    social_estimates_df$model[j] <- models[i]
    social_estimates_df$max.vif[j] <-  max(check_collinearity(mod)$VIF)
    
    j <- j + 1
    
  }
  
  if(i == 4){
    social_estimates_df$variable[j] <- mod_df$term[4]
    social_estimates_df$estimate[j] <- mod_df$estimate [4]
    social_estimates_df$conf.low[j] <- mod_df$conf.low[4]
    social_estimates_df$conf.high[j] <- mod_df$conf.high[4]
    social_estimates_df$polynomial[j] <- 1
    
    social_estimates_df$variable[j+1] <- mod_df$term[5]
    social_estimates_df$estimate[j+1] <- mod_df$estimate[5]
    social_estimates_df$conf.low[j+1] <- mod_df$conf.low[5]
    social_estimates_df$conf.high[j+1] <- mod_df$conf.high[5]
    social_estimates_df$polynomial[j+1] <- 2
    
    social_estimates_df$variable[j+2] <- mod_df$term[7]
    social_estimates_df$estimate[j+2] <- mod_df$estimate[7]
    social_estimates_df$conf.low[j+2] <- mod_df$conf.low[7]
    social_estimates_df$conf.high[j+2] <- mod_df$conf.high[7]
    social_estimates_df$polynomial[j+2] <- 1
    
    social_estimates_df$model[c(j:c(j+2))] <- models[i]
    social_estimates_df$max.vif[c(j:c(j+2))] <-  max(check_collinearity(mod)$VIF)
    
    j <- j + 3
    
  }
  rm(mod)
  rm(model_formula)
}

# Extended Data 2B
social_estimates_df$variable <- ifelse(grepl("poly", social_estimates_df$variable, fixed = TRUE), 
                                       "Age", "Longevity")
social_estimates_df$variable_2 <- paste(social_estimates_df$variable, 
                                        social_estimates_df$polynomial, sep = "_")
social_estimates_df$variable_2 <- factor(social_estimates_df$variable_2,
                                         levels = c("Age_3", "Age_2", "Age_1", "Longevity_1"))

sociality_mechanism_plot <- ggplot(subset(social_estimates_df, max.vif < 5),
                                   aes(x = variable_2, y = estimate, 
                                       colour = model, fill = model)) +
  geom_errorbar(data = subset(social_estimates_df, max.vif < 5), 
                aes(x = variable_2, 
                    ymin = conf.low, ymax = conf.high, colour = model),
                width = 0, position = position_dodge(0.5), linewidth = 1.5, alpha = 0.8) +
  geom_point(position = position_dodge(0.5), size = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  ylab("Effect estimate") +
  coord_flip() +
  scale_x_discrete(name = "Predictors",
                   labels = c("Age_1" = "Age - linear term",
                              "Age_2" = "Age - quadratic term",
                              "Longevity_1" = "Longevity")) +
  scale_colour_manual(name = "Model",
                      values = model_colours,
                      breaks = c("Age + Longevity + ID", "Age + ID", "Age"))+
  scale_fill_manual(values = model_colours,
                    breaks = c("Age + Longevity + ID", "Age + ID", "Age"))+
  theme_bw(base_size = 18)+
  theme(legend.position = "none",
        panel.grid = element_blank())
