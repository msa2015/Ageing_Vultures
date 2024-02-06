
################################################################################
## Sensitivity analyses 
################################################################################

# ---- Supplementary Figure S14 - Excluding individuals of unknown age ----
# Roost fidelity
fidelity_known_age <- glmmTMB(roost_fid ~ season * poly(age_scale, 3, raw=TRUE) + 
                                (1|ID) + (1|year),
                              REML = F,
                              family = binomial(link = "logit"),
                              data = subset(roost_df, real_age == 1))
summary(fidelity_known_age)

fidelity_known_age_plot <- plot_model(
  fidelity_known_age, title = "Roost fidelity", vline.color = "black",
  transform =  NULL, show.values = TRUE, value.offset = .3) + theme_bw(base_size = 16)


# Routine
routine_known_age <- glmmTMB(with_routine ~ age_scale + log(seq_length_scale) + 
                               (1|ID) + (1|ID),
                             REML = F,
                             family = binomial(link = "logit"),
                             data = subset(routine_df, real_age == 1))
summary(routine_known_age)

routine_known_age_plot <- plot_model(routine_known_age, title = "Routine", 
                                     vline.color = "black", transform = NULL,
                                     show.values = TRUE, value.offset = .3,
                                     axis.labels	= c("Sequence Length", "Age"), 
                                     axis.lim = c(-1, 5)) + 
  theme_bw(base_size = 16)


# Routine index
routine_index_known_age <- glmmTMB(R ~ expm1(age_scale) + (1|ID) + (1|year),
                                   family = beta_family(), 
                                   REML = F,
                                   data = subset(routine_index_df, real_age == 1))
summary(routine_index_known_age) 

routine_index_known_age_plot <- plot_model(routine_index_known_age, 
                                           title = "Routine index", vline.color = "black",
                                           show.values = TRUE, value.offset = .3,
                                           axis.labels = c("Age"), transform = NULL,
                                           axis.lim = c(-0.2, 0.2)) + 
  theme_bw(base_size = 16)


# Roost popularity
popularity_known_age <- glmmTMB(is_popular_roost ~ poly(age_scale, 3, raw = TRUE) * season + 
                                  (1|ID) + (1|year),
                                family = binomial(link="logit"),
                                data = subset(roost_df, real_age == 1))
summary(popularity_known_age)

popularity_known_age_plot <- plot_model(popularity_known_age, title = "Popular roosts", 
                                        vline.color = "black", transform = NULL, 
                                        show.values = TRUE, value.offset = .3) + 
  theme_bw(base_size = 16)


# Average strength
avg_strength_known_age <- glmmTMB(avg_strength ~ poly(age_scale,2,raw=TRUE) * season +
                                    (1|ID) + (1|year),
                                  family = gaussian(),
                                  REML = F,
                                  data = subset(social_df, real_age == 1))
summary(avg_strength_known_age)

avg_strength_known_age_plot <- plot_model(avg_strength_known_age, 
                                          title = "Average strength", vline.color = "black", 
                                          show.values = TRUE, value.offset = .3) + ylim(-0.1, 0.1) + 
  theme_bw(base_size = 16)


# ---- Supplementary Figure S15 - No LRFs and different threholds for roost fidelity -----
fidelity_no_lrf <- glmmTMB(roost_fid ~ season * poly(age_scale, 3, raw=TRUE) + 
                             (1|ID) + (1|year),
                           REML = F,
                           family = binomial(link = "logit"),
                           data = subset(roost_df, is_lrf == 0))
summary(fidelity_no_lrf)

fidelity_no_lrf_plot <- plot_model(
  fidelity_no_lrf, title = "Roost fidelity - no LRF", vline.color = "black",
  transform =  NULL, show.values = TRUE, value.offset = .3) + theme_bw(base_size = 16)


fidelity_1km <- glmmTMB(roost_fid_1km ~ season * poly(age_scale, 3, raw=TRUE) + 
                          (1|ID) + (1|year),
                        REML = F,
                        family = binomial(link = "logit"),
                        data = roost_df)
summary(fidelity_1km)

fidelity_1km_plot <- plot_model(
  fidelity_1km, title = "Roost fidelity - 1km", vline.color = "black",
  transform =  NULL, show.values = TRUE, value.offset = .3) + theme_bw(base_size = 16)


fidelity_20km <- glmmTMB(roost_fid_20km ~ season * poly(age_scale, 3, raw=TRUE) + 
                           (1|ID) + (1|year),
                         REML = F,
                         family = binomial(link = "logit"),
                         data = roost_df)
summary(fidelity_20km)

fidelity_20km_plot <- plot_model(
  fidelity_20km, title = "Roost fidelity - 20km", vline.color = "black",
  transform =  NULL, show.values = TRUE, value.offset = .3) + theme_bw(base_size = 16)


# ---- Supplementary Figure S17 - Different popular roosts thresholds and no LRF ----
popularity_25pc <- glmmTMB(popular_roost_25 ~ poly(age_scale, 3, raw = TRUE) * season + 
                             (1|ID) + (1|year),
                           family = binomial(link="logit"),
                           data = roost_df)
summary(popularity_25pc)

popularity_25pc_plot <- plot_model(popularity_25pc, title = "Roost popularity - 25% roosts", 
                                   vline.color = "black", transform = NULL, 
                                   show.values = TRUE, value.offset = .3) + 
  theme_bw(base_size = 16)


popularity_30pc <- glmmTMB(popular_roost_30 ~ poly(age_scale, 3, raw = TRUE) * season + 
                             (1|ID) + (1|year),
                           family = binomial(link="logit"),
                           data = roost_df)
summary(popularity_30pc)

popularity_30pc_plot <- plot_model(popularity_30pc, title = "Roost popularity - 30% roosts", 
                                   vline.color = "black", transform = NULL, 
                                   show.values = TRUE, value.offset = .3) + 
  theme_bw(base_size = 16)


popularity_no_lrf <- glmmTMB(is_popular_roost ~ season * poly(age_scale, 3, raw=TRUE) + 
                               (1|ID) + (1|year),
                             REML = F,
                             family = binomial(link = "logit"),
                             data = subset(roost_df, is_lrf == 0))
summary(popularity_no_lrf)

popularity_no_lrf_plot <- plot_model(popularity_no_lrf, title = "Roost popularity - no LRF", 
                                     vline.color = "black", transform = NULL, 
                                     show.values = TRUE, value.offset = .3) + 
  theme_bw(base_size = 16)


# ---- Supplementary Figure S18 - Are old individuals driving roost popularity? ---- 
# Randomly get 300 points for each age class (100 per season)
sub_sample <- roost_df %>%
  filter(in_study_area == 1) %>%
  group_by(age, season) %>%
  slice_sample(n = 100, replace = FALSE)

# Remove ages 0 and 24: they do not have all seasons, so they would be underrepresented
sub_sample <- filter(sub_sample, !(age %in% c(0,24)))

# Get the new popular roosts: 20% (or 6) of the total number of roosts within the study area
popular_roosts_ss <- sub_sample %>% 
  group_by(ID_roost_site) %>%
  summarise(n_points_ss = n()) %>%
  arrange(-n_points_ss) %>%
  head(6)

roost_df <- roost_df %>% 
  mutate(in_popular_roost_ss = ifelse(ID_roost_site %in% popular_roosts_ss$ID_roost_site, 1, 0))


popularity_sub_sample <- glmmTMB(in_popular_roost_ss ~ season * poly(age_scale, 3, raw=TRUE) + 
                                   (1|ID) + (1|year),
                                 REML = F,
                                 family = binomial(link = "logit"),
                                 data = roost_df)
summary(popularity_sub_sample)

popularity_sub_sample_plot <- plot_model(popularity_sub_sample, 
                                         vline.color = "black", transform = NULL, 
                                         show.values = TRUE, value.offset = .3) + 
  theme_bw(base_size = 16)


# ---- Supplementary Figure S20 - Longevity estimated as 2 years after last observation ----
roost_long_2y <- droplevels(subset(roost_df, !is.na(longevity_2y)))
roost_long_2y$age_scale <- scale(roost_long_2y$age, center = TRUE, scale = TRUE)
roost_long_2y$longevity_2y_scale <- scale(roost_long_2y$longevity_2y, center = TRUE, scale = TRUE)

social_long_2y <- droplevels(subset(social_df, !is.na(longevity_2y)))
social_long_2y$age_scale <- scale(social_long_2y$age, center = TRUE, scale = TRUE)
social_long_2y$longevity_2y_scale <- scale(social_long_2y$longevity_2y, center = TRUE, scale = TRUE)

models <- c("Age", "Age + ID", "Longevity + ID", "Age + Longevity + ID")
model_colours <- c("Age" = "#f3df6c", "Age + ID" = "#ceab07", 
                   "Longevity + ID" = "#c9c9c7", "Age + Longevity + ID" = "black")


### Roost fidelity and Roost popularity
columns_estimate = c('analysis', 'model', 'variable','polynomial', 'estimate', 'conf.low', 'conf.high')
movement_estimates_2y = data.frame(matrix(nrow = 22, 
                                          ncol = length(columns_estimate))) 
colnames(movement_estimates_2y) = columns_estimate
movement_estimates_2y$analysis[c(1:11)] <- "Roost fidelity"
movement_estimates_2y$analysis[c(12:22)] <- "Roost popularity"

variables <- list("roost_fid", "is_popular_roost")

movement_mechanisms <- list("* poly(age_scale, 3, raw = TRUE)",
                            "* poly(age_scale, 3, raw = TRUE) + (1|ID)", 
                            "+ longevity_2y_scale + (1|ID)", 
                            "* poly(age_scale, 3, raw = TRUE) + longevity_2y_scale + (1|ID)")

movement_mech_2y_models <- list()
movement_mech_2y_mod_summary <- list()

j <- 1

for(i in 1:length(variables)){
  for(k in 1:length(models)){
    model_formula <- paste(variables[i], " ~ season ", movement_mechanisms[[k]], 
                           " + (1|year)", sep = "")
    
    mod <- glmmTMB(formula = as.formula(model_formula), 
                   REML = F,                              # for model comparison
                   family = binomial(link="logit"),
                   data = roost_long_2y)
    
    mod_df <- broom.mixed::tidy(mod, conf.int = TRUE)
    
    movement_mech_2y_models[[k]] <- mod
    movement_mech_2y_mod_summary[[k]] <- summary(mod)
    
    if(k < 3){
      
      movement_estimates_2y$variable[j] <- mod_df$term[4]
      movement_estimates_2y$estimate[j] <- mod_df$estimate[4]
      movement_estimates_2y$conf.low[j] <- mod_df$conf.low[4]
      movement_estimates_2y$conf.high[j] <- mod_df$conf.high[4]
      movement_estimates_2y$polynomial[j] <- 1
      
      movement_estimates_2y$variable[j+1] <- mod_df$term[5]
      movement_estimates_2y$estimate[j+1] <- mod_df$estimate[5]
      movement_estimates_2y$conf.low[j+1] <- mod_df$conf.low[5]
      movement_estimates_2y$conf.high[j+1] <- mod_df$conf.high[5]
      movement_estimates_2y$polynomial[j+1] <- 2
      
      movement_estimates_2y$variable[j+2] <- mod_df$term[6]
      movement_estimates_2y$estimate[j+2] <- mod_df$estimate[6]
      movement_estimates_2y$conf.low[j+2] <- mod_df$conf.low[6]
      movement_estimates_2y$conf.high[j+2] <- mod_df$conf.high[6]
      movement_estimates_2y$polynomial[j+2] <- 3
      
      movement_estimates_2y$model[c(j:c(j+2))] <- models[k]
      
      j <- j + 3
    }
    
    if(k == 3){
      movement_estimates_2y$variable[j] <- mod_df$term[4]
      movement_estimates_2y$estimate[j] <- mod_df$estimate [4]
      movement_estimates_2y$conf.low[j] <- mod_df$conf.low[4]
      movement_estimates_2y$conf.high[j] <- mod_df$conf.high[4]
      movement_estimates_2y$polynomial[j] <- 1
      
      movement_estimates_2y$model[j] <- models[k]
      
      j <- j + 1
      
    }
    
    if(k == 4){
      movement_estimates_2y$variable[j] <- mod_df$term[4]
      movement_estimates_2y$estimate[j] <- mod_df$estimate [4]
      movement_estimates_2y$conf.low[j] <- mod_df$conf.low[4]
      movement_estimates_2y$conf.high[j] <- mod_df$conf.high[4]
      movement_estimates_2y$polynomial[j] <- 1
      
      movement_estimates_2y$variable[j+1] <- mod_df$term[5]
      movement_estimates_2y$estimate[j+1] <- mod_df$estimate[5]
      movement_estimates_2y$conf.low[j+1] <- mod_df$conf.low[5]
      movement_estimates_2y$conf.high[j+1] <- mod_df$conf.high[5]
      movement_estimates_2y$polynomial[j+1] <- 2
      
      movement_estimates_2y$variable[j+2] <- mod_df$term[6]
      movement_estimates_2y$estimate[j+2] <- mod_df$estimate[6]
      movement_estimates_2y$conf.low[j+2] <- mod_df$conf.low[6]
      movement_estimates_2y$conf.high[j+2] <- mod_df$conf.high[6]
      movement_estimates_2y$polynomial[j+2] <- 3
      
      movement_estimates_2y$variable[j+3] <- mod_df$term[7]
      movement_estimates_2y$estimate[j+3] <- mod_df$estimate[7]
      movement_estimates_2y$conf.low[j+3] <- mod_df$conf.low[7]
      movement_estimates_2y$conf.high[j+3] <- mod_df$conf.high[7]
      movement_estimates_2y$polynomial[j+3] <- 1
      
      movement_estimates_2y$model[c(j:c(j+3))] <- models[k]
      
      j <- j + 4
      
    }
  }
  rm(mod)
  rm(model_formula)
}

# Supplementary Figure S20 - Roost fidelity
movement_estimates_2y$variable <- ifelse(grepl("poly", movement_estimates_2y$variable, fixed = TRUE), 
                                         "Age", "Longevity")
movement_estimates_2y$variable_2 <- paste(movement_estimates_2y$variable, 
                                          movement_estimates_2y$polynomial, sep = "_")
movement_estimates_2y$variable_2 <- factor(movement_estimates_2y$variable_2,
                                           levels = c("Age_3", "Age_2", "Age_1", "Longevity_1"))

roost_mechanism_2y_plot <- ggplot(subset(movement_estimates_2y, analysis == "Roost fidelity"),
                                  aes(x = variable_2, y = estimate, 
                                      colour = model, fill = model)) +
  geom_errorbar(data = subset(movement_estimates_2y, analysis == "Roost fidelity"), 
                aes(x = variable_2, 
                    ymin = conf.low, ymax = conf.high, colour = model),
                width = 0, position = position_dodge(0.5), linewidth = 1.5, alpha = 0.8) +
  geom_point(position = position_dodge(0.5), size = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  ylab("Effect estimate") +
  coord_flip() +
  scale_x_discrete(name = "Predictors",
                   labels = c("Age_1" = expression(Age ~ 1^st ~ polynomial),
                              "Age_2" = expression(Age ~ 2^nd ~ polynomial),
                              "Age_3" = expression(Age ~ 3^rd ~ polynomial),
                              "Longevity_1" = "Longevity")) +
  scale_colour_manual(name = "Model",
                      values = model_colours,
                      breaks = c("Age + Longevity + ID", "Age + ID", "Age")) +
  scale_fill_manual(values = model_colours,
                    breaks = c("Age + Longevity + ID", "Age + ID", "Age")) +
  ggtitle("Roost fidelity - 2 years after last observation") +
  theme_bw(base_size = 18)+
  theme(legend.position = "none",
        panel.grid = element_blank())


# Supplementary Figure S20 - Roost popularity
popularity_mechanism_2y_plot <- ggplot(subset(movement_estimates_2y, analysis == "Roost popularity"),
                                       aes(x = variable_2, y = estimate, 
                                           colour = model, fill = model)) +
  geom_errorbar(data = subset(movement_estimates_2y, analysis == "Roost popularity"), 
                aes(x = variable_2, 
                    ymin = conf.low, ymax = conf.high, colour = model),
                width = 0, position = position_dodge(0.5), linewidth = 1.5, alpha = 0.8) +
  geom_point(position = position_dodge(0.5), size = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  ylab("Effect estimate") +
  coord_flip() +
  scale_x_discrete(name = "Predictors",
                   labels = c("Age_1" = expression(Age ~ 1^st ~ polynomial),
                              "Age_2" = expression(Age ~ 2^nd ~ polynomial),
                              "Age_3" = expression(Age ~ 3^rd ~ polynomial),
                              "Longevity_1" = "Longevity")) +
  scale_colour_manual(name = "Model",
                      values = model_colours,
                      breaks = c("Age + Longevity + ID", "Age + ID", "Age")) +
  scale_fill_manual(values = model_colours,
                    breaks = c("Age + Longevity + ID", "Age + ID", "Age")) +
  ggtitle("Roost popularity - 2 years after last observation") +
  theme_bw(base_size = 18)+
  theme(legend.position = "none",
        panel.grid = element_blank())


### Average strength
columns_estimate = c('analysis', 'model', 'variable','polynomial', 'estimate', 'conf.low', 'conf.high', 'max.vif')
social_estimates_2y = data.frame(matrix(nrow = 8,
                                        ncol = length(columns_estimate))) 
colnames(social_estimates_2y) = columns_estimate
social_estimates_2y$analysis <- "Average strength"

social_mechanisms <- list("+ poly(age_scale, 2, raw = TRUE)",
                          "+ poly(age_scale, 2, raw = TRUE) + (1|ID)", 
                          "+ longevity_2y_scale + (1|ID)", 
                          "+ poly(age_scale, 2, raw = TRUE) + longevity_2y_scale + (1|ID)")

social_mech_2y_models <- list()
social_mech_2y_mod_summary <- list()

j <- 1

for(i in 1:length(models)){
  model_formula <- paste("avg_strength ~ season ", social_mechanisms[[i]], 
                         " + (1|year)", sep = "")
  
  mod <- glmmTMB(formula = as.formula(model_formula), 
                 REML = F,                              # for model comparison
                 family = Gamma(link = "inverse"),
                 data = social_long_2y)
  
  mod_df <- broom.mixed::tidy(mod, conf.int = TRUE)
  
  social_mech_2y_models[[i]] <- mod
  social_mech_2y_mod_summary[[i]] <- summary(mod)
  
  if(i < 3){
    
    social_estimates_2y$variable[j] <- mod_df$term[4]
    social_estimates_2y$estimate[j] <- mod_df$estimate[4]
    social_estimates_2y$conf.low[j] <- mod_df$conf.low[4]
    social_estimates_2y$conf.high[j] <- mod_df$conf.high[4]
    social_estimates_2y$polynomial[j] <- 1
    
    social_estimates_2y$variable[j+1] <- mod_df$term[5]
    social_estimates_2y$estimate[j+1] <- mod_df$estimate[5]
    social_estimates_2y$conf.low[j+1] <- mod_df$conf.low[5]
    social_estimates_2y$conf.high[j+1] <- mod_df$conf.high[5]
    social_estimates_2y$polynomial[j+1] <- 2
    
    social_estimates_2y$model[c(j:c(j+1))] <- models[i]
    social_estimates_2y$max.vif[c(j:c(j+1))] <- max(check_collinearity(mod)$VIF)
    
    j <- j + 2
    
  }
  
  if(i == 3){
    social_estimates_2y$variable[j] <- mod_df$term[4]
    social_estimates_2y$estimate[j] <- mod_df$estimate [4]
    social_estimates_2y$conf.low[j] <- mod_df$conf.low[4]
    social_estimates_2y$conf.high[j] <- mod_df$conf.high[4]
    social_estimates_2y$polynomial[j] <- 1
    
    social_estimates_2y$model[j] <- models[i]
    social_estimates_2y$max.vif[j] <-  max(check_collinearity(mod)$VIF)
    
    j <- j + 1
    
  }
  
  if(i == 4){
    social_estimates_2y$variable[j] <- mod_df$term[4]
    social_estimates_2y$estimate[j] <- mod_df$estimate [4]
    social_estimates_2y$conf.low[j] <- mod_df$conf.low[4]
    social_estimates_2y$conf.high[j] <- mod_df$conf.high[4]
    social_estimates_2y$polynomial[j] <- 1
    
    social_estimates_2y$variable[j+1] <- mod_df$term[5]
    social_estimates_2y$estimate[j+1] <- mod_df$estimate[5]
    social_estimates_2y$conf.low[j+1] <- mod_df$conf.low[5]
    social_estimates_2y$conf.high[j+1] <- mod_df$conf.high[5]
    social_estimates_2y$polynomial[j+1] <- 2
    
    social_estimates_2y$variable[j+2] <- mod_df$term[7]
    social_estimates_2y$estimate[j+2] <- mod_df$estimate[7]
    social_estimates_2y$conf.low[j+2] <- mod_df$conf.low[7]
    social_estimates_2y$conf.high[j+2] <- mod_df$conf.high[7]
    social_estimates_2y$polynomial[j+2] <- 1
    
    social_estimates_2y$model[c(j:c(j+2))] <- models[i]
    social_estimates_2y$max.vif[c(j:c(j+2))] <-  max(check_collinearity(mod)$VIF)
    
    j <- j + 3
    
  }
  rm(mod)
  rm(model_formula)
}

# Supplementary Figure S20 - Average strength
social_estimates_2y$variable <- ifelse(grepl("poly", social_estimates_2y$variable, fixed = TRUE), 
                                       "Age", "Longevity")
social_estimates_2y$variable_2 <- paste(social_estimates_2y$variable, 
                                        social_estimates_2y$polynomial, sep = "_")
social_estimates_2y$variable_2 <- factor(social_estimates_2y$variable_2,
                                         levels = c("Age_3", "Age_2", "Age_1", "Longevity_1"))

sociality_mechanism_2y_plot <- ggplot(subset(social_estimates_2y, max.vif < 5),
                                      aes(x = variable_2, y = estimate, 
                                          colour = model, fill = model)) +
  geom_errorbar(data = subset(social_estimates_2y, max.vif < 5), 
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
                      breaks = c("Age + Longevity + ID", "Age + ID", "Age")) +
  scale_fill_manual(values = model_colours,
                    breaks = c("Age + Longevity + ID", "Age + ID", "Age")) +
  ggtitle("Average strength - 2 years after last observation") +
  theme_bw(base_size = 18)+
  theme(legend.position = "none",
        panel.grid = element_blank())


# ---- Supplementary Figure S21 - Longevity estimated as 10 x the usual observation rate ----
roost_long_10x <- droplevels(subset(roost_df, !is.na(longevity_10x)))
roost_long_10x$age_scale <- scale(roost_long_10x$age, center = TRUE, scale = TRUE)
roost_long_10x$longevity_10x_scale <- scale(roost_long_10x$longevity_10x, center = TRUE, scale = TRUE)

social_long_10x <- droplevels(subset(social_df, !is.na(longevity_10x)))
social_long_10x$age_scale <- scale(social_long_10x$age, center = TRUE, scale = TRUE)
social_long_10x$longevity_10x_scale <- scale(social_long_10x$longevity_10x, center = TRUE, scale = TRUE)

models <- c("Age", "Age + ID", "Longevity + ID", "Age + Longevity + ID")
model_colours <- c("Age" = "#f3df6c", "Age + ID" = "#ceab07", 
                   "Longevity + ID" = "#c9c9c7", "Age + Longevity + ID" = "black")

### Roost fidelity and Roost popularity
columns_estimate = c('analysis', 'model', 'variable','polynomial', 'estimate', 'conf.low', 'conf.high')
movement_estimates_10x = data.frame(matrix(nrow = 22, 
                                           ncol = length(columns_estimate))) 
colnames(movement_estimates_10x) = columns_estimate
movement_estimates_10x$analysis[c(1:11)] <- "Roost fidelity"
movement_estimates_10x$analysis[c(12:22)] <- "Roost popularity"

variables <- list("roost_fid", "is_popular_roost")

movement_mechanisms <- list("* poly(age_scale, 3, raw = TRUE)",
                            "* poly(age_scale, 3, raw = TRUE) + (1|ID)", 
                            "+ longevity_10x_scale + (1|ID)", 
                            "* poly(age_scale, 3, raw = TRUE) + longevity_10x_scale + (1|ID)")

movement_mech_10x_models <- list()
movement_mech_10x_mod_summary <- list()

j <- 1

for(i in 1:length(variables)){
  for(k in 1:length(models)){
    model_formula <- paste(variables[i], " ~ season ", movement_mechanisms[[k]], 
                           " + (1|year)", sep = "")
    
    mod <- glmmTMB(formula = as.formula(model_formula), 
                   REML = F,                              # for model comparison
                   family = binomial(link="logit"),
                   data = roost_long_10x)
    
    mod_df <- broom.mixed::tidy(mod, conf.int = TRUE)
    
    movement_mech_10x_models[[k]] <- mod
    movement_mech_10x_mod_summary[[k]] <- summary(mod)
    
    if(k < 3){
      
      movement_estimates_10x$variable[j] <- mod_df$term[4]
      movement_estimates_10x$estimate[j] <- mod_df$estimate[4]
      movement_estimates_10x$conf.low[j] <- mod_df$conf.low[4]
      movement_estimates_10x$conf.high[j] <- mod_df$conf.high[4]
      movement_estimates_10x$polynomial[j] <- 1
      
      movement_estimates_10x$variable[j+1] <- mod_df$term[5]
      movement_estimates_10x$estimate[j+1] <- mod_df$estimate[5]
      movement_estimates_10x$conf.low[j+1] <- mod_df$conf.low[5]
      movement_estimates_10x$conf.high[j+1] <- mod_df$conf.high[5]
      movement_estimates_10x$polynomial[j+1] <- 2
      
      movement_estimates_10x$variable[j+2] <- mod_df$term[6]
      movement_estimates_10x$estimate[j+2] <- mod_df$estimate[6]
      movement_estimates_10x$conf.low[j+2] <- mod_df$conf.low[6]
      movement_estimates_10x$conf.high[j+2] <- mod_df$conf.high[6]
      movement_estimates_10x$polynomial[j+2] <- 3
      
      movement_estimates_10x$model[c(j:c(j+2))] <- models[k]
      
      j <- j + 3
    }
    
    if(k == 3){
      movement_estimates_10x$variable[j] <- mod_df$term[4]
      movement_estimates_10x$estimate[j] <- mod_df$estimate [4]
      movement_estimates_10x$conf.low[j] <- mod_df$conf.low[4]
      movement_estimates_10x$conf.high[j] <- mod_df$conf.high[4]
      movement_estimates_10x$polynomial[j] <- 1
      
      movement_estimates_10x$model[j] <- models[k]
      
      j <- j + 1
      
    }
    
    if(k == 4){
      movement_estimates_10x$variable[j] <- mod_df$term[4]
      movement_estimates_10x$estimate[j] <- mod_df$estimate [4]
      movement_estimates_10x$conf.low[j] <- mod_df$conf.low[4]
      movement_estimates_10x$conf.high[j] <- mod_df$conf.high[4]
      movement_estimates_10x$polynomial[j] <- 1
      
      movement_estimates_10x$variable[j+1] <- mod_df$term[5]
      movement_estimates_10x$estimate[j+1] <- mod_df$estimate[5]
      movement_estimates_10x$conf.low[j+1] <- mod_df$conf.low[5]
      movement_estimates_10x$conf.high[j+1] <- mod_df$conf.high[5]
      movement_estimates_10x$polynomial[j+1] <- 2
      
      movement_estimates_10x$variable[j+2] <- mod_df$term[6]
      movement_estimates_10x$estimate[j+2] <- mod_df$estimate[6]
      movement_estimates_10x$conf.low[j+2] <- mod_df$conf.low[6]
      movement_estimates_10x$conf.high[j+2] <- mod_df$conf.high[6]
      movement_estimates_10x$polynomial[j+2] <- 3
      
      movement_estimates_10x$variable[j+3] <- mod_df$term[7]
      movement_estimates_10x$estimate[j+3] <- mod_df$estimate[7]
      movement_estimates_10x$conf.low[j+3] <- mod_df$conf.low[7]
      movement_estimates_10x$conf.high[j+3] <- mod_df$conf.high[7]
      movement_estimates_10x$polynomial[j+3] <- 1
      
      movement_estimates_10x$model[c(j:c(j+3))] <- models[k]
      
      j <- j + 4
      
    }
  }
  rm(mod)
  rm(model_formula)
}

# Supplementary Figure S21 - Roost fidelity 
movement_estimates_10x$variable <- ifelse(grepl("poly", movement_estimates_10x$variable, fixed = TRUE), 
                                          "Age", "Longevity")
movement_estimates_10x$variable_2 <- paste(movement_estimates_10x$variable, 
                                           movement_estimates_10x$polynomial, sep = "_")
movement_estimates_10x$variable_2 <- factor(movement_estimates_10x$variable_2,
                                            levels = c("Age_3", "Age_2", "Age_1", "Longevity_1"))

roost_mechanism_10x_plot <- ggplot(subset(movement_estimates_10x, analysis == "Roost fidelity"),
                                   aes(x = variable_2, y = estimate, 
                                       colour = model, fill = model)) +
  geom_errorbar(data = subset(movement_estimates_10x, analysis == "Roost fidelity"), 
                aes(x = variable_2, 
                    ymin = conf.low, ymax = conf.high, colour = model),
                width = 0, position = position_dodge(0.5), linewidth = 1.5, alpha = 0.8) +
  geom_point(position = position_dodge(0.5), size = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  ylab("Effect estimate") +
  coord_flip() +
  scale_x_discrete(name = "Predictors",
                   labels = c("Age_1" = expression(Age ~ 1^st ~ polynomial),
                              "Age_2" = expression(Age ~ 2^nd ~ polynomial),
                              "Age_3" = expression(Age ~ 3^rd ~ polynomial),
                              "Longevity_1" = "Longevity")) +
  scale_colour_manual(name = "Model",
                      values = model_colours,
                      breaks = c("Age + Longevity + ID", "Age + ID", "Age")) +
  scale_fill_manual(values = model_colours,
                    breaks = c("Age + Longevity + ID", "Age + ID", "Age")) +
  ggtitle("Roost fidelity - 10 times observation rate") +
  theme_bw(base_size = 18)+
  theme(legend.position = "none",
        panel.grid = element_blank())


# Supplementary Figure S21 - Roost popularity
popularity_mechanism_10x_plot <- ggplot(subset(movement_estimates_10x, analysis == "Roost popularity"),
                                        aes(x = variable_2, y = estimate, 
                                            colour = model, fill = model)) +
  geom_errorbar(data = subset(movement_estimates_10x, analysis == "Roost popularity"), 
                aes(x = variable_2, 
                    ymin = conf.low, ymax = conf.high, colour = model),
                width = 0, position = position_dodge(0.5), linewidth = 1.5, alpha = 0.8) +
  geom_point(position = position_dodge(0.5), size = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  ylab("Effect estimate") +
  coord_flip() +
  scale_x_discrete(name = "Predictors",
                   labels = c("Age_1" = expression(Age ~ 1^st ~ polynomial),
                              "Age_2" = expression(Age ~ 2^nd ~ polynomial),
                              "Age_3" = expression(Age ~ 3^rd ~ polynomial),
                              "Longevity_1" = "Longevity")) +
  scale_colour_manual(name = "Model",
                      values = model_colours,
                      breaks = c("Age + Longevity + ID", "Age + ID", "Age")) +
  scale_fill_manual(values = model_colours,
                    breaks = c("Age + Longevity + ID", "Age + ID", "Age")) +
  ggtitle("Roost popularity - 10 times observation rate") +
  theme_bw(base_size = 18)+
  theme(legend.position = "none",
        panel.grid = element_blank())


### Average strength
columns_estimate = c('analysis', 'model', 'variable','polynomial', 'estimate', 'conf.low', 'conf.high', 'max.vif')
social_estimates_10x = data.frame(matrix(nrow = 8,
                                         ncol = length(columns_estimate))) 
colnames(social_estimates_10x) = columns_estimate
social_estimates_10x$analysis <- "Average strength"

social_mechanisms <- list("+ poly(age_scale, 2, raw = TRUE)",
                          "+ poly(age_scale, 2, raw = TRUE) + (1|ID)", 
                          "+ longevity_10x_scale + (1|ID)", 
                          "+ poly(age_scale, 2, raw = TRUE) + longevity_10x_scale + (1|ID)")

social_mech_10x_models <- list()
social_mech_10x_mod_summary <- list()

j <- 1

for(i in 1:length(models)){
  model_formula <- paste("avg_strength ~ season ", social_mechanisms[[i]], 
                         " + (1|year)", sep = "")
  
  mod <- glmmTMB(formula = as.formula(model_formula), 
                 REML = F,                              # for model comparison
                 family = Gamma(link = "inverse"),
                 data = social_long_10x)
  
  mod_df <- broom.mixed::tidy(mod, conf.int = TRUE)
  
  social_mech_10x_models[[i]] <- mod
  social_mech_10x_mod_summary[[i]] <- summary(mod)
  
  if(i < 3){
    
    social_estimates_10x$variable[j] <- mod_df$term[4]
    social_estimates_10x$estimate[j] <- mod_df$estimate[4]
    social_estimates_10x$conf.low[j] <- mod_df$conf.low[4]
    social_estimates_10x$conf.high[j] <- mod_df$conf.high[4]
    social_estimates_10x$polynomial[j] <- 1
    
    social_estimates_10x$variable[j+1] <- mod_df$term[5]
    social_estimates_10x$estimate[j+1] <- mod_df$estimate[5]
    social_estimates_10x$conf.low[j+1] <- mod_df$conf.low[5]
    social_estimates_10x$conf.high[j+1] <- mod_df$conf.high[5]
    social_estimates_10x$polynomial[j+1] <- 2
    
    social_estimates_10x$model[c(j:c(j+1))] <- models[i]
    social_estimates_10x$max.vif[c(j:c(j+1))] <- max(check_collinearity(mod)$VIF)
    
    j <- j + 2
    
  }
  
  if(i == 3){
    social_estimates_10x$variable[j] <- mod_df$term[4]
    social_estimates_10x$estimate[j] <- mod_df$estimate [4]
    social_estimates_10x$conf.low[j] <- mod_df$conf.low[4]
    social_estimates_10x$conf.high[j] <- mod_df$conf.high[4]
    social_estimates_10x$polynomial[j] <- 1
    
    social_estimates_10x$model[j] <- models[i]
    social_estimates_10x$max.vif[j] <-  max(check_collinearity(mod)$VIF)
    
    j <- j + 1
    
  }
  
  if(i == 4){
    social_estimates_10x$variable[j] <- mod_df$term[4]
    social_estimates_10x$estimate[j] <- mod_df$estimate [4]
    social_estimates_10x$conf.low[j] <- mod_df$conf.low[4]
    social_estimates_10x$conf.high[j] <- mod_df$conf.high[4]
    social_estimates_10x$polynomial[j] <- 1
    
    social_estimates_10x$variable[j+1] <- mod_df$term[5]
    social_estimates_10x$estimate[j+1] <- mod_df$estimate[5]
    social_estimates_10x$conf.low[j+1] <- mod_df$conf.low[5]
    social_estimates_10x$conf.high[j+1] <- mod_df$conf.high[5]
    social_estimates_10x$polynomial[j+1] <- 2
    
    social_estimates_10x$variable[j+2] <- mod_df$term[7]
    social_estimates_10x$estimate[j+2] <- mod_df$estimate[7]
    social_estimates_10x$conf.low[j+2] <- mod_df$conf.low[7]
    social_estimates_10x$conf.high[j+2] <- mod_df$conf.high[7]
    social_estimates_10x$polynomial[j+2] <- 1
    
    social_estimates_10x$model[c(j:c(j+2))] <- models[i]
    social_estimates_10x$max.vif[c(j:c(j+2))] <-  max(check_collinearity(mod)$VIF)
    
    j <- j + 3
    
  }
  rm(mod)
  rm(model_formula)
}

# Extended Data 2B
social_estimates_10x$variable <- ifelse(grepl("poly", social_estimates_10x$variable, fixed = TRUE), 
                                        "Age", "Longevity")
social_estimates_10x$variable_2 <- paste(social_estimates_10x$variable, 
                                         social_estimates_10x$polynomial, sep = "_")
social_estimates_10x$variable_2 <- factor(social_estimates_10x$variable_2,
                                          levels = c("Age_3", "Age_2", "Age_1", "Longevity_1"))

sociality_mechanism_10x_plot <- ggplot(subset(social_estimates_10x, max.vif < 5),
                                       aes(x = variable_2, y = estimate, 
                                           colour = model, fill = model)) +
  geom_errorbar(data = subset(social_estimates_10x, max.vif < 5), 
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
                      breaks = c("Age + Longevity + ID", "Age + ID", "Age")) +
  scale_fill_manual(values = model_colours,
                    breaks = c("Age + Longevity + ID", "Age + ID", "Age")) +
  ggtitle("Average strength - 10 times observation rate") +
  theme_bw(base_size = 18)+
  theme(legend.position = "none",
        panel.grid = element_blank())

