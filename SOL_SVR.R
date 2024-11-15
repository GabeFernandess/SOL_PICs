#Delta F - SOL

# Libraries
library(readxl)
library(tidyverse)
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(sjstats)
library(car)
library(pwr)
library(viridis)
library(tidyr)
library(patchwork)
library(ggbeeswarm)


sol <-read_excel("C:/Users/gabep/OneDrive/Research projects/Study 2b - PIC/SOL_AT.xlsx")

colnames(sol)

#DATA MANIPULATION ----

names(sol)[names(sol) == "Brace Height_norm"] <- "bh"
names(sol)[names(sol) == "Attenuation"] <- "att"
names(sol)[names(sol) == "Peak_DR-Rec_DR"] <- "peak_dr_rec"  
names(sol)[names(sol) == "Peak_DRtime-Rec_time"] <- "peak_dr_time_rec"  
names(sol)[names(sol) == "Peak_DRtime-Derec_time"] <- "peak_dr_time_derec"
names(sol)[names(sol) == "deltaf"] <- "df"
names(sol)[names(sol) == "deltaf_k"] <- "dfk"
names(sol)[names(sol) == "peak_dr_SVR"] <- "peakdr"
names(sol)[names(sol) == "rec"] <- "rec"
names(sol)[names(sol) == "derec"] <- "derec"
names(sol)[names(sol) == "DR_rec"] <- "dr_rec"
names(sol)[names(sol) == "DR_drec"] <- "dr_derec"


# Create a data frame for the "sol" muscle

sol$group <- factor(sol$group, levels = c("AT", "Control"))

levels(as.factor(sol$group))

(refgrid <- list (group=c("AT","Control")))

sol <- sol %>%
  mutate(dfk_ratio = df / dfk)

sol <- sol %>%
  mutate(dr_ratio = peak_dr_time_rec / peak_dr_time_derec)

sol_clean <- sol[!is.na(sol$dr_ratio) & !is.nan(sol$dr_ratio) & !is.infinite(sol$dr_ratio), ]

#calculate mean
sol_mean <- sol %>%
  filter(!is.na(att)) %>%  
  group_by(participant, contraction, group) %>%
  summarise(
    df = mean(df, na.rm = TRUE),
    dfk = mean(dfk, na.rm = TRUE),
    dfk_ratio= mean(dfk_ratio, na.rm=TRUE),
    bh = mean(bh, na.rm = TRUE),
    att=mean(att, na.rm = TRUE),
    peakdr = mean(peakdr, na.rm = TRUE),
    rec = mean(rec, na.rm = TRUE),
    derec = mean(derec, na.rm = TRUE),
    dr_rec = mean(dr_rec, na.rm = TRUE),
    dr_derec = mean(dr_derec, na.rm = TRUE),
    peak_dr_time_rec = mean(peak_dr_time_rec, na.rm = TRUE),
    peak_dr_time_derec = mean(peak_dr_time_derec, na.rm = TRUE),
        dr_ratio = mean(dr_ratio, na.rm = TRUE)
      )
sol_mean

sol_mean2 <- sol_clean %>%
  group_by(participant, contraction, group) %>%
  summarise(dr_ratio = mean(dr_ratio, na.rm = TRUE))
    

#START ANALYSIS -----

# DELTA FK/ratio SOL----
    
fit_deltaf1 <- lmer(dfk_ratio ~ as.factor(group) + contraction+(1 | participant/mu_id), sol)
fit_deltaf2 <- lmer(dfk_ratio ~ as.factor(group) + age +(1| participant/mu_id), sol)

#Check for best model
anova(fit_deltaf1,fit_deltaf2,fit_deltaf3)

#Best model---
fit_deltafk <- fit_deltaf2 # no influence from age or rec - not used in model

summary(fit_deltafk)

anova(fit_deltafk) %>%
  knitr::kable()


fitdeltaf.emm.s <- emmeans(fit_deltafk, "group")
pairs(fitdeltaf.emm.s, adjust = "bonferroni")


#confint
anovafit_fit <- anova(fit_deltafk)
effectsize::omega_squared(anovafit_fit)

emm <- emmeans(fit_deltafk, pairwise ~ group)
confint(emm)


#confint
anovafit_fit <- anova(fit_deltafk)
effectsize::omega_squared(anovafit_fit)

#Refgrid
mar_df2 <- emmip(fit_deltafk, ~ as.factor(group), at = refgrid, CIs = T, plotit = F)
mar_df2


#FINAL PLOT delta f K  ----- 
ggplot(data = sol, aes(x = group, y = dfk_ratio)) +
  geom_quasirandom(data = subset(sol, group == "AT"), 
                   width = 0.12, size = 4.5, alpha = 0.3, 
                   shape = 23, colour ="#abdbe3" ,fill = "#76b5c5") +
  geom_quasirandom(data = subset(sol, group == "Control"), 
                   width = 0.12, size = 4.5, alpha = 0.3, 
                   shape = 23, colour = "#efa2b5", fill = "#e67391") +
  theme_light(base_size = 14) +
  guides(fill = "none", color = "none") +
  geom_jitter(data = subset(sol_mean, group == "AT"), 
              aes(x = group, y = dfk_ratio), 
              width = 0.02, size = 4.5, alpha = 1, 
              shape = 23, colour = "#154c79", fill="#127aab") +
  
  geom_jitter(data = subset(sol_mean, group == "Control"), 
              aes(x = group, y = dfk_ratio), 
              width = 0.02, size = 4.5, alpha = 1, 
              shape = 23, colour = "#730e27", fill = "#e31b62") +
  geom_errorbar(data = mar_df2, aes(ymin = LCL, ymax = UCL, y = yvar),
                position = position_nudge(x = -0.3), width = 0, size = 1, colour ="#242f36") +
  geom_point(data = mar_df2, aes(x = group, y = yvar), 
             size = 4.5,
             alpha = 1,
             position = position_nudge(x = -0.3), 
             shape = 16,
             colour = "#242f36") + 
    theme(
    axis.text = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  labs(y =  expression(Delta~"F/k (pps)"), x = "")-> plot_deltafk 

plot_deltafk <- plot_deltafk+ scale_x_discrete(labels = c("AT ", "Control"))

plot_deltafk

ggsave(file = "deltafK_ratio.tiff", units="in", width = 6.5, height = 7.5, dpi = 300, compression = "lzw")


#Deltaf-----

fit_deltaf1 <- lmer(df ~ as.factor(group) + rec +age +contraction+(1 | participant/mu_id), sol)
fit_deltaf2 <- lmer(df ~ as.factor(group) + age +(1| participant/mu_id), sol)
fit_deltaf3 <- lmer(df ~ as.factor(group) + (1 | participant/mu_id), sol)

#Check for best model
anova(fit_deltaf1,fit_deltaf2,fit_deltaf3)

#Best model---
fit_deltaf <- fit_deltaf2 # no diff in age or rec or contraction - not used in model

summary(fit_deltaf)

anova(fit_deltaf) %>%
  knitr::kable()


fitdeltaf.emm.s <- emmeans(fit_deltaf, "group")
pairs(fitdeltaf.emm.s, adjust = "bonferroni")


#confint
anovafit_fit <- anova(fit_deltaf)
effectsize::omega_squared(anovafit_fit)

emm <- emmeans(fit_deltaf, pairwise ~ group)
confint(emm)


#Refgrid
mar_df <- emmip(fit_deltaf, ~ as.factor(group), at = refgrid, CIs = T, plotit = F)
mar_df


#FINAL PLOT delta f  ----- 
ggplot(data = sol, aes(x = group, y = df)) +
  geom_quasirandom(data = subset(sol, group == "AT"), 
                   width = 0.12, size = 4.5, alpha = 0.25, 
                   shape = 23, colour ="#abdbe3" ,fill = "#76b5c5") +
  geom_quasirandom(data = subset(sol, group == "Control"), 
                   width = 0.12, size = 4.5, alpha = 0.25, 
                   shape = 23, colour = "#efa2b5", fill = "#e67391") +
  theme_light(base_size = 14) +
  guides(fill = "none", color = "none") +
  geom_jitter(data = subset(sol_mean, group == "AT"), 
              aes(x = group, y = df), 
              width = 0.02, size = 4.5, alpha = 1, 
              shape = 23, colour = "#154c79", fill="#127aab") +
  
  geom_jitter(data = subset(sol_mean, group == "Control"), 
              aes(x = group, y = df), 
              width = 0.02, size = 4.5, alpha = 1, 
              shape = 23, colour = "#730e27", fill = "#e31b62") +
  geom_errorbar(data = mar_df, aes(ymin = LCL, ymax = UCL, y = yvar),
                position = position_nudge(x = -0.3), width = 0, size = 1, colour ="#242f36") +
  geom_point(data = mar_df, aes(x = group, y = yvar), 
             size = 4.5,
             alpha = 1,
             position = position_nudge(x = -0.3), 
             shape = 16,
             colour = "#242f36") + 
  theme(
    axis.text = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  labs(y =  expression(Delta~"F (pps)"), x = "")-> plot_deltaf

plot_deltaf <- plot_deltaf+ scale_x_discrete(labels = c("AT ", "Control"))

plot_deltaf

ggsave(file = "deltaf.tiff", units="in", width = 6.5, height = 7.5, dpi = 300, compression = "lzw")



#Brace Height ----

# BH ---
  
fit_bh<- lmer(bh~ as.factor(group) + peak_dr_rec+(1 | participant/mu_id), sol)
summary(fit_bh)

anova(fit_bh) %>%
  knitr::kable()


fitbh.emm.s <- emmeans(fit_bh, "group", "peak_dr_rec")
pairs(fitbh.emm.s, adjust = "bonferroni")


#confint
anovafit_fit <- anova(fit_bh)
effectsize::omega_squared(anovafit_fit)

emm <- emmeans(fit_bh, pairwise ~ group)
confint(emm)


#Refgrid
mar_bh <- emmip(fit_bh, ~ as.factor(group), at = refgrid, CIs = T, plotit = F)
mar_bh


#FINAL Plot  BH_ v1-2-----
ggplot(data = sol, aes(x = group, y = bh)) +
  geom_quasirandom(data = subset(sol, group == "AT"), 
                   width = 0.12, size = 4.5, alpha = 0.3, 
                   shape = 23, colour ="#abdbe3" ,fill = "#76b5c5") +
  geom_quasirandom(data = subset(sol, group == "Control"), 
                   width = 0.12, size = 4.5, alpha = 0.3, 
                   shape = 23, colour = "#efa2b5", fill = "#e67391") +
  theme_light(base_size = 14) +
  guides(fill = "none", color = "none") +
  geom_jitter(data = subset(sol_mean, group == "AT"), 
              aes(x = group, y = bh), 
              width = 0.02, size = 4.5, alpha = 1, 
              shape = 23, colour = "#154c79", fill="#127aab") +
  
  geom_jitter(data = subset(sol_mean, group == "Control"), 
              aes(x = group, y = bh), 
              width = 0.02, size = 4.5, alpha = 1, 
              shape = 23, colour = "#730e27", fill = "#e31b62") +
  geom_errorbar(data = mar_bh, aes(ymin = LCL, ymax = UCL, y = yvar),
                position = position_nudge(x = -0.3), width = 0, size = 1, colour ="#242f36") +
  geom_point(data = mar_bh, aes(x = group, y = yvar), 
             size = 4.5,
             alpha = 1,
             position = position_nudge(x = -0.3), 
             shape = 16,
             colour = "#242f36") + 
  theme(
    axis.text = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  labs(y =  expression("Brace height (% rTri)"), x = "") -> plot_bh_ta

plot_bh_ta <- plot_bh_ta+ scale_x_discrete(labels = c("AT", "Control"))

plot_bh_ta

ggsave(file = "BH_TA_v1_2.tiff", units="in", width = 6.5, height = 7.5, dpi = 300, compression = "lzw")


#ATT ----

# ATTenuation ---


fit_att1<- lmer(att~ as.factor(group) + peak_dr_rec +(1 | participant/mu_id), sol)
fit_att2<- lmer(att~ as.factor(group) +(1 | participant/mu_id), sol)


#Check for best model
anova(fit_att1,fit_att2)

#Best model---
fit_att <- fit_att1

summary(fit_att)


anova(fit_att) %>%
  knitr::kable()



fitatt.emm.s <- emmeans(fit_att, "group")
pairs(fitatt.emm.s, adjust = "bonferroni")


#confint
anovafit_fit <- anova(fit_att)
effectsize::omega_squared(anovafit_fit)

emm <- emmeans(fit_att, pairwise ~ group)
confint(emm)

#Refgrid
mar_att <- emmip(fit_att, ~ as.factor(group), at = refgrid, CIs = T, plotit = F)
mar_att


#FINAL Plot  ATT v1-2-----
ggplot(data = sol, aes(x = group, y = att)) +
  geom_quasirandom(data = subset(sol, group == "AT"), 
                   width = 0.12, size = 4.5, alpha = 0.3, 
                   shape = 23, colour ="#abdbe3" ,fill = "#76b5c5") +
  geom_quasirandom(data = subset(sol, group == "Control"), 
                   width = 0.12, size = 4.5, alpha = 0.3, 
                   shape = 23, colour = "#efa2b5", fill = "#e67391") +
  theme_light(base_size = 14) +
  theme_light(base_size = 14) +
  guides(fill = "none", color = "none") +
  geom_jitter(data = subset(sol_mean, group == "AT"), 
              aes(x = group, y = att), 
              width = 0.02, size = 4.5, alpha = 1, 
              shape = 23, colour = "#154c79", fill="#127aab") +
  geom_jitter(data = subset(sol_mean, group == "Control"), 
              aes(x = group, y = att), 
              width = 0.02, size = 4.5, alpha = 1, 
              shape = 23, colour = "#730e27", fill = "#e31b62") +
  geom_errorbar(data = mar_att, aes(ymin = LCL, ymax = UCL, y = yvar),
                position = position_nudge(x = -0.3), width = 0, size = 1, colour ="#242f36") +
  geom_point(data = mar_att, aes(x = group, y = yvar), 
             size = 4.5,
             alpha = 1,
             position = position_nudge(x = -0.3), 
             shape = 16,
             colour = "#242f36") + 
  theme(
    axis.text = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  labs(y =  expression("Attenuation (pps/% peak torque)"), x = "") -> plot_att

plot_att <- plot_att+ scale_x_discrete(labels = c("AT", "Control"))

plot_att

ggsave(file = "ATT.tiff", units="in", width = 6.5, height = 7.5, dpi = 300, compression = "lzw")


# peak FR ----

fit_dr<- lmer(peakdr ~ as.factor(group) + contraction+(1 | participant/mu_id), sol) 


summary(fit_dr)


anova(fit_dr) %>%
  knitr::kable()


fitdr.emm.s3 <- emmeans(fit_dr, "group")
pairs(fitdr.emm.s3, adjust = "bonferroni")



#confint
anovafit_fit <- anova(fit_dr)
effectsize::omega_squared(anovafit_fit)

emm <- emmeans(fit_dr, pairwise ~ group)
confint(emm)

#Refgrid
mar_dr <- emmip(fit_dr, ~ as.factor(group), at = refgrid, CIs = T, plotit = F)
mar_dr


#Mean difference deltaf mean ------
#norm_RMS
emm4 <- emmeans(fit_dr, pairwise ~ group)
confint(emm4)


#FINAL Peak FR---plot dr_rec-----
ggplot(data = sol, aes(x = group, y = peakdr)) +
  geom_quasirandom(data = subset(sol, group == "AT"), 
                   width = 0.12, size = 4.5, alpha = 0.3, 
                   shape = 23, colour ="#abdbe3" ,fill = "#76b5c5") +
  geom_quasirandom(data = subset(sol, group == "Control"), 
                   width = 0.12, size = 4.5, alpha = 0.3, 
                   shape = 23, colour = "#efa2b5", fill = "#e67391") +
  theme_light(base_size = 14) +
  guides(fill = "none", color = "none") +
  geom_jitter(data = subset(sol_mean, group == "AT"), 
              aes(x = group, y = peakdr), 
              width = 0.02, size = 4.5, alpha = 1, 
              shape = 23, colour = "#154c79", fill="#127aab") +
  geom_jitter(data = subset(sol_mean, group == "Control"), 
              aes(x = group, y = peakdr), 
              width = 0.02, size = 4.5, alpha = 1, 
              shape = 23, colour = "#730e27", fill = "#e31b62") +
  geom_point(data = mar_dr, aes(x = group, y = yvar), 
             size = 4.5,
             alpha = 1,
             position = position_nudge(x = -0.3), 
             shape = 16,
             colour = "#242f36") + 
  geom_errorbar(data = mar_dr, aes(ymin = LCL, ymax = UCL, y = yvar),
                position = position_nudge(x = -0.3), width = 0, size = 1, colour ="#242f36") +
  theme(
    axis.text = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  labs(y =  expression("Peak FR (pps)"), x = "")-> plot_dr 

plot_dr <-plot_dr+ scale_x_discrete(labels = c("AT", "Control"))

plot_dr

ggsave(file = "fr_REC_sol.tiff", units="in", width = 13, height = 10, dpi = 300, compression = "lzw")



# FR at recruitment-----

fit_frrec<- lmer(dr_rec~ as.factor(group) + (1 | participant/mu_id), sol)
summary(fit_frrec)


anova(fit_frrec) %>%
  knitr::kable()


fitrec.emm.s <- emmeans(fit_frrec, "group")
pairs(fitrec.emm.s, adjust = "bonferroni")


#confint
anovafit_fit <- anova(fit_frrec)
effectsize::omega_squared(anovafit_fit)

emm <- emmeans(fit_frrec, pairwise ~ group)
confint(emm)


#Refgrid
mar_frrec <- emmip(fit_frrec, ~ as.factor(group), at = refgrid, CIs = T, plotit = F)
mar_frrec



#FINAL FR_REC---plot dr_rec-----
ggplot(data = sol, aes(x = group, y = dr_rec)) +
  geom_quasirandom(data = subset(sol, group == "AT"), 
                   width = 0.12, size = 4.5, alpha = 0.3, 
                   shape = 23, colour ="#abdbe3" ,fill = "#76b5c5") +
  geom_quasirandom(data = subset(sol, group == "Control"), 
                   width = 0.12, size = 4.5, alpha = 0.3, 
                   shape = 23, colour = "#efa2b5", fill = "#e67391") +
  theme_light(base_size = 14) +
  guides(fill = "none", color = "none") +
  geom_jitter(data = subset(sol_mean, group == "AT"), 
              aes(x = group, y = dr_rec), 
              width = 0.02, size = 4.5, alpha = 1, 
              shape = 23, colour = "#154c79", fill="#127aab") +
  geom_jitter(data = subset(sol_mean, group == "Control"), 
              aes(x = group, y = dr_rec), 
              width = 0.02, size = 4.5, alpha = 1, 
              shape = 23, colour = "#730e27", fill = "#e31b62") +
  geom_point(data = mar_frrec, aes(x = group, y = yvar), 
             size = 4.5,
             alpha = 1,
             position = position_nudge(x = -0.3), 
             shape = 16,
             colour = "#242f36") + 
  geom_errorbar(data = mar_frrec, aes(ymin = LCL, ymax = UCL, y = yvar),
                position = position_nudge(x = -0.3), width = 0, size = 1, colour ="#242f36") +
  theme(
    axis.text = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  labs(y =  expression("FR at recruitment (pps)"), x = "")-> plot_dr_rec 

plot_dr_rec <-plot_dr_rec+ scale_x_discrete(labels = c("AT", "Control"))

plot_dr_rec

ggsave(file = "fr_REC_sol.tiff", units="in", width = 13, height = 10, dpi = 300, compression = "lzw")



# FR at DErecruitment-----

fit_frderec<- lmer(dr_derec~ as.factor(group) +(1 | participant/mu_id), sol)
summary(fit_frderec)


anova(fit_frderec) %>%
  knitr::kable()



fitderec.emm.s <- emmeans(fit_frderec, "group")
pairs(fitrec.emm.s, adjust = "bonferroni")


#confint
anovafit_fit <- anova(fit_frderec)
effectsize::omega_squared(anovafit_fit)

emm <- emmeans(fit_frderec, pairwise ~ group)
confint(emm)


#Refgrid
mar_frderec <- emmip(fit_frderec, ~ as.factor(group), at = refgrid, CIs = T, plotit = F)
mar_frderec



#FINAL FR_DEREC---plot dr_rec-----
ggplot(data = sol, aes(x = group, y = dr_derec)) +
  geom_quasirandom(data = subset(sol, group == "AT"), 
                   width = 0.12, size = 4.5, alpha = 0.3, 
                   shape = 23, colour ="#abdbe3" ,fill = "#76b5c5") +
  geom_quasirandom(data = subset(sol, group == "Control"), 
                   width = 0.12, size = 4.5, alpha = 0.3, 
                   shape = 23, colour = "#efa2b5", fill = "#e67391") +
  theme_light(base_size = 14) +
  guides(fill = "none", color = "none") +
  geom_jitter(data = subset(sol_mean, group == "AT"), 
              aes(x = group, y = dr_derec), 
              width = 0.02, size = 4.5, alpha = 1, 
              shape = 23, colour = "#154c79", fill="#127aab") +
  geom_jitter(data = subset(sol_mean, group == "Control"), 
              aes(x = group, y = dr_rec), 
              width = 0.02, size = 4.5, alpha = 1, 
              shape = 23, colour = "#730e27", fill = "#e31b62") +
  geom_point(data = mar_frderec, aes(x = group, y = yvar), 
             size = 4.5,
             alpha = 1,
             position = position_nudge(x = -0.3), 
             shape = 16,
             colour = "#242f36") + 
  geom_errorbar(data = mar_frderec, aes(ymin = LCL, ymax = UCL, y = yvar),
                position = position_nudge(x = -0.3), width = 0, size = 1, colour ="#242f36") +
  theme(
    axis.text = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  labs(y =  expression("FR at derecruitment (pps)"), x = "")-> plot_dr_derec 

plot_dr_derec <-plot_dr_derec+ scale_x_discrete(labels = c("AT", "Control"))

plot_dr_derec

ggsave(file = "fr_DEREC_sol.tiff", units="in", width = 13, height = 10, dpi = 300, compression = "lzw")




# FR time recruitment to peak---

fit_timefrrec<- lmer(peak_dr_time_rec~ as.factor(group) +(1 | participant/mu_id), sol)
summary(fit_timefrrec)


anova(fit_timefrrec) %>%
  knitr::kable()



fitderec.emm.s <- emmeans(fit_timefrrec, "group")
pairs(fitrec.emm.s, adjust = "bonferroni")


#confint
anovafit_fit <- anova(fit_timefrrec)
effectsize::omega_squared(anovafit_fit)


emm <- emmeans(fit_timefrrec, pairwise ~ group)
confint(emm)

#Refgrid
mar_timefrrecc <- emmip(fit_timefrrec, ~ as.factor(group), at = refgrid, CIs = T, plotit = F)
mar_timefrrecc



#FINAL FR_DEREC---
ggplot(data = sol, aes(x = group, y = peak_dr_time_rec)) +
  geom_quasirandom(data = subset(sol, group == "AT"), 
                   width = 0.12, size = 4.5, alpha = 0.3, 
                   shape = 23, colour ="#abdbe3" ,fill = "#76b5c5") +
  geom_quasirandom(data = subset(sol, group == "Control"), 
                   width = 0.12, size = 4.5, alpha = 0.3, 
                   shape = 23, colour = "#efa2b5", fill = "#e67391") +
  theme_light(base_size = 14) +
  guides(fill = "none", color = "none") +
  geom_jitter(data = subset(sol_mean, group == "AT"), 
              aes(x = group, y = peak_dr_time_rec), 
              width = 0.02, size = 4.5, alpha = 1, 
              shape = 23, colour = "#154c79", fill="#127aab") +
  geom_jitter(data = subset(sol_mean, group == "Control"), 
              aes(x = group, y = peak_dr_time_rec), 
              width = 0.02, size = 4.5, alpha = 1, 
              shape = 23, colour = "#730e27", fill = "#e31b62") +
  geom_point(data = mar_timefrrecc, aes(x = group, y = yvar), 
             size = 4.5,
             alpha = 1,
             position = position_nudge(x = -0.3), 
             shape = 16,
             colour = "#242f36") + 
  geom_errorbar(data = mar_timefrrecc, aes(ymin = LCL, ymax = UCL, y = yvar),
                position = position_nudge(x = -0.3), width = 0, size = 1, colour ="#242f36") +
  theme(
    axis.text = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  labs(y =  expression("Ascending time to peak FR (s)"), x = "")-> plot_peak_dr_time_rec 

plot_peak_dr_time_rec <-plot_peak_dr_time_rec+ scale_x_discrete(labels = c("AT", "Control"))

plot_peak_dr_time_rec

ggsave(file = "time_to_peak_dr.tiff", units="in", width = 13, height = 10, dpi = 300, compression = "lzw")



# FR time from peak to derec---

fit_timefderrec<- lmer(peak_dr_time_derec~ as.factor(group) +(1 | participant/mu_id), sol)
summary(fit_timefderrec)


anova(fit_timefderrec) %>%
  knitr::kable()


fitderec.emm.s <- emmeans(fit_timefderrec, "group")
pairs(fitrec.emm.s, adjust = "bonferroni")


#confint
anovafit_fit <- anova(fit_timefderrec)
effectsize::omega_squared(anovafit_fit)

emm <- emmeans(fit_timefderrec, pairwise ~ group)
confint(emm)

#Refgrid
mar_timefrderecc <- emmip(fit_timefderrec, ~ as.factor(group), at = refgrid, CIs = T, plotit = F)
mar_timefrderecc



#FINAL Time_FR_deREC---plot dr_rec-----
ggplot(data = sol, aes(x = group, y = peak_dr_time_derec)) +
  geom_quasirandom(data = subset(sol, group == "AT"), 
                   width = 0.12, size = 4.5, alpha = 0.3, 
                   shape = 23, colour ="#abdbe3" ,fill = "#76b5c5") +
  geom_quasirandom(data = subset(sol, group == "Control"), 
                   width = 0.12, size = 4.5, alpha = 0.3, 
                   shape = 23, colour = "#efa2b5", fill = "#e67391") +
  theme_light(base_size = 14) +
  guides(fill = "none", color = "none") +
  geom_jitter(data = subset(sol_mean, group == "AT"), 
              aes(x = group, y = peak_dr_time_derec), 
              width = 0.02, size = 4.5, alpha = 1, 
              shape = 23, colour = "#154c79", fill="#127aab") +
  geom_jitter(data = subset(sol_mean, group == "Control"), 
              aes(x = group, y = peak_dr_time_derec), 
              width = 0.02, size = 4.5, alpha = 1, 
              shape = 23, colour = "#730e27", fill = "#e31b62") +
  geom_point(data = mar_timefrderecc, aes(x = group, y = yvar), 
             size = 4.5,
             alpha = 1,
             position = position_nudge(x = -0.3), 
             shape = 16,
             colour = "#242f36") + 
  geom_errorbar(data = mar_timefrderecc, aes(ymin = LCL, ymax = UCL, y = yvar),
                position = position_nudge(x = -0.3), width = 0, size = 1, colour ="#242f36") +
  theme(
    axis.text = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  labs(y =  expression("Descending time to peak FR(s)"), x = "")-> plot_peak_dr_time_derec 

plot_peak_dr_time_derec <-plot_peak_dr_time_derec+ scale_x_discrete(labels = c("AT", "Control"))

plot_peak_dr_time_derec

ggsave(file = "time_peak_dr_derect.tiff", units="in", width = 13, height = 10, dpi = 300, compression = "lzw")



 
# FR time from rec to derec---

fit_timedr_ratio <- lmer(dr_ratio~ as.factor(group) +(1 | participant/mu_id), sol_clean)

summary(fit_timedr_ratio)


anova(fit_timedr_ratio) %>%
  knitr::kable()



fitderectime.emm.s <- emmeans(fit_timedr_ratio, "group")
pairs(fitderectime.emm.s, adjust = "bonferroni")


#confint
anovafit_fit <- anova(fit_timedr_ratio)
effectsize::omega_squared(anovafit_fit)

emm <- emmeans(fit_timedr_ratio, pairwise ~ group)
confint(emm)

#Refgrid
mar_timedr_ratio <- emmip(fit_timedr_ratio, ~ as.factor(group), at = refgrid, CIs = T, plotit = F)
mar_timedr_ratio



#FINAL Time_FR_REC_to_Der---plot dr_rec-----
ggplot(data = sol_clean, aes(x = group, y = dr_ratio)) +
  geom_quasirandom(data = subset(sol_clean, group == "AT"), 
                   width = 0.12, size = 4.5, alpha = 0.3, 
                   shape = 23, colour ="#abdbe3" ,fill = "#76b5c5") +
  geom_quasirandom(data = subset(sol_clean, group == "Control"), 
                   width = 0.12, size = 4.5, alpha = 0.3, 
                   shape = 23, colour = "#efa2b5", fill = "#e67391") +
  theme_light(base_size = 14) +
  guides(fill = "none", color = "none") +
  geom_jitter(data = subset(sol_mean2, group == "AT"), 
              aes(x = group, y = dr_ratio), 
              width = 0.02, size = 4.5, alpha = 1, 
              shape = 23, colour = "#154c79", fill="#127aab") +
  geom_jitter(data = subset(sol_mean2, group == "Control"), 
              aes(x = group, y = dr_ratio), 
              width = 0.02, size = 4.5, alpha = 1, 
              shape = 23, colour = "#730e27", fill = "#e31b62") +
  geom_point(data = mar_timedr_ratio, aes(x = group, y = yvar), 
             size = 4.5,
             alpha = 1,
             position = position_nudge(x = -0.3), 
             shape = 16,
             colour = "#242f36") + 
  geom_errorbar(data = mar_timedr_ratio, aes(ymin = LCL, ymax = UCL, y = yvar),
                position = position_nudge(x = -0.3), width = 0, size = 1, colour ="#242f36") +
  theme(
    axis.text = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  labs(y =  expression(" Ascending / descending time (s)"), x = "")-> plot_dr_ratio 

plot_dr_ratio <-plot_dr_ratio+ scale_x_discrete(labels = c("AT", "Control"))

plot_dr_ratio

ggsave(file = "time_peak_dr_time.tiff", units="in", width = 8, height = 9
       , dpi = 300, compression = "lzw")



#combine Plots ------


p1 <- plot_deltafk
p2 <- plot_deltaf
p3 <- plot_bh_ta
p4 <- plot_att
  
  

p1b <- plot_dr 
p2b <- plot_dr_rec
p3b <- plot_dr_derec
p4b <- plot_peak_dr_time_rec
p5 <- plot_peak_dr_time_derec
p6 <- plot_dr_ratio



v <- p2 +p1 +p3 + p4

v
#ggsave(file = "combined_plot1c.tiff", units="in", width = 8.5, height = 8.5, dpi = 300, compression = "lzw")


v2 <- p1b +p2b +p3b + p4b + p5 +p6


v2

#ggsave(file = "combined_plot22.tiff", units="in", width = 8.5, height = 8.5, dpi = 300, compression = "lzw")


v +
  plot_annotation(
    tag_levels = 'a') & 
  theme(
    plot.tag = element_text(size = 15, face = "bold"),
    plot.title = element_text(size = 15, face = "bold") # Customize the title appearance
  )
ggsave(file = "combined_plot1c.tiff", units="in", width = 8.5, height = 9, dpi = 300, compression = "lzw")



v2 + 
  plot_annotation(
    tag_levels = 'a') & 
  theme(
    plot.tag = element_text(size = 15, face = "bold"),
    plot.title = element_text(size = 15, face = "bold") # Customize the title appearance
  )


ggsave(file = "combined_plott2c.tiff", units="in", width = 9.5, height = 9, dpi = 300, compression = "lzw")

