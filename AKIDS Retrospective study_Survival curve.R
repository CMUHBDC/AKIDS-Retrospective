library(riskRegression)
library(adjustedCurves)
library(survival)
library(prodlim)
library(dplyr)
library(patchwork)
library(survminer)

####################################################################################
# Fit Cause-specific Cox proportional hazard regression ====
## CKO
Adjust_model <- CSC(Hist(Lossfollow_CKO_FT_m, CKO_New)~ group + Age + gender + BMI + Smoking + 
                      dm_define + htn_define + Chf + Cvd + Dementia + Renal + NSAID_b90d + 
                      Contrast_b90d + ACEI_ARBs_drug_b90d + Antimicrobial_b90d + Diuretics_b90d +
                      Hb_b1y + min_SCr, data=Data, cause=1)
## Death
Adjust_model2 <- coxph(Surv(Lossfollow_death_FT_m, death_a1y)~ group + Age + gender + BMI + Smoking + 
                         dm_define + htn_define + Chf + Cvd + Dementia + Renal + NSAID_b90d + 
                         Contrast_b90d + ACEI_ARBs_drug_b90d + Antimicrobial_b90d + Diuretics_b90d + 
                         Hb_b1y + min_SCr, data = Data, x = T)

## De_Novo_CKD
Adjust_model3 <- CSC(Hist(Lossfollow_De_Novo_FT_m, De_Novo_New)~ group + Age + gender + BMI + Smoking +
                       dm_define + htn_define + Chf + Cvd + Dementia + Renal + NSAID_b90d + 
                       Contrast_b90d + ACEI_ARBs_drug_b90d + Antimicrobial_b90d + Diuretics_b90d + 
                       Hb_b1y + min_SCr, data=Data2, cause=1)


###############################################################################################
# Calculate Cause-Specific confounder-Adjusted Cumulative Incidence Functions ====
## CKO
Adjust_model_surv <- adjustedcif(data= Data, variable= "group", ev_time= "Lossfollow_CKO_FT_m",
                                 event= "CKO_New", method= "direct", outcome_model= Adjust_model,
                                 cause= 1, conf_int= TRUE)
## Death
Adjust_model2_surv <- adjustedsurv(data= Data, variable= "group", ev_time= "Lossfollow_death_FT_m",
                                   event= "death_a1y", method="direct", outcome_model= Adjust_model2,
                                   conf_int= TRUE)

## De_Novo_CKD
Adjust_model3_surv <- adjustedcif(data= Data2, variable= "group", ev_time= "Lossfollow_De_Novo_FT_m",
                                  event= "De_Novo_New", method= "direct", outcome_model= Adjust_model3,
                                  cause= 1, conf_int= TRUE)

###############################################################################################
# Plot Confounder-Adjusted Curve ====
## CKO
pval = format.pval(survdiff(Surv(Lossfollow_CKO_FT_m, CKO_a1y)~group, Data)$pvalue, digits = 4, eps = 0.01)
p = plot(Adjust_model_surv, 
     cif = T, 
     conf_int = F , 
     color = T, 
     linetype = T, 
     line_size = 1.5,
     xlab = 'Follow-up time (month)', 
     ylab = 'Probability of CKO',
     ylim = c(0.0,0.15),
     custom_colors = c('green4','red'),
     legend.title = paste('p-value', pval),
     legend.position = c(0.2,0.85),
     gg_theme = theme_survminer()) +
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20,vjust = -3),
        axis.text.x  = element_text(size = 18),
        axis.text.y = element_text(size = 18,angle = 90)
  ) +scale_x_continuous(breaks = seq(0,12,3))

pt = p + plot_layout(heights = c(6,1))



## Death
pval = format.pval(survdiff(Surv(Lossfollow_death_FT_m, death_a1y)~group,Data)$pvalue, digits = 4, eps = 0.01)
p2 = 
  plot(Adjust_model2_surv, 
     cif = T, 
     conf_int = F , 
     color = T, 
     linetype = T, 
     line_size = 1.5,
     xlab = 'Follow-up time (month)', 
     ylab = 'Probability of Mortality',
     ylim = c(0, .15),
     custom_colors = c('green4','red'),
     legend.title = paste('p-value', pval),
     legend.position = c(0.2,0.85),
     gg_theme = theme_survminer()) +
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20,vjust = -3),
        axis.text.x  = element_text(size = 18),
        axis.text.y = element_text(size = 18,angle = 90)
  ) + scale_x_continuous(breaks = seq(0,12,3))
pt2 = p2 + plot_layout(heights = c(6,1))


## De_Novo_CKD
pval = format.pval(survdiff(Surv(Lossfollow_De_Novo_FT_m, De_Novo_CKD_2)~group,Data2)$pvalue, digits = 4, eps = 0.01)
p3 = 
  plot(Adjust_model3_surv, 
     cif = T, 
     conf_int = F, 
     color = T, 
     linetype = T, 
     line_size = 1.5,
     xlab = 'Follow-up time (month)', 
     ylab = 'Probability of'~italic("De novo")~'CKD-ND',
     ylim = c(0.0,0.15),
     custom_colors = c('green4','red'),
     legend.title = paste('p-value', pval),
     legend.position = c(0.2,0.85),
     gg_theme = theme_survminer()) +
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20,vjust = -3),
        axis.text.x  = element_text(size = 18),
        axis.text.y = element_text(size = 18,angle = 90)
  ) + scale_x_continuous(breaks = seq(3,12,3), limits = c(3,12)) 
pt3 = p3 + plot_layout(heights = c(6,1))


cowplot::plot_grid(pt, pt2, pt3, nrow = 1, labels = c('(A)', '(B)', '(C)'), label_size = 24) 
ggsave("D:/(link...)/AdjustedCurve_two_groups_final.tiff", width=65, height=25, units='cm', dpi = 300)
