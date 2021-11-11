setwd("C:\\Users\\sirui.zhou\\work\\regeneron-wes")
setwd("C:\\Users\\Sirui\\Desktop\\WORKS\\regeneron")
dt <- data.table::fread("CD109_5cells_mineral_v2.txt")
dt[, Cell := factor(Cell, levels = unique(Cell))]
levels(dt$Cell)

dt_western <- data.table::fread("CD109_5cells_western.txt")

dt_main <- dplyr::left_join(dt,
                            dt_western[, .(med_scaled = median(Level_perc)), by = Cell])

# M2 <- lme4::lmer(Mineralization ~ Cell + (1|Exp),
#                  data=dt)
#
# M2_null <- lme4::lmer(Mineralization ~ 1 + (1|Exp),
#                       data=dt)
#
# # Test the exposure. I.e cell
# anova(M2, M2_null)

M2.cs <- nlme::gls(Mineralization ~ Cell, data = dt,
                   corr = nlme::corCompSymm(form = ~ 1 | Exp), method = "ML" )
M2.cs_null <- nlme::gls(Mineralization ~ 1, data = dt,
                     corr = nlme::corCompSymm(form = ~ 1 | Exp), method = "ML"  )
anova.res <- anova(M2.cs, M2.cs_null)
anova.res
anova.res$`p-value`[2]

# Effect of knock-down/median CD109 expression level on mineralization
M2.cs_level <- nlme::gls(Mineralization ~ med_scaled, data = dt_main,
                   corr = nlme::corCompSymm(form = ~ 1 | Exp) )
coef(summary(M2.cs_level))


###Plot Sirui###
library(tidyverse)
library(rstatix)
library(ggpubr)
require(gridExtra)
group_by(dt_main, Cell) %>% summarize(m = mean(Mineralization))

dt_western$Exp <- as.factor(dt_western$Exp)
dt_western %>%
  group_by(Cell) %>%
  get_summary_stats(Level_perc, type = "mean_sd")

dt_western$Cell <- factor(dt_western$Cell,
                       levels = c('control','70A116','72A144','72A124','72A123','70A146'),ordered = TRUE)

bxp_w <- ggboxplot(dt_western, x = "Cell", y = "Level_perc", add = "point", fill = "lightblue", alpha=0.7, color="black") +
  labs(y= "CD109 Level percentage to control") + ylim(0,1)


res.aov_w <- anova_test(data = dt_western, dv = Level_perc, wid = Exp, within = Cell)
get_anova_table(res.aov_w)


A=bxp_w + 
  labs(
    subtitle = get_test_label(res.aov_w, detailed = F)
  ) + xlab("")


dt_main$Exp <- as.factor(dt_main$Exp)
dt_main$Cell <- factor(dt_main$Cell,
                     levels = c('control','70A116','72A144','72A124','72A123','70A146'), ordered = TRUE)
dt_main %>%
  group_by(Cell) %>%
  get_summary_stats(Mineralization, type = "mean_sd")
bxp_m <- ggboxplot(dt_main, x = "Cell", y = "Mineralization", fill = "#f03b20", color="black", alpha = 0.7, add = "point")


B=bxp_m + 
  labs(
    subtitle = expression(paste("Change in CD109 expression level on mineralization: Beta = -1.71, p = 1.8x10"^{"-7"})),
    y= "Mineralization per ug CD109") + 
  theme(text=element_text(size=13, 
                          family="Sans")) + xlab("")

####
s_m2 <- summary(M2.cs)
coef_dt <- data.table::data.table(param = row.names(coef(s_m2)), coef(s_m2), confint(M2.cs))
# Remove intercept
coef_dt <- coef_dt[grep("intercept", param, ignore.case = T, invert = T)]

# Order by median cd109 level


coef_dt[, Cell := gsub("Cell", "", param)]
coef_dt[, Cell := factor(Cell, levels = c('70A116','72A144','72A124','72A123','70A146'))]




C <- ggplot(coef_dt, aes(x = Cell, y = Value, ymin = `2.5 %`, ymax = `97.5 %`)) +
  geom_pointrange(size = 0.6, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.5)) +
  xlab("Cell") +
  ylab("Effect on CD109 mineralization\ncompared to control (95% CI)") +
  labs(
    subtitle = expression(paste("Anova, p = 6.9x10"^{"-10"}))) + theme_classic() + 
  theme(text=element_text(size=13, 
                          family="Sans"),
        axis.text.x=element_text(size=13, color = "black", hjust=0.5),
        axis.text.y=element_text(size=13, color = "black"))

grid.arrange(A, B, C)

###end###



# Plot
library(ggplot2)
ggplot(coef_dt, aes(x = Cell, y = Value, ymin = `2.5 %`, ymax = `97.5 %`)) +
  geom_pointrange(size = 1.5, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.5)) +
  xlab("Cell line") +
  ylab("Effect on CD109 mineralization\ncompared to control (95% CI)") +
  theme_bw() +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.y=element_text(size = 10),
        axis.text.x=element_text(face = "bold"),
        legend.direction = "horizontal",
        legend.position = c(0.2,0.9),
        legend.title = element_text(face="bold"))

ggsave(file.path(wd, "cell_x_mineralization.pdf"), dpi = "retina", width = 10, height = 5)



plot_dt <- data.table::data.table(dt_main)
order_var <- dt_western[, .(med_scaled = median(Level_perc)), by = Cell][order(med_scaled)][, .(Cell)]
plot_dt[, Cell_ord := factor(Cell, levels = order_var$Cell)]
med_control <- plot_dt[Cell == "control", median(Mineralization)]
ggplot(plot_dt[Cell != "control"], aes(x = Cell_ord, y = Mineralization)) +
  geom_boxplot() +
  geom_hline(yintercept = med_control, linetype = "dashed", col = "blue") +
  xlab("Cell line") +
  ylab("CD109 Mineralization") +
  theme_bw() +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.y=element_text(size = 10),
        axis.text.x=element_text(face = "bold"),
        legend.direction = "horizontal",
        legend.position = c(0.2,0.9),
        legend.title = element_text(face="bold"))

ggsave(file.path(wd, "level_x_mineralization.pdf"), dpi = "retina", width = 10, height = 5)
