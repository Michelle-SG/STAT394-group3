##################
## NEW IMPROVED ##
##################
library(tidyverse)
library(corrplot)
library(GGally)
library(knitr)
library(pander)
library(ggrepel)
library(scales)

df <- readr::read_csv("orichalcum_ingots_dataset.csv", skip = 1)
colnames(df) <- c("Sample","Cu","Zn","Pb","Fe","Ni","Ag","Sb","As","Co",
                  "Cd","Sn","Te","Bi","Mn","Li","Al","V","Cr","Rb","Sr","Ba","Cu_Zn","w_w")
df <- df %>% dplyr::select(-w_w, -Cu_Zn)


##########################
## Concentration in PPM ##
##########################

ppm <- df %>%
  dplyr::select(Ni, Ag, Sb, As, Co, Cd, Sn, Te, Bi, Mn,
                Li, Al, V, Cr, Rb, Sr, Ba) %>%
  mutate(across(everything(), as.numeric))
ppm
element_order <- c("Ni","Ag","Sb","As","Co","Cd","Sn","Te","Bi","Mn",
                   "Li","Al","V","Cr","Rb","Sr","Ba")



# Summary stats histograms plot
ppm %>% 
  pivot_longer(cols = everything(), names_to = "Element", values_to = "Concentration") %>%
  mutate(
    Element = factor(Element, levels = element_order),
    Concentration = as.numeric(Concentration)   # make sure it’s numeric
  ) %>%
  ggplot(aes(x = Concentration)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  facet_wrap(~Element, scales = "free") +
  theme_minimal()



#correlation plot
corrplot(cor(ppm, use = "pairwise.complete.obs"), method = "color", type   = "upper")

# Cr vs V by Ni
ggplot(df, aes(x = Cr, y = V, color = Ni)) +
  geom_point(size = 3) +
  labs(title = "Cr vs V coloured by Ni")

# Ni vs Co by Ag
ggplot(df, aes(x=Ni, y=Co, color=Ag)) +
  geom_point(size=3) +
  labs(title = "Ni vs Co coloured by Ag")

# As vs Sb by Cd
ggplot(df, aes(x=As, y=Sb, color=Cd)) +
  geom_point(size=3) +
  labs(title = "As vs Sb coloured by Cd")

# Bi vs Ag by Co
ggplot(df, aes(x=Bi, y=Ag, color=Co)) +
  geom_point(size=3) +
  labs(title = "Bi vs Ag coloured by Co")



# boxplots for each element
ppm %>%
  pivot_longer(cols = everything(),
               names_to = "Element",
               values_to = "Concentration") %>%
  mutate(Concentration = as.numeric(Concentration)) %>%
  ggplot(aes(x = Element, y = Concentration)) +
  geom_boxplot() +
  coord_flip() +
  theme_minimal() +
  labs(title = "Boxplots for each element")


# log-scale boxplots for each element
ppm %>%
  pivot_longer(cols = everything(),
               names_to = "Element",
               values_to = "Concentration") %>%
  mutate(Concentration = as.numeric(Concentration)) %>%
  ggplot(aes(x = Element, y = Concentration)) +
  geom_boxplot() +
  coord_flip() +
  scale_y_continuous(trans = "log") +
  theme_minimal() +
  labs(title = "Natural log-scale boxplots for each element")



# Q-Q plots for each element
ppm %>%
  pivot_longer(cols = everything(), names_to = "Element", values_to = "Concentration") %>%
  mutate(Concentration = as.numeric(Concentration)) %>%
  ggplot(aes(sample = Concentration)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  facet_wrap(~Element, scales = "free") +
  labs(title = "Q–Q plots for each element")



# Mahalanobis distance
ppm_cc <- ppm %>%
  mutate(across(everything(), as.numeric)) %>%
  tidyr::drop_na()

mu.hat <- colMeans(ppm_cc)
Sigma.hat <- cov(ppm_cc)
dM <- mahalanobis(ppm_cc, center = mu.hat, cov = Sigma.hat)

df_chi <- ncol(ppm_cc)                               # e.g., 17 for Ni..Ba
upper.quantiles <- qchisq(c(.90, .95, .99), df = df_chi)
density.at.quantiles <- dchisq(x = upper.quantiles, df = df_chi)
cut.points <- data.frame(upper.quantiles, density.at.quantiles)

ggplot(data.frame(dM), aes(x = dM)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = nclass.FD(dM), fill = "white", col = "black") +
  geom_rug() +
  stat_function(fun = dchisq, args = list(df = df_chi),
                col = "red", size = 2, alpha = .7,
                xlim = range(c(0, dM, upper.quantiles))) +
  geom_segment(data = cut.points,
               aes(x = upper.quantiles, xend = upper.quantiles,
                   y = 0, yend = density.at.quantiles),
               col = "blue", size = 2) +
  labs(title = "Mahalanobis Distance",
       x = "Mahalanobis distances and cut points",
       y = "Histogram and density")



# Table of 'level of surprise' for the 38 ingots
ppm$dM <- dM
ppm$surprise <- cut(ppm$dM,
                    breaks= c(0, upper.quantiles, Inf),
                    labels=c("Typical", "Somewhat", "Surprising", "Very"))
pander(table(ppm$surprise))



# pairs plot
ggpairs(ppm, columns=1:17,
        ggplot2::aes(col=surprise, alpha=.5),
        upper = list(continuous = "density", combo = "box_no_facet")) +
  ggplot2::scale_color_manual(values=c("lightgray", "green", "blue", "red")) +
  ggplot2::theme(axis.text.x = element_text(angle=90, hjust=1))+
  labs(title = "Pairs plot")



# summary table (mean, med, LQ, UP, Min, Max)
summary_table <- ppm %>%
  dplyr::select(any_of(element_order)) %>%   # force dplyr::select
  pivot_longer(cols = everything(),
               names_to = "Element", values_to = "Concentration") %>%
  mutate(Element = factor(Element, levels = element_order),
         Concentration = as.numeric(Concentration)) %>%
  group_by(Element) %>%
  summarise(
    Mean   = mean(Concentration, na.rm = TRUE),
    Median = median(Concentration, na.rm = TRUE),
    LQ     = quantile(Concentration, 0.25, na.rm = TRUE),
    UQ     = quantile(Concentration, 0.75, na.rm = TRUE),
    Min    = min(Concentration, na.rm = TRUE),
    Max    = max(Concentration, na.rm = TRUE),
    .groups = "drop"
  )

kable(summary_table, digits = 3)



# PCA
X <- ppm_cc %>%
  mutate(across(everything(), log1p)) %>%  # natural log; use log() if no zeros
  scale() %>%
  as.matrix()

## 1) PCA
pca <- prcomp(X, center = FALSE, scale. = FALSE)  # already centered/scaled above

# Variance explained (for scree)
var_explained <- (pca$sdev^2)
pev <- var_explained / sum(var_explained)
scree_df <- tibble(PC = paste0("PC", seq_along(pev)),
                   Proportion = pev,
                   Cumulative = cumsum(pev))

# Scores (PC coordinates of samples)
scores <- as_tibble(pca$x[, 1:4, drop = FALSE]) %>%
  rename(PC1 = 1, PC2 = 2, PC3 = 3, PC4 = 4)

# Loadings (element contributions)
loadings <- as_tibble(pca$rotation[, 1:2, drop = FALSE], rownames = "Element") %>%
  rename(PC1 = 2, PC2 = 3)

col_var <- if ("surprise" %in% names(ppm)) {
  factor(ppm$surprise[complete.cases(ppm)], 
         levels = c("Typical","Somewhat","Surprising","Very"))
} else {
  factor(rep("All", nrow(scores)))
}

# Scree plot
ggplot(scree_df, aes(x = factor(PC, levels = paste0("PC", seq_along(pev))),
                     y = Proportion)) +
  geom_col(fill = "grey70") +
  geom_line(aes(y = Cumulative, group = 1)) +
  geom_point(aes(y = Cumulative)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "PCA Scree Plot",
       x = NULL, y = "Variance explained (and cumulative)") +
  theme_minimal(base_size = 12)

# PC1 vs PC2 scatter
pc12_df <- scores %>% mutate(Color = col_var)

ggplot(pc12_df, aes(x = PC1, y = PC2, color = Color)) +
  geom_point(size = 2.8, alpha = 0.9) +
  stat_ellipse(type = "norm", level = 0.68, linetype = 2, linewidth = 0.5,
               show.legend = FALSE) +
  labs(title = "PCA Scores: PC1 vs PC2",
       x = paste0("PC1 (", percent(pev[1], accuracy = 0.1), ")"),
       y = paste0("PC2 (", percent(pev[2], accuracy = 0.1), ")"),
       color = "Group") +
  theme_minimal(base_size = 12)

# Biplot
arrow_scale <- 0.9 * max(abs(scores$PC1), abs(scores$PC2)) /
  max(sqrt(loadings$PC1^2 + loadings$PC2^2))
loadings_scaled <- loadings %>%
  mutate(PC1 = PC1 * arrow_scale,
         PC2 = PC2 * arrow_scale)

ggplot() +
  geom_point(data = pc12_df, aes(PC1, PC2, color = Color), size = 2.5, alpha = 0.9) +
  geom_segment(data = loadings_scaled,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.012, "npc")), linewidth = 0.6) +
  geom_text_repel(data = loadings_scaled, aes(PC1, PC2, label = Element),
                  size = 3.3, max.overlaps = Inf) +
  labs(title = "PCA Biplot: PC1 vs PC2",
       x = paste0("PC1 (", percent(pev[1], accuracy = 0.1), ")"),
       y = paste0("PC2 (", percent(pev[2], accuracy = 0.1), ")"),
       color = "Group") +
  theme_minimal(base_size = 12)
