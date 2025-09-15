####################
## STAT394 GP EDA ##
####################

df <- read_csv("orichalcum_ingots_dataset.csv", skip = 1)

colnames(df) <- c("Sample","Cu","Zn","Pb","Fe","Ni","Ag","Sb","As","Co",
                  "Cd","Sn","Te","Bi","Mn","Li","Al","V","Cr","Rb","Sr","Ba","Cu_Zn","w_w")
df <- df %>% select(-w_w)

# puts elements in same order as dataset otherwise it goes in alphabetical order
element_order <- c("Cu","Zn","Pb","Fe","Ni","Ag","Sb","As","Co","Cd","Sn",
                   "Te","Bi","Mn","Li","Al","V","Cr","Rb","Sr","Ba","Cu_Zn")


# Summary stats histograms plot
df %>% 
  pivot_longer(-Sample, names_to="Element", values_to="Concentration") %>%
  mutate(Element = factor(Element, levels = element_order)) %>%
  ggplot(aes(x=Concentration)) +
  geom_histogram(bins=30, fill="skyblue") +
  facet_wrap(~Element, scales="free")

#correlation plot
corrplot(cor(df[,-1]), method="color")


ggplot(df, aes(x=Cr, y=V, color=Ni)) +
  geom_point(size=3)

ggplot(df, aes(x=Zn, y=Cu_Zn, color=Pb)) +
  geom_point(size=3)

ggplot(df, aes(x=Pb, y=Sn, color=Sb)) +
  geom_point(size=3)

ggplot(df, aes(x=Ni, y=Co, color=Fe)) +
  geom_point(size=3)

df %>% 
  pivot_longer(-Sample, names_to="Element", values_to="Concentration") %>%
  ggplot(aes(x=Element, y=Concentration)) +
  geom_boxplot() +
  coord_flip()

df %>% 
  pivot_longer(-Sample, names_to = "Element", values_to = "Concentration") %>%
  ggplot(aes(x = Element, y = Concentration)) +
  geom_boxplot() +
  coord_flip() +
  scale_y_log10()
