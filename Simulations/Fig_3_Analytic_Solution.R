# Revised Figure 3 Solution

require(tidyverse)

pp = seq(from = 0.1, to = 1, by = 0.1)
moi = 10^(seq(from = -6, to = 3, by = 0.05))
n = 1:1000

df1 = matrix(nrow = length(n) * length(moi) * length(pp),
             ncol = 3,
             data = 0) %>% data.frame

colnames(df1) = c("PP","MOI","n")

df1$PP = rep(pp, times = length(n),each = length(moi))
df1$MOI = rep(moi, times = length(n) * length(pp))
df1$V = rep(n,each = length(moi) * length(pp))

df2 = df1 %>%
  mutate(p_V = dpois(x = V, lambda = MOI),
         p8_V_Comp = (1 - (1 - PP) ^ V) ^ 8,
         p8_V_NoComp = 1 - (1 - PP^8)^V,
         p8_NoComp = p_V * p8_V_NoComp,
         p8_Comp = p_V * p8_V_Comp) %>%
  group_by(PP,MOI) %>%
  summarise(p8_NoComp = sum(p8_NoComp),
            p8_Comp = sum(p8_Comp)) %>%
  mutate(Comp_Benefit = p8_Comp / p8_NoComp)

df3 = df2 %>%
  ungroup() %>%
  mutate(p_Infected_Comp = pmin(p8_Comp * 1e6,1),
         p_Infected_NoComp = pmin(p8_NoComp * 1e6,1))

# Fig. 3A

ggplot() +
  geom_line(data = df3,
            aes(x = MOI %>% log10,
                y = p_Infected_Comp,
                group = PP,
                color = PP)) +
  geom_line(data = df3,
            aes(x = MOI %>% log10,
                y = p_Infected_NoComp,
                group = PP,
                color = PP),
            lty = 2) +
  scale_color_viridis_c()

df4 = df3 %>% 
  filter(p_Infected_Comp >= 0.5) %>%
  group_by(PP) %>%
  summarise(min_Comp_MOI = min(MOI))
df5 = df3 %>% 
  filter(p_Infected_NoComp >= 0.5) %>%
  group_by(PP) %>%
  summarise(min_NoComp_MOI = min(MOI))

#Fig. 3B
ggplot() +
  geom_line(data = df4,
            aes(x = PP,
                y = min_Comp_MOI %>% log10)) +
  geom_line(data = df5,
            aes(x = PP,
                y = min_NoComp_MOI %>% log10),
            lty = 2)