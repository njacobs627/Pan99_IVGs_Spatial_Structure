# Revised Figure 3 Solution

require(tidyverse)

pp = c(1:10/10)
moi = 10^(seq(from = -4, to = 1, by = 0.05))
n = 1:5000

df1 = matrix(nrow = length(n) * length(moi) * length(pp),
             ncol = 3,
             data = 0) %>% data.frame

colnames(df1) = c("PP","MOI","n")

df1$PP = rep(pp, times = length(n),each = length(moi))
df1$MOI = rep(moi, times = length(n) * length(pp))
df1$N = rep(n,each = length(moi) * length(pp))

df2 = df1 %>%
  mutate(p_N = dpois(x = N, lambda = MOI),
         p_Complete_from_N = (1 - (1 - PP) ^ N) ^ 8,
         p_Complete = p_N * p_Complete_from_N) %>%
  group_by(PP,MOI) %>%
  summarise(p_Complete = sum(p_Complete))

df3 = df2 %>%
  ungroup() %>%
  mutate(p_Infected_Comp = pmin(p_Complete * 1e4,1),
         p_Infected_NoComp = pmin(1 / PP^8 * MOI * 1e4,1))

# Fig. 3A

ggplot() +
  geom_line(data = df3,
            aes(x = MOI %>% log10,
                y = p_Infected_Comp,
                group = PP,
                color = PP)) +
  scale_color_viridis_c()

df4 = df3 %>% 
  filter(p_Infected_Comp >= 0.5) %>%
  group_by(PP) %>%
  summarise(min_Comp_MOI = min(MOI))
df5 = df3 %>% 
  filter(p_Infected_NoComp >= 0.5) %>%
  group_by(PP) %>%
  summarise(min_NoComp_MOI = min(MOI))

df5 = matrix(nrow = 10,
             ncol = 2,
             data = 0) %>% data.frame
colnames(df5) = c("PP","NoComp_MOI")
df5$PP = seq(from = 0.1, to = 1, by = 0.1)
df5$min_NoComp_MOI = 1 / 1e4 / PP^8

#Fig. 3B
ggplot() +
  geom_line(data = df4,
            aes(x = PP,
                y = min_Comp_MOI %>% log10)) +
  geom_line(data = df5,
            aes(x = PP,
                y = min_NoComp_MOI %>% log10),
            lty = 2)