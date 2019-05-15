require(tidyverse)
require(dplyr)
Proj.Home = "/Users/Nate/Box Sync/Lowen_Lab/Data_Notebook/Documentation" #For Nate's Macbook Pro

Data.Path = file.path(Proj.Home,"Simulation_Test_Results")
setwd(Data.Path)

# Analyze ----
Fig4 = read.csv(file = file.path(Data.Path,"Fig4_Data_SHORT.csv"))

Fig4.Peak = Fig4 %>%
  group_by(Sim,Spread_Frac,Pp,Virus_D) %>%
  dplyr::summarise(Max_Productive = max(Productive),
            Max_Virions = max(Virions),
            Max_Semi = max(Infected - Productive))

Fig4.Peak

Fig4.Hours = Fig4 %>%
  filter((t * 4) %in% (1:(96 * 4)))
write.csv(Fig4.Peak,file=file.path(Data.Path,"Fig4_Peak_Data.csv"),row.names=FALSE)
write.csv(Fig4.Hours,file=file.path(Data.Path,"Fig4_Hours_Data.csv"),row.names=FALSE)
# Visualize ----

Fig4.Peak = read.csv(file = file.path(Data.Path,"Fig4_Peak_Data.csv")) %>%
  mutate(Virus_D = log10(Virus_D / 3600)) %>%
  filter(Virus_D < 5,
         Spread_Frac %in% c(0,0.25)) %>%
  mutate(Spread_Frac = recode(Spread_Frac, "0.25" = "25% Local Spread","0" = "Diffusion Only"))

Fig4.Hours = read.csv(file=file.path(Data.Path,"Fig4_Hours_Data.csv")) %>%
  mutate(Virus_D = log10(Virus_D / 3600)) %>%
  filter(Virus_D < 5,
         Spread_Frac %in% c(0,0.25)) %>%
  mutate(Spread_Frac = recode(Spread_Frac, "0.25" = "25% Local Spread","0" = "Diffusion Only"))
  
  

# Plot Stock ----
Plot.Stock = ggplot() +
  scale_color_manual(values = c("Diffusion Only" = "midnightblue","25% Local Spread" = "steelblue1"), guide = FALSE) +
  scale_fill_manual(values = c("Diffusion Only" = "midnightblue","25% Local Spread" = "steelblue1"), guide = FALSE) +
  scale_y_continuous(limits = c(0,105)) +
  geom_vline(xintercept = log10((5.825)), lty = 2) +
  coord_cartesian(xlim = c(-0.5,4)) +
  theme(text=element_text(size = 13,face="bold"),
        strip.text.x=element_text(size=rel(1.5),margin=margin(0,0,3,0)),
        strip.text.y=element_text(size=rel(1.5),margin=margin(0,0,0,0),angle=0),
        strip.background = element_blank(),
        plot.title = element_text(size=rel(1.5)),
        axis.text.x=element_text(angle = 0,vjust=0.75,size=rel(1.5),color = "black"),
        axis.text.y=element_text(size=rel(1.5),color = "black"),
        axis.line.x = element_line(size=0.5),
        axis.line.y = element_line(size=0.5),
        axis.ticks.x = element_line(size=0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.title.y = element_text(size=rel(1.2),color = "black"),
        axis.title.x = element_text(size=rel(1.2)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

# Raw r0 ----
Early.Time = 12
df4 = Fig4.Hours %>%
  filter(t == Early.Time) %>%
  mutate(Slope = pmax(log10(Productive),0)) %>%
  dplyr::select(Pp,Virus_D,Spread_Frac,Sim,Slope) %>%
  group_by(Pp,Virus_D,Spread_Frac) %>%
  dplyr::summarise(r0 = mean(Slope)) %>%
  ungroup %>%
  mutate(Pp = Pp + 0.005) %>%
  mutate(Pp = round(Pp,2))
df4

Plot.Stock +
  geom_point(data = df4,
             aes(x = Virus_D,
                 y = r0,
                 color = Spread_Frac)) +
  geom_smooth(data = df4,
              aes(x = Virus_D,
                  y = r0,
                  linetype = factor(Pp),
                  color = Spread_Frac,
                  fill = Spread_Frac,
                  group = interaction(Spread_Frac,Pp)),
              lwd = 1.2,
              alpha = 0.3) +
  labs(y = TeX("\\textbf{$\\log_1_0$ Cells Infected in 12 Hours}"),
       x = TeX("\\textbf{Diffusion Coefficient ($\\log_1_0$ $\\mu$m^2/s)}"),
       linetype = TeX("\\textbf{$\\P_P$}"),
       size = TeX("\\textbf{$\\P_P$}"),
       color = NULL,
       fill = NULL) +
  scale_linetype_manual(values = c("1" = 1,"0.58" = 1212)) +
  scale_y_continuous(limits = c(0,4)) +
  scale_color_manual(values = c("Diffusion Only" = "midnightblue","25% Local Spread" = "steelblue1")) +
  scale_fill_manual(values = c("Diffusion Only" = "midnightblue","25% Local Spread" = "steelblue1"))

ggsave('Figures/4A_r0_Raw.pdf',
       width = 7,
       height = 5,
       unit = "in")

#r0 Cost (log10 Productive cells in 12 hours) ----

Early.Time = 12

df1 = Fig4.Hours %>%
  filter(t == Early.Time,
         Pp == 1) %>%
  mutate(Slope = log10(Productive) / t) %>%
  group_by(Virus_D,Spread_Frac) %>%
  dplyr::summarise(Intact_r0 = mean(Slope))

df2 = Fig4.Hours %>%  
  filter(t == Early.Time,
         Pp < 1) %>%
  mutate(Slope = pmax(log10(Productive) / t,0)) %>%
  dplyr::select(Pp,Virus_D,Spread_Frac,Sim,Slope)

df3 = right_join(df1,df2) %>%
  mutate(Cost = (1 - (Slope / Intact_r0)) * 100) %>%
  group_by(Virus_D,Spread_Frac) %>%
  dplyr::summarise(Cost = mean(Cost))

Plot.Stock +
  geom_point(data = df3,
             aes(x = Virus_D,
                 y = Cost,
                 color = Spread_Frac,
                 group = Spread_Frac)) +
  geom_smooth(data = df3,
              aes(x = Virus_D,
                  y = Cost,
                  group = Spread_Frac,
                  color = Spread_Frac,
                  fill = Spread_Frac),
              lty = 1212,
              alpha = 0.3) +
  labs(y = "% Reduction in Initial Growth Rate",
       x = TeX("\\textbf{Diffusion Coefficient ($\\log_1_0$ $\\mu$m^2/s)}"))

ggsave('Figures/4B_r0_Reduction.pdf',
       width = 6,
       height = 5,
       unit = "in")

#Time to 100 Cells ----
Target_Cells = 100

df1 = Fig4.Hours %>%
  filter(Pp == 1,
         Productive >= Target_Cells) %>%
  group_by(Spread_Frac,Virus_D,Sim) %>%
  dplyr::summarise(To_Target_Cells = min(t)) %>%
  ungroup() %>%
  group_by(Spread_Frac,Virus_D) %>%
  dplyr::summarise(Intact_To_Target_Cells = mean(To_Target_Cells))

df2 = Fig4.Hours %>%
  filter(Pp < 1,
         Productive >= Target_Cells) %>%
  group_by(Spread_Frac,Virus_D,Sim) %>%
  dplyr::summarise(To_Target_Cells = min(t)) %>%
  dplyr::select(Spread_Frac,Virus_D,Sim,To_Target_Cells)

df3 = right_join(df1,df2) %>%
  mutate(Target_Cell_Cost = (To_Target_Cells - Intact_To_Target_Cells) / Intact_To_Target_Cells * 100) %>%
  group_by(Virus_D,Spread_Frac) %>%
  dplyr::summarise(Target_Cell_Cost = mean(Target_Cell_Cost))

Plot.Stock +
  geom_point(data = df3,
             aes(x = Virus_D,
                 y = Target_Cell_Cost,
                 color = Spread_Frac,
                 group = Spread_Frac)) + 
  geom_smooth(data = df3,
              aes(x = Virus_D,
                  y = Target_Cell_Cost,
                  group = Spread_Frac,
                  color = Spread_Frac,
                  fill = Spread_Frac),
              lty = 1212,
              alpha = 0.3) +
  scale_y_continuous(limits = c(0,600)) +
  labs(y = "% Increase in Time to Infect 100 Cells",
       x = TeX("\\textbf{Diffusion Coefficient ($\\log_1_0$ $\\mu$m^2/s)}"))

ggsave('Figures/4C_Increase_in_Time_to_100_Cells.pdf',
       width = 6,
       height = 5,
       unit = "in")

#Time to 1e5 Virions ----

Target_Virions = 1e5
df1 = Fig4.Hours %>%
  filter(Pp == 1,
         Virions >= Target_Virions) %>%
  group_by(Spread_Frac,Virus_D,Sim) %>%
  dplyr::summarise(To_Target_Virions = min(t)) %>%
  ungroup() %>%
  group_by(Spread_Frac,Virus_D) %>%
  dplyr::summarise(Intact_To_Target_Virions = mean(To_Target_Virions))


df2 = Fig4.Hours %>%
  filter(Pp < 1,
         Virions >= Target_Virions) %>%
  group_by(Spread_Frac,Virus_D,Sim) %>%
  dplyr::summarise(To_Target_Virions = min(t)) %>%
  dplyr::select(Spread_Frac,Virus_D,Sim,To_Target_Virions)

df3 = right_join(df1,df2) %>%
  mutate(Target_Virion_Cost = (To_Target_Virions - Intact_To_Target_Virions) / Intact_To_Target_Virions * 100) %>%
  group_by(Virus_D,Spread_Frac) %>%
  dplyr::summarise(Target_Virion_Cost = mean(Target_Virion_Cost))

Plot.Stock +
  geom_point(data = df3,
             aes(x = Virus_D,
                 y = Target_Virion_Cost,
                 color = Spread_Frac,
                 group = Spread_Frac)) + 
  geom_smooth(data = df3,
              aes(x = Virus_D,
                  y = Target_Virion_Cost,
                  group = Spread_Frac,
                  color = Spread_Frac,
                  fill = Spread_Frac),
              lty = 1212,
              alpha = 0.3) +
  scale_y_continuous(limits = c(0,400)) +
  labs(y = TeX("\\textbf{% Increase in Time to Yield 10^5 Virions}"),
       x = TeX("\\textbf{Diffusion Coefficient ($\\log_1_0$ $\\mu$m^2/s)}"))

ggsave('Figures/4D_Increase_in_Time_to_1e5_Virions.pdf',
       width = 6,
       height = 5,
       unit = "in")

# Supplement ----
# Representative Time Course ----
ggplot() +
  geom_line(data = Fig4.Hours %>% filter(Pp == 1,Sim == 1,Spread_Frac == "25% Local Spread"),
            aes(x = t,
                y = Productive,
                color = Virus_D,
                group = Virus_D),
            lwd = 1.2) +
  scale_color_viridis(guide = FALSE) +
  scale_y_continuous(limits = c(0,10000),labels = comma) +
  scale_x_continuous(limits = c(0,96),breaks = c(0,24,48,72)) +
  labs(x = TeX("\\textbf{Time Post-Inoculation (Hours)}"),
       y = TeX("\\textbf{Productively Infected Cells}")) +
  theme(text=element_text(size = 13,face="bold"),
        strip.text.x=element_text(size=rel(1.5),margin=margin(0,0,3,0)),
        strip.text.y=element_text(size=rel(1.5),margin=margin(0,0,0,0),angle=0),
        strip.background = element_blank(),
        plot.title = element_text(size=rel(1.5)),
        axis.text.x=element_text(angle = 0,vjust=0.75,size=rel(1.5),color = "black"),
        axis.text.y=element_text(size=rel(1.5),color = "black"),
        axis.line.x = element_line(size=0.5),
        axis.line.y = element_line(size=0.5),
        axis.ticks.x = element_line(size=0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.title.y = element_text(size=rel(1.2),color = "black"),
        axis.title.x = element_text(size=rel(1.2)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

ggsave('Figures/Supp2C_Time_Course_PP1_Local_Spread.pdf',
       width = 5,
       height = 5,
       unit = "in")

ggplot() +
  geom_line(data = Fig4.Hours %>% filter(Pp == 0.575,Sim == 1,Spread_Frac == "25% Local Spread"),
            aes(x = t,
                y = Productive,
                color = Virus_D,
                group = Virus_D),
            lwd = 1.2) +
  scale_color_viridis(guide = FALSE) +
  scale_y_continuous(limits = c(0,10000),labels = comma) +
  scale_x_continuous(limits = c(0,96),breaks = c(0,24,48,72)) +
  labs(x = TeX("\\textbf{Time Post-Inoculation (Hours)}"),
       y = TeX("\\textbf{Productively Infected Cells}")) +
  theme(text=element_text(size = 13,face="bold"),
        strip.text.x=element_text(size=rel(1.5),margin=margin(0,0,3,0)),
        strip.text.y=element_text(size=rel(1.5),margin=margin(0,0,0,0),angle=0),
        strip.background = element_blank(),
        plot.title = element_text(size=rel(1.5)),
        axis.text.x=element_text(angle = 0,vjust=0.75,size=rel(1.5),color = "black"),
        axis.text.y=element_text(size=rel(1.5),color = "black"),
        axis.line.x = element_line(size=0.5),
        axis.line.y = element_line(size=0.5),
        axis.ticks.x = element_line(size=0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.title.y = element_text(size=rel(1.2),color = "black"),
        axis.title.x = element_text(size=rel(1.2)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

ggsave('Figures/Supp2G_Time_Course_PP_058_Local_Spread.pdf',
       width = 6,
       height = 5,
       unit = "in")

ggplot() +
  geom_line(data = Fig4.Hours %>% filter(Pp == 1,Sim == 1,Spread_Frac == "Diffusion Only"),
            aes(x = t,
                y = Productive,
                color = Virus_D,
                group = Virus_D),
            lwd = 1.2) +
  scale_color_viridis(guide = FALSE) +
  scale_y_continuous(limits = c(0,10000),labels = comma) +
  scale_x_continuous(limits = c(0,96),breaks = c(0,24,48,72)) +
  labs(x = TeX("\\textbf{Time Post-Inoculation (Hours)}"),
       y = TeX("\\textbf{Productively Infected Cells}")) +
  theme(text=element_text(size = 13,face="bold"),
        strip.text.x=element_text(size=rel(1.5),margin=margin(0,0,3,0)),
        strip.text.y=element_text(size=rel(1.5),margin=margin(0,0,0,0),angle=0),
        strip.background = element_blank(),
        plot.title = element_text(size=rel(1.5)),
        axis.text.x=element_text(angle = 0,vjust=0.75,size=rel(1.5),color = "black"),
        axis.text.y=element_text(size=rel(1.5),color = "black"),
        axis.line.x = element_line(size=0.5),
        axis.line.y = element_line(size=0.5),
        axis.ticks.x = element_line(size=0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.title.y = element_text(size=rel(1.2),color = "black"),
        axis.title.x = element_text(size=rel(1.2)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

ggsave('Figures/Supp2A_Time_Course_PP1_No_Spread.pdf',
       width = 5,
       height = 5,
       unit = "in")

ggplot() +
  geom_line(data = Fig4.Hours %>% filter(Pp == 0.575,Sim == 1,Spread_Frac == "Diffusion Only"),
            aes(x = t,
                y = Productive,
                color = Virus_D,
                group = Virus_D),
            lwd = 1.2) +
  scale_color_viridis(guide = FALSE) +
  scale_y_continuous(limits = c(0,10000),labels = comma) +
  scale_x_continuous(limits = c(0,96),breaks = c(0,24,48,72)) +
  labs(x = TeX("\\textbf{Time Post-Inoculation (Hours)}"),
       y = TeX("\\textbf{Productively Infected Cells}")) +
  theme(text=element_text(size = 13,face="bold"),
        strip.text.x=element_text(size=rel(1.5),margin=margin(0,0,3,0)),
        strip.text.y=element_text(size=rel(1.5),margin=margin(0,0,0,0),angle=0),
        strip.background = element_blank(),
        plot.title = element_text(size=rel(1.5)),
        axis.text.x=element_text(angle = 0,vjust=0.75,size=rel(1.5),color = "black"),
        axis.text.y=element_text(size=rel(1.5),color = "black"),
        axis.line.x = element_line(size=0.5),
        axis.line.y = element_line(size=0.5),
        axis.ticks.x = element_line(size=0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.title.y = element_text(size=rel(1.2),color = "black"),
        axis.title.x = element_text(size=rel(1.2)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

ggsave('Figures/Supp2E_Time_Course_PP_058_No_Spread.pdf',
       width = 6,
       height = 5,
       unit = "in")

# Virion Time Course ----
ggplot() +
  geom_line(data = Fig4.Hours %>% filter(Pp == 1,Sim == 1,Spread_Frac == "25% Local Spread"),
            aes(x = t,
                y = Virions %>% log10,
                color = Virus_D,
                group = Virus_D),
            lwd = 1.2) +
  scale_color_viridis(guide = FALSE) +
  scale_x_continuous(limits = c(0,96),breaks = c(0,24,48,72)) +
  labs(x = TeX("\\textbf{Time Post-Inoculation (Hours)}"),
       y = TeX("\\textbf{$\\log_1_0$ Virions}")) +
  theme(text=element_text(size = 13,face="bold"),
        strip.text.x=element_text(size=rel(1.5),margin=margin(0,0,3,0)),
        strip.text.y=element_text(size=rel(1.5),margin=margin(0,0,0,0),angle=0),
        strip.background = element_blank(),
        plot.title = element_text(size=rel(1.5)),
        axis.text.x=element_text(angle = 0,vjust=0.75,size=rel(1.5),color = "black"),
        axis.text.y=element_text(size=rel(1.5),color = "black"),
        axis.line.x = element_line(size=0.5),
        axis.line.y = element_line(size=0.5),
        axis.ticks.x = element_line(size=0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.title.y = element_text(size=rel(1.2),color = "black"),
        axis.title.x = element_text(size=rel(1.2)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

ggsave('Figures/Supp2D_Virion_Time_Course_PP1_Local_Spread.pdf',
       width = 5,
       height = 5,
       unit = "in")

ggplot() +
  geom_line(data = Fig4.Hours %>% filter(Pp == 0.575,Sim == 1,Spread_Frac == "25% Local Spread"),
            aes(x = t,
                y = Virions %>% log10,
                color = Virus_D,
                group = Virus_D),
            lwd = 1.2) +
  scale_color_viridis(guide = FALSE) +
  scale_x_continuous(limits = c(0,96),breaks = c(0,24,48,72)) +
  labs(x = TeX("\\textbf{Time Post-Inoculation (Hours)}"),
       y = TeX("\\textbf{$\\log_1_0$ Virions}")) +
  theme(text=element_text(size = 13,face="bold"),
        strip.text.x=element_text(size=rel(1.5),margin=margin(0,0,3,0)),
        strip.text.y=element_text(size=rel(1.5),margin=margin(0,0,0,0),angle=0),
        strip.background = element_blank(),
        plot.title = element_text(size=rel(1.5)),
        axis.text.x=element_text(angle = 0,vjust=0.75,size=rel(1.5),color = "black"),
        axis.text.y=element_text(size=rel(1.5),color = "black"),
        axis.line.x = element_line(size=0.5),
        axis.line.y = element_line(size=0.5),
        axis.ticks.x = element_line(size=0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.title.y = element_text(size=rel(1.2),color = "black"),
        axis.title.x = element_text(size=rel(1.2)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

ggsave('Figures/Supp2H_Virion_Time_Course_PP_058_Local_Spread.pdf',
       width = 6,
       height = 5,
       unit = "in")

ggplot() +
  geom_line(data = Fig4.Hours %>% filter(Pp == 1,Sim == 1,Spread_Frac == "Diffusion Only"),
            aes(x = t,
                y = Virions %>% log10,
                color = Virus_D,
                group = Virus_D),
            lwd = 1.2) +
  scale_color_viridis(guide = FALSE) +
  scale_x_continuous(limits = c(0,96),breaks = c(0,24,48,72)) +
  labs(x = TeX("\\textbf{Time Post-Inoculation (Hours)}"),
       y = TeX("\\textbf{$\\log_1_0$ Virions}")) +
  theme(text=element_text(size = 13,face="bold"),
        strip.text.x=element_text(size=rel(1.5),margin=margin(0,0,3,0)),
        strip.text.y=element_text(size=rel(1.5),margin=margin(0,0,0,0),angle=0),
        strip.background = element_blank(),
        plot.title = element_text(size=rel(1.5)),
        axis.text.x=element_text(angle = 0,vjust=0.75,size=rel(1.5),color = "black"),
        axis.text.y=element_text(size=rel(1.5),color = "black"),
        axis.line.x = element_line(size=0.5),
        axis.line.y = element_line(size=0.5),
        axis.ticks.x = element_line(size=0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.title.y = element_text(size=rel(1.2),color = "black"),
        axis.title.x = element_text(size=rel(1.2)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

ggsave('Figures/Supp2B_Virion_Time_Course_PP1_No_Spread.pdf',
       width = 5,
       height = 5,
       unit = "in")

ggplot() +
  geom_line(data = Fig4.Hours %>% filter(Pp == 0.575,Sim == 1,Spread_Frac == "Diffusion Only"),
            aes(x = t,
                y = Virions %>% log10,
                color = Virus_D,
                group = Virus_D),
            lwd = 1.2) +
  scale_color_viridis(guide = FALSE) +
  scale_x_continuous(limits = c(0,96),breaks = c(0,24,48,72)) +
  labs(x = TeX("\\textbf{Time Post-Inoculation (Hours)}"),
       y = TeX("\\textbf{$\\log_1_0$ Virions}")) +
  theme(text=element_text(size = 13,face="bold"),
        strip.text.x=element_text(size=rel(1.5),margin=margin(0,0,3,0)),
        strip.text.y=element_text(size=rel(1.5),margin=margin(0,0,0,0),angle=0),
        strip.background = element_blank(),
        plot.title = element_text(size=rel(1.5)),
        axis.text.x=element_text(angle = 0,vjust=0.75,size=rel(1.5),color = "black"),
        axis.text.y=element_text(size=rel(1.5),color = "black"),
        axis.line.x = element_line(size=0.5),
        axis.line.y = element_line(size=0.5),
        axis.ticks.x = element_line(size=0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.title.y = element_text(size=rel(1.2),color = "black"),
        axis.title.x = element_text(size=rel(1.2)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

ggsave('Figures/Supp2F_Virion_Time_Course_PP_058_No_Spread.pdf',
       width = 6,
       height = 5,
       unit = "in")

D_Legend = ggplot() +
  geom_line(data = Fig4.Hours %>% filter(Pp == 0.575,Sim == 1,Spread_Frac == "Diffusion Only"),
            aes(x = t,
                y = Virions %>% log10,
                color = Virus_D,
                group = Virus_D),
            lwd = 1.2) +
  scale_color_viridis() +
  scale_x_continuous(limits = c(0,96),breaks = c(0,24,48,72)) +
  labs(color = TeX("\\textbf{Diffusion Coefficient ($\\log_1_0$ $\\mu$m^2/s)}")) +
  theme(text=element_text(size = 13,face="bold"),
        strip.text.x=element_text(size=rel(1.5),margin=margin(0,0,3,0)),
        strip.text.y=element_text(size=rel(1.5),margin=margin(0,0,0,0),angle=0),
        strip.background = element_blank(),
        plot.title = element_text(size=rel(1.5)),
        axis.text.x=element_text(angle = 0,vjust=0.75,size=rel(1.5),color = "black"),
        axis.text.y=element_text(size=rel(1.5),color = "black"),
        axis.line.x = element_line(size=0.5),
        axis.line.y = element_line(size=0.5),
        axis.ticks.x = element_line(size=0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.title.y = element_text(size=rel(1.2),color = "black"),
        axis.title.x = element_text(size=rel(1.2)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

plot(g_legend(D_Legend))
ggsave('Figures/Supp2_D_Legend.pdf',
       width = 4,
       height = 4,
       unit = "in")

# Max # Semi-Infected Cells ----

df1 = Fig4.Peak %>%
  group_by(Virus_D,Spread_Frac) %>%
  dplyr::summarise(Max_Semi = mean(Max_Semi))

Plot.Stock +
  geom_point(data = df1,
             aes(x = Virus_D,
                 y = Max_Semi,
                 color = Spread_Frac)) +
  geom_smooth(data = df1,
              aes(x = Virus_D,
                  y = Max_Semi,
                  group = Spread_Frac,
                  color = Spread_Frac,
                  fill = Spread_Frac),
              lty = 1212,
              alpha = 0.3) +
  scale_y_continuous(limits = c(-100,3200), breaks = c(0,1000,2000,3000)) +
  labs(y = TeX("\\textbf{Peak # Cells with IVGs}"),
       x = TeX("\\textbf{Diffusion Coefficient ($\\log_1_0$ $\\mu$m^2/s)}")) 

ggsave('Figures/6B_Max_Semi.pdf',
       width = 6,
       height = 5,
       unit = "in")