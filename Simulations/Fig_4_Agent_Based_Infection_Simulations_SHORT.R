#### Load Packages ----

require(foreach)
require(doMC)
require(Rmisc)
require(ggplot2)

### Set File Paths ----

Proj.Home = "/Users/Nate/Box Sync/Lowen_Lab/Data_Notebook/Documentation" #For Nate's Macbook Pro

Data.Path = file.path(Proj.Home,"Simulation_Test_Results")

#### Define Parameters ----

# Static parameters
Grid.X = 100
Grid.Y = Grid.X
Cell_Width = 30 # Assume each cell is a 30 um x 30 um square

# Rates from Handel et al. J R Soc Interface 2014
Death_Rate = 2 / 24
Birth_Rate = 2 / 24

Decay_Rate = -log(0.5) / 4 # Assume 1/2 life of 4 hours
Exclude_Rate = 1 / 8 # From Nicolle's thesis
Revert_Rate = 1 / 12

Infect_Rate = 1 / 6 # 6 hour eclipse phase
# Incorporates the 6-hour lag before segments are "useful"
# So "Segment #" is the number of useful segments (cells become productive instantly once 8th segment is delivered)

# Varied parameters: first, set defaults
Virus.D = 5.825 * 60 * 60 # 5.825 um^2 / s based on Einsten-Stokes equation
Attach_Rate = 1 / (10 / 60) # 20 minutes
Detach_Rate = 10 / 24 # 10/day = 10/24 hr
Burst_Rate = 962 / 24
Start_Radius = 0

Sim_Num = 2
Pp = c(0.575,1)
Spread_Frac = c(0,0.25)
Virus.D = (5.825 * 60 * 60) *
  c(0.1, 1, 10, 100)
Spread_Frac = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)
ts = 1/20

Duration = 48

# Scaling to ts
Infect_Rate = Infect_Rate * ts
Attach_Rate = Attach_Rate * ts
Detach_Rate = Detach_Rate * ts
Burst_Rate = Burst_Rate * ts
Death_Rate = Death_Rate * ts
Birth_Rate = Birth_Rate * ts
Exclude_Rate = Exclude_Rate * ts
Decay_Rate = Decay_Rate * ts
Revert_Rate = Revert_Rate * ts
Virus.D = Virus.D * ts

Infect_Rate_Num = length(Infect_Rate)
Attach_Rate_Num = length(Attach_Rate)
Detach_Rate_Num = length(Detach_Rate)
Burst_Rate_Num = length(Burst_Rate)
Death_Rate_Num = length(Death_Rate)
Birth_Rate_Num = length(Birth_Rate)
Decay_Rate_Num = length(Decay_Rate)
Revert_Rate_Num = length(Revert_Rate)
Exclude_Rate_Num = length(Exclude_Rate)
Virus.D_Num = length(Virus.D)
Pp_Num = length(Pp)
Spread_Frac_Num = length(Spread_Frac)

Start_Radius_Num = length(Start_Radius)

#### Populations -----
# Initialize Cells ----

Cell.Headers = c("Cell.ID","X","Y","Distance","Living","Segments","Virion.Hits","Infected","Complete","Productive","Susceptible","Naive")
Cells = matrix(data = 0,
               nrow = Grid.X * Grid.Y,
               ncol = length(Cell.Headers))
colnames(Cells) = Cell.Headers
Cells[,"Cell.ID"] = 1:(Grid.X * Grid.Y)
Cells[,"X"] = rep(1:100, each = 100)  #1,1,1,1,1... 2,2,2,2,2
Cells[,"Y"] = rep(1:100, times = 100) #1,2,3,4,5...1,2,3,4,5
Cells[,"Distance"] = sqrt((Cells[,"X"] - 50) ^ 2 + (Cells[,"Y"] - 50)^ 2)
Cells[,"Living"] = 1
Cells[,"Susceptible"] = 1
Cells[,"Naive"] = 1
Init.Cells = Cells

# Initialize Virions ----
Virion.Headers = c("Virion.ID","Cell.ID","X","Y","Cell.X","Cell.Y","Bound","Diffuse","Movable")

Virions = matrix(data = 0,
                 nrow = 0,
                 ncol = length(Virion.Headers))
colnames(Virions) = Virion.Headers
Init.Virions = Virions

#### Executive Function ----

Sim.Exec = function(Attach_Rate,
                    Detach_Rate,
                    Infect_Rate,
                    Birth_Rate,
                    Death_Rate,
                    Burst_Rate,
                    Exclude_Rate,
                    Revert_Rate,
                    Decay_Rate,
                    Start_Radius,
                    Virus.D,
                    Pp,
                    Spread_Frac,
                    Cells,
                    Virions,
                    Duration,
                    ts,
                    Sim,
                    Varied.Par) {
  set.seed(Sim)
  
  Cell.Summary.Headers = c("Sim","t",
                           "Attach_Rate","Detach_Rate","Birth_Rate","Death_Rate","Burst_Rate","Infect_Rate","Exclude_Rate","Revert_Rate","Decay_Rate","Virus_D","Pp","N0","Spread_Frac", # Parameters
                           "Virions","Free.Virions","Bound.Virions","New.Virions","Segments","Cell.MOI","Inf.MOI","Semi.MOI","Prod.MOI","Productive","Semi.Cells","Dead.Cells","Refractory.Cells","Wasted.Cells","Infected","Susceptible","Naive",
                           "Radius.25","Radius.50","Radius.75","Radius.90","Radius.95","Radius.975","Radius.Max",
                           "Varied.Par")
  
  Cell.Summary = matrix(data = 0,
                        nrow = Duration/ts,
                        ncol = length(Cell.Summary.Headers))
  colnames(Cell.Summary) = Cell.Summary.Headers
  
  Cell.Summary[,"Sim"] = Sim
  Cell.Summary[,"Varied.Par"] = Varied.Par
  Cell.Summary[,"t"] = 1:(Duration / ts) * ts
  Cell.Summary[,"Virus_D"] = Virus.D / ts
  Cell.Summary[,"Pp"] = Pp
  Cell.Summary[,"Spread_Frac"] = Spread_Frac
  Cell.Summary[,"Attach_Rate"] = Attach_Rate / ts
  Cell.Summary[,"Detach_Rate"] = Detach_Rate / ts
  Cell.Summary[,"Birth_Rate"] = Birth_Rate / ts
  Cell.Summary[,"Death_Rate"] = Death_Rate / ts
  Cell.Summary[,"Burst_Rate"] = Burst_Rate / ts
  Cell.Summary[,"Infect_Rate"] = Infect_Rate / ts
  Cell.Summary[,"Exclude_Rate"] = Exclude_Rate / ts
  Cell.Summary[,"Revert_Rate"] = Revert_Rate / ts
  Cell.Summary[,"Decay_Rate"] = Decay_Rate / ts
  Cell.Summary[,"N0"] = (Start_Radius * 2 + 1) ^ 2
  
  Cell.Summary[,"New.Virions"] = 0
  
  # Initialize Infected Cells ----
  Start.X = 50 - Start_Radius
  End.X = 50 + Start_Radius
  Start.Y = 50 - Start_Radius
  End.Y = 50 + Start_Radius
  Start.Index = Cells[which((Cells[,"X"] %in% (Start.X:End.X)) & (Cells[,"Y"] %in% (Start.Y:End.Y))),"Cell.ID"]
  
  Cells[Start.Index,"Segments"] = 8
  Cells[Start.Index,"Virion.Hits"] = 1
  Cells[Start.Index,"Infected"] = 1
  Cells[Start.Index,"Productive"] = 1
  
  Grid.X = max(Cells[,"X"])
  Grid.Y = max(Cells[,"Y"])
  Cell_Width = 30
  Grid_Width_X = Cell_Width * Grid.X
  Grid_Width_Y = Cell_Width * Grid.Y
  
  Inf.Cells = Cells[Cells[,"Infected"] == 1,,drop = FALSE]
  Prod.Cells = Inf.Cells[Inf.Cells[,"Productive"] == 1,,drop = FALSE]
  
  #Translocation probability for locally-dispersing virions
  Trans.Prob = matrix(data = 0,
                      ncol = 3,
                      nrow = 4)
  colnames(Trans.Prob) = c("Prob","X.Adjust","Y.Adjust")
  Trans.Prob[,"Prob"] = rep(1,4)
  Trans.Prob[,"X.Adjust"] = c(-1,0,1,0)
  Trans.Prob[,"Y.Adjust"] = c(0,-1,0,1)
  
  # Iterate Time Series ----
  for (t in 1:(Duration/ts)) { # Start Time
    
    # Exit if no more cells can be infected (Productive cells have died, and no more cells can be infected)
    if ((sum(Cells[,"Susceptible"]) + sum(Cells[,"Productive"])) == 0) {  
      Cell.Summary = Cell.Summary[1:(t - 1),]
      break
    }
    
    # Characterize Populations
    free = which(Virions[,"Bound"] == 0)
    diffuse.free = which(((Virions[,"Diffuse"]) * (1 - Virions[,"Bound"])) == 1)
    adherent = which(Virions[,"Bound"] == 1)
    diffuse.adherent = which(((Virions[,"Diffuse"]) * (Virions[,"Bound"])) == 1)
    
    # Characterize Cells ----
    Cells[,"Infected"] = (Cells[,"Segments"] > 0) * (Cells[,"Living"])
    Cells[,"Productive"] = (Cells[,"Segments"] == 8) * (Cells[,"Living"])
    
    # Obtain Infected, Productive, Semi-Infected, and Dead cell subpopulations
    Inf.Cells = Cells[Cells[,"Infected"] == 1,,drop = FALSE]
    Prod.Cells = Inf.Cells[Inf.Cells[,"Productive"] == 1,,drop = FALSE]
    
    Semi.Cells = Inf.Cells[Inf.Cells[,"Productive"] == 0,,drop = FALSE]
    Dead.Cells = Cells[Cells[,"Living"] == 0,,drop = FALSE]
    
    Refract.Cells = Cells[((1 - Cells[,"Living"]) + Cells[,"Susceptible"]) == 0,,drop = FALSE]
    Wasted.Cells = Refract.Cells[Refract.Cells[,"Segments"] > 0,,drop = FALSE]
    
    # Free virions diffuse ----
    if (length(diffuse.free) > 0) {
      
      Virion.Num = length(diffuse.free)
      
      # Add diffusion distance
      Distance = rnorm(mean = 0, sd = sqrt(2 * Virus.D), n = Virion.Num)
      Theta = runif(min = 0 , max = 2*pi, n = Virion.Num)
      X.Distance = cos(Theta) * Distance
      Y.Distance = sin(Theta) * Distance
      
      Virions[diffuse.free,"X"] = Virions[diffuse.free,"X"] + X.Distance
      Virions[diffuse.free,"Y"] = Virions[diffuse.free,"Y"] + Y.Distance
      
      # Ensure virions stay in grid -- if anything is outside, it comes around from the other side
      
      Virions[diffuse.free,"X"] = Virions[diffuse.free,"X"] %% Grid_Width_X
      Virions[diffuse.free,"Y"] = Virions[diffuse.free,"Y"] %% Grid_Width_Y
      
      # Update Cell locations of virions
      Virions[diffuse.free,"Cell.X"] = floor(Virions[diffuse.free,"X"] / Cell_Width) + 1
      Virions[diffuse.free,"Cell.Y"] = floor(Virions[diffuse.free,"Y"] / Cell_Width) + 1
      #Virions[diffuse.free,"Cell.X"] = pmin(Virions[diffuse.free,"Cell.X"],Grid.X)
      #Virions[diffuse.free,"Cell.Y"] = pmin(Virions[diffuse.free,"Cell.Y"],Grid.Y)
      Virions[diffuse.free,"Cell.ID"] = (Virions[diffuse.free,"Cell.X"] - 1) * Grid.X + Virions[diffuse.free,"Cell.Y"]
      
    }
    
    # Free virions attach ----
    attach.num = round(Attach_Rate * length(diffuse.free), 0)
    attach.num = min(length(diffuse.free),
                     rpois(n = 1, lambda = Attach_Rate * length(diffuse.free)))
    if (attach.num > 0) {
      attach = sample(x = diffuse.free, size = attach.num, replace = FALSE)
      Virions[attach,"Bound"] = 1
    }
    
    # Bound virions detach ----
    detach.num = round(Detach_Rate * length(diffuse.adherent), 0)
    detach.num = min(length(diffuse.adherent),
                     rpois(n = 1, lambda = Detach_Rate * length(diffuse.adherent)))
    if (detach.num > 0) {
      detach = sample(x = diffuse.adherent, size = detach.num, replace = FALSE)
      Virions[detach,"Bound"] = 0
    }
    # Bound virions infect ----
    
    infect.num = round(Infect_Rate * length(adherent), 0)
    infect.num = min(length(adherent),
                     rpois(n = 1, lambda = Infect_Rate * length(adherent)))
    
    if (infect.num > 0) {
      infect = sample(x = adherent, size = infect.num, replace = FALSE)
      
      target = table(Virions[infect,"Cell.ID"])
      
      # Add segments to cell
      Cell.Target = as.numeric(names(target))
      Cell.Hit = as.numeric(target)
      
      Cells[Cell.Target,"Virion.Hits"] = Cells[Cell.Target,"Virion.Hits"] + Cell.Hit * Cells[Cell.Target,"Susceptible"]
      Cells[Cell.Target,"Naive"] = 0
      Delivered = rbinom(n = length(Cell.Hit),
                         prob = 1 - ((1 - Pp) ^ Cell.Hit),
                         size = 8 - Cells[Cell.Target,"Segments"])
      Cells[Cell.Target,"Segments"] = (Cells[Cell.Target,"Segments"] + Delivered) * Cells[Cell.Target,"Susceptible"]
      
      Virions[infect,"Bound"] = NA
      Virions = Virions[-infect,]
    }
    
    # Productive cells produce virions ----
    Burst.Size = rpois(lambda = Burst_Rate, n = nrow(Prod.Cells))
    Burst.Sum = sum(Burst.Size)
    
    if (Burst.Sum > 0) {
      New.Virions = matrix(data = 1,
                           nrow = Burst.Sum,
                           ncol = length(Virion.Headers))
      
      colnames(New.Virions) = Virion.Headers
      New.Virions[,"Bound"] = 1
      
      # Determine Location ----
      rows = rep(1:nrow(Prod.Cells), Burst.Size)
      New.Virions[,"Cell.ID"] = Prod.Cells[rows,"Cell.ID"]
      New.Virions[,"Cell.X"] = Prod.Cells[rows,"X"]
      New.Virions[,"Cell.Y"] = Prod.Cells[rows,"Y"]
      New.Virions[,"X"] = (New.Virions[,"Cell.X"] - 1) * Cell_Width + runif(min = 0, max = Cell_Width, n = Burst.Sum)
      New.Virions[,"Y"] = (New.Virions[,"Cell.Y"] - 1) * Cell_Width + runif(min = 0, max = Cell_Width, n = Burst.Sum)
      
      # Local Spread ----
      
      local.num = round(Spread_Frac * Burst.Sum,0)
      
      if (local.num > 0) {
        local = sample(x = 1:Burst.Sum, size = local.num, replace = FALSE)
        New.Virions[local,"Bound"] = 1
        New.Virions[local,"Diffuse"] = 0
        
        # Sample 1 - 4 for each to determine translocation direction of each locally spreading virion
        adjust = sample(x = 1:nrow(Trans.Prob),
                        size = local.num,
                        prob = Trans.Prob[,"Prob"],
                        replace = TRUE)
        
        # Move to neighboring cell
        New.Virions[local,"Cell.X"] = New.Virions[local,"Cell.X"] + Trans.Prob[adjust,"X.Adjust"]
        New.Virions[local,"Cell.Y"] = New.Virions[local,"Cell.Y"] + Trans.Prob[adjust,"Y.Adjust"]
        
        # Ensure that virions are inside the grid
        New.Virions[local,"Cell.X"] = pmax(1, New.Virions[local,"Cell.X"])
        New.Virions[local,"Cell.Y"] = pmax(1, New.Virions[local,"Cell.Y"])
        
        New.Virions[local,"Cell.X"] = pmin(Grid.X, New.Virions[local,"Cell.X"])
        New.Virions[local,"Cell.Y"] = pmin(Grid.Y, New.Virions[local,"Cell.Y"])
        
        # Update Cell.ID
        New.Virions[local,"Cell.ID"] = (New.Virions[local,"Cell.X"] - 1) * 100 + New.Virions[local,"Cell.Y"]
        
      }
      Virions = rbind(Virions,New.Virions)
    }
    
    # Productive cells become refractory to super-infection ----
    prod = which(Cells[,"Productive"] == 1)
    
    exclude.num = round(Exclude_Rate * length(prod), 0)
    exclude.num = min(length(prod),
                      rpois(n = 1, lambda = Exclude_Rate * length(prod)))
    
    if (exclude.num > 0) {
      exclude = sample(x = prod, size = exclude.num, replace = FALSE)
      Cells[exclude,"Susceptible"] = 0
    }
    # Productive cells die ----
    death.num = round(Death_Rate * length(prod), 0)
    death.num = min(length(prod),
                    rpois(n = 1, lambda = Death_Rate * length(prod)))
    if (death.num > 0) {
      death = sample(x = prod, size = death.num, replace = FALSE)
      Cells[death,"Living"] = 0
      Cells[death,"Productive"] = 0
      Cells[death,"Segments"] = 0
      Cells[death,"Susceptible"] = 0
      Cells[death,"Infected"] = 0
    }
    
    #Summarize system ----
    Cell.Summary[t,"Productive"] = sum(Cells[,"Productive"])
    Cell.Summary[t,"Dead.Cells"] = nrow(Cells) - sum(Cells[,"Living"])
    Cell.Summary[t,"Semi.Cells"] = sum(Cells[,"Infected"]) - sum(Cells[,"Productive"])
    Cell.Summary[t,"Refractory.Cells"] = nrow(Refract.Cells)
    Cell.Summary[t,"Wasted.Cells"] = nrow(Wasted.Cells)
    Cell.Summary[t,"Infected"] = sum(Cells[,"Infected"])
    Cell.Summary[t,"Susceptible"] = sum(Cells[,"Susceptible"])
    Cell.Summary[t,"Naive"] = sum(Cells[,"Naive"])
    # 
    Cell.Summary[t,"Virions"] = nrow(Virions)
    Cell.Summary[t,"Bound.Virions"] = sum(Virions[,"Bound"])
    Cell.Summary[t,"Free.Virions"] = nrow(Virions) - sum(Virions[,"Bound"])
    Cell.Summary[t,"New.Virions"] = Burst.Sum
    Cell.Summary[t,"Segments"] = mean(Inf.Cells[,"Segments"])
    Cell.Summary[t,"Inf.MOI"] = mean(Inf.Cells[,"Virion.Hits"])
    Cell.Summary[t,"Radius.25"] = quantile(Prod.Cells[,"Distance"],0.25)
    Cell.Summary[t,"Radius.50"] = quantile(Prod.Cells[,"Distance"],0.50)
    Cell.Summary[t,"Radius.75"] = quantile(Prod.Cells[,"Distance"],0.75)
    Cell.Summary[t,"Radius.90"] = quantile(Prod.Cells[,"Distance"],0.90)
    Cell.Summary[t,"Radius.95"] = quantile(Prod.Cells[,"Distance"],0.95)
    Cell.Summary[t,"Radius.975"] = quantile(Prod.Cells[,"Distance"],0.975)
    Cell.Summary[t,"Radius.Max"] = max(Prod.Cells[,"Distance"])
  } # Next time step (t)
  return(Cell.Summary)
}

# Run Simulation ----

registerDoMC(cores = 2) # Nate Macbook

Sim.Start = Sys.time()
Cell.Summary = foreach(pp = 1:Pp_Num, .combine = rbind) %:%
  foreach(d = 1:Virus.D_Num, .combine = rbind) %:%
  foreach(spread = 1:Spread_Frac_Num, .combine = rbind) %:%
  foreach(sim = 1:Sim_Num, .combine = rbind) %dopar%
  try(Sim.Exec(Attach_Rate = Attach_Rate[1],
           Detach_Rate = Detach_Rate[1],
           Death_Rate = Death_Rate[1],
           Birth_Rate = Birth_Rate[1],
           Infect_Rate = Infect_Rate[1],
           Burst_Rate = Burst_Rate[1],
           Exclude_Rate = Exclude_Rate[1],
           Revert_Rate = Revert_Rate[1],
           Decay_Rate = Decay_Rate[1],
           Start_Radius = Start_Radius[1],
           Virus.D = Virus.D[d],
           Pp = Pp[pp],
           Spread_Frac = Spread_Frac[spread],
           Cells = Init.Cells,
           Virions = Init.Virions,
           Duration = Duration,
           ts = ts,
           Sim = sim,
           Varied.Par = "All"))

# QC ----
Sim.End = Sys.time()

Sim.Time = Sim.End - Sim.Start
# Sim.Time = 1.3 hours on a 2-core Macbook Pro (Late 2013)

#### Export Data ----
NTJ21.Results = Cell.Summary
write.csv(NTJ21.Results,file=file.path(Data.Path,"Fig4_Data_SHORT.csv"),row.names=FALSE)
NTJ21.Results = read.csv(file=file.path(Data.Path,"Fig4_Data_SHORT.csv"))
x = subset(NTJ21.Results, as.numeric(Sim) <= 10)
write.csv(x,file=file.path(Data.Path,"Fig4_Data_SHORT.csv"),row.names=FALSE)

# Import and Analyze Data ----

require(tidyverse)
require(dplyr)

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

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