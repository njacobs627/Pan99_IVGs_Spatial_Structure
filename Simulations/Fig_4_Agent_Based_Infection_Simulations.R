#### Load Packages ----

require(foreach)
require(doMC)
require(Rmisc)
require(ggplot2)
set.seed(666)

### Set File Paths ----

Proj.Home="/Users/Nate/Dropbox/Modeling" # Nate Macbook
Proj.Home="/home/rstudio/Dropbox/Modeling" # Rstudio Amazon
Data.Path = file.path(Proj.Home,"Amazon_Output")

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

# Local Simulation Parameters
Pp = c(0.575,1)
Sim_Num = 2
Spread_Frac = c(0.5)
Virus.D = 5.825 * c(0.001, 0.01, 0.1, 1, 10, 100)
Virus.D = 5.825 * c(0.001, 0.1, 1, 100)
Virus.D = 5.825 * 60 * 60
Spread_Frac = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)
ts = 1/20

Duration = 1

#Amazon Simulation Parameters
Sim_Num = 10
Pp = c(0.575, 1.0)
Spread_Frac = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)
Virus.D = 5.825 * 60 * 60
Virus.D = Virus.D * c(rep(c(1, 0.1, 10, 100, 1e3),5) * rep(c(0.2, 0.4, 0.6, 0.8,1),each = 5),1e5, 1e8, 1e11)
Duration = 96
ts = 1/20
sort(Virus.D / 3600)

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
      
      # Virions[diffuse.free,"X"] = abs(Virions[diffuse.free,"X"])
      # Virions[diffuse.free,"X"] = Virions[diffuse.free,"X"] %% (2 * Grid_Width_X)
      # Virions[diffuse.free,"X"] = Virions[diffuse.free,"X"] - 2 * pmax(0,
      #                                                                  Virions[diffuse.free,"X"] - (Grid_Width_X))
      # 
      # Virions[diffuse.free,"Y"] = abs(Virions[diffuse.free,"Y"])
      # Virions[diffuse.free,"Y"] = Virions[diffuse.free,"Y"] %% (2 * Grid_Width_Y)
      # Virions[diffuse.free,"Y"] = Virions[diffuse.free,"Y"] - 2 * pmax(0,
      #                                                                  Virions[diffuse.free,"Y"] - (Grid_Width_Y))
      
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
  write.table(Cell.Summary[1,],file="/home/rstudio/Dropbox/Modeling/Amazon_Output/NTJ_Sim_Progress.txt",row.names=F)
  return(Cell.Summary)
}

# Run Simulation ----

Sim.Start = Sys.time()

registerDoMC(cores = 2) # Nate Macbook
registerDoMC(cores = 64) # Amazon Server

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
Sim.End - Sim.Start

#### Export Data ----
NTJ21.Results = Cell.Summary
write.csv(NTJ21.Results,file=file.path(Data.Path,"NTJ_Sim21_Data_Amazon.csv"),row.names=FALSE)
NTJ21.Results = read.csv(file=file.path(Data.Path,"NTJ_Sim21_Data_Amazon.csv"))
write.table(c(paste("Start ",Sim.Start),paste("End ",Sim.End),paste("Run Time ",difftime(Sim.End,Sim.Start,units="hours")," hours")),
            file=file.path(Data.Path,"NTJ_Sim21_Log.txt"))
x = subset(NTJ21.Results, as.numeric(Sim) <= 10)
write.csv(x,file=file.path(Data.Path,"NTJ_Sim21_Data_Amazon.csv"),row.names=FALSE)
NTJ21.Results = read.csv(file=file.path(Data.Path,"NTJ_Sim21_Data_Amazon.csv"))
# Simulation causing the error has the following parameters
# Sim = 10, 
# Pp = 0.575, 
# Spread_Frac = 0.1, 
# Virus_D = 2.097e+12 (5.825 * 60 * 60 * 1e8)
table(NTJ21.Results$Pp)
x1 = subset(NTJ21.Results, t == 1)
table(x1$Sim)
