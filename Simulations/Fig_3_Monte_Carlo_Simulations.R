# Simulation 3

#### Documentation ----

# Objective: Test analytic solution through simulation, and estimate how the infectious dose changes with Pp (using MOI)
# Approach: Iterating over MOI range of only virus A, infect a population of cells,
#           and record how many become productively infected.
# Output: 
#     % Productive vs. MOI -- how many is enough/what is the point of diminishing returns?
#     Segments/Cell vs. MOI -- slower accumulation of segments at lower Pp
#     % Populations infected vs. MOI -- What is the infectious dose (MOI) at the population level?
# Requirements:
# 10 minutes/1 Sim
# 1 Day = 150 Simulations
# 1 Week = 1000 Simulations

### Define Simulation Function

Sim.3.Exec = function(s,
                      m,
                      p,
                      Cell.Num,
                      Summary.Names) {
  set.seed(s)
  # Define Summary for Export
  Export.Summary=matrix(data = 0,
                        ncol = length(Summary.Names),
                        nrow = 1)
  colnames(Export.Summary) = Summary.Names
  
  # Simulate Infection
  Infected.Cells = Single_Infection(MOI = m,
                                    Cell.Num = Cell.Num,
                                    Yewdell = rep(p,8))
  
  Export.Summary[1,"Sim"] = s
  Export.Summary[1,"Pp"] = p
  Export.Summary[1,"MOI"] = m
  Export.Summary[1, "Segments"] = mean(Infected.Cells[,"Segments"])
  Export.Summary[1, "Productive"] = mean(Infected.Cells[,"Productive"]) * 100
  Export.Summary[1, "Infected_Cells"] = mean(Infected.Cells[,"Infected"]) * 100
  Export.Summary[1, "Semi"] = (mean(Infected.Cells[,"Infected"]) - mean(Infected.Cells[,"Productive"])) * 100
  
  #if (p==1) write.table(Export.Summary,file="/home/rstudio/Dropbox/Modeling/Simulation_Data/NTJ_Sim_3_Progress.txt",row.names=F)
  
  return(Export.Summary)
}

#### Housekeeping
require(foreach)
require(doMC)
require(Rmisc)

#Proj.Home="/Users/Nate/Box Sync/Lowen_Lab/Dropbox/Modeling"
Proj.Home="/home/rstudio/Dropbox/Modeling"
Function.File=file.path(Proj.Home,"Functions","Flu_Functions.R")
Data.Path=file.path(Proj.Home,"Simulation_Data")
Data.Path = file.path(Proj.Home,"Amazon_Output")
source(Function.File)

##### Parameters ----

# Segment Order --- PB2, PB1, PA, HA, NP, NA, M, NS



Sim.Num = 1000

Pp=matrix(data = c(rep(1:10,each=8)/10,rep(0.58,8)),
          byrow = TRUE,ncol=8)

Pp.Num=nrow(Pp)

MOI=c(1e-5,2e-5,3e-5,4e-5,5e-5,6e-5,7e-5,8e-5,9e-5,
      1e-4,2e-4,3e-4,4e-4,5e-4,6e-4,7e-4,8e-4,9e-4,
      1e-3,2.5e-3,5e-3,7.5e-3,
      1e-2,2.5e-2,5e-2,7.5e-2,
      1e-1,2.5e-1,5e-1,7.5e-1,
      1,2.5,5,10,
      25,50,100)
MOI.Num=length(MOI)

Flu.Segments=c("PB2","PB1","PA","HA","NP","NA","M","NS")

Cell.Num = 100000 # Usually inoculate 1,000,000 cells, but similar results obtained by simulating 100,000 cells

Cell.Summary.Names=c(
  # Model Parameters
  "Sim","MOI","Pp",
  # Outcomes
  "Segments","Productive","Infected","Semi","Infected_Cells"
)

# Loop ----

Sim.Start = Sys.time()

#registerDoMC(cores = 2) # Nate Macbook
registerDoMC(cores = 64) # Amazon Server

Cell.Summary = foreach(sim = 1:Sim.Num, .combine=rbind) %:%
               foreach(pp = 1:Pp.Num, .combine=rbind) %:%
               foreach(moi = 1:MOI.Num, .combine=rbind) %dopar%
                  Sim.3.Exec(s = sim,
                             p = Pp[pp,1],
                             m = MOI[moi],
                             Cell.Num = Cell.Num,
                             Summary.Names = Cell.Summary.Names)

Cell.Summary[, "Infected"] = Cell.Summary[,"Productive"] > 0
Sim.End = Sys.time()

#### Export Data ----
write.csv(Cell.Summary,file=file.path(Data.Path,"NTJ_Sim3_Summary_Data.csv"),row.names=FALSE)

#### Log
write.table(c(paste("Start ",Sim.Start),paste("End ",Sim.End),paste("Run Time ",difftime(Sim.End,Sim.Start,units="hours")," hours")),
            file=file.path(Data.Path,"NTJ_Sim3_Log.txt"))
