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

#### Housekeeping
require(foreach)
require(doMC)
require(Rmisc)
require(viridis)
### Define Infection Function ----
Single_Infection=function(Cell.Num,MOI,Pp) {
  Cell.Names = c(
    "MOI","Virions","X","Y",
    "Cell","Infected","Productive","Poly","Surface.HA","New.Segments","Segments","Convert",
    "PB2","PB1","PA","HA","NP","NA","M","NS",
    "Copy.PB2","Copy.PB1","Copy.PA","Copy.HA","Copy.NP","Copy.NA","Copy.M","Copy.NS"
  )
  Cells = matrix(nrow = Cell.Num,
                 ncol = length(Cell.Names),
                 data = 0)
  
  Cells=as.data.frame(Cells)
  
  colnames(Cells) = Cell.Names
  
  Cells[, "Cell"] = seq(1, Cell.Num)
  Cells[, "MOI"] = MOI
  
  withA = sort(ceiling(runif(n = MOI * Cell.Num) * Cell.Num))
  
  x = table(withA)
  x_index = as.numeric(names(x))
  x_values = as.numeric(x)
  
  Cells[x_index, "Infected"] = x_values > 0
  Cells[x_index, "Virions"] = x_values
  
  PB2.Index = which((colnames(Cells)) == "PB2")
  Copy.PB2.Index = which((colnames(Cells)) == "Copy.PB2")
  
  for (segment in 1:8) {
    x = table(sort(sample(
      withA,
      size = length(withA) * Pp[segment],
      replace = FALSE
    )))
    
    x_index = as.numeric(names(x))
    x_values = as.numeric(x)
    
    Cells[x_index, PB2.Index + segment - 1] = x_values > 0     #Boolean Presence/Absence of Segment
    Cells[x_index, Copy.PB2.Index + segment - 1] = x_values                  #Copy Number of each segment
  }
  
  Cells[, "Segments"] = rowSums(Cells[,(PB2.Index:(PB2.Index + 8 - 1))])
  Cells[, "Poly"] = Cells[, "PB2"] * Cells[, "PB1"] * Cells[, "PA"] * Cells[, "NP"]
  Cells[, "Surface.HA"] = Cells[, "Poly"] * Cells[, "HA"]
  Cells[, "Productive"] = Cells[, "Surface.HA"] * Cells[, "M"] * Cells[, "NS"] * Cells[, "NA"]
  
  Cells
}
### Define Simulation Function ----

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
                                    Pp = rep(p,8))
  
  Export.Summary[1,"Sim"] = s
  Export.Summary[1,"Pp"] = p
  Export.Summary[1,"MOI"] = m
  Export.Summary[1, "Segments"] = mean(Infected.Cells[,"Segments"])
  Export.Summary[1, "Productive"] = mean(Infected.Cells[,"Productive"]) * 100
  Export.Summary[1, "Infected_Cells"] = mean(Infected.Cells[,"Infected"]) * 100
  Export.Summary[1, "Semi"] = (mean(Infected.Cells[,"Infected"]) - mean(Infected.Cells[,"Productive"])) * 100
  
  return(Export.Summary)
}


### Set File Paths ----

Proj.Home = "/Users/Nate/Box Sync/Lowen_Lab/Data_Notebook/Documentation" #For Nate's Macbook Pro

Data.Path = file.path(Proj.Home,"Simulation_Test_Results")

### Define Simulation Parameters ----

Sim.Num = 40

Pp=matrix(data = c(rep(1:10,each=8)/10,rep(0.58,8)),
          byrow = TRUE,ncol=8)

Pp.Num=nrow(Pp)

MOI=c(1e-5,3.16e-5,
      1e-4,3.16e-4,
      1e-3,3.16e-3,
      1e-2,3.16e-2,
      1e-1,3.16e-1,
      1,3.16,
      10)

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

registerDoMC(cores = 2) # Nate Macbook

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
Sim.Time = Sim.End - Sim.Start

Sim.Time
# Sim.Time = 2.5 minutes per sim

#### Export Data ----
write.csv(Cell.Summary,file=file.path(Data.Path,"Fig3_Data_SHORT.csv"),row.names=FALSE)

Sim3 = read.csv(file.path(Data.Path,"Fig3_Data_SHORT.csv"),sep=",",header=TRUE)

Sim3 = subset(Sim3,MOI>0)
Sim3$Infected = as.numeric(Sim3$Infected)
#Segments vs. MOI
Sum.Sim3 = Sim3 %>%
  mutate(MOI = (MOI)) %>%
  group_by(MOI,Pp) %>%
  dplyr::summarise(Segments = mean(Segments),
            Productive = mean(Productive),
            Infected = mean(Infected) * 100,
            Semi = mean(Semi),
            Infected_Cells = mean(Infected_Cells),
            Semi_Prop = mean(Semi) / mean(Infected_Cells))
# Populations Infected ----

ggplot() + 
  geom_line(data=Sum.Sim3 %>% filter(Pp != 0.58),
            aes(x = log10(MOI),
                y = Infected,
                color = Pp,
                group = Pp),
            lwd = 1) +
  scale_color_viridis(guide = FALSE) +
  labs(x = TeX("\\textbf{MOI (log_1_0 Virions / Cell)}"),
       y = TeX("\\textbf{% Populations Infected}"),
       color = NULL) +
  theme(text=element_text(size=12,face="bold"),
        strip.text.x=element_text(size=rel(1.5),margin=margin(0,0,3,0)),
        strip.text.y=element_text(size=rel(1.5),margin=margin(0,0,0,0),angle=0),
        strip.background = element_blank(),
        plot.title = element_text(size=rel(1.5)),
        axis.text.x=element_text(angle=0,vjust=0,size=rel(1.5),color = "black"),
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

ggsave('Figures/3A_Populations_Infected.pdf',
       width = 4,
       height = 4,
       unit = "in")

# Population ID50 ----

Pop.ID50 = matrix(nrow = 10,
                  ncol = 2,
                  data = 0)

colnames(Pop.ID50) = c("Pp","ID50")
Pop.ID50[,"Pp"] = 1:10 / 10

for (p in 1:9) {
  model = drm(data = Sum.Sim3 %>% filter(Pp == (p / 10), Infected < 100),
              formula = Infected ~ MOI,
              fct = LL.5()) %>%
    ED(50)
  Pop.ID50[p,"ID50"] = model[,1]
}
Pop.ID50[10,2] = 1e-5
Pop.ID50 = data.frame(Pop.ID50)
Pop.ID50$ID50 = log10(Pop.ID50$ID50)

ggplot() +
  geom_line(data = Pop.ID50[1:9,],
            aes(x = Pp,
                y = ID50),
            lwd = 1) +
  scale_x_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1)) +
  labs(x = TeX("\\textbf{P_P}"),
       y = TeX("\\textbf{Population ID_5_0 (log_1_0 Virions / Cell)}")) +
  theme(text=element_text(size=12,face="bold"),
        strip.text.x=element_text(size=rel(1.5),margin=margin(0,0,3,0)),
        strip.text.y=element_text(size=rel(1.5),margin=margin(0,0,0,0),angle=0),
        strip.background = element_blank(),
        plot.title = element_text(size=rel(1.5)),
        axis.text.x=element_text(angle=0,vjust=0,size=rel(1.5),color = "black"),
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


ggsave('Figures/3B_Population_ID50.pdf',
       width = 4,
       height = 4,
       unit = "in")

ggplot() +
  geom_line(data = Sum.Sim3 %>% filter(Pp %in% c(0.58)),
            aes(x = log10(MOI),
                y = Semi_Prop * 100,
                group = Pp), lwd = 1.2, color = "black") +
  scale_color_viridis(guide = FALSE, begin = 0.1, end = 0.9) +
  labs(x = TeX("\\textbf{MOI (log_1_0 Virions / Cell)}"),
       y = TeX("\\textbf{% Infected Cells with IVGs}")) +
  scale_x_continuous(limits = c(-2,2)) +
  theme(text=element_text(size=12,face="bold"),
        strip.text.x=element_text(size=rel(1.5),margin=margin(0,0,3,0)),
        strip.text.y=element_text(size=rel(1.5),margin=margin(0,0,0,0),angle=0),
        strip.background = element_blank(),
        plot.title = element_text(size=rel(1.5)),
        axis.text.x=element_text(angle=0,vjust=0,size=rel(1.5),color = "black"),
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

ggsave('Figures/6A_Proportion_Infected_Cells_with_IVGs.pdf',
       width = 4,
       height = 4,
       unit = "in")
