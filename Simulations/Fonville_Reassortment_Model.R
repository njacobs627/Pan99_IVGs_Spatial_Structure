require(ggplot2)
require(knitr)
require(Rmisc)
require(ggthemes)
require(scales)
require(Rmisc)
require(stringr)
require(gridExtra)
require(viridis)
require(dplyr)
require(latex2exp)
theme_set(theme_grey())

Proj.Home = "/Users/Nate/Box Sync/Lowen_Lab/Data_Notebook/Documentation" #For Nate's Macbook Pro
Data.Path = file.path(Proj.Home,"Data")
CellInfection=function(Cell.Num,MOI,Cell.Type,Pp) {
  
  set.seed(627) #
  
  Cell.Names = c(
    "MOI","Cell.Type",
    "Cell","Infected","Co.Infected","Productive","Productive.Co.Infected","Poly","Surface.HA","Inf.A","Inf.B",
    "pParent.A","pParent.B","pReassort","pHN.Reassort",
    "PB2","PB1","PA","HA","NP","NA","M","NS",
    "A.PB2","A.PB1","A.PA","A.HA","A.NP","A.NA","A.M","A.NS",
    "B.PB2","B.PB1","B.PA","B.HA","B.NP","B.NA","B.M","B.NS",
    "Copy.A.PB2","Copy.A.PB1","Copy.A.PA","Copy.A.HA","Copy.A.NP","Copy.A.NA","Copy.A.M","Copy.A.NS",
    "Copy.B.PB2","Copy.B.PB1","Copy.B.PA","Copy.B.HA","Copy.B.NP","Copy.B.NA","Copy.B.M","Copy.B.NS",
    "Copy.PB2","Copy.PB1","Copy.PA","Copy.HA","Copy.NP","Copy.NA","Copy.M","Copy.NS",
    "Six","HA.Neg","NS.Neg"
  )
  Cells = matrix(nrow = Cell.Num,
                 ncol = length(Cell.Names),
                 data = 0)
  
  colnames(Cells) = Cell.Names
  
  Cells[, "Cell"] = seq(1, Cell.Num)
  Cells[, "MOI"] = MOI
  
  withA = sort(ceiling(runif(n = MOI * Cell.Num) * Cell.Num))
  withB = sort(ceiling(runif(n = MOI * Cell.Num) * Cell.Num))
  
  x = table(withA)
  x_index = as.numeric(names(x))
  x_values = as.numeric(x)
  Cells[x_index, "Inf.A"] = x_values
  
  x = table(withB)
  x_index = as.numeric(names(x))
  x_values = as.numeric(x)
  Cells[x_index, "Inf.B"] = x_values
  
  Cells[, "Co.Infected"] = (Cells[, "Inf.A"] * Cells[, "Inf.B"] > 0)
  Cells[, "Infected"] = (Cells[, "Inf.A"] + Cells[, "Inf.B"] > 0)
  Cells[, "pParent.A"] = 1
  Cells[, "pParent.B"] = 1
  Cells[, "pHN.Reassort"] = 1
  
  A.PB2.Index = which((colnames(Cells)) == "A.PB2")
  PB2.Index = which((colnames(Cells)) == "PB2")
  
  for (segment in 1:8) {
    x = table(sort(sample(
      withA,
      size = length(withA) * Pp[segment],
      replace = FALSE
    )))
    x_index = as.numeric(names(x))
    x_values = as.numeric(x)
    
    Cells[x_index, A.PB2.Index + segment - 1] = x_values > 0     #Boolean Presence/Absence of Segment, A strain
    Cells[x_index, A.PB2.Index + segment - 1 + 16] = x_values     #Copy Number of each segment, A strain
    
    x = table(sort(sample(
      withB,
      size = length(withB) * Pp[segment],
      replace = FALSE
    )))
    x_index = as.numeric(names(x))
    x_values = as.numeric(x)
    
    Cells[x_index, A.PB2.Index + segment - 1 + 8] = x_values > 0     #Boolean Presence/Absence of Segment, B strain
    Cells[x_index, A.PB2.Index + segment - 1 + 16 + 8] = x_values     #Copy Number of each segment, B strain
    
    Cells[, A.PB2.Index + segment - 1 + 16 + 8 + 8] = Cells[, A.PB2.Index + segment - 1 + 16] + Cells[, A.PB2.Index + segment - 1 + 16 + 8] #Copy number of each segment, irrespective of source
    Cells[, PB2.Index + segment - 1] = (Cells[, A.PB2.Index + segment - 1] +
                                          Cells[, A.PB2.Index + segment - 1 + 8]) > 0  #Presence/Absence of Gene, irrespective of source
    
    Cells[, "pParent.A"] = Cells[, "pParent.A"] * (Cells[, A.PB2.Index + segment -
                                                           1 + 16] / Cells[, A.PB2.Index + segment - 1 + 16 + 8 + 8])
    Cells[, "pParent.B"] = Cells[, "pParent.B"] * (Cells[, A.PB2.Index + segment -
                                                           1 + 16 + 8] / Cells[, A.PB2.Index + segment - 1 + 16 + 8 + 8])
  }
  
  Cells[, "Poly"] = Cells[, "PB2"] * Cells[, "PB1"] * Cells[, "PA"] * Cells[, "NP"]
  Cells[, "Surface.HA"] = Cells[, "Poly"] * Cells[, "HA"]
  Cells[, "Productive"] = Cells[, "Surface.HA"] * Cells[, "M"] * Cells[, "NS"] *
    Cells[, "NA"]
  Cells[, "Productive.Co.Infected"] = Cells[, "Productive"] * Cells[, "Co.Infected"]
  Cells[, "pReassort"] = 1 -
    Cells[, "pParent.A"] -
    Cells[, "pParent.B"]
  Cells[, "pHN.Reassort"] = 1 -
    (Cells[, "Copy.A.HA"] / Cells[, "Copy.HA"]) * (Cells[, "Copy.A.NA"] / Cells[, "Copy.NA"]) -
    (Cells[, "Copy.B.HA"] / Cells[, "Copy.HA"]) * (Cells[, "Copy.B.NA"] / Cells[, "Copy.NA"])
  
  Cells
}

### Parameters ----

Exp.Pp.2 = read.csv(file.path(Data.Path,"Pp_Summary.csv"),na.strings=c("","-"))  %>% filter(Delay == 0)

PB2.Index = which(colnames(Exp.Pp.2) == "PB2")

Pp = Exp.Pp.2[,PB2.Index:(PB2.Index + 7)]

Rep.Num=nrow(Pp)

rownames(Pp)=as.character(1:Rep.Num)

MOI=matrix(data=rep(c(.001,.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1,1.5,2,5,10),Rep.Num), #Range of MOIs
          byrow=TRUE,nrow=Rep.Num)

rownames(MOI)=rownames(Pp)

Flu.Segments=c("PB2","PB1","PA","HA","NP","NA","M","NS")

MOI.Num=ncol(MOI)

Cell.Num = 100000 # Usually inoculate 1,000,000 cells, but similar results obtained by simulating 100,000 cells
Burst.Size = 986
Rep.Num=nrow(Pp)

### Initialize Summary Data Frame ----

Cell.Summary.Names=c(
  # Model Parameters
  "Rep","MOI",
  
  # Outcomes
  "Productive", "Productive.Co.Infected","Expressing.HA",
  "Virion.Number","Parent.A.Percent","Parent.B.Percent","Reassortant.Percent","Reassortant.HN.Percent",
  "Parent.A.Number","Parent.B.Number","Reassortant.Number","Reassortant.HN.Number"
)

Cell.Summary=matrix(nrow=(MOI.Num*Rep.Num),
                    ncol=length(Cell.Summary.Names),
                    data=0)

colnames(Cell.Summary)=Cell.Summary.Names
Cell.Summary=as.data.frame(Cell.Summary)
Cell.Summary$Rep=rownames(Pp)

Complete.Summary = Cell.Summary[1:MOI.Num,]

Summary.Index = 1

Pp = as.matrix(Pp)


# Simulation ----
Sim.Start = Sys.time()

for (moi in 1:MOI.Num) {
  for (cell in 1:Rep.Num) {
    
    Infected.Cells = CellInfection(MOI = MOI[cell,moi],
                                   Cell.Num = Cell.Num,
                                   Pp = Pp[cell,])
    
    Cell.Summary[Summary.Index,"Productive"] =  mean(Infected.Cells[,"Productive"]) * 100
    Cell.Summary[Summary.Index,"Productive.Co.Infected"] = mean(Infected.Cells[,"Productive.Co.Infected"]) * 100
    Cell.Summary[Summary.Index,"Expressing.HA"] = mean(Infected.Cells[,"Surface.HA"]) * 100
    Cell.Summary[Summary.Index,"Rep"] = rownames(Pp)[cell]
    Cell.Summary[Summary.Index,"MOI"] = MOI[cell,moi]
    Cell.Summary[Summary.Index,"Virion.Number"] = Burst.Size * Cell.Summary[Summary.Index,"Productive"]
    
    Infected.Cells=Infected.Cells[Infected.Cells[,"Productive"]==1,,drop=FALSE] 
    Cell.Summary[Summary.Index,"Parent.A.Percent"] = mean(Infected.Cells[,"pParent.A"]) * 100
    Cell.Summary[Summary.Index,"Parent.B.Percent"] = mean(Infected.Cells[,"pParent.B"]) * 100
    Cell.Summary[Summary.Index,"Reassortant.Percent"] = mean(Infected.Cells[,"pReassort"]) * 100
    Cell.Summary[Summary.Index,"Reassortant.HN.Percent"] = mean(Infected.Cells[,"pHN.Reassort"]) * 100
    
    Cell.Summary[Summary.Index,"Parent.A.Number"] = Cell.Summary[Summary.Index,"Parent.A.Percent"] * Cell.Summary[Summary.Index,"Virion.Number"]
    Cell.Summary[Summary.Index,"Parent.B.Number"] = Cell.Summary[Summary.Index,"Parent.B.Percent"] * Cell.Summary[Summary.Index,"Virion.Number"]
    Cell.Summary[Summary.Index,"Reassortant.Number"] = Cell.Summary[Summary.Index,"Reassortant.Percent"] * Cell.Summary[Summary.Index,"Virion.Number"]
    Cell.Summary[Summary.Index,"Reassortant.HN.Number"] = Cell.Summary[Summary.Index,"Reassortant.HN.Percent"] * Cell.Summary[Summary.Index,"Virion.Number"]
    
    Summary.Index=Summary.Index+1
  }
  
  Infected.Cells = CellInfection(MOI = MOI[1,moi],
                                 Cell.Num = Cell.Num,
                                 Pp = rep(1,8))
  
  Complete.Summary[moi,"Productive"] =  mean(Infected.Cells[,"Productive"]) * 100
  Complete.Summary[moi,"Productive.Co.Infected"] = mean(Infected.Cells[,"Productive.Co.Infected"]) * 100
  Complete.Summary[moi,"Expressing.HA"] = mean(Infected.Cells[,"Surface.HA"]) * 100
  Complete.Summary[moi,"Rep"] = 1
  Complete.Summary[moi,"MOI"] = MOI[cell,moi]
  Complete.Summary[moi,"Virion.Number"] = Burst.Size * Cell.Summary[Summary.Index,"Productive"]
  
  Infected.Cells=Infected.Cells[Infected.Cells[,"Productive"]==1,,drop=FALSE] 
  Complete.Summary[moi,"Parent.A.Percent"] = mean(Infected.Cells[,"pParent.A"]) * 100
  Complete.Summary[moi,"Parent.B.Percent"] = mean(Infected.Cells[,"pParent.B"]) * 100
  Complete.Summary[moi,"Reassortant.Percent"] = mean(Infected.Cells[,"pReassort"]) * 100
  Complete.Summary[moi,"Reassortant.HN.Percent"] = mean(Infected.Cells[,"pHN.Reassort"]) * 100
  
  Complete.Summary[moi,"Parent.A.Number"] = Cell.Summary[Summary.Index,"Parent.A.Percent"] * Cell.Summary[Summary.Index,"Virion.Number"]
  Complete.Summary[moi,"Parent.B.Number"] = Cell.Summary[Summary.Index,"Parent.B.Percent"] * Cell.Summary[Summary.Index,"Virion.Number"]
  Complete.Summary[moi,"Reassortant.Number"] = Cell.Summary[Summary.Index,"Reassortant.Percent"] * Cell.Summary[Summary.Index,"Virion.Number"]
  Complete.Summary[moi,"Reassortant.HN.Number"] = Cell.Summary[Summary.Index,"Reassortant.HN.Percent"] * Cell.Summary[Summary.Index,"Virion.Number"]
}

Sim.End = Sys.time()
Sim.Time = Sim.End - Sim.Start


# Analysis ----
X="Expressing.HA"
Y="Reassortant.Percent"

FM.Data=read.csv(file = file.path(Data.Path,"Fonville_Marshall_Reassortment_Data.csv"),header=TRUE,sep=",")

Cell.Summary$Rep = as.numeric(Cell.Summary$Rep)

ggplot() +
  geom_line(data=Cell.Summary,aes(y=Reassortant.Percent,x=Expressing.HA,group=as.factor(Rep),color=as.factor(Rep)),lwd=1) +
  geom_line(data = Complete.Summary, aes(x = Expressing.HA, y = Reassortant.Percent), lty = 2, lwd = 1) +
  geom_point(data=FM.Data,aes_string(x=X,y=Y),size=2,col="black") +
  annotate("text",x = 75, y = 15,size = 6,label = TeX("\\textbf{$\\P_P$ = 1}")) +
  labs(x=expression(bold("% Cells HA"^"+")),
       y="% Reassortment",
       color="Rep") +
  guides(color = FALSE) +
  scale_x_continuous(limits=c(0,100)) +
  scale_y_continuous(limits=c(0,100)) +
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

ggsave(file.path(Proj.Home,"Simulation_Test_Results/Figures/1B_Predicted_Reassortment.pdf"),
       dpi = 300,
       width = 4,
       height = 4,
       unit = "in")

write.csv(Cell.Summary, file = file.path(Proj.Home,"Simulation_Test_Results/Fig1_Exp_Pp_Reassortment.csv"))
write.csv(Complete.Summary, file = file.path(Proj.Home,"Simulation_Test_Results/Fig1_Pp1_Reassortment.csv"))
