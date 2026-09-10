###########################################################################################################
##                                                                                                       ##
## Title: Medusa cell speed and movement characteristics analysis                                        ##
## Author: Fernando J. Bascon                                                                            ##
## Affiliations: Multicellgenome Lab, Institut de Biologia Evolutiva, CSIC-UPF                           ##
## Date: 2023-2026                                                                                       ##
##                                                                                                       ##
###########################################################################################################


## For this analysis, 37 different Medusa cells where filmed while moving with the objective of characterize
## this behavior. The videos were analysed using Fiji (1), specifically the Manual Tracking plugin (2). 

## For this analysis, all the data was compiled in a .csv file, having their ID Number changed for MC (Medusa cell) or
## C (Control) and the corresponding number. For controls, both coenocytes and cells with only one nuclei where used. 
## At first each sample was processed and analysed independently. The comparison of both samples can be found at the end.
## The information regardless each value belonged to which type of cell was added manually to the .csv file using R. 

## For this analysis the R packages used are: 
## - ggplot2 (3) for the graphics generation
## - ggsignif (4) for adding the statistical significance to the plots
## - celltrackR (5), a package designed to analyse trajectories and movement characteristics of objects

library(ggplot2)
library(ggsignif)
library(celltrackR)

###########################################################################################################
## 1. Medusa Cells Analysis                                                                              ##
###########################################################################################################

## The data from Medusa cells is located in the same folder, with the name: MCRawData.csv

MCData <- read.table("MCRawData.csv", header=T, sep=",")

####################
## Speed Analysis ##
####################

## The first thing we do is plot the data and calculate the speed mean of Medusa cells:

MCData$Track <- factor(MCData$Track, levels = c("MC1","MC2","MC3","MC4","MC5","MC6","MC7","MC8","MC9","MC10",
                                            "MC11","MC12","MC13","MC14","MC15","MC16","MC17","MC18","MC19","MC20",
                                            "MC21","MC22","MC23","MC24","MC25","MC26","MC27","MC28","MC29","MC30",
                                            "MC31","MC32","MC33","MC34","MC35","MC36","MC37"), ordered = T)

ggplot(MCData, aes(x=Type, y=Velocity)) +
  geom_boxplot() + 
  labs(title = "Average Velocity of Taxi Cells while moving", x="Cell", y="Velocity (µm/min)") +
  theme_light()

summary(MCData$Velocity) ## Mean: 7.484µm/min

## To illustrate the variation in speed between each Medusa cell, we also plot each individual cell:

ggplot(MCData, aes(x=Track, y=Velocity, fill=Type)) +
  geom_boxplot() + 
  labs(title = "Speed of individual medusa cells over time", x="Individual cells", y="Speed (µm/min)") +
  theme_light() + 
  theme(axis.title = element_blank()) +
  scale_fill_manual(values=c("#4D12DE33")) +
  theme(axis.title = element_text(face = "bold.italic", size = 11)) +
  theme(plot.title = element_text(face = "bold.italic", size = 13)) +
  theme(axis.text = element_text(face = "bold")) +
  theme(legend.position = "none")

####################
## Track Analysis ##
####################

## With celltrackR (5) we can plot the trajectories of each cell. When reading the data from Fiji (1) with R the 
## coordinates of the Y axis get inverted. Also, celtrackR (5) needs to load a new file to apply the format it needs
## to work. For this reason, after applying the following functions (also accounting for the pixel to µm ratio so the
## plot is shown in µm) the data was saved in a new file: MCTrackData.csv, located in the same folder as the raw data.

MCTrackData <- MCData
MCTrackData$Y <- ((MCTrackData$Y)-1921)*(-1)
MCTrackData$Y <- MCTrackData$Y*0.339
MCTrackData$X <- MCTrackData$X*0.339

## After these adjustments and save it into the separate file we can continue plotting directly with celtrackR (5).

MCTrackData <- read.tracks.csv( "MCTrackData.csv", sep=",",
                         header = TRUE, 
                         id.column = 2, time.column = 3, pos.columns = 4:5 )

plot(MCTrackData, ylab = "Y Axis (µm)", xlab = "X Axis (µm)")


#####################
## Speed Over Time ##
#####################

## The objective of this plot is to represent how the speed of each Medusa cell varies over time
## For this purpose we decided to do a HeatMap. 

MCDataHeatMap <- MCData

## The time is recollected in slides, not in any time measure, so before plotting we need to adjust 
## the scale. Each image was taken every 1.5min:

MCDataHeatMap$Time = (MCDataHeatMap$Time)*1.5

## After that, we can proceed to generating the plot.
MCDataHeatMap$Track <- factor(MCDataHeatMap$Track, levels = c("MC37","MC36","MC35","MC34","MC33","MC32","MC31","MC30","MC29","MC28",
                                            "MC27","MC26","MC25","MC24","MC23","MC22","MC21","MC20","MC19","MC18",
                                            "MC17","MC16","MC15","MC14","MC13","MC12","MC11","MC10","MC9","MC8",
                                            "MC7","MC6","MC5","MC4","MC3","MC2","MC1"), ordered = T)

ggplot(MCDataHeatMap, aes(Time, Track, fill= Velocity)) + 
  geom_tile() + scale_fill_gradientn(colours = c("#9A031EBF", "#F4845FBF", "#E6AF2EBF")) +
  ggtitle("Changes in speed of individual medusa cells over time") +
  xlab("Time (min)") + ylab("Medusa cell") + labs(fill = "Speed (µm/min)") + 
  theme(axis.title = element_text(face = "bold.italic", size = 11)) +
  theme(plot.title = element_text(face = "bold.italic", size = 13)) +
  theme(axis.text = element_text(face = "bold")) +
  theme(legend.title = element_text(face = "bold.italic", size = 11))

###########################################################################################################
## 2. Control Analysis                                                                                   ##
###########################################################################################################

## As before, al data taken with Fiji (1) was joined together into a .csv file. The tag "control" and the
## clarification of each one belonging to a coenocyte or a mononucleated cell was aded by hand using R. 
## The data is located in the same folder under the name CtrRawData.csv. Everything that was done with the
## Medusa cell data was repeated for the controls. 

CtrData <- read.table("CtrRawData.csv", header=T, sep=",")

####################
## Speed Analysis ##
####################

ggplot(CtrData, aes(x=Type, y=Velocity)) +
  geom_boxplot() + 
  labs(title = "Average Velocity of Control Cells", x="Cell", y="Velocity (µm/min)") +
  theme_light()

summary(CtrData$Velocity) ## Mean: 0.7103µm/min

## Here, we decided to add an extra analysis to check if there was some difference between the coenocytes and
## the mononuclear cells:

ggplot(CtrData, aes(x=Subtype, y=Velocity)) +
  geom_boxplot() + 
  labs(title = "Average Velocity of Coenocytes vs Mononuclear Controls", x="Cell", y="Velocity (µm/min)") +
  theme_light()

## To decide which statistical test to use to calculate the p-value between the speeds of both controls, wee need
## to test first the normal distribution of the data:

shapiro.test(CtrData$Velocity) # p-value < 2.2e-16

## Since it does not follow a normal distribution, the statistical test needed is the wilcoxon test:

wilcox.test(CtrData[(CtrData$Subtype == "Coenocyte"),]$Velocity, CtrData[(CtrData$Subtype == "Uninuclear"),]$Velocity) 
# p-value = 0.08483, not significant

## As we can see, there is no statistical difference between the speeds of both controls, so from this point on,
## all the sample was counted as a unique control sample. 

#####################
## Tracks Analysis ##
#####################

## As before, the data needed to be edited to adjust for the inversion of the Y axis and applying the pixel to µm
## ratio for the plots. The code used was the same one as for the Medusa cells, so we will proceed directly with 
## the plotting. The file generated after the editing is located in the same file as everything else with the name
## CtrTrackData.csv.

CtrTrackData <- read.tracks.csv( "CtrTrackData.csv", sep=",",
                              header = TRUE, 
                              id.column = 2, time.column = 3, pos.columns = 4:5 )

plot(CtrTrackData, ylab = "Y Axis (µm)", xlab = "X Axis (µm)")


###########################################################################################################
## 3. Medusa cells vs Controls                                                                           ##
###########################################################################################################

## Now we follow with the comparison between Medusa cells and controls, following the same workflow as before. 
## Both Raw Data files were combined into one, located in the same folder as the rest and called CombinedRawData.csv.

CombinedData <- read.table("CombinedRawData.csv", header=T, sep=",")

####################
## Speed Analysis ##
####################

ggplot(CombinedData, aes(x=Type, y=Velocity, fill=Type)) +
  geom_violin() + 
  labs(title = "Average speed of Medusa cells", x="Sample", y="Speed (µm/min)") +
  theme_light() + 
  theme(axis.title = element_blank()) +
  scale_fill_manual(values=c("#E6AF2E66","#F3A68D99")) +
  theme(axis.title = element_text(face = "bold.italic", size = 11)) +
  theme(plot.title = element_text(face = "bold.italic", size = 13)) +
  theme(axis.text = element_text(face = "bold")) + 
  theme(legend.position = "none") +
  geom_signif(data=CombinedData,comparisons=list(c("Control","Medusa Cell")),map_signif_level = TRUE,annotations="****") +
  theme(legend.position = "none")

shapiro.test(CombinedData$Velocity) # p-value < 2.2e-16

wilcox.test(CombinedData[(CombinedData$Type == "Medusa Cell"),]$Velocity, CombinedData[(CombinedData$Type == "Control"),]$Velocity) 
# p-value < 2.2e-16, completely significant

##########################################
## Mean Square Displacement (MSD) Plots ##
##########################################

## This time, instead of analyzing the tracks we generated a MSD Plot.

## MSD plots are a kind of visual representation that allows us to see if any given movement is caused by
## diffusion or if it has any kind of directionality. For a range of time, the average square displacement
## of all independent movements (frame by frame) is calculated, resulting in a linear plot. The presence or
## not of a slope in the plot and the size of that slope is proportional to the presence of a directional
## force driving the movement. 

## AThis time, since we already have all the tracks loaded and the analysis is donde independently, we do not
## need to combine all in one file:

MedusaCell.msd <- aggregate( MCTrackData, squareDisplacement, FUN = "mean.se" )
MedusaCell.msd$cells <- "Medusa cells"
MedusaCell.msd$dt <- MedusaCell.msd$i * timeStep( MCTrackData )
str( MedusaCell.msd )

Control.msd <- aggregate( CtrTrackData, squareDisplacement, FUN = "mean.se" )
Control.msd$cells <- "Control"
Control.msd$dt <- Control.msd$i * timeStep( CtrTrackData )
str( Control.msd )

msddata <- rbind( MedusaCell.msd, Control.msd )
names(msddata)[5] <- paste("Sample")
head(msddata)


ggplot( msddata, aes( x = dt , y = mean, color = Sample, fill = Sample ) ) +
  geom_ribbon( aes( ymin = lower, ymax = upper) , alpha = 0.2 ,color=NA ) +
  geom_line( ) +
  scale_color_manual(values=c("#E6AF2E","#F3A68D")) +
  scale_fill_manual(values=c("#E6AF2E66","#F3A68D99")) +
  labs( x = expression( paste(Delta,"t (minutes)") ),
        y = "Mean square displacement", title = "Medusa cells' movement directionality") +
  theme_light() + theme(legend.position = "top") + xlim(0,100) +
  theme(axis.title = element_text(face = "bold.italic", size = 11)) +
  theme(plot.title = element_text(face = "bold.italic", size = 13)) +
  theme(axis.text = element_text(face = "bold")) +
  theme(legend.title = element_text(face="bold.italic", size =11)) +
  theme(legend.text = element_text(face="bold")) +
  theme(legend.position = "right")


###########################################################################################################

## Citations:
## 1. Schindelin, J., Arganda-Carreras, I., Frise, E. et al. Fiji: an open-source platform for 
##    biological-image analysis. Nat Methods 9, 676-682 (2012). https://doi.org/10.1038/nmeth.2019
## 2. Cordeli?res, F.P. Manual Tracking Plugin, Fiji, Available online: 
##    https://imagej.net/ij/plugins/track/track.html
## 3. H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
## 4. Ahlmann-Eltze, C., & Patil, I. (2021). ggsignif: R Package for Displaying Significance Brackets for 'ggplot2'. 
##    PsyArxiv. doi:10.31234/osf.io/7awm6
## 5. Inge M.N. Wortel, Annie Y. Liu, Katharina Dannenberg, Jeffrey C. Berry, Mark J. Miller, Johannes Textor (2021). 
##    CelltrackR: an R package for fast and flexible analysis of immune cell migration data. ImmunoInformatics, 1-2:p100003. 
##    DOI: https://doi.org/10.1016/j.immuno.2021.100003 URL: https://ingewortel.github.io/celltrackR/.
