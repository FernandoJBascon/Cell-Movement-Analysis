###########################################################################################################
##                                                                                                       ##
## Title: Medusa cell speed and movement characteristics analysis before and after head/tail separation  ##
## Author: Fernando J. Bascon                                                                            ##
## Affiliations: Multicellgenome Lab, Institut de Biologia Evolutiva, CSIC-UPF                           ##
## Date: 2023-2026                                                                                       ##
##                                                                                                       ##
###########################################################################################################

## For this experiment 14 different Medusa cells were filmed from the moment of their birth until time after
## head/tail separates. The aim of this experiment is to analyse the speed and movement characteristics of
## Medusa cells vs those of head and tail independently after the moment of separation. 

## The videos were analysed using Fiji (1), specifically the Manual Tracking plugin (2). For comparing 
## the movement of both parts before and after, both Head and Tail were considered independent objects 
## in both scenarios. 

## The way Manual Tracking (2) works, it produce an independent raw data file for each sample (these sample
## being Head Before, Head After, Tail Before and Tail After), named numerically, For the purpose of making
## the subsequent analysis easier, all files were fused together, all samples name based on Head or Tail, 
## Before or After and a new column was added with the Track Number. The independent files for each Medusa cell
## are located in the folder "Independent Raw Data" and the fused file used for this analysis is called: RawData.csv.

## For this analysis the R packages used are: 
## - ggplot2 (3) for the graphics generation
## - ggsignif (4) for adding the statistical significance to the plots
## - celltrackR (5), a package designed to analyse trajectories and movement characteristics of objects

DATA <- read.table("RawData.csv", header=T, sep=",")

library(ggplot2)
library(ggsignif)
library(celltrackR)


###########################################################################################################
## 1. Analysis of speed                                                                                  ##
###########################################################################################################

## For this we decided to plot the data using a Violin Plot, since it represents better the distribution
## of our data. Statistical analysis done for this plot is written after the plot itself. 

DATA$Track <- factor(DATA$Track, levels = c("Head Before","Head After","Tail Before","Tail After"), 
                     ordered = TRUE) 

ggplot(DATA, aes(x=Track, y=Velocity, fill=Track)) +
  geom_violin() + 
  labs(title = "Speed of head and tail before and after separation during Medusa cell movement", x="Medusa cell part", y="Speed (µm/min)") +
  theme_light() + 
  scale_fill_manual(values=c("#4D12DE33","#4D12DE33","#F9B27AB3","#F9B27AB3")) +
  theme(axis.title = element_blank()) +
  theme(axis.title = element_text(face = "bold.italic", size = 11)) +
  theme(plot.title = element_text(face = "bold.italic", size = 13)) +
  theme(axis.text = element_text(face = "bold"))  +
  geom_signif(data=DATA,comparisons=list(c("Head Before","Head After")),map_signif_level = TRUE,annotations="****", y_position = 40) +
  geom_signif(data=DATA,comparisons=list(c("Tail Before","Tail After")),map_signif_level = TRUE,annotations="****", y_position = 40) +
  geom_signif(data=DATA,comparisons=list(c("Head After","Tail After")),map_signif_level = TRUE,annotations="****", y_position = 35) +
  theme(legend.position = "none")

## The mean speed for each sample is: 
summary(DATA$Velocity[DATA$Track == "Head Before"]) ## 6.47µm/min
summary(DATA$Velocity[DATA$Track == "Head After"]) ## 1.70µm/min
summary(DATA$Velocity[DATA$Track == "Tail Before"]) ## 5.06µm/min
summary(DATA$Velocity[DATA$Track == "Tail After"]) ## 0.92µm/min

## In order to do the statistical analysis we need to decide which is the most appropiate test for our data:
## t-test or wilcoxon test. For this we check first if our data follow a normal distribution or not. 

shapiro.test(DATA$Velocity) ## p-value < 2.2e-16. The data is not normally distributed

## To further confirm that this is the case, we also perform the normality analysis sample by sample:

HeadBefore <- DATA[DATA$Track == "Head Before",]
TailBefore <- DATA[DATA$Track == "Tail Before",]
HeadAfter <- DATA[DATA$Track == "Head After",]
TailAfter <- DATA[DATA$Track == "Tail After",]

shapiro.test(HeadBefore$Velocity) ## p-value < 2.2e-16. The data is not normally distributed
shapiro.test(TailBefore$Velocity) ## p-value < 2.2e-16. The data is not normally distributed
shapiro.test(HeadAfter$Velocity) ## p-value < 2.2e-16. The data is not normally distributed
shapiro.test(TailAfter$Velocity) ## p-value < 2.2e-16. The data is not normally distributed

## Since our data does not follow a normal distribution, the most appropiate test for significance calculation
## is the wilcoxon test:

wilcox.test(HeadBefore$Velocity, HeadAfter$Velocity, alternative = "greater") # p-value < 2.2e-16
wilcox.test(TailBefore$Velocity, TailAfter$Velocity, alternative = "greater") # p-value < 2.2e-16
wilcox.test(HeadBefore$Velocity, TailBefore$Velocity, alternative = "greater") # p-value < 2.2e-16
wilcox.test(HeadAfter$Velocity, TailAfter$Velocity, alternative = "greater") # p-value < 2.2e-16

###########################################################################################################
## 2. Analysis of Tracks                                                                                 ##
###########################################################################################################

## This analysis consist on generating plots in which the path taken by each Medusa cell can be seen as they
## moved throughout the plate.The program used for generating the raw data is not read properly,
## since R seems to read the Y axis values inverted. In order to generate these track plots we need to 
## solve this axis inversion first and save it into a new data file. This needs to be done since
## celltrackR (5) needs to read an external file to apply the necessary format. 

DATATrack <- DATA
DATATrack$Y <- ((DATATrack$Y)-1921)*(-1)
DATATrack$Concatenated <- paste(DATATrack$Track, DATATrack$TrackNumber)

## After doing this changes, the data was saved into TrackData.csv, which we need to read directly with 
## celltrackR (5):

TrackData <- read.tracks.csv( "TrackData.csv", sep=",",
                         header = TRUE, 
                         id.column = 12, time.column = 5, pos.columns = 6:7 )
plot(TrackData)

## For plotting each independent Medusa cell, this data was also separated in 14 different Tracks. At the 
## moment this was done, both, the X and Y axis were multiplied by 0.339 to account for the pixel value.
## This way, and for future plots, the X and Y axis will show micrometers instead of number of pixels. 

## This is an example of how each independent track was processed. All of them where then saved into
## independent files, located in the folder "Independent Tracks Data":

Track512 <- DATATrack[DATATrack$TrackNumber == 512,]
Track512$Y <- Track512$Y*0.339
Track512$X <- Track512$X*0.339

## With all the new independent files, we proceed to load them and plot them all together:

Track221 <- read.tracks.csv( "Independent Tracks Data/Track221.csv", sep=",",
                             header = TRUE, 
                             id.column = 4, time.column = 5, pos.columns = 6:7 )
Track222 <- read.tracks.csv( "Independent Tracks Data/Track222.csv", sep=",",
                             header = TRUE, 
                             id.column = 4, time.column = 5, pos.columns = 6:7 )
Track223 <- read.tracks.csv( "Independent Tracks Data/Track223.csv", sep=",",
                             header = TRUE, 
                             id.column = 4, time.column = 5, pos.columns = 6:7 )
Track224 <- read.tracks.csv( "Independent Tracks Data/Track224.csv", sep=",",
                             header = TRUE, 
                             id.column = 4, time.column = 5, pos.columns = 6:7 )
Track225 <- read.tracks.csv( "Independent Tracks Data/Track225.csv", sep=",",
                             header = TRUE, 
                             id.column = 4, time.column = 5, pos.columns = 6:7 )
Track311 <- read.tracks.csv( "Independent Tracks Data/Track311.csv", sep=",",
                             header = TRUE, 
                             id.column = 4, time.column = 5, pos.columns = 6:7 )
Track312 <- read.tracks.csv( "Independent Tracks Data/Track312.csv", sep=",",
                             header = TRUE, 
                             id.column = 4, time.column = 5, pos.columns = 6:7 )
Track321 <- read.tracks.csv( "Independent Tracks Data/Track321.csv", sep=",",
                             header = TRUE, 
                             id.column = 4, time.column = 5, pos.columns = 6:7 )
Track411 <- read.tracks.csv( "Independent Tracks Data/Track411.csv", sep=",",
                             header = TRUE, 
                             id.column = 4, time.column = 5, pos.columns = 6:7 )
Track421 <- read.tracks.csv( "Independent Tracks Data/Track421.csv", sep=",",
                             header = TRUE, 
                             id.column = 4, time.column = 5, pos.columns = 6:7 )
Track422 <- read.tracks.csv( "Independent Tracks Data/Track422.csv", sep=",",
                             header = TRUE, 
                             id.column = 4, time.column = 5, pos.columns = 6:7 )
Track423 <- read.tracks.csv( "Independent Tracks Data/Track423.csv", sep=",",
                             header = TRUE, 
                             id.column = 4, time.column = 5, pos.columns = 6:7 )
Track511 <- read.tracks.csv( "Independent Tracks Data/Track511.csv", sep=",",
                             header = TRUE, 
                             id.column = 4, time.column = 5, pos.columns = 6:7 )
Track512 <- read.tracks.csv( "Independent Tracks Data/Track512.csv", sep=",",
                             header = TRUE, 
                             id.column = 4, time.column = 5, pos.columns = 6:7 )

par(mfrow = c(3, 5))
plot(Track221, col = c("#F9B27AFF","#F35B04FF","#4D12DE8C","#4D12DEE6"), ylab = "Y Axis (µm)", xlab = "X Axis (µm)")
plot(Track222, col = c("#F9B27AFF","#F35B04FF","#4D12DE8C","#4D12DEE6"), ylab = "Y Axis (µm)", xlab = "X Axis (µm)")
plot(Track223, col = c("#F9B27AFF","#F35B04FF","#4D12DE8C","#4D12DEE6"), ylab = "Y Axis (µm)", xlab = "X Axis (µm)")
plot(Track224, col = c("#F9B27AFF","#F35B04FF","#4D12DE8C","#4D12DEE6"), ylab = "Y Axis (µm)", xlab = "X Axis (µm)")
plot(Track225, col = c("#F9B27AFF","#F35B04FF","#4D12DE8C","#4D12DEE6"), ylab = "Y Axis (µm)", xlab = "X Axis (µm)")
plot(Track311, col = c("#4D12DEE6","#4D12DE8C","#F35B04FF","#F9B27AFF"), ylab = "Y Axis (µm)", xlab = "X Axis (µm)")
plot(Track312, col = c("#F9B27AFF","#F35B04FF","#4D12DE8C","#4D12DEE6"), ylab = "Y Axis (µm)", xlab = "X Axis (µm)")
plot(Track321, col = c("#F9B27AFF","#F35B04FF","#4D12DE8C","#4D12DEE6"), ylab = "Y Axis (µm)", xlab = "X Axis (µm)")
plot(Track411, col = c("#F9B27AFF","#F35B04FF","#4D12DE8C","#4D12DEE6"), ylab = "Y Axis (µm)", xlab = "X Axis (µm)")
plot(Track421, col = c("#F9B27AFF","#F35B04FF","#4D12DE8C","#4D12DEE6"), ylab = "Y Axis (µm)", xlab = "X Axis (µm)")
plot(Track422, col = c("#F9B27AFF","#F35B04FF","#4D12DE8C","#4D12DEE6"), ylab = "Y Axis (µm)", xlab = "X Axis (µm)")
plot(Track423, col = c("#F9B27AFF","#F35B04FF","#4D12DE8C","#4D12DEE6"), ylab = "Y Axis (µm)", xlab = "X Axis (µm)")
plot(Track511, col = c("#F9B27AFF","#F35B04FF","#4D12DE8C","#4D12DEE6"), ylab = "Y Axis (µm)", xlab = "X Axis (µm)")
plot(Track512, col = c("#F9B27AFF","#F35B04FF","#4D12DE8C","#4D12DEE6"), ylab = "Y Axis (µm)", xlab = "X Axis (µm)")
plot(NULL,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", legend=c("Tail Before","Tail After","Head Before", "Head After"),
       col=c("#F9B27AFF","#F35B04FF","#4D12DE8C","#4D12DEE6"), lty=1, cex=1)
mtext("Legend", cex=1)
par(mfrow = c(1, 1))

## Here is the code for ploting only one of the Medusa cells as an example: 
plot(Track223, col = c("#F9B27AFF","#F35B04FF","#4D12DE8C","#4D12DEE6"), ylab ="Y Axis (µm)", xlab ="X Axis (µm)", main = "Tracks of one medusa cell before and after head/tail separation")

###########################################################################################################
## 3. Mean Square Displacement (MSD) plots                                                               ##
###########################################################################################################

## MSD plots are a kind of visual representation that allows us to see if any given movement is caused by
## diffusion or if it has any kind of directionality. For a range of time, the average square displacement
## of all independent movements (frame by frame) is calculated, resulting in a linear plot. The presence or
## not of a slope in the plot and the size of that slope is proportional to the presence of a directional
## force driving the movement. 

## For this analysis, we need the data to be classified in independent files for each sample (Head Before and 
## After, Tail Before and After) instead of individual tracks for each Medusa cell. The following is an example of how 
## one of the files was processed before being saved into a file:

HeadBeforeData <- DATA[DATA$Track == "Head Before",]
HeadBeforeData$Y <- HeadBeforeData$Y*0.339 
HeadBeforeData$X <- HeadBeforeData$X*0.339

## The 4 files created are located in the folder "MSD Data". With this, we can read and plot them using celltrackR (5):

HeadBeforeTrack <- read.tracks.csv( "MSD Data/HeadBeforeMSD.csv", sep=",",
                             header = TRUE, 
                             id.column = 11, time.column = 5, pos.columns = 6:7 )
HeadAfterTrack <- read.tracks.csv( "MSD Data/HeadAfterMSD.csv", sep=",",
                                  header = TRUE, 
                                  id.column = 11, time.column = 5, pos.columns = 6:7 )
TailBeforeTrack <- read.tracks.csv( "MSD Data/TailBeforeMSD.csv", sep=",",
                                  header = TRUE, 
                                  id.column = 11, time.column = 5, pos.columns = 6:7 )
TailAfterTrack <- read.tracks.csv( "MSD Data/TailAfterMSD.csv", sep=",",
                                  header = TRUE, 
                                  id.column = 11, time.column = 5, pos.columns = 6:7 )

HeadBefore.msd <- aggregate( HeadBeforeTrack, squareDisplacement, FUN = "mean.se" )
HeadBefore.msd$cells <- "Head Before"
HeadBefore.msd$dt <- HeadBefore.msd$i * timeStep( HeadBeforeTrack )
str( HeadBefore.msd )

HeadAfter.msd <- aggregate( HeadAfterTrack, squareDisplacement, FUN = "mean.se" )
HeadAfter.msd$cells <- "Head After"
HeadAfter.msd$dt <- HeadAfter.msd$i * timeStep( HeadAfterTrack )
str( HeadAfter.msd )

TailBefore.msd <- aggregate( TailBeforeTrack, squareDisplacement, FUN = "mean.se" )
TailBefore.msd$cells <- "Tail Before"
TailBefore.msd$dt <- TailBefore.msd$i * timeStep( TailBeforeTrack )
str( TailBefore.msd )

TailAfter.msd <- aggregate( TailAfterTrack, squareDisplacement, FUN = "mean.se" )
TailAfter.msd$cells <- "Tail After"
TailAfter.msd$dt <- TailAfter.msd$i * timeStep( TailAfterTrack )
str( TailAfter.msd )

msddata <- rbind( HeadBefore.msd, HeadAfter.msd, TailBefore.msd, TailAfter.msd )
head(msddata)

ggplot( msddata, aes( x = dt , y = mean, color = cells, fill = cells ) ) +
  geom_ribbon( aes( ymin = lower, ymax = upper) , alpha = 0.2 ,color=NA ) +
  geom_line( ) +
  labs( x = expression( paste(Delta,"t (minutes)") ),
        y = "mean square displacement") +
  theme_light() + theme(legend.position = "top")

###########################################################################################################
## 4. Speed Over Time                                                                                    ##
###########################################################################################################

## The objective of this plot is to represent how the speed of each Medusa cell varies over the spam of
## life before the moment head and tail separate. For this purpose we decided to do a HeatMap. Before
## generating the plot we needed to extract the independend data for each Medusa cell only from the moment
## before separation happens: 

DataHeatMap <- DATATrack[(DATATrack$Track == "Head Before" | DATATrack$Track == "Tail Before"),]

## During the filming of the original video, each image was taken every 90sec. The file data takes each
## time point as if it were one unit, so in order to get a better representation we adjusted this, converting
## this unit into minutes:
DataHeatMap$Slice = (DataHeatMap$Slice)*1.5

## The samples names were renamed for better understanding once plotted:
DataHeatMap$Concatenated[DataHeatMap$Concatenated == "Head Before 221"] <- "H1"
DataHeatMap$Concatenated[DataHeatMap$Concatenated == "Tail Before 221"] <- "T1"
DataHeatMap$Concatenated[DataHeatMap$Concatenated == "Head Before 222"] <- "H2"
DataHeatMap$Concatenated[DataHeatMap$Concatenated == "Tail Before 222"] <- "T2"
DataHeatMap$Concatenated[DataHeatMap$Concatenated == "Head Before 223"] <- "H3"
DataHeatMap$Concatenated[DataHeatMap$Concatenated == "Tail Before 223"] <- "T3"
DataHeatMap$Concatenated[DataHeatMap$Concatenated == "Head Before 224"] <- "H4"
DataHeatMap$Concatenated[DataHeatMap$Concatenated == "Tail Before 224"] <- "T4"
DataHeatMap$Concatenated[DataHeatMap$Concatenated == "Head Before 225"] <- "H5"
DataHeatMap$Concatenated[DataHeatMap$Concatenated == "Tail Before 225"] <- "T5"
DataHeatMap$Concatenated[DataHeatMap$Concatenated == "Head Before 311"] <- "H6"
DataHeatMap$Concatenated[DataHeatMap$Concatenated == "Tail Before 311"] <- "T6"
DataHeatMap$Concatenated[DataHeatMap$Concatenated == "Head Before 312"] <- "H7"
DataHeatMap$Concatenated[DataHeatMap$Concatenated == "Tail Before 312"] <- "T7"
DataHeatMap$Concatenated[DataHeatMap$Concatenated == "Head Before 321"] <- "H8"
DataHeatMap$Concatenated[DataHeatMap$Concatenated == "Tail Before 321"] <- "T8"
DataHeatMap$Concatenated[DataHeatMap$Concatenated == "Head Before 411"] <- "H9"
DataHeatMap$Concatenated[DataHeatMap$Concatenated == "Tail Before 411"] <- "T9"
DataHeatMap$Concatenated[DataHeatMap$Concatenated == "Head Before 421"] <- "H10"
DataHeatMap$Concatenated[DataHeatMap$Concatenated == "Tail Before 421"] <- "T10"
DataHeatMap$Concatenated[DataHeatMap$Concatenated == "Head Before 422"] <- "H11"
DataHeatMap$Concatenated[DataHeatMap$Concatenated == "Tail Before 422"] <- "T11"
DataHeatMap$Concatenated[DataHeatMap$Concatenated == "Head Before 423"] <- "H12"
DataHeatMap$Concatenated[DataHeatMap$Concatenated == "Tail Before 423"] <- "T12"
DataHeatMap$Concatenated[DataHeatMap$Concatenated == "Head Before 511"] <- "H13"
DataHeatMap$Concatenated[DataHeatMap$Concatenated == "Tail Before 511"] <- "T13"
DataHeatMap$Concatenated[DataHeatMap$Concatenated == "Head Before 512"] <- "H14"
DataHeatMap$Concatenated[DataHeatMap$Concatenated == "Tail Before 512"] <- "T14"

DataHeatMap$Concatenated <- factor(DataHeatMap$Concatenated, levels = c("T14","H14","T13","H13","T12","H12","T11","H11",
                                                                        "T10","H10","T9","H9","T8","H8","T7","H7","T6","H6",
                                                                        "T5","H5","T4","H4","T3","H3","T2","H2",
                                                                        "T1","H1"), ordered = TRUE)

## Finally we generate the plot:
ggplot(DataHeatMap, aes(Slice, Concatenated, fill= Velocity)) + 
  geom_tile() + scale_fill_gradientn(colours = c("blue", "green", "yellow")) +
  ggtitle("Changes in speed over time") +
  xlab("Time (min)") + ylab("Cell") + labs(fill = "Velocity (µm/min)")

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
