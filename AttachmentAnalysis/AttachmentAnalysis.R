###########################################################################################################
##                                                                                                       ##
## Title: Medusa cell speed and movement in environments with different adherence to the substrate       ##
## Author: Fernando J. Bascon                                                                            ##
## Affiliations: Multicellgenome Lab, Institut de Biologia Evolutiva, CSIC-UPF                           ##
## Date: 2023-2026                                                                                       ##
##                                                                                                       ##
###########################################################################################################

## The objective of this experiment was to analyse which was the role of the attachment of the Medusa cells
## to the substrate during their movement. For this purpose, we treated the plates where the Medusa cells were
## filmed moving with 3 different conditions: 
##    - Control: untreated surface, 37 Medusa cells, same ones used in MedusaMovementAnalysis folder
##    - Fibronectin: increased surface attachment, 45 Medusa cells
##    - FaCellitate: decreased surface attachment, 31 Medusa cells

## The videos were then analysed with Fiji (1), using the Manual Tracking Pluggin (2). The data obtained from
## this was compiled in a file called AttachmentRawData.csv located in the same folder as this document. 
## This data was then analysed using R. The R packages used are: 
## - ggplot2 (3) for the graphics generation
## - ggsignif (4) for adding the statistical significance to the plots

library(ggplot2)
library(ggsignif)

DATA <- read.table("AttachmentRawData.csv", header=T, sep=";")

ggplot(DATA, aes(x=Sample, y=Velocity, fill=Sample)) +
  geom_violin() + 
  labs(title = "Average speed of medusa cells with different adherence", x="Adherence Sample", y="Speed (µm/min)") +
  theme_light() +
  theme(axis.title = element_blank()) +
  scale_fill_manual(values=c("#E6AF2E66","#F3A68D99","#F3665F99")) +
  theme(axis.title = element_text(face = "bold.italic", size = 11)) +
  theme(plot.title = element_text(face = "bold.italic", size = 13)) +
  theme(axis.text = element_text(face = "bold"))  +
  geom_signif(data=DATA,comparisons=list(c("Control","FaCellitate")),map_signif_level = TRUE,annotations="****") +
  geom_signif(data=DATA,comparisons=list(c("Control","Fibronectin")),map_signif_level = TRUE,annotations="****", y_position = 45) +
  theme(legend.position = "none")

summary(DATA[(DATA$Sample == "Control"),]) ## Mean: 7.484µm/min
summary(DATA[(DATA$Sample == "Fibronectin"),]) ## Mean: 1.724µm/min
summary(DATA[(DATA$Sample == "FaCellitate"),]) ## Mean: 2.290µm/min

## Before dping the statistical analysis comparing each condition, we need to test the normal distribution
## of the data to choose the best test possible. Our dataset has more than 5000 observations, making it 
## imposible to use the Shapiro Test. We decided to plot the data in histograms and chech visually:

ggplot(DATA, aes(x=Velocity)) +
  geom_histogram() + 
  theme_light() +
  theme(axis.title = element_blank()) +
  scale_fill_manual(values=c("#E6AF2E66","#F3A68D99","#F3665F99")) +
  theme(axis.title = element_text(face = "bold.italic", size = 11)) +
  theme(plot.title = element_text(face = "bold.italic", size = 13)) +
  theme(axis.text = element_text(face = "bold"))  +
  theme(legend.position = "none")

## Since the data does not have a normal distribution, we decided to use Wilcoxon Test for the statistical analysis:
wilcox.test(DATA[(DATA$Sample == "Control"),]$Velocity, DATA[(DATA$Sample == "Fibronectin"),]$Velocity) # p-value < 2.2e-16
wilcox.test(DATA[(DATA$Sample == "Control"),]$Velocity, DATA[(DATA$Sample == "FaCellitate"),]$Velocity) # p-value < 2.2e-16

###########################################################################################################

## Citations:
## 1. Schindelin, J., Arganda-Carreras, I., Frise, E. et al. Fiji: an open-source platform for 
##    biological-image analysis. Nat Methods 9, 676-682 (2012). https://doi.org/10.1038/nmeth.2019
## 2. Cordeli?res, F.P. Manual Tracking Plugin, Fiji, Available online: 
##    https://imagej.net/ij/plugins/track/track.html
## 3. H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
## 4. Ahlmann-Eltze, C., & Patil, I. (2021). ggsignif: R Package for Displaying Significance Brackets for 'ggplot2'. 
##    PsyArxiv. doi:10.31234/osf.io/7awm6

