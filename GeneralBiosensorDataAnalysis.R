##########NOTE##########
#If this code doesn't work when downloaded, check that R hasn't automatically changed the "Temp. [°C]" label in lines 37 and 54.
########################

#Code written by Anya Lindstrom Battle in March 2020 for the analysis of data obtained during rotation project with Prof Wei Huang.

#Before running this code make a plate plan in the correct format.
#Also separate the OD600 and Luminescence data and remove metadata.
#Please see example files in GitHub Repository; all files need to be saved as csv files.

#Installing packages needed to run this code - only do this once; it will take a while
install.packages(c("reshape2", "dplyr", "ggplot2"))

###########Do not touch this bit - this is a function which is used below for data analysis#############
# Create a function that calculates 95% confidence intervals for the given
# data vector using a t-distribution
conf_int95 <- function(data) {
  n <- length(data)
  error <- qt(0.95, df=n-1) * sd(data)/sqrt(n)
  return(error)
}
########################################################################################################

#Set the working directory (where R imports all the data from)
#Example: setwd("C:/Users/Anya/Dropbox/PhD Oxford 2/Huang rotation/Exp005_Results")

#Manipulate OD600 data to be in right format
rawdata <- read.csv("ExampleBiosensorGrowthData.csv", header = FALSE)#Change this line to reflect your own file name 
rawdata_transpose <- as.data.frame(t(rawdata))

rawdata_transpose[] <- lapply(rawdata_transpose, as.character)
colnames(rawdata_transpose) <- rawdata_transpose[1, ]
rawdata_transpose <- rawdata_transpose[-1 ,]
rawdata_transpose = rawdata_transpose[, -1]

#Re-organise the data by well
library(reshape2)
reshaped <- (melt(rawdata_transpose, id = c("Time [s]", "Temp. [°C]"), variable.name = "Well", value.name = "OD600"))

#View your reshaped data
summary(reshaped)
head(reshaped)

#Manipulate Lum Data to be in right format - this is essentially exact same as done above for OD600

rawdata_Lum <- read.csv("ExampleBiosensorLumData.csv", header = FALSE) #Again, change this file name
rawdata_Lum_transpose <- as.data.frame(t(rawdata_Lum))
rawdata_Lum_transpose[] <- lapply(rawdata_Lum_transpose, as.character)
colnames(rawdata_Lum_transpose) <- rawdata_Lum_transpose[1, ]
rawdata_Lum_transpose <- rawdata_Lum_transpose[-1 ,]
rawdata_Lum_transpose <- rawdata_Lum_transpose[, -1]
rawdata_Lum_transpose <- rawdata_Lum_transpose[-c(100:103) ]

library(reshape2)
reshaped_Lum <- (melt(rawdata_Lum_transpose, id = c("Time [s]", "Temp. [°C]"), variable.name = "Well", value.name = "Lum"))
head(reshaped_Lum)

#Import the platemap
platemap <- read.csv("ExamplePlatePlan.csv") #Change this to reflect the name of your file
library(dplyr)

#Organise the data by the information provided in the plate map
annotated <- inner_join(reshaped, platemap, by="Well")
annotated_Lum <- inner_join(reshaped_Lum, platemap, by="Well")
write.csv(annotated, "BiosensorGrowthDataAnnotated.csv") #Saves your annotated file in the specified directory
write.csv(annotated_Lum, "BiosensorLumDataAnnotated.csv")

#Calculate the mean, standard deviation, and 95 % confidence interval for your OD600 and Lum data
#OD600
stats <- annotated %>%
  group_by(Environment,Strain,`Time [s]`) %>%
  summarise(N = length(as.numeric(OD600)), Average = mean(as.numeric(OD600)), StDev = sd(as.numeric(OD600)), CI95=conf_int95(OD600)) %>%
  filter(!is.na(Strain))
#Lum
stats_Lum <- annotated_Lum %>%
  group_by(Environment,Strain,`Time [s]`) %>%
  summarise(N = length(as.numeric(Lum)), Average = mean(as.numeric(Lum)), StDev = sd(as.numeric(Lum)), CI95=conf_int95(Lum)) %>%
  filter(!is.na(Strain))

####################Plotting results########################

#Creating custom labels for the graph
strain.labs <- c("BL21\n Alk-", "BL21\n Alk+ 1", "BL21\n Alk+ 2", "BL21\n Alk+ 3", "C43\n Alk-", "C43\n Alk+ 1", "C43\n Alk+ 2", "C43\n Alk+ 3")
names(strain.labs) <- c("1","2","3","4", "5", "6", "7", "8")
env.labs <- c("400 uM IPTG", "25 uM IPTG", "0 uM IPTG")
names(env.labs) <- c("1","2","3")

#Plotting OD over time
library(ggplot2)
OD_biosensor <- ggplot(stats, aes(x=as.numeric(`Time [s]`), y=Average)) + 
  geom_ribbon(aes(ymin=Average-CI95, ymax=Average+CI95),
             color=NA, alpha=0.3) +
  geom_line() +
  scale_y_log10() +
  facet_grid(
    Environment ~ Strain,
  labeller = labeller(Environment = env.labs, Strain = strain.labs))+
  labs(x="Time [s]", y=expression("OD"[600]))

#Plotting Lum over time
Lum_biosensor <- ggplot(stats_Lum, aes(x=as.numeric(`Time [s]`), y=Average)) + 
  geom_ribbon(aes(ymin=Average-CI95, ymax=Average+CI95),
              color=NA, alpha=0.3) +
  geom_line() +
  facet_grid(
    Environment ~ Strain,
  labeller = labeller(Strain = strain.labs, Environment = env.labs))+
  labs(x="Time [s]", y="Luminescence")

#To plot Lum/OD we first need to calculate this and add this value as a column to our data
#Making Lum/OD dataset
OD_Lum <- reshaped
OD_Lum$Lum <- as.numeric(reshaped_Lum$Lum)
OD_Lum$OD600 <- as.numeric(OD_Lum$OD600)
OD_Lum$Lum_OD <- OD_Lum$Lum/OD_Lum$OD600

annotated_OD_Lum <- inner_join(OD_Lum, platemap, by="Well")
write.csv(annotated_OD_Lum, "BiosensorODLumDataAnnotated.csv")#This saves the final dataset as a new file in the working directory

#Calculate same stats as above for the Lum/OD data
stats_OD_Lum <- annotated_OD_Lum %>%
  group_by(Environment,Strain,`Time [s]`) %>%
  summarise(N = length(as.numeric(Lum_OD)), Average = mean(as.numeric(Lum_OD)), StDev = sd(as.numeric(Lum_OD)), CI95=conf_int95(Lum_OD)) %>%
  filter(!is.na(Strain))

#Plot the Lum/OD Data
Lum_OD_biosensor <- ggplot(stats_OD_Lum, aes(x=as.numeric(`Time [s]`), y=Average)) + 
  geom_ribbon(aes(x=as.numeric(`Time [s]`), ymin=Average-CI95, ymax=Average+CI95),inherit.aes = FALSE,
              color=NA, alpha=0.3) +
  geom_line(colour = "black", lwd = 1) +
  facet_grid(
    Environment ~ Strain,
    labeller = labeller(Strain = strain.labs, Environment = env.labs))+
  labs(x="Time [s]", y=expression("Luminescence/OD"[600]))



#############################Plotting the positive and negative controls#########################
#This bit you will need to change as I am specifying well numbers based on my own experimental design.

#Extract positive control data
#Extract it for OD600 Data and plot it
rawdata_transpose_posctrl <- rawdata_transpose[c(1:2, 12,17,22, 27,32,37,42)] # first 3 are 0.1% oil, remaining 4 are 0.01 % oil.

#Plotting one of the positive control traces to check it looks OK
postctrl <- ggplot(rawdata_transpose_posctrl, aes(x=as.numeric(`Time [s]`), y=as.numeric(B5))) + 
  geom_line() +
  labs(x="Time [s]", y="OD600")

#Now do the same for the Lum data
rawdata_Lum_posctrl <- rawdata_Lum_transpose[c(1:2, 12,17,22, 27,32,37,42)]
postctrl_Lum <- ggplot(rawdata_Lum_posctrl, aes(x=as.numeric(`Time [s]`), y=as.numeric(B5))) + 
  geom_line() +
  labs(x="Time [s]", y="Luminescence")

#Calculate Lum/OD values for positive control
rawdata_ODLum_posctrl <- rawdata_transpose_posctrl #Create the rawdata_ODLUm_posctrl dataframe from the OD600-containing dataframe

#Add Lum data to this dataframe 
rawdata_ODLum_posctrl <- cbind(rawdata_ODLum_posctrl, Lum = rawdata_Lum_posctrl[,3:9])

#Calculate OD/Lum data and add it to the dataframe
#Define a for loop to calculate Lum/OD

for(i in 3:9){
  Lum_OD <- as.numeric(rawdata_ODLum_posctrl[,i+7])/as.numeric(rawdata_ODLum_posctrl[,i])
  df <- data.frame(Lum_OD)
  rawdata_ODLum_posctrl <- cbind(rawdata_ODLum_posctrl, df)
}

names(rawdata_ODLum_posctrl)[17:23] <- c("Lum_OD.B5", "Lum_OD.C5", "Lum_OD.D5", "Lum_OD.E5", "Lum_OD.F5", "Lum_OD.G5", "Lum_OD.H5")

#Calculate the mean and 95 % confidence interval for your data
rawdata_ODLum_posctrl$AverageLumOD001 <- rowMeans(rawdata_ODLum_posctrl[,17:19])
rawdata_ODLum_posctrl$conf_int_LumOD001 <- apply(rawdata_ODLum_posctrl[,17:19],1,conf_int95)
rawdata_ODLum_posctrl$AverageLumOD01 <- rowMeans(rawdata_ODLum_posctrl[,20:23])
rawdata_ODLum_posctrl$conf_int_LumOD01 <- apply(rawdata_ODLum_posctrl[,20:23],1,conf_int95)

#Plotting the positive control data
postctrl_ODLum <- ggplot(rawdata_ODLum_posctrl, aes(x=as.numeric(`Time [s]`), y=as.numeric(AverageLumOD01)), inherit.aes = FALSE) + 
  geom_ribbon(aes(ymin=AverageLumOD01-conf_int_LumOD01, ymax=AverageLumOD01+conf_int_LumOD01),
             color=NA, alpha=0.3) +
  geom_line(colour = "chartreuse", lwd = 1) +
  geom_line(data = rawdata_ODLum_posctrl, aes(x=as.numeric(`Time [s]`), y=as.numeric(AverageLumOD001)), inherit.aes = FALSE, colour = "chartreuse3", lwd = 1) + 
  geom_ribbon(aes(ymin=AverageLumOD001-conf_int_LumOD001, ymax=AverageLumOD001+conf_int_LumOD001),
              color=NA, alpha=0.3)+
  labs(x="Time [s]", y=expression("Luminescence/OD"[600]))

#Negative control
#Lum/OD Data
#Extract the negative control data
rawdata_transpose_negctrl <- rawdata_transpose[c(1:2, 6,11,16,21,26,31,36,41)]
rawdata_Lum_negctrl <- rawdata_Lum_transpose[c(1:2, 6,11,16,21,26,31,36,41)]

#Calculate Lum/OD values for negative control
rawdata_ODLum_negctrl <- rawdata_transpose_negctrl#Create the rawdata_ODLUm_negctrl dataframe from the OD600-containing dataframe

#Add Lum data to this dataframe 
rawdata_ODLum_negctrl <- cbind(rawdata_ODLum_negctrl, Lum = rawdata_Lum_negctrl[,3:10])

#Calculate OD/Lum data and add it to the dataframe
#Define a for loop to calculate Lum/OD

for(i in 3:10){
  Lum_ODneg <- as.numeric(rawdata_ODLum_negctrl[,i+8])/as.numeric(rawdata_ODLum_negctrl[,i])
  dfneg <- data.frame(Lum_ODneg)
  rawdata_ODLum_negctrl <- cbind(rawdata_ODLum_negctrl, Lum_ODneg)
}
names(rawdata_ODLum_negctrl)[19:26] <- c("Lum_OD.A4", "Lum_OD.B4", "Lum_OD.C4", "Lum_OD.D4", "Lum_OD.E4", "Lum_OD.F4", "Lum_OD.G4", "Lum_OD.H4")

#Calculate the mean and 95 % confidence interval for your data
rawdata_ODLum_negctrl$AverageLumOD <- rowMeans(rawdata_ODLum_negctrl[,19:26])
rawdata_ODLum_negctrl$conf_int_LumOD <- apply(rawdata_ODLum_negctrl[,19:26],1,conf_int95)

negctrl_ODLum <- ggplot(rawdata_ODLum_negctrl, aes(x=as.numeric(`Time [s]`), y=as.numeric(AverageLumOD))) + 
  geom_ribbon(aes(ymin=AverageLumOD-conf_int_LumOD, ymax=AverageLumOD+conf_int_LumOD),
             color=NA, alpha=0.3) +
  geom_line() +
  labs(x="Time (minutes)", y="Luminescence/OD")

###################Make summary plot with positive and negative controls##########################
Lum_OD_biosensor_ctrls <-  Lum_OD_biosensor +
 geom_line(data = rawdata_ODLum_posctrl, mapping = aes(x=as.numeric(`Time [s]`), y=as.numeric(AverageLumOD01)), colour = "chartreuse", lwd = 1)+
  geom_ribbon(data = rawdata_ODLum_posctrl, aes(x=as.numeric(`Time [s]`), ymin=AverageLumOD01-conf_int_LumOD01, ymax=AverageLumOD01+conf_int_LumOD01),
              inherit.aes = FALSE, color=NA, alpha=0.3) +
  geom_line(data = rawdata_ODLum_posctrl, mapping = aes(x=as.numeric(`Time [s]`), y=as.numeric(AverageLumOD001)), colour = "chartreuse3", lwd = 1)+
  geom_ribbon(data = rawdata_ODLum_posctrl, aes(x=as.numeric(`Time [s]`), ymin=AverageLumOD001-conf_int_LumOD001, ymax=AverageLumOD001+conf_int_LumOD001),
              inherit.aes = FALSE, color=NA, alpha=0.3) +
  geom_line(data = rawdata_ODLum_negctrl, mapping = aes(x=as.numeric(`Time [s]`), y=AverageLumOD), colour = "red", lwd = 1)+
  geom_ribbon(data = rawdata_ODLum_negctrl, aes(x=as.numeric(`Time [s]`), ymin=AverageLumOD-conf_int_LumOD, ymax=AverageLumOD+conf_int_LumOD),
             inherit.aes = FALSE, color=NA, alpha=0.3)+
  theme(text = element_text(size = 16))+
  theme(axis.text.x = element_text(angle=90, size = 10, vjust = 0.5, hjust = 1))+
  theme(axis.text.y = element_text(size=10))+
  scale_y_continuous(labels = scales::scientific)
 
############################################################## 
