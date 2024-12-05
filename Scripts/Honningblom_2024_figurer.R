#### Script to reproduce figures in earlier reports ####

####Packages####
library(sciplot)

####Import data####
#setwd("C:/Users/marianne.evju/OneDrive - NINA/R/honningblom")# spesifiserer hvor filer skal hentes fra og lagres
getwd()
dta <- rutedata <- import("Data/Honningblom_2024.xlsx",sheet="Honningblom_data",na=c("NA"), stringsasfactors= TRUE)
#dta <- read.delim("honningblom_2021.txt")#laster inn data
names(dta)#viser navn p? hver av kolonnene i r?datafilen
str(dta)#sjekker om data er kontinuerlige, faktorvariabler mm. 
dta$Lokalitet <- as.factor(dta$Lokalitet)#Lokalitet som faktor
levels(dta$Lokalitet)#Check if lokalitet is correct

## antall honningblomindivider per rute
## lager to figurer per site: en med de rutene som inng?r alle ?r
## en med alle rutene (inkl. nye 2021-ruter)

# opprinnelige ruter
tapply(dta$Ant.pr.0.25m2[dta$All.years==1],list(dta$Year[dta$All.years==1],dta$Lokalitet[dta$All.years==1]),mean,na.rm=T)

# alle ruter
tapply(dta$Ant.pr.0.25m2,list(dta$Year,dta$Lokalitet),mean,na.rm=T)

## Skipstadsand

jpeg(filename="Skipstadsand_antall_2024.jpg",width=35,height=20,units="cm",pointsize=20,res=300)
par(mfrow=c(1,2))
lineplot.CI(Year,Ant.pr.0.25m2,data=dta[dta$Lokalitet=="Skipstadsand" & dta$All.years==1,],xlab="Aar",ylab="Ant. individer per rute",cex.lab=1.3,las=1,ylim=c(0,100))
title("Opprinnelige ruter")
lineplot.CI(Year,Ant.pr.0.25m2,data=dta[dta$Lokalitet=="Skipstadsand" & !dta$All.years=="NA",],xlab="Aar",ylab="Ant. individer per rute",cex.lab=1.3,las=1,ylim=c(0,100))
title("Alle ruter")
par(mfrow=c(1,1))
dev.off()

## Skjellvik

jpeg(filename="Skjellvik_antall_2024.jpg",width=35,height=20,units="cm",pointsize=20,res=300)
par(mfrow=c(1,2))
lineplot.CI(Year,Ant.pr.0.25m2,data=dta[dta$Lokalitet=="Skjellvik" & dta$All.years==1,],xlab="Aar",ylab="Ant. individer per rute",cex.lab=1.3,las=1,ylim=c(0,40))
title("Opprinnelige ruter")
lineplot.CI(Year,Ant.pr.0.25m2,data=dta[dta$Lokalitet=="Skjellvik"& !dta$All.years=="NA",],xlab="Aar",ylab="Ant. individer per rute",cex.lab=1.3,las=1,ylim=c(0,40))
title("Alle ruter")
par(mfrow=c(1,1))
dev.off()

#Teneskjaer

jpeg(filename="Teneskjaer_antall_2024.jpg",width=35,height=20,units="cm",pointsize=20,res=300)
par(mfrow=c(1,2))
lineplot.CI(Year,Ant.pr.0.25m2,data=dta[dta$Lokalitet=="Teneskjaer" & dta$All.years==1,],xlab="Aar",ylab="Ant. individer per rute",cex.lab=1.3,las=1,ylim=c(0,350))
title("Opprinnelige ruter")
lineplot.CI(Year,Ant.pr.0.25m2,data=dta[dta$Lokalitet=="Teneskjaer"& !dta$All.years=="NA",],xlab="Aar",ylab="Ant. individer per rute",cex.lab=1.3,las=1,ylim=c(0,350))
title("Alle ruter")
par(mfrow=c(1,1))
dev.off()

#Filletassen

jpeg(filename="Filletassen_antall_2024.jpg",width=35,height=20,units="cm",pointsize=20,res=300)
par(mfrow=c(1,2))
lineplot.CI(Year,Ant.pr.0.25m2,data=dta[dta$Lokalitet=="Filletassen",],xlab="Aar",ylab="Ant. individer per rute",cex.lab=1.3,las=1,ylim=c(0,10))
lineplot.CI(Year,And.fertile,data=dta[dta$Lokalitet=="Filletassen",],xlab="Aar",ylab="Andel fertile",cex.lab=1.3,las=1,ylim=c(0,0.6))
par(mfrow=c(1,1))
dev.off()

## Ekstra figur: Andel fertile i rutene på Teneskjær
tapply(dta$And.fertile,list(dta$Year,dta$Lokalitet),mean,na.rm=T)


jpeg(filename="Teneskjaer_fertile_2024.jpg",width=35,height=20,units="cm",pointsize=20,res=300)
lineplot.CI(Ar,And.fertile,data=dta[dta$Lokalitet=="Teneskjaer",],xlab="Aar",ylab="Andel fertile",cex.lab=1.3,las=1,ylim=c(0,0.4))
title("Teneskjaer")
dev.off()

jpeg(filename = "Figures/fertile_2024.jpg", width = 35, height = 50, units = "cm", pointsize = 20, res = 300)
par(mfrow = c(4,1))
lineplot.CI(Year, And.fertile, data = dta[dta$Lokalitet=="Skipstadsand",], xlab = "Aar", ylab = "Andel fertile individer", cex.lab = 1.3, las = 1, ylim = c(0,0.5),main="Skipstadsand")
lineplot.CI(Year, And.fertile, data = dta[dta$Lokalitet=="Skjellvik",], xlab = "Aar", ylab = "Andel fertile individer", cex.lab = 1.3, las = 1, ylim = c(0,0.5),main="Skjellvik")
lineplot.CI(Year, And.fertile, data = dta[dta$Lokalitet=="Teneskjaer",], xlab = "Aar", ylab = "Andel fertile individer", cex.lab = 1.3, las = 1, ylim = c(0,0.5),main="Teneskjaer")
lineplot.CI(Year, And.fertile, data = dta[dta$Lokalitet=="Filletassen",], xlab = "Aar", ylab = "Andel fertile individer", cex.lab = 1.3, las = 1, ylim = c(0,0.5),main="Filletassen")
par(mfrow = c(1,1))
dev.off()

# ekstra figur - jorddybde
jpeg(filename = "Figures/soil_2023.jpg", width = 35, height = 20, units = "cm", pointsize = 20, res = 300)
bargraph.CI(Lokalitet,Jorddybde,data=dta[dta$Year==2023,],cex.lab=1.3,las=1,ylim=c(0,30),ylab="Jorddybde (cm)")
abline(h=0)
dev.off()




