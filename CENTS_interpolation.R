#NO3 interpolation for a dataset similar to CENTS results
#Basic structure: Blocks -> Plots -> Subplot -> Subsets -> Ugly NO3 time series
#Need to have EVACROP simulated precolation (see other files in repository)

library(dplyr)
library(tidyr)
library(lubridate)

#We begin with a data frame "nit_all" containing ALL the NO3 measurements from CENTS Foulumgaard. We trim the dates and subset it to a single Tillage treatment in a single Rotation treatment
nit_all <- read.csv(file="") #Read the nit_all.csv
nit_inter <- nit_all[order(nit_all$Date),] %>% filter(Date>="2002-01-01",Rotation=="R5",Tillage=="P") %>% group_by(Date,Block,Cover,Depth) %>% summarise(nitN=mean(NitrateN,na.rm=T)) %>% data.frame() %>% na.omit() #Produces an ordered data frame with mean NitrateN measurements in R5
row.names(nit_inter) <- c(1:nrow(nit_inter))
nit_inter$Cover <- factor(nit_inter$Cover)
nit_inter$Depth <- factor(nit_inter$Depth)
nit_inter$Block <- factor(nit_inter$Block)

#Difference between consecutive measurements for the interpolation loop
#Let's see if we can run all covers and both depths in one go
#Let's see if we can do Blocks separately in one go
bruh <- data.frame()
for(l in levels(nit_inter$Block)){
  for(k in levels(nit_inter$Depth)){
    for(j in levels(nit_inter$Cover)){
      huh <- filter(nit_inter,Depth==k,Cover==j,Block==l)
      huh$Cdiff <- NA
      for (i in 1:nrow(huh)) {
        huh$Cdiff[i] <- huh$nitN[i+1]-huh$nitN[i]
      }
      huh$Cdiff[i] <- 0
      bruh <- rbind(bruh,huh)
    }
  }
}
nit_inter <- full_join(nit_inter,select(bruh,-nitN),by=c("Date","Cover","Block","Depth"))

#Load the EVACROP percolation data
EVA_Foulum_02_12 <- read.table(file="",sep="")
EVA_Foulum_12_22 <- read.table(file="",sep="")
EVA_Foulum <- rbind(EVA_Foulum_02_12,EVA_Foulum_12_22)
names(EVA_Foulum) <- c("Date","mTemp","Prec","pEvap","Irrigation","aEvap","Bottom_Drainage")
EVA_Foulum$Date <- date(ymd(EVA_Foulum$Date))
rm(EVA_Foulum_02_12,EVA_Foulum_12_22)

Interpol <- data.frame()
bruh <- data.frame()

for(l in levels(nit_inter$Block)){
  for(k in levels(nit_inter$Depth)){
    for(j in levels(nit_inter$Cover)){
      wha <- filter(nit_inter,Depth==k,Cover==j,Block==l)
      huh <- filter(EVA_Foulum,Date>=min(wha$Date)) %>% select(Date,"Bottom_Drainage") #Grab drainage from EVACROP
      bruh <- left_join(huh,wha,by="Date")
      bruh$Depth <- k
      bruh$Cover <- j
      bruh$Block <- l
      Interpol <- rbind(Interpol,bruh)
    }
  }
}
Interpol$check <- TRUE #A check mark for when there is a measurement
Interpol$check[is.na(Interpol$nitN)==T] <- FALSE #Interpolation periods (i.e. with no measurement) are maked FALSE
Interpol$Depth <- factor(Interpol$Depth)
Interpol$Cover <- factor(Interpol$Cover)
Interpol$Block <- factor(Interpol$Block)
rm(nit_inter,EVA_Foulum)

#Adding the sum of simulated drainage for inter_measured period for the interpolation loop. Ugliest code ever written, but yolo.
D_index <- 0 #First we make an index of inter-measurement periods
Interpol$D_index <- 0 #Index keeping track of discrete interpolation periods (i.e. gaps between consective measurements)
bruh <- data.frame()

for(l in levels(Interpol$Block)){
  for(k in levels(Interpol$Depth)){
    for(j in levels(Interpol$Cover)){
      huh <- filter(Interpol,Block==l,Depth==k,Cover==j)
      for (i in 1:nrow(huh)){
        if (huh$check[i]==TRUE) {
          huh$D_index[i] <- NA #Remember not to include the sim drainage in the days there is nitrate measurements
          D_index <- D_index + 1 #If there is a measurement, count one more gap
        }
        else huh$D_index[i] <- D_index
      }
      bruh <- rbind(bruh,huh)
      D_index <- 0
    }
  }
}

Interpol$D_index <- bruh$D_index

Sample_tots <- Interpol %>% group_by(Block,Depth,Cover,D_index) %>% summarize(D_tot=sum(Bottom_Drainage)) %>% data.frame() %>% na.omit()
Interpol$D_tot <- 0
Interpol <- rows_update(Interpol,Sample_tots,by=c("Block","Depth","Cover","D_index")) #Ugly, but it works
rm(Sample_tots,D_index)

#Alrighty. More looooops! Yeah yeah, it's inefficient, I don't care.
Interpol$fd <- 0 #Infamous weighing factor
Interpol$Ci <- 0 #Interpolated daily nitrate concentration

bruh <- data.frame()
for(l in levels(Interpol$Block)){

for(k in levels(Interpol$Depth)){
  for(j in levels(Interpol$Cover)){
    huh <- filter(Interpol,Block==l,Depth==k,Cover==j)
    Dbs <- huh$Bottom_Drainage[huh$check==TRUE] #A vector of all the final sim drainages of interpolation periods. Please kill me.
    #These keep track of the last measured value through the loop
    Da <- 0 #Initial sim drainage in interpolation period
    Db <- 0 #Final sim drainage in interpolation period
    Ca <- 0 #Measured concentration at the start of interpolation period
    Cd <- 0 #Difference between consecutive measurements for the interpolation loop
    fd <- 1 #Infamous weighing factor
    Dsum <- 0 #This one is the partial sum of simulated drainage used in the loop
    for (i in 1:nrow(huh)) {
      if(is.na(huh$nitN[i])==F){ #Keeping track of the last measured value and last difference between measured values
        Ca <- huh$nitN[i]
        Cd <- huh$Cdiff[i]
        Da <- huh$Bottom_Drainage[i] #Initial sim drainage in interpolation period
        Dsum <- 0 #Reset the partial sum of sim drainage if we have a measurement
        huh$Ci[i] <- huh$nitN[i]
      } else{
      Dsum <- Dsum+huh$Bottom_Drainage[i] #Partial sum of sim drainage. Resets above if there is a nitrate measurement.
      Db <- Dbs[huh$D_index[i]+1]
      fd <- (0.5*Da + Dsum + 0.5*huh$Bottom_Drainage[i])/(0.5*Da + huh$D_tot[i] + 0.5*Db) #Calculate the infamous weighing factor
      huh$fd[i] <- fd #Let's keep track of it for now
      if(is.nan(huh$fd[i])==T){huh$Ci[i] <- Ca} else huh$Ci[i] <- Ca+Cd*fd #Yes, we're nesting ifs. I beg for the sweet release of death.
      }
    }
    bruh <- rbind(bruh,huh)
  }
}

}

Interpol <- bruh

#And finally. Interpolated N-leaching
Interpol$Nleach <- Interpol$Bottom_Drainage*Interpol$Ci
#Remember: 1 mm = 1 L/m1 = 10,000 L/ha. An nitrate conc. is in units of mg N/L = 10^-6 kg/L. 
#So Leaching is in units of 10^4*10^-6 kg N/ha = 1/100 kg N/ha.
#Therefore:
Interpol$Nleach <- Interpol$Nleach/100 #Is in units of kg N/ha
Interpol <- filter(Interpol,D_index!=0|is.na(D_index)==T) #Trim the edges. There seem to be no edges tho.

#Finally. Let's write it down, lest we forget it.
write.csv(select(Interpol,Date,Depth,Cover,Block,nitN,Bottom_Drainage,Ci,Nleach),file="../BUP-KLIMON/WWheat_1v2_m/Interpol_R5_P.csv",quote=F,row.names=F)

#Look at what you've done
library(ggplot2)
ggplot(Interpol,aes(x=Date,y=Ci))+
  geom_line(aes(linetype=Block,col=Depth))+
  geom_point(aes(y=nitN,shape=Block,col=Depth))+
  facet_wrap(~Cover,nrow=3)
