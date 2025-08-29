#setwd()

####################################################
### 1. Write function to reformat morphological data
reformat<-function(data){
  
  ## Get first individual per species
  ind1<-data.frame(t(data[c(3:8,12:17,153),])) #pull out specimen id, measurements, and measurer
  ind1<-ind1[-1,] #drop column and row names
  colnames(ind1)<-cols #rename columns
  rownames(ind1)<-NULL
  ind1<-cbind(Species,ind1) #add species as a new column
  ind1<-subset(ind1,ind1$SpecimenID!="") # drop any species without 1 specimen
  
  ## Get second individual per species
  ind2<-data.frame(t(data[c(18:23,27:32,153),])) #pull out specimen id, measurements, and measurer
  ind2<-ind2[-1,] #drop column and row names
  colnames(ind2)<-cols #rename columns
  rownames(ind2)<-NULL
  ind2<-cbind(Species,ind2) #add species as a new column
  ind2<-subset(ind2,ind2$SpecimenID!="") #drop any species without 2 specimens
  #View(ind2)
  
  ## Get third individual per species
  ind3<-data.frame(t(data[c(33:38,42:47,153),])) #pull out specimen id, measurements, and measurer
  ind3<-ind3[-1,] #drop column and row names
  colnames(ind3)<-cols #rename columns
  rownames(ind3)<-NULL
  ind3<-cbind(Species,ind3) #add species as a new column
  ind3<-subset(ind3,ind3$SpecimenID!="") #drop any species without 3 specimens
  #View(ind3)
  
  ## Get fourth individual per species
  ind4<-data.frame(t(data[c(48:53,57:62,153),])) #pull out specimen id, measurements, and measurer
  ind4<-ind4[-1,] #drop column and row names
  colnames(ind4)<-cols #rename columns
  rownames(ind4)<-NULL
  ind4<-cbind(Species,ind4) #add species as a new column
  ind4<-subset(ind4,ind4$SpecimenID!="") #drop any species without 4 specimens
  #View(ind4)
  
  ## Get fifth individual per species
  ind5<-data.frame(t(data[c(63:68,72:77,153),])) #pull out specimen id, measurements, and measurer
  ind5<-ind5[-1,] #drop column and row names
  colnames(ind5)<-cols #rename columns
  rownames(ind5)<-NULL
  ind5<-cbind(Species,ind5) #add species as a new column
  ind5<-subset(ind5,ind5$SpecimenID!="") #drop any species without 5 specimens
  #View(ind5)
  
  ## Get sixth individual per species
  ind6<-data.frame(t(data[c(78:83,87:92,153),])) #pull out specimen id, measurements, and measurer
  ind6<-ind6[-1,] #drop column and row names
  colnames(ind6)<-cols #rename columns
  rownames(ind6)<-NULL
  ind6<-cbind(Species,ind6) #add species as a new column
  ind6<-subset(ind6,ind6$SpecimenID!="") #drop any species without 6 specimens
  #View(ind6)
  
  ## Get seventh individual per species
  ind7<-data.frame(t(data[c(93:98,102:107,153),])) #pull out specimen id, measurements, and measurer
  ind7<-ind7[-1,] #drop column and row names
  colnames(ind7)<-cols #rename columns
  rownames(ind7)<-NULL
  ind7<-cbind(Species,ind7) #add species as a new column
  ind7<-subset(ind7,ind7$SpecimenID!="") #drop any species without 7 specimens
  #View(ind7)
  
  ## Get eighth individual per species
  ind8<-data.frame(t(data[c(108:113,117:122,153),])) #pull out specimen id, measurements, and measurer
  ind8<-ind8[-1,] #drop column and row names
  colnames(ind8)<-cols #rename columns
  rownames(ind8)<-NULL
  ind8<-cbind(Species,ind8) #add species as a new column
  ind8<-subset(ind8,ind8$SpecimenID!="") #drop any species without 8 specimens
  #View(ind8)
  
  ## Get ninth individual per species
  ind9<-data.frame(t(data[c(123:128,132:137,153),])) #pull out specimen id, measurements, and measurer
  ind9<-ind9[-1,] #drop column and row names
  colnames(ind9)<-cols #rename columns
  rownames(ind9)<-NULL
  ind9<-cbind(Species,ind9) #add species as a new column
  ind9<-subset(ind9,ind9$SpecimenID!="") #drop any species without 9 specimens
  #View(ind9)
  
  ## Get tenth individual per species
  ind10<-data.frame(t(data[c(138:143,147:152,153),])) #pull out specimen id, measurements, and measurer
  ind10<-ind10[-1,] #drop column and row names
  colnames(ind10)<-cols #rename columns
  rownames(ind10)<-NULL
  ind10<-cbind(Species,ind10) #add species as a new column
  ind10<-subset(ind10,ind10$SpecimenID!="") #drop any species without 10 specimens
  #View(ind10)
  
  ## Reformat to tidy dataset
  new.data<-rbind(ind1,ind2,ind3,ind4,ind5,
                  ind6,ind7,ind8,ind9,ind10) #combine all data together
  new.data[,4:13]<-apply(new.data[,4:13],2,as.numeric) #make measurements numeric
  new.data
}

######################################
### 2. Load and reformat hindlimb data
data <- read.csv("Blackbird Hindlimb Dataset.csv")
data <- subset(data, select = IOC.Name:Xanthopsar.flavus)

cols<-c("SpecimenID",as.character(data$IOC.Name[c(4:8,12:17,153)])) # get new column names from "IOC.Name" column
Species<-colnames(data)[2:51] # change this number with # of species
Species<-gsub(".","_",Species,fixed=T) # replace periods with underscores in names to match phylogeny tip labels

tidy <- reformat(data) # reformat specimen dataset!

tidy <- tidy[order(tidy$Species,tidy$Sex),] # reorder dataset by species then sex

write.csv(tidy,file="Blackbirdhindlimb_tidy.csv",row.names=F)





###########################################
### 3. Add ecological data to hindlimb data
eco <- read.csv("Blackbird Ecology Data.csv")

eco$Species<-gsub(" ","_",eco$Species,fixed=T) # replace spaces with underscores in names to match phylogeny tip labels

eco2 <- eco[,c(1,3,5,6,7,9,11)] # subset data to not have descriptive text columns

eco2$Parasite <- eco2$Nesting.Score # create new column to represent parasite, and replace values in nesting score
eco2$Parasite <- replace(eco2$Parasite, eco2$Parasite != "Parasite", 0)
eco2$Parasite <- replace(eco2$Parasite, eco2$Parasite == "Parasite", 1)
eco2$Nesting.Score <- replace(eco2$Nesting.Score, eco2$Nesting.Score == "Parasite", NA)


all <- merge(tidy,eco2) # merge morpho and ecological data

write.csv(all,file = "AllBlackbirdData.csv",row.names=F)
