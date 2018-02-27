#Packages
library(reshape2);library(raster);library(sp);library(ggplot2);library(gridExtra);library(data.table);library(plyr);library(dplyr);library(POT)
########################################################################################################################################################################################################################################################################################################################################################################################
dev.off()
remove (list=ls())



                                                  #####################################################
                                                  #    Data ROCKYFOR3D V - Simulations JANUARY 2018   #
                                                  #                   FINAL FILE                      #
                                                  #####################################################


# First inputs - See Crolles-IV for more details
#######################################################################################################################################################################################################################################################################################################################################################################################
##Cliffs
rS                = raster("C:/Users/manon.farvacque/Documents/_CROLLES-V/_Simu&Results/_Urban-Risk-Map/Final-InputsRF3D/blshape.asc") 
dataS             = rasterToPoints(rS); remove(rS) 
##DEM Resolution
Resolution_DEM    = 5; NbSimu_Cell = 10000; NbTot_CliffCells = length(dataS[,1]); Cliff_TOT_Surface = 257518.35 #(m2)
##Volume classes 
VolumesClass      = seq(from = 1.5, to = 19.5, by = 1); breaks = 19
##Price per square metre [e/m2]
Price_m2          = 3000
##Frequency determination - See Crolles-IV for more details
DataField         = fread("C:/Users/manon.farvacque/Documents/CROLLES-IIBIS/FrequencyDetermination_Crolles/CubageArdillais.txt"); colnames(DataField) = c("Volume")
DataField_Vol     = DataField$Volume[c(1:29)]
Field_Cliff       = 11.5; Window = 100 
threshold         = 1; mle = fitgpd(DataField_Vol, threshold = threshold, est ="mle")
Scale             = as.numeric(mle$param["scale"]); Shape = as.numeric(mle$param["shape"]); Lambda = (length(which(DataField_Vol>threshold))/Window)*(1/Field_Cliff)  
########################################################################################################################################################################################################################################################################################################################################################################################


####
########
################
########################################################################################################################################################################################################################################################################################################################################################################################
# URBAN 


##RF3D simu results
dataRF3D_URBAN           = fread("C:/Users/manon.farvacque/Documents/_CROLLES-V/_Simu&Results/_Urban-Risk-Map/_Results/resultatPy")
colnames(dataRF3D_URBAN) = c("NumSimu","BatNumber","RowBat","ColumnBat","RowRock","ColumnRock","Ebrute","Volume")
dataRF3D_URBAN           = as.data.frame(dataRF3D_URBAN)

##District 
rB    = raster("C:/Users/manon.farvacque/Documents/_CROLLES-V/_Simu&Results/_Urban-Risk-Map/Final-InputsRF3D/net_number.asc")
dataB = rasterToPoints(rB); colnames(dataB) = c("XBat","YBat","NetNumber")

##Nb of impacts
ImpactNumber  = length(dataRF3D_URBAN$NumSimu)

##List of the houses reached by rockfalls (ID number)  
HousesReached = data.frame(unique(dataRF3D_URBAN$BatNumber)); colnames(HousesReached) = c("NetNumber")

##Quantitative Risk Assessment QRA - MEAN MODEL 
# Houses reached loop
QRA_Table = matrix(0, nrow=length(HousesReached$NetNumber), ncol=20)
for (k in 1:length(HousesReached$NetNumber)){
     QRA_Table[k,1] = HousesReached$NetNumber[k]
     # Rockyfor3D results extract for the specifics k (one building, all introduced with the loop)
     dataRF3D_URBAN_specific_k = subset(dataRF3D_URBAN, dataRF3D_URBAN$BatNumber == HousesReached$NetNumber[k])
     # Degree of loss - According to Agliardi F. et al (2009) 
     dataRF3D_URBAN_specific_k[,9] = 1-(1.358/(1 + exp((((dataRF3D_URBAN_specific_k$Ebrute*10^3)-129000)/120300))))
     dataRF3D_URBAN_specific_k[,9] = replace(dataRF3D_URBAN_specific_k[,9], dataRF3D_URBAN_specific_k[,9] < 0, 0)
     colnames(dataRF3D_URBAN_specific_k) = c("NumSimu","BatNumber","RowBat","ColumnBat","RowRock","ColumnRock","Ebrute","Volume","Damage_Agliardi")
     # Bulding k price estimation
     PriceBat = Price_m2*length(subset(dataB,dataB[,3]==HousesReached$NetNumber[k])[,3])*Resolution_DEM^2
     # Volume classification
     dataRF3D_URBAN_specific_k[,8]    = replace(dataRF3D_URBAN_specific_k[,8],dataRF3D_URBAN_specific_k[,8] == 1, 1.01) # Ceiling not efficient if vol == 1 
     dataRF3D_URBAN_specific_k$Volume = ceiling(dataRF3D_URBAN_specific_k$Volume)-0.5 
  
  for (i in 1:length(VolumesClass)){
       # GPD 
       GPDLaw           = Lambda*(((1+(Shape*(((VolumesClass[i]-0.5)-threshold)/Scale)))^(-1/Shape))-((1+(Shape*(((VolumesClass[i]+0.5)-threshold)/Scale)))^(-1/Shape)))
       # Risk by volumes class
       Int5             = subset(dataRF3D_URBAN_specific_k, dataRF3D_URBAN_specific_k$Volume == VolumesClass[i])
       # Number of impacts in building k for class volume i 
       NbImpact         = length(Int5[,1])
       # Damage mean
       DamageMean       = mean(Int5[,9]); DamageMean[is.na(DamageMean)] = 0 
       # Number of events per year in class volume i and for Crolles area (hm2)
       Frequency        = GPDLaw*Cliff_TOT_Surface*10^-4
       # Impact probability 
       Pz               = NbImpact/((NbSimu_Cell/breaks)*NbTot_CliffCells)
       # QRA calculation
       QRA              = Frequency*Pz*DamageMean*PriceBat
       QRA_Table[k,i+1] = QRA 
       # Result (compilation on all volume classes)
       remove(GPDLaw,Int5,NbImpact,DamageMean,Frequency,Pz,QRA)}
  
remove(dataRF3D_URBAN_specific_k,PriceBat)}
colnames(QRA_Table) = c("Bat","1.5","2.5","3.5","4.5","5.5","6.5","7.5","8.5","9.5","10.5","11.5","12.5","13.5","14.5","15.5","16.5","17.5","18.5","19.5")

##Results
ImpactNumber
length(HousesReached$NetNumber)
sum(QRA_Table[,c(2:20)])
########################################################################################################################################################################################################################################################################################################################################################################################


####
########
################
########################################################################################################################################################################################################################################################################################################################################################################################
# INDIVIDUAL

##Area 
rB              = raster("C:/Users/manon.farvacque/Documents/_CROLLES-V/_Simu&Results/_Individual-Risk-Map/Final-InputsRF3D/net_number.asc")
dataB           = rasterToPoints(rB); remove(rB)
colnames(dataB) = c("XBat","YBat","NetNumber"); dataB = as.data.frame(dataB)

##RF3D simu results - ALL
     
      #dataINDIVIDUAL           = fread("C:/Users/manon.farvacque/Documents/_CROLLES-V/_Simu&Results/_Individual-Risk-Map/_ResultsRF3D/resultatPy")
      #colnames(dataINDIVIDUAL) = c("IDcell","Energy","Volume"); A = nrow(dataINDIVIDUAL)

##Split results table to small tables 

      #z = 11; i = 1; j = A/z
      #  for (x in 1:z){
      #       spt = dataINDIVIDUAL[c(i:j),]
      #       write.table(spt,file = paste0("dataINDIVIDUAL", x, ".csv"))
      #       i = j + 1; j = i + ((A/z)-1)
      #  }

##3D Matrix QRA construction - PART 1 / Performed individually on each subset, as:

      #RF3D results - specific to subset x 
      #A              = "C:/Users/manon.farvacque/Documents/_CROLLES-V/_Simu&Results/_Individual-Risk-Map/dataINDIVIDUALx.csv"
      #dataINDIVIDUAL = fread(A, select = c(2:4))

      ##QRA MATRIX 3D - Construction part 1
      #QRA_3Dmatrix = array(0,dim=c(nrow(dataB),breaks,3))                # dimension/layer 1 : will used for frequency
                                                                          # dimension/layer 2 : Nr passage
                                                                          # dimension/layer 3 : sum(Damage)
      #pb <- winProgressBar(title="progress bar",min = 0,max=nrow(dataINDIVIDUAL),width = 300)
      #  for (i in 1:nrow(dataINDIVIDUAL)){
      #       row    = as.numeric(dataINDIVIDUAL[i,1])                    # row corresponds to IDcell 
      #       column = floor(as.numeric(dataINDIVIDUAL[i,3]))             # column corresponds to the volume class
      #       if(column == 20){column = 19}
      #       QRA_3Dmatrix[row,column,2] = QRA_3Dmatrix[row,column,2] + 1 # Nr passage, initially zero, then +1 
      #       QRA_3Dmatrix[row,column,3] = QRA_3Dmatrix[row,column,3] + as.numeric(1-(1.358/(1 + exp((((dataINDIVIDUAL[i,2]*10^3)-129000)/120300)))))
                                                                          # sum damage, beforehand energy to damage by Agliardi vulnerability curve
      #       remove(row,column)
      #       Sys.sleep(0.1)
      #       setWinProgressBar(pb, i, title=paste(round(i/nrow(dataINDIVIDUAL)*100, 0),"% done"))
      #  }
      #saveRDS(QRA_3Dmatrix, file="QRA-Partx.rds"); close(pb)

##Energy distribution extraction 

      #int00      = "C:/Users/manon.farvacque/Documents/_CROLLES-V/_Simu&Results/_Individual-Risk-Map/_ResultsRF3D/ResultatPy-VolumeCL/ResultsVCL-01.csv"
      #dataVCL    = fread(int00); dataVCL = as.data.frame(dataVCL)
       
      ###Energy distribution/cell 
      #DistEnergy = list()
      #int1       = unique(dataVCL[,1])
      #   for (k in 1:length(int1)){
      #        DistEnergy[[k]] = subset(dataVCL[,2], dataVCL[,1] == int1[k])
      #        print(k)
      #   }
      #remove(dataVCL)
      
      ##Compute maximum length
      #maxL     <- max(sapply(DistEnergy, length))
      ##Add NA values to list elements
      #A        <- lapply(DistEnergy, function(v) { c(v, rep(NA, maxL-length(v)))})
      ##Cbind
      #B        <- do.call(cbind, A)     # Distribution in a same table 
      #EDTable  <- rbind(int1, B)        # First line : ID Cell, then Energy distribution
      ##Saving
      #saveRDS(EDTable, file="DistributionEnergy-VolumeCL-01.rds")
     
##3D Matrix QRA construction - Mean Damage
      
# Merging matrix 
Part01 = readRDS("QRA-Part1.rds"); Part02 = readRDS("QRA-Part2.rds"); Part03 = readRDS("QRA-Part3.rds");Part04 = readRDS("QRA-Part4.rds")
Part05 = readRDS("QRA-Part5.rds"); Part06 = readRDS("QRA-Part6.rds"); Part07 = readRDS("QRA-Part7.rds");Part08 = readRDS("QRA-Part8.rds")
Part09 = readRDS("QRA-Part9.rds"); Part10 = readRDS("QRA-Part10.rds"); Part11 = readRDS("QRA-Part11.rds")

# QRA construction 
QRA = array(0,dim=c(nrow(dataB),breaks,7))  

        # Frequency layer 
for (i in 1:19){
     QRA[,i,1] = Cliff_TOT_Surface*10^-4*Lambda*(((1+(Shape*(((i)-threshold)/Scale)))^(-1/Shape))-((1+(Shape*(((i+1)-threshold)/Scale)))^(-1/Shape)))
     } 
        # Nr passage layer - ALL
QRA[,,2] = Part01[,,2]+Part02[,,2]+Part03[,,2]+Part04[,,2]+Part05[,,2]+Part06[,,2]+Part07[,,2]+Part08[,,2]+Part09[,,2]+Part10[,,2]+Part11[,,2]
        # Pz IDcell/V(CL) 
QRA[,,3] = QRA[,,2]/((NbSimu_Cell/breaks)*NbTot_CliffCells)
        # Sum(Damage)      - ALL
QRA[,,4] = Part01[,,3]+Part02[,,3]+Part03[,,3]+Part04[,,3]+Part05[,,3]+Part06[,,3]+Part07[,,3]+Part08[,,3]+Part09[,,3]+Part10[,,3]+Part11[,,3]
        # Damage Mean
QRA[,,5] = QRA[,,4]/QRA[,,2]
        # ID value - Surface 25 m2 to have m2 detroyed/yr (the future of urban planning)
QRA[,,6] = Resolution_DEM^2
        # QRA Results 
QRA[,,7] = QRA[,,1]*QRA[,,3]*QRA[,,5]*QRA[,,6]

Risk = data.frame()
for (j in 1:nrow(QRA)){
     Risk[j,1] = j
     Risk[j,2] = sum(QRA[j,,7], na.rm = TRUE)}
colnames(Risk) = c("NetNumber","Risk")

#ALL = full_join(dataB, Risk, by ="NetNumber")
#r01 = rasterFromXYZ(ALL[,c(1,2,4)])
#writeRaster(r01,'RISK.asc',format='ascii'); remove(ALL)

##3D Matrix QRA construction - Q90 Damage

A = readRDS("C:/Users/manon.farvacque/Documents/_CROLLES-V/_Simu&Results/_Individual-Risk-Map/DistributionEnergy/DistributionEnergy-VolumeCL-08.rds")
B = as.data.frame(A[,1])
C = A[-1,1890]
D = C[!is.na(C)]
D = 1-(1.358/(1 + exp((((D*10^3)-129000)/120300))))
quantile(D, 0.95)
