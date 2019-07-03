library(MASS)
bt_frost <- read.delim('bt_frost2.txt')
bt_frost$Warm <- bt_frost$bt0510
bt_frost$Cold <- bt_frost$tmin
bt_frost$frostbin <- log(1/(365.25/bt_frost$Frost50 - 1))
frostmodelmax <- lm(formula = Frost50 ~ Warm + I(Warm^2) + poly(Cold, 2) + Warm:Cold, 
                    data = bt_frost)
frostmodelmin <- lm(Frost50 ~  1, data = bt_frost)

dontdo<-function(x){
  step(frostmodelmax, scope=list(upper = frostmodelmax,
                                 lower= frostmodelmin), direction="both", trace=0)
  
  
  step(frostmodelmin, scope=list(upper = frostmodelmax,
                                 lower= frostmodelmin), direction="forward")
stepmodelfrost <- lm(formula = Frost50 ~ I(pmax(Cold,(Warm-18)/2)), data = bt_frost)
summary(stepmodelfrost)
}

stepmodelfrost <- lm(formula = Frost50 ~ Warm + poly(Cold, 2) + Warm:Cold, data = bt_frost)
summary(stepmodelfrost)


si <- read.delim('SimpleExport.txt')
si <- subset(si, !is.na(Cold) & !is.na(pH50) & AWC150 >0 & WOOD_SPGR_GREENVOL_DRYWT >0 & elev < 3600)
si$frost50 <- predict(stepmodelfrost, si)
formodel <- si[,6:ncol(si)]
######################################
dontdo<-function(x){
  newmodel<- stepAIC(simodel, scope=list(upper= ~ lnppt+Cold+Warm+AWC150+WaterTable+om150+Wet+pH50+Carya+Fagus+Juglans+Liriodendron+Populus+Prunus+Tilia+Ulmus+ABBA+ACRU+ACSA3+FRAM2+FRPE+JUVI+PIST+PIRE+PIBA2+QUAL+QURU+QUVE+PODE+THOC2, lower= ~1))
}
################################################

vartab <- cor(subset(formodel, !is.na(Cold)&!is.na(pH50)))
#newmodel<- stepAIC(simodel, scope=list(upper= ~ ., lower= ~1))
dontdo<-function(x){  #test dataset for affects of a truly random variable
  randomnumbers <- rnorm(nrow(formodel), 1, 1)
  
  formodel2 <- cbind(formodel[,c("LNSI","lnppt","Cold","Warm","AWC150",'Bhs',"WaterTable",'Inceptisols','Liriodendron')], randomnumbers)
}
dontdo<-function(x){ #ways to stepwise add or subtract variables to optimize model
  simodelmax <- lm(LNSI ~ lnppt + Cold + Warm + I(Warm^0.5) + AWC150 + poly(WaterTable,2) + om150 + Wet + poly(pH50,2) + clay150 + sand150 + carbdepth + Wet*WaterTable + Carya + Fagus + Juglans + Liriodendron + Populus + Prunus + Tilia + Ulmus + ABBA + ACRU + ACSA3 + FRAM2 + FRPE + JUVI + PIST + PIRE + PIBA2 + QUAL + QURU + QUVE + PODE + THOC2 + Pinus*WaterTable + Pinus*AWC150 + Pinus*pH50 + Quercus*WaterTable + Quercus*AWC150 + Quercus*pH50 + WaterTable*AWC150 + WaterTable*lnppt + lnppt*Warm + Bhs + Bhs*AWC150 + Bhs*WaterTable + WOOD_SPGR_GREENVOL_DRYWT + Shade + Ultisols + Alfisols + Inceptisols + Mollisols + Vertisols + Andisols + Pinus + Quercus + Acer + Fraxinus + Picea + Liquidambar + Abies + Betula + Juniperus + Larix + Celtis + Thuja + Platanus + Robinia + Nyssa + Pseudotsuga + Tsuga + Sequoia + frost50, data = formodel)
  
  summary(simodelmax)
  simodelmin <- lm(LNSI ~  1, data = formodel)
  step(simodelmin, scope=list(upper = simodelmax,
                              lower= simodelmin), direction="forward")
  step(simodelmax, scope=list(upper = simodelmax,
                              lower= simodelmin), direction="backward")
  step(simodelmax, scope=list(upper = simodelmax,
                              lower= simodelmin), direction="both", trace=0)
}




##### trained the model with a frostfree day model

stepforward <-   
  lm(formula = lm(formula = LNSI ~ lnppt + frost50 + Pseudotsuga + Liriodendron + 
                    AWC150 + Juniperus + Populus + om150 + Warm + I(Warm^0.5) + 
                    Wet + Cold + poly(pH50, 2) + QUAL + Pinus + THOC2 + clay150 + 
                    FRAM2 + Shade + Picea + ABBA + WOOD_SPGR_GREENVOL_DRYWT + 
                    Prunus + QUVE + QURU + sand150 + carbdepth + Acer + PIRE + 
                    PIBA2 + Ultisols + poly(WaterTable, 2) + PIST + Betula + 
                    Quercus + Robinia + Liquidambar + Juglans + JUVI + Tilia + 
                    Abies + Larix + Platanus + Nyssa + Carya + Vertisols + Mollisols + 
                    Alfisols + Fagus + Fraxinus + PODE + Inceptisols + ACRU + 
                    ACSA3 + AWC150:Pinus + lnppt:Warm, data = formodel), data = formodel)
summary(stepforward)
##### trained the model with a frostfree day model
stepbackward <-    
  lm(formula = lm(formula = LNSI ~ lnppt + Cold + Warm + I(Warm^0.5) + AWC150 + 
                    poly(WaterTable, 2) + om150 + Wet + poly(pH50, 2) + clay150 + 
                    sand150 + carbdepth + Carya + Fagus + Juglans + 
                    Liriodendron + Populus + Prunus + Tilia + ABBA + ACRU + ACSA3 + 
                    FRAM2 + FRPE + JUVI + PIRE + PIBA2 + QUAL + QURU + QUVE + 
                    PODE + THOC2 + Pinus + Quercus + Bhs + WOOD_SPGR_GREENVOL_DRYWT + 
                    Shade + Ultisols + Mollisols + Vertisols + Andisols + Acer + 
                    Picea + Liquidambar + Abies + Betula + Juniperus + Larix + 
                    Platanus + Robinia + Nyssa + Pseudotsuga + frost50 + Wet:WaterTable + 
                    WaterTable:Pinus + AWC150:Pinus + Pinus:pH50 + WaterTable:Quercus + 
                    AWC150:Quercus + pH50:Quercus + AWC150:WaterTable + lnppt:WaterTable + 
                    lnppt:Warm + AWC150:Bhs + WaterTable:Bhs, data = formodel), data = formodel)

summary(stepbackward)

species <- unique(si[,c("SYMBOL","Wet","WOOD_SPGR_GREENVOL_DRYWT", "Shade",
                        "Carya", "Fagus", "Juglans", "Liriodendron", "Populus", "Prunus", "Tilia", "Ulmus", 
                        "ABBA", "ACRU", "ACSA3", "FRAM2", "FRPE", "JUVI", "PIST", "PIRE", "PIBA2", "QUAL", 
                        "QURU", "QUVE", "PODE", "THOC2",  "Pinus", "Quercus", "Acer", "Fraxinus", "Picea", 
                        "Liquidambar", "Abies", "Betula", "Juniperus", "Larix", "Celtis", "Thuja", "Platanus", 
                        "Robinia", "Nyssa", "Pseudotsuga", "Tsuga", "Sequoia")]
)

FIPS <- aggregate(si[,c("lnppt","Cold", "Warm","frost50")], by = list(si$fips, si$elev), FUN = 'mean')
colnames(FIPS) <- c("fips","elev","lnppt", "Cold", "Warm","frost50")

Series <- unique(si[,c("series", "AWC150", "om150", "pH50", "clay150", "sand150","WaterTable",
                       "carbdepth", "Bhs", "Ultisols", "Alfisols", "Inceptisols", "Mollisols", "Vertisols", "Andisols")])


################### graph the frost free index -------
library(ggplot2)
library(RColorBrewer)
forplot1 <- aggregate(bt_frost[,c("Frost50")], by=list(bt_frost$Cold, bt_frost$Warm), FUN='mean')
colnames(forplot1) <- c('Cold','Warm','Frost50')
forplot1 <- unique(forplot1)

ggplot(aes(x=Cold, y=Warm), data = forplot1)+
  geom_point(aes(color=Frost50))+
  scale_colour_gradient2(low='green', high='red', mid = 'white', midpoint = 150)  

forplot2 <- aggregate(si[,c("frost50")], by=list(si$Cold, si$Warm), FUN='mean')
colnames(forplot2) <- c('Cold','Warm','frost50')
forplot2 <- unique(forplot2)

ggplot(aes(x=Cold, y=Warm), data = forplot2)+
  geom_point(aes(color=frost50))+
  scale_colour_gradient2(low='green', high='red', mid = 'white', midpoint = 150)  

ggplot(aes(x=(pH50), y=exp(LNSI)/0.3048), data = si[si$SYMBOL %in% 'PIST',])+#, data = si[si$Wet == 4,]
  geom_point(aes(), size= 0.1)+
  geom_smooth(method = 'loess', color = 'red')+
  geom_smooth(method = 'lm', color = 'blue')+
  scale_y_continuous(breaks = c(0:5)*25)+
  coord_cartesian(ylim = c(0, 150))
  

ggplot(aes(x=SYMBOL, y=exp(LNSI)/0.3048), data = si)+#, data = si[si$Wet == 4,]
  geom_boxplot()+
  scale_y_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,200))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  
###### select param

SelectSeries <- c('TAWAS','LUPTON', 'KALKASKA', 'GRAYLING', 'HOUGHTON', 'GRAYCALM', 'AU GRES', 'RUBICON', 'CROSWELL', 'KINROSS', 'Iosco', 'DEFORD', 'MORGANLAKE', 'EMMET', 'ALLENDALE', 'LEAFRIVER', 'BATTLEFIELD', 'WHEATLEY', 'EAST LAKE', 'BLUE LAKE', 'HALFADAY',  'ALGONQUIN', 'SPRINGPORT', 'WAKELEY','CHINWHISKER','MOSSBACK','LEELANAU', 'ISLANDLAKE','SOUTHWELLS','SPRINGLAKE','ORTISCO','MENOMINEE','ONAWAY','CHESTONIA','CHARLEVOIX','ENSLEY')

SelectSeries <- c('CATHRO', 'TAWAS', 'WAKELEY', 'ALGONQUIN', 'CHESTONIA', 'SPRINGLAKE', 'LUPTON', 'ANGELICA', 'IOSCO', 'ENSLEY', 'EAST LAKE', 'CHARLEVOIX', 'AU GRES', 'RUBICON', 'KAWKAWLIN', 'SOUTHWELLS', 'OSSINEKE', 'BLUE LAKE', 'ISLANDLAKE', 'KALKASKA', 'KALKASKA, BURNED', 'EMMET', 'ONAWAY', 'MOSSBACK', 'MANCELONA', 'DEFORD', 'DEER PARK', 'CROSWELL', 'HALFADAY', 'KINROSS', 'OTISCO', 'CHINWHISKER', 'DAWSON', 'MENOMINEE', 'MORGANLAKE', 'LEELANAU', 'LEAFRIVER', 'ALLENDALE', 'BATTLEFIELD')


SelectFIPS <- c(26009)
SelectSpecies <- c('ABBA','PIMA', 'PIBA2','FAGR','ACSA3','QUAL','POGR4', 'PIST', 'PIRE', 'QURU', 'POTR5', 'ACRU','FRNI','THOC2','LALA')

SelectSpecies <- c('POGR4')
SelectSeries <- c('RUBICON', 'OTISCO', 'CHINWHISKER', 'LEELANAU', 'CROSWELL', 'HALFADAY', 'KINROSS', 'EMMET', 'MENOMINEE', 'MORGANLAKE', 'KALKASKA')
ggplot(aes(x=paste(SYMBOL,series), y=exp(LNSI)/0.3048), data = si[si$SYMBOL %in% SelectSpecies & si$series %in% SelectSeries,])+
  geom_boxplot()+
  scale_y_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,200))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


inputSpecies <- species[species$SYMBOL %in% SelectSpecies,]
inputFIPS <- FIPS[FIPS$fips %in% SelectFIPS,]
inputSeries <- Series[Series$series %in% SelectSeries,]

inputsall <- merge(inputFIPS, inputSeries)
inputsall <- merge(inputSpecies, inputsall)

inputsall$predSI <- predict(stepbackward,inputsall)
inputsall$SI_m <- exp(inputsall$predSI)
inputsall$SI_ft <- inputsall$SI_m/0.3048
inputsall$ppt <- exp(inputsall$lnppt)
summarysi <- aggregate(inputsall[,c('Warm','Cold','ppt', 'SI_m','SI_ft')], 
                                 by= list(inputsall$fips, inputsall$SYMBOL,inputsall$series), FUN='mean')
colnames(summarysi) <- c('fips', 'SYMBOL', 'series', 'Warm','Cold','ppt', 'SI_m','SI_ft')
showtable <- (summarysi[order(summarysi$SYMBOL, summarysi$SI_m),c('fips', 'Warm','Cold','ppt', 'series', 'SYMBOL','SI_m','SI_ft')])


################# weighted model
formodel <- si
formodel$wts <- 1
formodel[formodel$SYMBOL %in% SelectSpecies,]$wts <- 100
formodel[formodel$series %in% SelectSeries,]$wts <- 1000
simodelmax <- lm(LNSI ~  I(1/(WaterTable+25)) + lnppt + Cold + Warm + I(Warm^0.5) + AWC150 + poly(WaterTable,2) + om150 + Wet + poly(pH50,2) + clay150 + sand150 + carbdepth + Wet*WaterTable + Carya + Fagus + Juglans + Liriodendron + Populus + Prunus + Tilia + Ulmus + ABBA + ACRU + ACSA3 + FRAM2 + FRPE + JUVI + PIST + PIRE + PIBA2 + QUAL + QURU + QUVE + PODE + THOC2 + Pinus*WaterTable + Pinus*AWC150 + Pinus*pH50 + Quercus*WaterTable + Quercus*AWC150 + Quercus*pH50 + WaterTable*AWC150 + WaterTable*lnppt + lnppt*Warm + Bhs + Bhs*AWC150 + Bhs*WaterTable + WOOD_SPGR_GREENVOL_DRYWT + Shade + Ultisols + Alfisols + Inceptisols + Mollisols + Vertisols + Andisols + Pinus + Quercus + Acer + Fraxinus + Picea + Liquidambar + Abies + Betula + Juniperus + Larix + Celtis + Thuja + Platanus + Robinia + Nyssa + Pseudotsuga + Tsuga + Sequoia + frost50, data = formodel, weights = formodel$wts)

summary(simodelmax)
simodelmin <- lm(LNSI ~  1, data = formodel, weights = formodel$wts)
step(simodelmin, scope=list(upper = simodelmax,
                            lower= simodelmin), direction="forward", trace=0)
wtmodel <-
  lm(formula = LNSI ~ om150 + frost50 + THOC2 + Populus + Picea + 
       Pinus + lnppt + AWC150 + ABBA + poly(pH50, 2) + Larix + Cold + 
       Quercus + poly(WaterTable, 2) + Warm + I(Warm^0.5) + ACSA3 + 
       Wet + Juniperus + Shade + carbdepth + QUAL + PIRE + Fagus + 
       Inceptisols + Tilia + I(1/(WaterTable + 25)) + QURU + clay150 + 
       Liriodendron + sand150 + Bhs + Ultisols + Alfisols + PIBA2 + 
       PIST + WOOD_SPGR_GREENVOL_DRYWT + Sequoia + PODE + Pseudotsuga + 
       Tsuga + FRAM2 + Prunus + Fraxinus + Carya + JUVI + Andisols + 
       Liquidambar + Platanus + Vertisols + Mollisols + Pinus:AWC150 + 
       lnppt:Warm + AWC150:Quercus + AWC150:Bhs, data = formodel, 
     weights = formodel$wts)
summary(wtmodel)

inputSpecies <- species[species$SYMBOL %in% SelectSpecies,]
inputFIPS <- FIPS[FIPS$fips %in% SelectFIPS,]
inputSeries <- Series[Series$series %in% SelectSeries,]

inputsall <- merge(inputFIPS, inputSeries)
inputsall <- merge(inputSpecies, inputsall)

inputsall$predSI <- predict(wtmodel,inputsall)
inputsall$SI_m <- exp(inputsall$predSI)
inputsall$SI_ft <- inputsall$SI_m/0.3048
inputsall$ppt <- exp(inputsall$lnppt)
summarysi <- aggregate(inputsall[,c('Warm','Cold','ppt', 'SI_m','SI_ft')], 
                       by= list(inputsall$fips, inputsall$SYMBOL,inputsall$series), FUN='mean')
colnames(summarysi) <- c('fips', 'SYMBOL', 'series', 'Warm','Cold','ppt', 'SI_m','SI_ft')
showtable <- (summarysi[order(summarysi$SYMBOL, summarysi$SI_m),c('fips', 'Warm','Cold','ppt', 'series', 'SYMBOL','SI_m','SI_ft')])

write.csv(showtable, 'showtable.csv')
#########weight standard formula with species of interest

SelectSpecies <- c('ABBA', 'ACRU', 'ACSA3', 'BEPA', 'FRNI', 'LALA', 'PIBA2', 'PIGL', 'PIMA', 'PIRE', 'PIST', 'POGR4', 'POTR5', 'QURU', 'THOC2', 'TIAM')
SelectSpecies <- c('POGR4')
SelectSeries <- c('KALKASKA')
formodel <- si
formodel$wts <- 1
formodel[formodel$SYMBOL %in% SelectSpecies,]$wts <- 1000
formodel[formodel$series %in% SelectSeries,]$wts <- 1000

stepforward <-   
  lm(formula = LNSI ~ lnppt + frost50 + Pseudotsuga + Liriodendron + 
       AWC150 + Juniperus + Populus + om150 + Warm + I(Warm^0.5) + 
       Wet + Cold + poly(pH50, 2) + QUAL + Pinus + THOC2 + clay150 + 
       FRAM2 + Shade + Picea + ABBA + WOOD_SPGR_GREENVOL_DRYWT + 
       Prunus + QUVE + QURU + sand150 + carbdepth + Acer + PIRE + 
       PIBA2 + Ultisols + poly(WaterTable, 2) + PIST + Betula + 
       Quercus + Robinia + Liquidambar + Juglans + JUVI + Tilia + 
       Abies + Larix + Platanus + Nyssa + Carya + Vertisols + Mollisols + 
       Alfisols + Fagus + Fraxinus + PODE + Inceptisols + ACRU + 
       ACSA3 + AWC150:Pinus + lnppt:Warm,
     data = formodel, weights = formodel$wts)
summary(stepforward)


inputSpecies <- species[species$SYMBOL %in% SelectSpecies,]
inputFIPS <- FIPS[FIPS$fips %in% SelectFIPS,]
inputSeries <- Series[Series$series %in% SelectSeries,]

inputsall <- merge(inputFIPS, inputSeries)
inputsall <- merge(inputSpecies, inputsall)

inputsall$predSI <- predict(stepforward,inputsall)
inputsall$SI_m <- exp(inputsall$predSI)
inputsall$SI_ft <- inputsall$SI_m/0.3048
inputsall$ppt <- exp(inputsall$lnppt)
summarysi <- aggregate(inputsall[,c('Warm','Cold','ppt', 'SI_m','SI_ft')], 
                       by= list(inputsall$fips, inputsall$SYMBOL,inputsall$series), FUN='mean')
colnames(summarysi) <- c('fips', 'SYMBOL', 'series', 'Warm','Cold','ppt', 'SI_m','SI_ft')
showtable <- (summarysi[order(summarysi$SYMBOL, summarysi$SI_m),c('fips', 'Warm','Cold','ppt', 'series', 'SYMBOL','SI_m','SI_ft')])

########### find residuals by series
stepbackward <-    
  lm(formula = lm(formula = LNSI ~ lnppt + Cold + Warm + I(Warm^0.5) + AWC150 + 
                    poly(WaterTable, 2) + om150 + Wet + poly(pH50, 2) + clay150 + 
                    sand150 + carbdepth + Carya + Fagus + Juglans + 
                    Liriodendron + Populus + Prunus + Tilia + ABBA + ACRU + ACSA3 + 
                    FRAM2 + FRPE + JUVI + PIRE + PIBA2 + QUAL + QURU + QUVE + 
                    PODE + THOC2 + Pinus + Quercus + Bhs + WOOD_SPGR_GREENVOL_DRYWT + 
                    Shade + Ultisols + Mollisols + Vertisols + Andisols + Acer + 
                    Picea + Liquidambar + Abies + Betula + Juniperus + Larix + 
                    Platanus + Robinia + Nyssa + Pseudotsuga + frost50 + Wet:WaterTable + 
                    WaterTable:Pinus + AWC150:Pinus + Pinus:pH50 + WaterTable:Quercus + 
                    AWC150:Quercus + pH50:Quercus + AWC150:WaterTable + lnppt:WaterTable + 
                    lnppt:Warm + AWC150:Bhs + WaterTable:Bhs, data = formodel), data = formodel)

summary(stepbackward)

stepforward <-   
  lm(formula = LNSI ~ lnppt + frost50 + Pseudotsuga + Liriodendron + 
       AWC150 + Juniperus + Populus + om150 + Warm + I(Warm^0.5) + 
       Wet + Cold + poly(pH50, 2) + QUAL + Pinus + THOC2 + clay150 + 
       FRAM2 + Shade + Picea + ABBA + WOOD_SPGR_GREENVOL_DRYWT + 
       Prunus + QUVE + QURU + sand150 + carbdepth + Acer + PIRE + 
       PIBA2 + Ultisols + poly(WaterTable, 2) + PIST + Betula + 
       Quercus + Robinia + Liquidambar + Juglans + JUVI + Tilia + 
       Abies + Larix + Platanus + Nyssa + Carya + Vertisols + Mollisols + 
       Alfisols + Fagus + Fraxinus + PODE + Inceptisols + ACRU + 
       ACSA3 + AWC150:Pinus + lnppt:Warm,
     data = formodel, weights = formodel$wts)
summary(stepforward)
inputsall <- si

inputsall$predSI <- predict(stepforward,inputsall)
inputsall$resids <- inputsall$LNSI - inputsall$predSI

SelectSeries <- c('MYAKKA', 'ELVE', 'SNOWDON', 'GILPIN', 'HOUGHTON','DAWSON','TAWAS', 'MORROCO','BREMS', 'RUBICON', 'GRAYLING','GRAYCALM','KALEVA','GRATTAN', 'PLAINFIELD', 'NESTER', 'PERRINTON', 'CROSWELL', 'PIPESTONE', 'KALKASKA')
ggplot(aes(x=series, y=resids), data = inputsall[si$series %in% SelectSeries,])+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
