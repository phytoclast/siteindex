library(ggplot2)
bt_frost <- read.delim('bt_frost2.txt')
bt_frost$Warm <- bt_frost$bt0510
bt_frost$Cold <- bt_frost$tmin
bt_frost$frostbin <- log(1/(365.25/bt_frost$Frost50 - 1))
frostmodelmax <- lm(formula = Frost50 ~ Warm + I(Warm^2) + poly(Cold, 2) + Warm:Cold, 
                    data = bt_frost)
frostmodelmin <- lm(Frost50 ~  1, data = bt_frost)

stepmodelfrost <- lm(formula = Frost50 ~ Warm + poly(Cold, 2) + Warm:Cold, data = bt_frost)
summary(stepmodelfrost)


si <- read.delim('SimpleExport.txt')
si <- subset(si, !is.na(Cold) & !is.na(pH50) & AWC150 >0 & WOOD_SPGR_GREENVOL_DRYWT >0 & elev < 3600)
si$frost50 <- predict(stepmodelfrost, si)
formodel <- si[,6:ncol(si)]

# create test set

SelectSeries <- c('CATHRO', 'TAWAS', 'WAKELEY', 'ALGONQUIN', 'CHESTONIA', 'SPRINGLAKE', 'LUPTON', 'ANGELICA', 'IOSCO', 'ENSLEY', 'EAST LAKE', 'CHARLEVOIX', 'AU GRES', 'RUBICON', 'KAWKAWLIN', 'SOUTHWELLS', 'OSSINEKE', 'BLUE LAKE', 'ISLANDLAKE', 'KALKASKA', 'KALKASKA, BURNED', 'EMMET', 'ONAWAY', 'MOSSBACK', 'MANCELONA', 'DEFORD', 'DEER PARK', 'CROSWELL', 'HALFADAY', 'KINROSS', 'OTISCO', 'CHINWHISKER', 'DAWSON', 'MENOMINEE', 'MORGANLAKE', 'LEELANAU', 'LEAFRIVER', 'ALLENDALE', 'BATTLEFIELD')


SelectFIPS <- c(26009)
SelectSpecies <- c('BEPA', 'ABBA','PIMA', 'PIBA2','FAGR','ACSA3','QUAL','POGR4', 'PIST', 'PIRE', 'QURU', 'POTR5', 'ACRU','FRNI','THOC2','LALA','TIAM')



inputSpecies <- species[species$SYMBOL %in% SelectSpecies,]
inputFIPS <- FIPS[FIPS$fips %in% SelectFIPS,]
inputSeries <- Series[Series$series %in% SelectSeries,]

inputsall <- merge(inputFIPS, inputSeries)
inputsall <- merge(inputSpecies, inputsall)



########### find residuals by series
################# weighted model
formodel <- si
formodel$wts <- 1
formodel[formodel$SYMBOL %in% SelectSpecies,]$wts <- 100
formodel[formodel$series %in% SelectSeries,]$wts <- 1000
wtmodel <-    
  lm(formula = LNSI ~ om150 + frost50 + Thuja + Populus + QUAL + 
       Picea + Pinus + lnppt + AWC150 + ABBA + Larix + Cold + Warm + 
       I(Warm^0.5) + poly(pH50, 2) + poly(WaterTable, 2) + Quercus + 
       carbdepth + Juniperus + ACSA3 + PIBA2 + Inceptisols + Wet + 
       Fagus + Betula + QURU + clay150 + sand150 + WOOD_SPGR_GREENVOL_DRYWT + 
       Liriodendron + Shade + PIRE + PIST + Tilia + I(1/(WaterTable + 
                                                           25)) + Bhs + Alfisols + Ultisols + FRAM2 + Pseudotsuga + 
       Sequoia + Tsuga + PODE + Carya + JUVI + ACRU + Acer + Liquidambar + 
       Prunus + Platanus + Vertisols + Fraxinus + Andisols + Mollisols + 
       lnppt:Warm + Pinus:AWC150 + AWC150:Quercus + AWC150:Bhs, 
     data = formodel, weights = formodel$wts)
summary(wtmodel)

inputself <- si
#determine residual by series
inputself$predSI <- predict(wtmodel,inputself)
inputself$resids <- inputself$LNSI - inputself$predSI

resids <- aggregate(inputself[,c('resids')], by= list(inputself$series), FUN='mean')
colnames(resids) <- c('series','resids')

inputsall <- merge(inputsall, resids, by='series')
inputsall$predSI <- predict(stepbackward,inputsall)+inputsall$resids
inputsall$SI_m <- exp(inputsall$predSI)
inputsall$SI_ft <- inputsall$SI_m/0.3048
inputsall$ppt <- exp(inputsall$lnppt)
summarysi <- aggregate(inputsall[,c('Warm','Cold','ppt', 'SI_m','SI_ft')], 
                       by= list(inputsall$fips, inputsall$SYMBOL,inputsall$series), FUN='mean')
colnames(summarysi) <- c('fips', 'SYMBOL', 'series', 'Warm','Cold','ppt', 'SI_m','SI_ft')
showtable <- (summarysi[order(summarysi$SYMBOL, summarysi$SI_m),c('fips', 'Warm','Cold','ppt', 'series', 'SYMBOL','SI_m','SI_ft')])


SelectSeries2 <- c('MYAKKA', 'ELVE', 'SNOWDON', 'GILPIN', 'HOUGHTON','DAWSON','TAWAS', 'MORROCO','BREMS', 'RUBICON', 'GRAYLING','GRAYCALM','KALEVA','GRATTAN', 'PLAINFIELD', 'NESTER', 'PERRINTON', 'CROSWELL', 'PIPESTONE', 'KALKASKA')
ggplot(aes(x=series, y=SI_ft), data = inputsall[si$series %in% SelectSeries2,])+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

write.csv(showtable, 'showtable.csv')