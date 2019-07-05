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



#---- Simple filter
library(plyr)
SelectSeries2 <- c('MYAKKA', 'ELVE', 'SNOWDON', 'GILPIN', 'HOUGHTON','DAWSON','TAWAS', 'MORROCO','BREMS', 'RUBICON', 'GRAYLING','GRAYCALM','KALEVA','GRATTAN', 'PLAINFIELD', 'NESTER', 'PERRINTON', 'CROSWELL', 'PIPESTONE', 'KALKASKA')

Selected_si <- subset(si, Warm <=16 & Warm >= 13 & Cold <= -9 & Cold >= -15 &  exp(lnppt) >= 700 & exp(lnppt) <= 900 )

Selected_si_agg <-    ddply(Selected_si, c("series", 'SYMBOL'), summarise, Nrows = length(series),
                                          siMin = exp(quantile(LNSI, 0))/0.3048, siMax = exp(quantile(LNSI, 1))/0.3048, siMean = exp(mean(LNSI))/0.3048)
write.csv(Selected_si_agg, 'Selected_si_agg.csv')



ggplot(aes(x=series, y=exp(LNSI)), data = Selected_si[Selected_si$series %in% SelectSeries2,])+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
