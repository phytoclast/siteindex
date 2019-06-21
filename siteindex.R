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
si <- subset(si, !is.na(Cold) & !is.na(pH50) & AWC150 >0 & WOOD_SPGR_GREENVOL_DRYWT >0)
si$frost50 <- predict(stepmodelfrost, si)
formodel <- si[,5:ncol(si)]
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
  simodelmax <- lm(LNSI ~  lnppt + Cold + Warm + I(Warm^0.5) + AWC150 + poly(WaterTable,2) + om150 + Wet + poly(pH50,2) + clay150 + sand150 + carbdepth + Wet*WaterTable + Carya + Fagus + Juglans + Liriodendron + Populus + Prunus + Tilia + Ulmus + ABBA + ACRU + ACSA3 + FRAM2 + FRPE + JUVI + PIST + PIRE + PIBA2 + QUAL + QURU + QUVE + PODE + THOC2 + Pinus*WaterTable + Pinus*AWC150 + Pinus*pH50 + Quercus*WaterTable + Quercus*AWC150 + Quercus*pH50 + WaterTable*AWC150 + WaterTable*lnppt + lnppt*Warm + Bhs + Bhs*AWC150 + Bhs*WaterTable + WOOD_SPGR_GREENVOL_DRYWT + Shade + Ultisols + Alfisols + Inceptisols + Mollisols + Vertisols + Andisols + Pinus + Quercus + Acer + Fraxinus + Picea + Liquidambar + Abies + Betula + Juniperus + Larix + Celtis + Thuja + Platanus + Robinia + Nyssa + Pseudotsuga + Tsuga + Sequoia + frost50, data = formodel)
  
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
  lm(formula = LNSI ~ lnppt + frost50 + Pseudotsuga + Liriodendron + 
       AWC150 + Juniperus + Populus + om150 + Warm + I(Warm^0.5) + 
       Wet + Cold + poly(pH50, 2) + QUAL + Pinus + THOC2 + clay150 + 
       FRAM2 + Shade + Picea + ABBA + WOOD_SPGR_GREENVOL_DRYWT + 
       Prunus + QUVE + QURU + sand150 + carbdepth + Acer + PIRE + 
       PIBA2 + Ultisols + poly(WaterTable, 2) + PIST + Betula + 
       Quercus + Robinia + Liquidambar + Juglans + JUVI + Tilia + 
       Abies + Larix + Platanus + Nyssa + Carya + Vertisols + Mollisols + 
       Alfisols + Fagus + Fraxinus + PODE + Inceptisols + ACRU + 
       ACSA3 + AWC150:Pinus + lnppt:Warm, data = formodel)
summary(stepforward)
##### trained the model with a frostfree day model
stepbackward <-    
  lm(formula = LNSI ~ lnppt + Cold + Warm + AWC150 + poly(WaterTable, 
                                                          2) + om150 + Wet + poly(pH50, 2) +clay150 + sand150 + carbdepth + WaterTable + Carya + Fagus + Juglans + Liriodendron + Populus + Prunus + Tilia + ABBA + ACRU + ACSA3 + FRAM2 + JUVI + PIST + PIRE + PIBA2 + QUAL + QURU + QUVE + PODE + THOC2 + Pinus + pH50 + Quercus + Bhs + WOOD_SPGR_GREENVOL_DRYWT + Shade + Ultisols + Inceptisols + Mollisols + Acer + Fraxinus + Picea + Liquidambar + Betula + Juniperus + Larix + Platanus + Robinia + Nyssa + Pseudotsuga + Tsuga + Sequoia + frost50 + Wet:WaterTable + AWC150:Pinus + Pinus:pH50 + WaterTable:Quercus + AWC150:Quercus + AWC150:WaterTable + lnppt:WaterTable + lnppt:Warm + AWC150:Bhs + WaterTable:Bhs, data = formodel)

summary(stepbackward)

species <- unique(si[,c("SYMBOL","Wet","WOOD_SPGR_GREENVOL_DRYWT", "Shade",
                        "Carya", "Fagus", "Juglans", "Liriodendron", "Populus", "Prunus", "Tilia", "Ulmus", 
                        "ABBA", "ACRU", "ACSA3", "FRAM2", "FRPE", "JUVI", "PIST", "PIRE", "PIBA2", "QUAL", 
                        "QURU", "QUVE", "PODE", "THOC2",  "Pinus", "Quercus", "Acer", "Fraxinus", "Picea", 
                        "Liquidambar", "Abies", "Betula", "Juniperus", "Larix", "Celtis", "Thuja", "Platanus", 
                        "Robinia", "Nyssa", "Pseudotsuga", "Tsuga", "Sequoia")]
)

FIPS <- unique(si[,c("fips","lnppt", "Cold", "Warm","frost50")])

Series <- unique(si[,c("series", "om150", "pH50", "clay150", "sand150","WaterTable",
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


