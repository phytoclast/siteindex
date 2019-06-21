library(MASS)
bt_frost <- read.delim('bt_frost2.txt')
bt_frost$Warm <- bt_frost$bt0510
bt_frost$Cold <- bt_frost$tmin
frostmodelmax <- lm(Frost50 ~  poly(Warm,2) + poly(Cold,2) + Warm:Cold, data = bt_frost)
frostmodelmin <- lm(Frost50 ~  1, data = bt_frost)

step(frostmodelmax, scope=list(upper = frostmodelmax,
                               lower= frostmodelmin), direction="both", trace=0)


step(frostmodelmin, scope=list(upper = frostmodelmax,
                               lower= frostmodelmin), direction="forward")

stepmodelfrost <- lm(formula = Frost50 ~ poly(Warm, 2) + poly(Cold, 2) + Warm:Cold, 
   data = bt_frost)
summary(stepmodelfrost)


si <- read.delim('SimpleExport.txt')
si <- subset(si, !is.na(Cold) & !is.na(pH50) & AWC150 >0 & WOOD_SPGR_GREENVOL_DRYWT >0)
si$frost50 <- predict(stepmodelfrost, si)
formodel <- si[,5:ncol(si)]
######################################

newmodel<- stepAIC(simodel, scope=list(upper= ~ lnppt+Cold+Warm+AWC150+WaterTable+om150+Wet+pH50+Carya+Fagus+Juglans+Liriodendron+Populus+Prunus+Tilia+Ulmus+ABBA+ACRU+ACSA3+FRAM2+FRPE+JUVI+PIST+PIRE+PIBA2+QUAL+QURU+QUVE+PODE+THOC2, lower= ~1))
################################################

simodel <- lm(LNSI ~  ., data = formodel)

summary(simodel)
vartab <- cor(subset(formodel, !is.na(Cold)&!is.na(pH50)))
#newmodel<- stepAIC(simodel, scope=list(upper= ~ ., lower= ~1))

formodel$rand <- sample(1:100,1)


formodel <- si[,5:ncol(si)]
formodel <- subset(formodel, !is.na(Cold)&!is.na(pH50))
  randomnumbers <- rnorm(nrow(formodel), 1, 1)
  #test the affects of a truly random variable
  formodel2 <- cbind(formodel[,c("LNSI","lnppt","Cold","Warm","AWC150",'Bhs',"WaterTable",'Inceptisols','Liriodendron')], randomnumbers)
 
  
    simodelmax <- lm(LNSI ~  lnppt + Cold + Warm + AWC150 + poly(WaterTable,2) + om150 + Wet + poly(pH50,2) + clay150 + sand150 + carbdepth + Wet*WaterTable + Carya + Fagus + Juglans + Liriodendron + Populus + Prunus + Tilia + Ulmus + ABBA + ACRU + ACSA3 + FRAM2 + FRPE + JUVI + PIST + PIRE + PIBA2 + QUAL + QURU + QUVE + PODE + THOC2 + Pinus*WaterTable + Pinus*AWC150 + Pinus*pH50 + Quercus*WaterTable + Quercus*AWC150 + Quercus*pH50 + WaterTable*AWC150 + WaterTable*lnppt + lnppt*Warm + Bhs + Bhs*AWC150 + Bhs*WaterTable + WOOD_SPGR_GREENVOL_DRYWT + Shade + Ultisols + Alfisols + Inceptisols + Mollisols + Vertisols + Andisols + Pinus + Quercus + Acer + Fraxinus + Picea + Liquidambar + Abies + Betula + Juniperus + Larix + Celtis + Thuja + Platanus + Robinia + Nyssa + Pseudotsuga + Tsuga + Sequoia + frost50, data = formodel)
  
    summary(simodelmax)
  simodelmin <- lm(LNSI ~  1, data = formodel)
  step(simodelmin, scope=list(upper = simodelmax,
                              lower= simodelmin), direction="forward")
  step(simodelmax, scope=list(upper = simodelmax,
                              lower= simodelmin), direction="backward")
  step(simodelmax, scope=list(upper = simodelmax,
                              lower= simodelmin), direction="both", trace=0)


  
  
  
  ##### trained the model with a frostfree day model
  
forwardfrost <-   
  lm(formula = LNSI ~ lnppt + frost50 + Pseudotsuga + Liriodendron + 
       AWC150 + Juniperus + Populus + om150 + Warm + Pinus + poly(pH50, 
                                                                  2) + Picea + Wet + Cold + QUAL + PIST + THOC2 + FRAM2 + Prunus + 
       Acer + QURU + PIRE + Ultisols + clay150 + sand150 + PIBA2 + 
       Tilia + QUVE + Betula + WOOD_SPGR_GREENVOL_DRYWT + carbdepth + 
       Quercus + Liquidambar + Juglans + Shade + Robinia + JUVI + 
       ABBA + Platanus + Tsuga + Larix + poly(WaterTable, 2) + Nyssa + 
       Carya + Sequoia + Mollisols + Inceptisols + ACSA3 + ACRU + 
       PODE + Fagus + Fraxinus + Abies + lnppt:Warm + AWC150:Pinus + 
       AWC150:Quercus, data = formodel)
  
  summary(forwardfrost)
  ##### trained the model with a frostfree day model
  backwardfrost <-    
  lm(formula = LNSI ~ lnppt + Cold + Warm + AWC150 + poly(WaterTable, 
                                                          2) + om150 + Wet + poly(pH50, 2) +clay150 + sand150 + carbdepth + WaterTable + Carya + Fagus + Juglans + Liriodendron + Populus + Prunus + Tilia + ABBA + ACRU + ACSA3 + FRAM2 + JUVI + PIST + PIRE + PIBA2 + QUAL + QURU + QUVE + PODE + THOC2 + Pinus + pH50 + Quercus + Bhs + WOOD_SPGR_GREENVOL_DRYWT + Shade + Ultisols + Inceptisols + Mollisols + Acer + Fraxinus + Picea + Liquidambar + Betula + Juniperus + Larix + Platanus + Robinia + Nyssa + Pseudotsuga + Tsuga + Sequoia + frost50 + Wet:WaterTable + AWC150:Pinus + Pinus:pH50 + WaterTable:Quercus + AWC150:Quercus + AWC150:WaterTable + lnppt:WaterTable + lnppt:Warm + AWC150:Bhs + WaterTable:Bhs, data = formodel)
  
  summary(backwardfrost)
 
  
allsi <- si[,c("lnppt", "Cold", "Warm", "AWC150", "poly(WaterTable,2)", "om150", "Wet", "poly(pH50,2)", "clay150", "sand150", "carbdepth", "Wet*WaterTable", "Carya", "Fagus", "Juglans", "Liriodendron", "Populus", "Prunus", "Tilia", "Ulmus", "ABBA", "ACRU", "ACSA3", "FRAM2", "FRPE", "JUVI", "PIST", "PIRE", "PIBA2", "QUAL", "QURU", "QUVE", "PODE", "THOC2", "Pinus*WaterTable", "Pinus*AWC150", "Pinus*pH50", "Quercus*WaterTable", "Quercus*AWC150", "Quercus*pH50", "WaterTable*AWC150", "WaterTable*lnppt", "lnppt*Warm", "Bhs", "Bhs*AWC150", "Bhs*WaterTable", "WOOD_SPGR_GREENVOL_DRYWT", "Shade", "Ultisols", "Alfisols", "Inceptisols", "Mollisols", "Vertisols", "Andisols", "Pinus", "Quercus", "Acer", "Fraxinus", "Picea", "Liquidambar", "Abies", "Betula", "Juniperus", "Larix", "Celtis", "Thuja", "Platanus", "Robinia", "Nyssa", "Pseudotsuga", "Tsuga", "Sequoia", "frost50")]

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


   ###################
  library(ggplot2)
  library(RColorBrewer)
  forplot <- aggregate(si[,c("frost50")], by=list(si$Cold, si$Warm), FUN='mean')
  colnames(forplot) <- c('Cold','Warm','frost50')
  forplot <- unique(forplot)
  
  ggplot(aes(x=Cold, y=Warm), data = forplot)+
    geom_point(aes(color=frost50))+
    scale_colour_gradient2(low='green', high='red', mid = 'white', midpoint = 150)  
  
  
  