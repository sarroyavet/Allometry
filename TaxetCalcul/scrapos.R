#### Kilbourne et al 2013 for the mass in gr #### 
#### ...moment of inertia of the forelimb #### 
Data_Animal_de_pd['MoI'] <- 
  10*(10^(-2.02)*(1/1000^1.78)*Data_Animal_de_pd['massAvg']^(1.78))  # g/cm2 *10 to kg/m2
#### ...mass of the forelimb  #### 
Data_Animal_de_pd['Limb_mass'] <- 
  (10^(-1.15)*(1/1000^1.01)*Data_Animal_de_pd['massAvg']^(1.01))/1000 # g /1000 to kg
#### ...forelimb center of mass #### 
Data_Animal_de_pd['CoM'] <- 
  (10^(-0.40)*(1/1000^0.37)*Data_Animal_de_pd['massAvg']^(0.37))/100 # cm /100 to m
#### Deduced equation for MoI at the center of gravity (kilbourne data) ####
Data_Animal_de_pd['MoI_G'] <- 
  Data_Animal_de_pd['MoI']- Data_Animal_de_pd['Limb_mass']*
  Data_Animal_de_pd['CoM']^2

#### Moment of inertia in the at the almost elbow ####
Data_Animal_de_pd['MoI_H'] <- 
  Data_Animal_de_pd['MoI_G']+ 
  Data_Animal_de_pd['Limb_mass']*
  (Data_Animal_de_pd['Long_m']-Data_Animal_de_pd['CoM'])^2

#### Angular acceleration of the elbow ####
Data_Animal_de_pd['Alpha'] <- 
  Data_Animal_de_pd['Mscl_Force']*Data_Animal_de_pd['MoA_m']/
  (Data_Animal_de_pd['MoI_H']*0.6646727) 

#### Angular velocity of the elbow ####
Data_Animal_de_pd['Omg'] <- Data_Animal_de_pd['Alpha']*0.254/
  (Data_Animal_de_pd['Sf']) # Angular velocity of the joint for only 0.254 of 
# the period (assuming it is when we have the 
# max speed)

#### Sliding velocity  ####
Data_Animal_de_pd['Vsup'] <- Data_Animal_de_pd['Omg']*
  Data_Animal_de_pd['Dmed_m']/2


# ALL and MASS ####
Long.SMA <- getsma(Long, massAvg, Data_Animal_de_pd, 1/3, TRUE)
summary(Long.SMA)
plot(Long.SMA)

Long.SMA.stance <- getsmaGroup(Long, massAvg, Stance, Data_Animal_de_pd, FALSE)
summary(Long.SMA.stance)
plotgrps(Long.SMA.stance)

# 
#### Long and MASS ####
Long.SMA <- sma(Long~massAvg, 
                log = "xy",
                data = Data_Animal_de_pd, 
                slope.test = 1/3, 
                robust=T
)
print(Long.SMA)
summary(Long.SMA)
plot(Long.SMA)

# Stance
Long.SMA.stance <- sma(Long~massAvg*Stance, 
                       log = "xy",
                       data = Data_Animal_de_pd, 
                       slope.test = 1/3, 
                       robust = T
)
summary(Long.SMA.stance)
plot(Long.SMA.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Long.SMA.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#Morphotype
Long.SMA.Morphotype <- sma(Long~massAvg*Morphotype, 
                           log = "xy",
                           data = Data_Animal_de_pd, 
                           slope.test = 1/3, 
                           robust = T
)
summary(Long.SMA.Morphotype)
plot(Long.SMA.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Long.SMA.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Long.SMA.order <- sma(Long~massAvg*order, 
                      log = "xy",
                      data = Data_Animal_de_pd, 
                      slope.test = 1/3,
                      robust = T
)
summary(Long.SMA.order)
plot(Long.SMA.order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Long.SMA.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Family
Long.SMA.family <- sma(Long~massAvg*family, 
                       log = "xy",
                       data = Data_Animal_de_pd, 
                       slope.test = 1/3,
                       robust = T
)
summary(Long.SMA.family)
plot(Long.SMA.family, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Long.SMA.family[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
#Diet
Long.SMA.diet <- sma(Long~massAvg*Diet, 
                     log = "xy", 
                     data = Data_Animal_de_pd,
                     slope.test = 1/3,
                     robust = T
)
summary(Long.SMA.diet)
plot(Long.SMA.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Long.SMA.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#

#### La and MASS ####
La.SMA <- sma(La~massAvg, 
              log = "xy",
              data = Data_Animal_de_pd,
              slope.test = 1/3,
              #robust = T
)
summary(La.SMA)
plot(La.SMA)

# Stance
La.SMA.stance <- sma(La~massAvg*Stance,
                     log = "xy",
                     data = Data_Animal_de_pd, 
                     slope.test = 1/3,
                     robust = T
)
summary(La.SMA.stance)
plot(La.SMA.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = La.SMA.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype 
La.SMA.Morphotype <- sma(La~massAvg*Morphotype, 
                         log = "xy",
                         data = Data_Animal_de_pd, 
                         slope.test = 1/3, 
                         robust = T
)
summary(La.SMA.Morphotype)
plot(La.SMA.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = La.SMA.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
La.SMA.order <- sma(La~massAvg*order,
                    log = "xy",
                    data = Data_Animal_de_pd,
                    slope.test = 1/3,
                    robust = T
)
summary(La.SMA.order)
plot(La.SMA.order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = La.SMA.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Family
La.SMA.family <- sma(La~massAvg*family, 
                     log = "xy",
                     data = Data_Animal_de_pd, 
                     slope.test = 1/3,
                     robust = T
)
summary(La.SMA.family)
plot(La.SMA.family, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = La.SMA.family[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
La_SMA_diet <- sma(La~massAvg*Diet, 
                   log = "xy",
                   data = Data_Animal_de_pd,
                   slope.test = 1/3, 
                   robust = T
)
summary(La_SMA_diet)
plot(La_SMA_diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = La_SMA_diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#
#### Rmax and MASS ####
Rmax.SMA <- sma(Rmax~massAvg, 
                log = "xy",
                data = Data_Animal_de_pd, 
                slope.test = 1/3, 
                robust = T
)
summary(Rmax.SMA)
plot(Rmax.SMA)
# Stance
Rmax.SMA.stance <- sma(Rmax~massAvg*Stance, 
                       log = "xy",
                       data = Data_Animal_de_pd,
                       slope.test = 1/3, 
                       robust = T
)
summary(Rmax.SMA.stance)
plot(Rmax.SMA.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Rmax.SMA.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Rmax.SMA.Morphotype <- sma(Rmax~massAvg*Morphotype, 
                           log = "xy",
                           data = Data_Animal_de_pd, 
                           slope.test = 1/3, 
                           #robust = T
)
summary(Rmax.SMA.Morphotype)
plot(Rmax.SMA.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Rmax.SMA.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Rmax.SMA.order <- sma(Rmax~massAvg*order,
                      log = "xy",
                      data = Data_Animal_de_pd,
                      slope.test = 1/3,
                      robust = T
)
summary(Rmax_SMA_order)
plot(Rmax_SMA_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Rmax_SMA_order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Family
Rmax.SMA.family <- sma(Rmax~massAvg*family, 
                       log = "xy",
                       data = Data_Animal_de_pd, 
                       slope.test = 1/3,
                       robust = T
)
summary(Rmax.SMA.family)
plot(Rmax.SMA.family, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Rmax.SMA.family[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Rmax_SMA_diet <- sma(Rmax~massAvg*Diet, log = "xy" ,data = Data_Animal_de_pd, slope.test = 1/3)
summary(Rmax_SMA_diet)
plot(Rmax_SMA_diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Rmax_SMA_diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#
#### Rmin and MASS ####
Rmin.SMA <- sma(Rmin~massAvg, 
                log = "xy",
                data = Data_Animal_de_pd,
                slope.test = 1/3, 
                robust = T
)
summary(Rmin.SMA)
plot(Rmin.SMA)
# Stance
Rmin.SMA.stance <- sma(Rmin~massAvg*Stance, 
                       log = "xy",
                       data = Data_Animal_de_pd,
                       slope.test = 1/3, 
                       robust = T
)
summary(Rmin.SMA.stance)
plot(Rmin.SMA.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Rmin.SMA.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Rmin.SMA.Morphotype <- sma(Rmin~massAvg*Morphotype, 
                           log = "xy",
                           data = Data_Animal_de_pd, 
                           slope.test = 1/3, 
                           robust = T
)
summary(Rmin.SMA.Morphotype)
plot(Rmin.SMA.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Rmin.SMA.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Rmin.SMA.order <- sma(Rmin~massAvg*order, 
                      log = "xy",
                      data = Data_Animal_de_pd, 
                      slope.test = 1/3,
                      robust = T
)
summary(Rmin.SMA.order)
plot(Rmin_SMA_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Rmin.SMA.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Family
Rmin.SMA.family <- sma(Rmin~massAvg*family, 
                       log = "xy",
                       data = Data_Animal_de_pd, 
                       slope.test = 1/3,
                       robust = T
)
summary(Rmin.SMA.family)
plot(Rmin.SMA.family, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Rmin.SMA.family[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Rmin.SMA.diet <- sma(Rmin~massAvg*Diet, 
                     log = "xy",
                     data = Data_Animal_de_pd,
                     slope.test = 1/3,
                     robust = T
)
summary(Rmin.SMA.diet)
plot(Rmin.SMA.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Rmin.SMA.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#
#### DeltaR and MASS ####
DeltaR.SMA <- sma(DeltaR~massAvg, 
                  log = "xy",
                  data = Data_Animal_de_pd,
                  slope.test = 1/3, 
                  robust = T
)
summary(DeltaR.SMA)
plot(DeltaR.SMA)
# Stance
DeltaR.SMA.stance <- sma(DeltaR~massAvg*Stance, 
                         log = "xy",
                         data = Data_Animal_de_pd,
                         slope.test = 1/3, 
                         robust = T
)
summary(DeltaR.SMA.stance)
plot(DeltaR.SMA.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = DeltaR.SMA.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
DeltaR.SMA.Morphotype <- sma(DeltaR~massAvg*Morphotype, 
                             log = "xy",
                             data = Data_Animal_de_pd, 
                             slope.test = 1/3, 
                             robust = T
)
summary(DeltaR.SMA.Morphotype)
plot(DeltaR.SMA.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = DeltaR.SMA.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
DeltaR.SMA.order <- sma(DeltaR~massAvg*order, 
                        log = "xy",
                        data = Data_Animal_de_pd, 
                        slope.test = 1/3,
                        robust = T
)
summary(DeltaR.SMA.order)
plot(DeltaR.SMA.order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = DeltaR.SMA.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Family
DeltaR.SMA.family <- sma(DeltaR~massAvg*family, 
                         log = "xy",
                         data = Data_Animal_de_pd, 
                         slope.test = 1/3,
                         robust = T
)
summary(DeltaR.SMA.family)
plot(DeltaR.SMA.family, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = DeltaR.SMA.family[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
DeltaR.SMA.diet <- sma(DeltaR~massAvg*Diet, 
                       log = "xy",
                       data = Data_Animal_de_pd,
                       slope.test = 1/3,
                       robust = T
)
summary(DeltaR.SMA.diet)
plot(DeltaR.SMA.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = DeltaR.SMA.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### Dmed and MASS ####
Dmed.SMA <- sma(Dmed~massAvg, 
                log = "xy",
                data = Data_Animal_de_pd, 
                slope.test = 1/3, 
                robust = T
)
summary(Dmed.SMA)
coef(Dmed.SMA)
plot(Dmed.SMA)
# Stance
Dmed.SMA.stance <- sma(Dmed~massAvg*Stance, 
                       log = "xy",
                       data = Data_Animal_de_pd,
                       slope.test = 1/3, 
                       robust = T
)
summary(Dmed.SMA.stance)
coef(Dmed.SMA.stance)
plot(Dmed.SMA.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Dmed.SMA.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Dmed.SMA.Morphotype <- sma(Dmed~massAvg*Morphotype, 
                           log = "xy",
                           data = Data_Animal_de_pd, 
                           slope.test = 1/3, 
                           robust = T
)
summary(Dmed.SMA.Morphotype)
plot(Dmed.SMA.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Dmed.SMA.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Dmed.SMA.order <- sma(Dmed~massAvg*order, 
                      log = "xy",
                      data = Data_Animal_de_pd,
                      slope.test = 1/3,
                      robust = T
)
summary(Dmed.SMA.order)
plot(Dmed.SMA.order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Dmed.SMA.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Family
Dmed.SMA.family <- sma(Dmed~massAvg*family, 
                       log = "xy",
                       data = Data_Animal_de_pd, 
                       slope.test = 1/3,
                       robust = T
)
summary(Dmed.SMA.family)
plot(Dmed.SMA.family, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Dmed.SMA.family[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Dmed_SMA_diet <- sma(Dmed~massAvg*Diet, 
                     log = "xy",
                     data = Data_Animal_de_pd,
                     slope.test = 1/3,
                     robust = T
)
summary(Dmed_SMA_diet)
plot(Dmed_SMA_diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Dmed_SMA_diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#
#### angeqv and MASS ####
angeqv.SMA <- sma(angeqv~massAvg, 
                  log = "xy",
                  data = Data_Animal_de_pd,
                  slope.test = 0.1, 
                  #robust = T
)
summary(angeqv.SMA)
plot(angeqv.SMA)
# Stance
angeqv.SMA.stance <- sma(angeqv~massAvg*Stance, 
                         log = "xy",
                         data = Data_Animal_de_pd,
                         slope.test = 1/3, 
                         robust = T
)
summary(angeqv.SMA.stance)
plot(angeqv.SMA.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = angeqv.SMA.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
angeqv.SMA.Morphotype <- sma(angeqv~massAvg*Morphotype, 
                             log = "xy",
                             data = Data_Animal_de_pd, 
                             slope.test = 1/3, 
                             robust = T
)
summary(angeqv.SMA.Morphotype)
plot(angeqv.SMA.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = angeqv.SMA.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
angeqv.SMA.order <- sma(angeqv~massAvg*order, 
                        log = "xy",
                        data = Data_Animal_de_pd, 
                        slope.test = 1/3,
                        robust = T
)
summary(angeqv.SMA.order)
plot(angeqv_SMA_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = angeqv.SMA.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Family
angeqv.SMA.family <- sma(angeqv~massAvg*family, 
                         log = "xy",
                         data = Data_Animal_de_pd, 
                         slope.test = 1/3,
                         robust = T
)
summary(angeqv.SMA.family)
plot(angeqv.SMA.family, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = angeqv.SMA.family[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
angeqv.SMA.diet <- sma(angeqv~massAvg*Diet, 
                       log = "xy",
                       data = Data_Animal_de_pd,
                       slope.test = 1/3,
                       robust = T
)
summary(angeqv.SMA.diet)
plot(angeqv.SMA.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = angeqv.SMA.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
#### Pry_Area and MASS ####
Pry_Area.SMA <- sma(Pry_Area~massAvg, 
                    log = "xy",
                    data = Data_Animal_de_pd,
                    slope.test = 1/3, 
                    robust = T
)
summary(Pry_Area.SMA)
plot(Pry_Area.SMA)
# Stance
Pry_Area.SMA.stance <- sma(Pry_Area~massAvg*Stance, 
                           log = "xy",
                           data = Data_Animal_de_pd,
                           slope.test = 1/3, 
                           robust = T
)
summary(Pry_Area.SMA.stance)
plot(Pry_Area.SMA.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Pry_Area.SMA.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Pry_Area.SMA.Morphotype <- sma(Pry_Area~massAvg*Morphotype, 
                               log = "xy",
                               data = Data_Animal_de_pd, 
                               slope.test = 1/3, 
                               robust = T
)
summary(Pry_Area.SMA.Morphotype)
plot(Pry_Area.SMA.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Pry_Area.SMA.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Pry_Area.SMA.order <- sma(Pry_Area~massAvg*order, 
                          log = "xy",
                          data = Data_Animal_de_pd, 
                          slope.test = 1/3,
                          robust = T
)
summary(Pry_Area.SMA.order)
plot(Pry_Area_SMA_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Pry_Area.SMA.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Family
Pry_Area.SMA.family <- sma(Pry_Area~massAvg*family, 
                           log = "xy",
                           data = Data_Animal_de_pd, 
                           slope.test = 1/3,
                           robust = T
)
summary(Pry_Area.SMA.family)
plot(Pry_Area.SMA.family, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Pry_Area.SMA.family[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Pry_Area.SMA.diet <- sma(Pry_Area~massAvg*Diet, 
                         log = "xy",
                         data = Data_Animal_de_pd,
                         slope.test = 1/3,
                         robust = T
)
summary(Pry_Area.SMA.diet)
plot(Pry_Area.SMA.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Pry_Area.SMA.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### Sf and err MASS ####
Sf.SMA <- sma(Sf~massAvg, 
              log = "xy",
              data = Data_Animal_de_pd,
              slope.test = 1/3, 
              robust = T
)
summary(Sf.SMA)
plot(Sf.SMA)
# Stance
Sf.SMA.stance <- sma(Sf~massAvg*Stance, 
                     log = "xy",
                     data = Data_Animal_de_pd,
                     slope.test = 1/3, 
                     robust = T
)
summary(Sf.SMA.stance)
plot(Sf.SMA.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Sf.SMA.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Sf.SMA.Morphotype <- sma(Sf~massAvg*Morphotype, 
                         log = "xy",
                         data = Data_Animal_de_pd, 
                         slope.test = 1/3, 
                         robust = T
)
summary(Sf.SMA.Morphotype)
plot(Sf.SMA.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Sf.SMA.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Sf.SMA.order <- sma(Sf~massAvg*order, 
                    log = "xy",
                    data = Data_Animal_de_pd, 
                    slope.test = 1/3,
                    robust = T
)
summary(Sf.SMA.order)
plot(Sf_SMA_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Sf.SMA.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Family
Sf.SMA.family <- sma(Sf~massAvg*family, 
                     log = "xy",
                     data = Data_Animal_de_pd, 
                     slope.test = 1/3,
                     robust = T
)
summary(Sf.SMA.family)
plot(Sf.SMA.family, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Sf.SMA.family[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Sf.SMA.diet <- sma(Sf~massAvg*Diet, 
                   log = "xy",
                   data = Data_Animal_de_pd,
                   slope.test = 1/3,
                   robust = T
)
summary(Sf.SMA.diet)
plot(Sf.SMA.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Sf.SMA.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### Mscl_mass_kg err and MASS ####
Mscl_mass_kg.SMA <- sma(Mscl_mass_kg~massAvg, 
                        log = "xy",
                        data = Data_Animal_de_pd,
                        slope.test = 1/3, 
                        robust = T
)
summary(Mscl_mass_kg.SMA)
plot(Mscl_mass_kg.SMA)
# Stance
Mscl_mass_kg.SMA.stance <- sma(Mscl_mass_kg~massAvg*Stance, 
                               log = "xy",
                               data = Data_Animal_de_pd,
                               slope.test = 1/3, 
                               robust = T
)
summary(Mscl_mass_kg.SMA.stance)
plot(Mscl_mass_kg.SMA.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Mscl_mass_kg.SMA.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Mscl_mass_kg.SMA.Morphotype <- sma(Mscl_mass_kg~massAvg*Morphotype, 
                                   log = "xy",
                                   data = Data_Animal_de_pd, 
                                   slope.test = 1/3, 
                                   robust = T
)
summary(Mscl_mass_kg.SMA.Morphotype)
plot(Mscl_mass_kg.SMA.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Mscl_mass_kg.SMA.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Mscl_mass_kg.SMA.order <- sma(Mscl_mass_kg~massAvg*order, 
                              log = "xy",
                              data = Data_Animal_de_pd, 
                              slope.test = 1/3,
                              robust = T
)
summary(Mscl_mass_kg.SMA.order)
plot(Mscl_mass_kg_SMA_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Mscl_mass_kg.SMA.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Family
Mscl_mass_kg.SMA.family <- sma(Mscl_mass_kg~massAvg*family, 
                               log = "xy",
                               data = Data_Animal_de_pd, 
                               slope.test = 1/3,
                               robust = T
)
summary(Mscl_mass_kg.SMA.family)
plot(Mscl_mass_kg.SMA.family, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Mscl_mass_kg.SMA.family[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Mscl_mass_kg.SMA.diet <- sma(Mscl_mass_kg~massAvg*Diet, 
                             log = "xy",
                             data = Data_Animal_de_pd,
                             slope.test = 1/3,
                             robust = T
)
summary(Mscl_mass_kg.SMA.diet)
plot(Mscl_mass_kg.SMA.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Mscl_mass_kg.SMA.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### Mscl_length_m err and MASS ####
Mscl_length_m.SMA <- sma(Mscl_length_m~massAvg, 
                         log = "xy",
                         data = Data_Animal_de_pd,
                         slope.test = 1/3, 
                         robust = T
)
summary(Mscl_length_m.SMA)
plot(Mscl_length_m.SMA)
# Stance
Mscl_length_m.SMA.stance <- sma(Mscl_length_m~massAvg*Stance, 
                                log = "xy",
                                data = Data_Animal_de_pd,
                                slope.test = 1/3, 
                                robust = T
)
summary(Mscl_length_m.SMA.stance)
plot(Mscl_length_m.SMA.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Mscl_length_m.SMA.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Mscl_length_m.SMA.Morphotype <- sma(Mscl_length_m~massAvg*Morphotype, 
                                    log = "xy",
                                    data = Data_Animal_de_pd, 
                                    slope.test = 1/3, 
                                    robust = T
)
summary(Mscl_length_m.SMA.Morphotype)
plot(Mscl_length_m.SMA.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Mscl_length_m.SMA.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Mscl_length_m.SMA.order <- sma(Mscl_length_m~massAvg*order, 
                               log = "xy",
                               data = Data_Animal_de_pd, 
                               slope.test = 1/3,
                               robust = T
)
summary(Mscl_length_m.SMA.order)
plot(Mscl_length_m_SMA_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Mscl_length_m.SMA.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Family
Mscl_length_m.SMA.family <- sma(Mscl_length_m~massAvg*family, 
                                log = "xy",
                                data = Data_Animal_de_pd, 
                                slope.test = 1/3,
                                robust = T
)
summary(Mscl_length_m.SMA.family)
plot(Mscl_length_m.SMA.family, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Mscl_length_m.SMA.family[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Mscl_length_m.SMA.diet <- sma(Mscl_length_m~massAvg*Diet, 
                              log = "xy",
                              data = Data_Animal_de_pd,
                              slope.test = 1/3,
                              robust = T
)
summary(Mscl_length_m.SMA.diet)
plot(Mscl_length_m.SMA.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Mscl_length_m.SMA.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### MoI and err MASS ####
MoI.SMA <- sma(MoI~massAvg, 
               log = "xy",
               data = Data_Animal_de_pd,
               slope.test = 1/3, 
               robust = T
)
summary(MoI.SMA)
plot(MoI.SMA)
# Stance
MoI.SMA.stance <- sma(MoI~massAvg*Stance, 
                      log = "xy",
                      data = Data_Animal_de_pd,
                      slope.test = 1/3, 
                      robust = T
)
summary(MoI.SMA.stance)
plot(MoI.SMA.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI.SMA.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
MoI.SMA.Morphotype <- sma(MoI~massAvg*Morphotype, 
                          log = "xy",
                          data = Data_Animal_de_pd, 
                          slope.test = 1/3, 
                          robust = T
)
summary(MoI.SMA.Morphotype)
plot(MoI.SMA.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI.SMA.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
MoI.SMA.order <- sma(MoI~massAvg*order, 
                     log = "xy",
                     data = Data_Animal_de_pd, 
                     slope.test = 1/3,
                     robust = T
)
summary(MoI.SMA.order)
plot(MoI_SMA_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI.SMA.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Family
MoI.SMA.family <- sma(MoI~massAvg*family, 
                      log = "xy",
                      data = Data_Animal_de_pd, 
                      slope.test = 1/3,
                      robust = T
)
summary(MoI.SMA.family)
plot(MoI.SMA.family, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI.SMA.family[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
MoI.SMA.diet <- sma(MoI~massAvg*Diet, 
                    log = "xy",
                    data = Data_Animal_de_pd,
                    slope.test = 1/3,
                    robust = T
)
summary(MoI.SMA.diet)
plot(MoI.SMA.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI.SMA.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### Limb_mass err and MASS ####
Limb_mass.SMA <- sma(Limb_mass~massAvg, 
                     log = "xy",
                     data = Data_Animal_de_pd,
                     slope.test = 1/3, 
                     robust = T
)
summary(Limb_mass.SMA)
plot(Limb_mass.SMA)
# Stance
Limb_mass.SMA.stance <- sma(Limb_mass~massAvg*Stance, 
                            log = "xy",
                            data = Data_Animal_de_pd,
                            slope.test = 1/3, 
                            robust = T
)
summary(Limb_mass.SMA.stance)
plot(Limb_mass.SMA.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Limb_mass.SMA.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Limb_mass.SMA.Morphotype <- sma(Limb_mass~massAvg*Morphotype, 
                                log = "xy",
                                data = Data_Animal_de_pd, 
                                slope.test = 1/3, 
                                robust = T
)
summary(Limb_mass.SMA.Morphotype)
plot(Limb_mass.SMA.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Limb_mass.SMA.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Limb_mass.SMA.order <- sma(Limb_mass~massAvg*order, 
                           log = "xy",
                           data = Data_Animal_de_pd, 
                           slope.test = 1/3,
                           robust = T
)
summary(Limb_mass.SMA.order)
plot(Limb_mass_SMA_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Limb_mass.SMA.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Family
Limb_mass.SMA.family <- sma(Limb_mass~massAvg*family, 
                            log = "xy",
                            data = Data_Animal_de_pd, 
                            slope.test = 1/3,
                            robust = T
)
summary(Limb_mass.SMA.family)
plot(Limb_mass.SMA.family, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Limb_mass.SMA.family[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Limb_mass.SMA.diet <- sma(Limb_mass~massAvg*Diet, 
                          log = "xy",
                          data = Data_Animal_de_pd,
                          slope.test = 1/3,
                          robust = T
)
summary(Limb_mass.SMA.diet)
plot(Limb_mass.SMA.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Limb_mass.SMA.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### CoM and err MASS ####
CoM.SMA <- sma(CoM~massAvg, 
               log = "xy",
               data = Data_Animal_de_pd,
               slope.test = 1/3, 
               robust = T
)
summary(CoM.SMA)
plot(CoM.SMA)
# Stance
CoM.SMA.stance <- sma(CoM~massAvg*Stance, 
                      log = "xy",
                      data = Data_Animal_de_pd,
                      slope.test = 1/3, 
                      robust = T
)
summary(CoM.SMA.stance)
plot(CoM.SMA.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = CoM.SMA.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
CoM.SMA.Morphotype <- sma(CoM~massAvg*Morphotype, 
                          log = "xy",
                          data = Data_Animal_de_pd, 
                          slope.test = 1/3, 
                          robust = T
)
summary(CoM.SMA.Morphotype)
plot(CoM.SMA.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = CoM.SMA.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
CoM.SMA.order <- sma(CoM~massAvg*order, 
                     log = "xy",
                     data = Data_Animal_de_pd, 
                     slope.test = 1/3,
                     robust = T
)
summary(CoM.SMA.order)
plot(CoM_SMA_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = CoM.SMA.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Family
CoM.SMA.family <- sma(CoM~massAvg*family, 
                      log = "xy",
                      data = Data_Animal_de_pd, 
                      slope.test = 1/3,
                      robust = T
)
summary(CoM.SMA.family)
plot(CoM.SMA.family, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = CoM.SMA.family[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
CoM.SMA.diet <- sma(CoM~massAvg*Diet, 
                    log = "xy",
                    data = Data_Animal_de_pd,
                    slope.test = 1/3,
                    robust = T
)
summary(CoM.SMA.diet)
plot(CoM.SMA.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = CoM.SMA.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### MoA_m and err MASS ####
MoA_m.SMA <- sma(MoA_m~massAvg, 
                 log = "xy",
                 data = Data_Animal_de_pd,
                 slope.test = 1/3, 
                 robust = T
)
summary(MoA_m.SMA)
plot(MoA_m.SMA)
# Stance
MoA_m.SMA.stance <- sma(MoA_m~massAvg*Stance, 
                        log = "xy",
                        data = Data_Animal_de_pd,
                        slope.test = 1/3, 
                        robust = T
)
summary(MoA_m.SMA.stance)
plot(MoA_m.SMA.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoA_m.SMA.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
MoA_m.SMA.Morphotype <- sma(MoA_m~massAvg*Morphotype, 
                            log = "xy",
                            data = Data_Animal_de_pd, 
                            slope.test = 1/3, 
                            robust = T
)
summary(MoA_m.SMA.Morphotype)
plot(MoA_m.SMA.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoA_m.SMA.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
MoA_m.SMA.order <- sma(MoA_m~massAvg*order, 
                       log = "xy",
                       data = Data_Animal_de_pd, 
                       slope.test = 1/3,
                       robust = T
)
summary(MoA_m.SMA.order)
plot(MoA_m_SMA_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoA_m.SMA.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Family
MoA_m.SMA.family <- sma(MoA_m~massAvg*family, 
                        log = "xy",
                        data = Data_Animal_de_pd, 
                        slope.test = 1/3,
                        robust = T
)
summary(MoA_m.SMA.family)
plot(MoA_m.SMA.family, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoA_m.SMA.family[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
MoA_m.SMA.diet <- sma(MoA_m~massAvg*Diet, 
                      log = "xy",
                      data = Data_Animal_de_pd,
                      slope.test = 1/3,
                      robust = T
)
summary(MoA_m.SMA.diet)
plot(MoA_m.SMA.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoA_m.SMA.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### P_eq and MASS ####
P_eq.SMA <- sma(P_eq~massAvg, 
                log = "xy",
                data = Data_Animal_de_pd,
                slope.test = 1/3, 
                robust = T
)
summary(P_eq.SMA)
plot(P_eq.SMA)
# Stance
P_eq.SMA.stance <- sma(P_eq~massAvg*Stance, 
                       log = "xy",
                       data = Data_Animal_de_pd,
                       slope.test = 1/3, 
                       robust = T
)
summary(P_eq.SMA.stance)
plot(P_eq.SMA.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = P_eq.SMA.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
P_eq.SMA.Morphotype <- sma(P_eq~massAvg*Morphotype, 
                           log = "xy",
                           data = Data_Animal_de_pd, 
                           slope.test = 1/3, 
                           robust = T
)
summary(P_eq.SMA.Morphotype)
plot(P_eq.SMA.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = P_eq.SMA.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
P_eq.SMA.order <- sma(P_eq~massAvg*order, 
                      log = "xy",
                      data = Data_Animal_de_pd, 
                      slope.test = 1/3,
                      robust = T
)
summary(P_eq.SMA.order)
plot(P_eq_SMA_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = P_eq.SMA.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Family
P_eq.SMA.family <- sma(P_eq~massAvg*family, 
                       log = "xy",
                       data = Data_Animal_de_pd, 
                       slope.test = 1/3,
                       robust = T
)
summary(P_eq.SMA.family)
plot(P_eq.SMA.family, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = P_eq.SMA.family[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
P_eq.SMA.diet <- sma(P_eq~massAvg*Diet, 
                     log = "xy",
                     data = Data_Animal_de_pd,
                     slope.test = 1/3,
                     robust = T
)
summary(P_eq.SMA.diet)
plot(P_eq.SMA.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = P_eq.SMA.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### MoI_G and err MASS ####
MoI_G.SMA <- sma(MoI_G~massAvg, 
                 log = "xy",
                 data = Data_Animal_de_pd,
                 slope.test = 1/3, 
                 robust = T
)
summary(MoI_G.SMA)
plot(MoI_G.SMA)
# Stance
MoI_G.SMA.stance <- sma(MoI_G~massAvg*Stance, 
                        log = "xy",
                        data = Data_Animal_de_pd,
                        slope.test = 1/3, 
                        robust = T
)
summary(MoI_G.SMA.stance)
plot(MoI_G.SMA.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_G.SMA.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
MoI_G.SMA.Morphotype <- sma(MoI_G~massAvg*Morphotype, 
                            log = "xy",
                            data = Data_Animal_de_pd, 
                            slope.test = 1/3, 
                            robust = T
)
summary(MoI_G.SMA.Morphotype)
plot(MoI_G.SMA.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_G.SMA.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
MoI_G.SMA.order <- sma(MoI_G~massAvg*order, 
                       log = "xy",
                       data = Data_Animal_de_pd, 
                       slope.test = 1/3,
                       robust = T
)
summary(MoI_G.SMA.order)
plot(MoI_G_SMA_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_G.SMA.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Family
MoI_G.SMA.family <- sma(MoI_G~massAvg*family, 
                        log = "xy",
                        data = Data_Animal_de_pd, 
                        slope.test = 1/3,
                        robust = T
)
summary(MoI_G.SMA.family)
plot(MoI_G.SMA.family, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_G.SMA.family[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
MoI_G.SMA.diet <- sma(MoI_G~massAvg*Diet, 
                      log = "xy",
                      data = Data_Animal_de_pd,
                      slope.test = 1/3,
                      robust = T
)
summary(MoI_G.SMA.diet)
plot(MoI_G.SMA.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_G.SMA.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### MoI_H and MASS ####
MoI_H.SMA <- sma(MoI_H~massAvg, 
                 log = "xy",
                 data = Data_Animal_de_pd,
                 slope.test = 1/3, 
                 robust = T
)
summary(MoI_H.SMA)
plot(MoI_H.SMA)
# Stance
MoI_H.SMA.stance <- sma(MoI_H~massAvg*Stance, 
                        log = "xy",
                        data = Data_Animal_de_pd,
                        slope.test = 1/3, 
                        robust = T
)
summary(MoI_H.SMA.stance)
plot(MoI_H.SMA.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_H.SMA.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
MoI_H.SMA.Morphotype <- sma(MoI_H~massAvg*Morphotype, 
                            log = "xy",
                            data = Data_Animal_de_pd, 
                            slope.test = 1/3, 
                            robust = T
)
summary(MoI_H.SMA.Morphotype)
plot(MoI_H.SMA.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_H.SMA.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
MoI_H.SMA.order <- sma(MoI_H~massAvg*order, 
                       log = "xy",
                       data = Data_Animal_de_pd, 
                       slope.test = 1/3,
                       robust = T
)
summary(MoI_H.SMA.order)
plot(MoI_H_SMA_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_H.SMA.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Family
MoI_H.SMA.family <- sma(MoI_H~massAvg*family, 
                        log = "xy",
                        data = Data_Animal_de_pd, 
                        slope.test = 1/3,
                        robust = T
)
summary(MoI_H.SMA.family)
plot(MoI_H.SMA.family, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_H.SMA.family[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
MoI_H.SMA.diet <- sma(MoI_H~massAvg*Diet, 
                      log = "xy",
                      data = Data_Animal_de_pd,
                      slope.test = 1/3,
                      robust = T
)
summary(MoI_H.SMA.diet)
plot(MoI_H.SMA.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_H.SMA.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
#### Alpha and Mass ####
Alpha.SMA <- sma(Alpha~massAvg, 
                 log = "xy",
                 data = Data_Animal_de_pd,
                 slope.test = 1/3, 
                 robust = T
)
summary(Alpha.SMA)
plot(Alpha.SMA)
# Stance
Alpha.SMA.stance <- sma(Alpha~massAvg*Stance, 
                        log = "xy",
                        data = Data_Animal_de_pd,
                        slope.test = 1/3, 
                        robust = T
)
summary(Alpha.SMA.stance)
plot(Alpha.SMA.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Alpha.SMA.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Alpha.SMA.Morphotype <- sma(Alpha~massAvg*Morphotype, 
                            log = "xy",
                            data = Data_Animal_de_pd, 
                            slope.test = 1/3, 
                            robust = T
)
summary(Alpha.SMA.Morphotype)
plot(Alpha.SMA.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Alpha.SMA.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Alpha.SMA.order <- sma(Alpha~massAvg*order, 
                       log = "xy",
                       data = Data_Animal_de_pd, 
                       slope.test = 1/3,
                       robust = T
)
summary(Alpha.SMA.order)
plot(Alpha_SMA_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Alpha.SMA.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Family
Alpha.SMA.family <- sma(Alpha~massAvg*family, 
                        log = "xy",
                        data = Data_Animal_de_pd, 
                        slope.test = 1/3,
                        robust = T
)
summary(Alpha.SMA.family)
plot(Alpha.SMA.family, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Alpha.SMA.family[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Alpha.SMA.diet <- sma(Alpha~massAvg*Diet, 
                      log = "xy",
                      data = Data_Animal_de_pd,
                      slope.test = 1/3,
                      robust = T
)
summary(Alpha.SMA.diet)
plot(Alpha.SMA.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Alpha.SMA.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
#### Omg and Mass ####
Omg.SMA <- sma(Omg~massAvg, 
               log = "xy",
               data = Data_Animal_de_pd,
               slope.test = 1/3, 
               robust = T
)
summary(Omg.SMA)
plot(Omg.SMA)
# Stance
Omg.SMA.stance <- sma(Omg~massAvg*Stance, 
                      log = "xy",
                      data = Data_Animal_de_pd,
                      slope.test = 1/3, 
                      robust = T
)
summary(Omg.SMA.stance)
plot(Omg.SMA.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Omg.SMA.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Omg.SMA.Morphotype <- sma(Omg~massAvg*Morphotype, 
                          log = "xy",
                          data = Data_Animal_de_pd, 
                          slope.test = 1/3, 
                          robust = T
)
summary(Omg.SMA.Morphotype)
plot(Omg.SMA.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Omg.SMA.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Omg.SMA.order <- sma(Omg~massAvg*order, 
                     log = "xy",
                     data = Data_Animal_de_pd, 
                     slope.test = 1/3,
                     robust = T
)
summary(Omg.SMA.order)
plot(Omg_SMA_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Omg.SMA.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Family
Omg.SMA.family <- sma(Omg~massAvg*family, 
                      log = "xy",
                      data = Data_Animal_de_pd, 
                      slope.test = 1/3,
                      robust = T
)
summary(Omg.SMA.family)
plot(Omg.SMA.family, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Omg.SMA.family[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Omg.SMA.diet <- sma(Omg~massAvg*Diet, 
                    log = "xy",
                    data = Data_Animal_de_pd,
                    slope.test = 1/3,
                    robust = T
)
summary(Omg.SMA.diet)
plot(Omg.SMA.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Omg.SMA.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
#### Vsup and Mass ####
Vsup.SMA <- sma(Vsup~massAvg, 
                log = "xy",
                data = Data_Animal_de_pd,
                slope.test = 1/3, 
                robust = T
)
summary(Vsup.SMA)
plot(Vsup.SMA)
# Stance
Vsup.SMA.stance <- sma(Vsup~massAvg*Stance, 
                       log = "xy",
                       data = Data_Animal_de_pd,
                       slope.test = 1/3, 
                       robust = T
)
summary(Vsup.SMA.stance)
plot(Vsup.SMA.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Vsup.SMA.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Vsup.SMA.Morphotype <- sma(Vsup~massAvg*Morphotype, 
                           log = "xy",
                           data = Data_Animal_de_pd, 
                           slope.test = 1/3, 
                           robust = T
)
summary(Vsup.SMA.Morphotype)
plot(Vsup.SMA.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Vsup.SMA.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Vsup.SMA.order <- sma(Vsup~massAvg*order, 
                      log = "xy",
                      data = Data_Animal_de_pd, 
                      slope.test = 1/3,
                      robust = T
)
summary(Vsup.SMA.order)
plot(Vsup_SMA_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Vsup.SMA.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Family
Vsup.SMA.family <- sma(Vsup~massAvg*family, 
                       log = "xy",
                       data = Data_Animal_de_pd, 
                       slope.test = 1/3,
                       robust = T
)
summary(Vsup.SMA.family)
plot(Vsup.SMA.family, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Vsup.SMA.family[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Vsup.SMA.diet <- sma(Vsup~massAvg*Diet, 
                     log = "xy",
                     data = Data_Animal_de_pd,
                     slope.test = 1/3,
                     robust = T
)
summary(Vsup.SMA.diet)
plot(Vsup.SMA.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Vsup.SMA.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### MoI_HC and err MASS ####
MoI_HC.SMA <- sma(MoI_HC~massAvg, 
                  log = "xy",
                  data = Data_Animal_de_pd,
                  slope.test = 1/3, 
                  robust = T
)
summary(MoI_HC.SMA)
plot(MoI_HC.SMA)
# Stance
MoI_HC.SMA.stance <- sma(MoI_HC~massAvg*Stance, 
                         log = "xy",
                         data = Data_Animal_de_pd,
                         slope.test = 1/3, 
                         robust = T
)
summary(MoI_HC.SMA.stance)
plot(MoI_HC.SMA.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_HC.SMA.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
MoI_HC.SMA.Morphotype <- sma(MoI_HC~massAvg*Morphotype, 
                             log = "xy",
                             data = Data_Animal_de_pd, 
                             slope.test = 1/3, 
                             robust = T
)
summary(MoI_HC.SMA.Morphotype)
plot(MoI_HC.SMA.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_HC.SMA.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
MoI_HC.SMA.order <- sma(MoI_HC~massAvg*order, 
                        log = "xy",
                        data = Data_Animal_de_pd, 
                        slope.test = 1/3,
                        robust = T
)
summary(MoI_HC.SMA.order)
plot(MoI_HC_SMA_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_HC.SMA.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Family
MoI_HC.SMA.family <- sma(MoI_HC~massAvg*family, 
                         log = "xy",
                         data = Data_Animal_de_pd, 
                         slope.test = 1/3,
                         robust = T
)
summary(MoI_HC.SMA.family)
plot(MoI_HC.SMA.family, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_HC.SMA.family[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
MoI_HC.SMA.diet <- sma(MoI_HC~massAvg*Diet, 
                       log = "xy",
                       data = Data_Animal_de_pd,
                       slope.test = 1/3,
                       robust = T
)
summary(MoI_HC.SMA.diet)
plot(MoI_HC.SMA.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_HC.SMA.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### Alpha_C and err MASS ####
Alpha_C.SMA <- sma(Alpha_C~massAvg, 
                   log = "xy",
                   data = Data_Animal_de_pd,
                   slope.test = 1/3, 
                   robust = T
)
summary(Alpha_C.SMA)
plot(Alpha_C.SMA)
# Stance
Alpha_C.SMA.stance <- sma(Alpha_C~massAvg*Stance, 
                          log = "xy",
                          data = Data_Animal_de_pd,
                          slope.test = 1/3, 
                          robust = T
)
summary(Alpha_C.SMA.stance)
plot(Alpha_C.SMA.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Alpha_C.SMA.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Alpha_C.SMA.Morphotype <- sma(Alpha_C~massAvg*Morphotype, 
                              log = "xy",
                              data = Data_Animal_de_pd, 
                              slope.test = 1/3, 
                              robust = T
)
summary(Alpha_C.SMA.Morphotype)
plot(Alpha_C.SMA.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Alpha_C.SMA.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Alpha_C.SMA.order <- sma(Alpha_C~massAvg*order, 
                         log = "xy",
                         data = Data_Animal_de_pd, 
                         slope.test = 1/3,
                         robust = T
)
summary(Alpha_C.SMA.order)
plot(Alpha_C_SMA_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Alpha_C.SMA.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Family
Alpha_C.SMA.family <- sma(Alpha_C~massAvg*family, 
                          log = "xy",
                          data = Data_Animal_de_pd, 
                          slope.test = 1/3,
                          robust = T
)
summary(Alpha_C.SMA.family)
plot(Alpha_C.SMA.family, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Alpha_C.SMA.family[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Alpha_C.SMA.diet <- sma(Alpha_C~massAvg*Diet, 
                        log = "xy",
                        data = Data_Animal_de_pd,
                        slope.test = 1/3,
                        robust = T
)
summary(Alpha_C.SMA.diet)
plot(Alpha_C.SMA.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Alpha_C.SMA.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### OmgC and err MASS ####
OmgC.SMA <- sma(OmgC~massAvg, 
                log = "xy",
                data = Data_Animal_de_pd,
                slope.test = 1/3, 
                robust = T
)
summary(OmgC.SMA)
plot(OmgC.SMA)
# Stance
OmgC.SMA.stance <- sma(OmgC~massAvg*Stance, 
                       log = "xy",
                       data = Data_Animal_de_pd,
                       slope.test = 1/3, 
                       robust = T
)
summary(OmgC.SMA.stance)
plot(OmgC.SMA.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = OmgC.SMA.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
OmgC.SMA.Morphotype <- sma(OmgC~massAvg*Morphotype, 
                           log = "xy",
                           data = Data_Animal_de_pd, 
                           slope.test = 1/3, 
                           robust = T
)
summary(OmgC.SMA.Morphotype)
plot(OmgC.SMA.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = OmgC.SMA.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
OmgC.SMA.order <- sma(OmgC~massAvg*order, 
                      log = "xy",
                      data = Data_Animal_de_pd, 
                      slope.test = 1/3,
                      robust = T
)
summary(OmgC.SMA.order)
plot(OmgC_SMA_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = OmgC.SMA.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Family
OmgC.SMA.family <- sma(OmgC~massAvg*family, 
                       log = "xy",
                       data = Data_Animal_de_pd, 
                       slope.test = 1/3,
                       robust = T
)
summary(OmgC.SMA.family)
plot(OmgC.SMA.family, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = OmgC.SMA.family[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
OmgC.SMA.diet <- sma(OmgC~massAvg*Diet, 
                     log = "xy",
                     data = Data_Animal_de_pd,
                     slope.test = 1/3,
                     robust = T
)
summary(OmgC.SMA.diet)
plot(OmgC.SMA.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = OmgC.SMA.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### VsupC and Mass ####
VsupC.SMA <- sma(VsupC~massAvg, 
                 log = "xy",
                 data = Data_Animal_de_pd,
                 slope.test = 1/3, 
                 robust = T
)
summary(VsupC.SMA)
plot(VsupC.SMA)
# Stance
VsupC.SMA.stance <- sma(VsupC~massAvg*Stance, 
                        log = "xy",
                        data = Data_Animal_de_pd,
                        slope.test = 1/3, 
                        robust = T
)
summary(VsupC.SMA.stance)
plot(VsupC.SMA.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = VsupC.SMA.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
VsupC.SMA.Morphotype <- sma(VsupC~massAvg*Morphotype, 
                            log = "xy",
                            data = Data_Animal_de_pd, 
                            slope.test = 1/3, 
                            robust = T
)
summary(VsupC.SMA.Morphotype)
plot(VsupC.SMA.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = VsupC.SMA.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
VsupC.SMA.order <- sma(VsupC~massAvg*order, 
                       log = "xy",
                       data = Data_Animal_de_pd, 
                       slope.test = 1/3,
                       robust = T
)
summary(VsupC.SMA.order)
plot(VsupC_SMA_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = VsupC.SMA.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Family
VsupC.SMA.family <- sma(VsupC~massAvg*family, 
                        log = "xy",
                        data = Data_Animal_de_pd, 
                        slope.test = 1/3,
                        robust = T
)
summary(VsupC.SMA.family)
plot(VsupC.SMA.family, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = VsupC.SMA.family[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
VsupC.SMA.diet <- sma(VsupC~massAvg*Diet, 
                      log = "xy",
                      data = Data_Animal_de_pd,
                      slope.test = 1/3,
                      robust = T
)
summary(VsupC.SMA.diet)
plot(VsupC.SMA.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = VsupC.SMA.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
#### Maxspeed_ms and MASS ####
Maxspeed_ms.SMA <- sma(Maxspeed_ms~massAvg, 
                       log = "xy",
                       data = Data_Animal_de_pd,
                       slope.test = 1/3, 
                       robust = T
)
summary(Maxspeed_ms.SMA)
plot(Maxspeed_ms.SMA)
# Stance
Maxspeed_ms.SMA.stance <- sma(Maxspeed_ms~massAvg*Stance, 
                              log = "xy",
                              data = Data_Animal_de_pd,
                              slope.test = 1/3, 
                              robust = T
)
summary(Maxspeed_ms.SMA.stance)
plot(Maxspeed_ms.SMA.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Maxspeed_ms.SMA.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Maxspeed_ms.SMA.Morphotype <- sma(Maxspeed_ms~massAvg*Morphotype, 
                                  log = "xy",
                                  data = Data_Animal_de_pd, 
                                  slope.test = 1/3, 
                                  robust = T
)
summary(Maxspeed_ms.SMA.Morphotype)
plot(Maxspeed_ms.SMA.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Maxspeed_ms.SMA.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Maxspeed_ms.SMA.order <- sma(Maxspeed_ms~massAvg*order, 
                             log = "xy",
                             data = Data_Animal_de_pd, 
                             slope.test = 1/3,
                             robust = T
)
summary(Maxspeed_ms.SMA.order)
plot(Maxspeed_ms.SMA.order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Maxspeed_ms.SMA.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Family
Maxspeed_ms.SMA.family <- sma(Maxspeed_ms~massAvg*family, 
                              log = "xy",
                              data = Data_Animal_de_pd, 
                              slope.test = 1/3,
                              robust = T
)
summary(Maxspeed_ms.SMA.family)
plot(Maxspeed_ms.SMA.family, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Maxspeed_ms.SMA.family[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Maxspeed_ms.SMA.diet <- sma(Maxspeed_ms~massAvg*Diet, 
                            log = "xy",
                            data = Data_Animal_de_pd,
                            slope.test = 1/3,
                            robust = T
)
summary(Maxspeed_ms.SMA.diet)
plot(Maxspeed_ms.SMA.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Maxspeed_ms.SMA.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)


################################################################################
# ALL and LONG #################################
#### La and LONG ####
La.SMAlong <- sma(La~Long, 
                  log = "xy",
                  data = Data_Animal_de_pd,
                  slope.test = 1/3,
                  robust = T
)
summary(La.SMAlong)
plot(La.SMAlong)

# Stance
La.SMAlong.stance <- sma(La~Long*Stance,
                         log = "xy",
                         data = Data_Animal_de_pd, 
                         slope.test = 1/3,
                         robust = T
)
summary(La.SMAlong.stance)
plot(La.SMAlong.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = La.SMAlong.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype 
La.SMAlong.Morphotype <- sma(La~Long*Morphotype, 
                             log = "xy",
                             data = Data_Animal_de_pd, 
                             slope.test = 1/3, 
                             robust = T
)
summary(La.SMAlong.Morphotype)
plot(La.SMAlong.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = La.SMAlong.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
La.SMAlong.order <- sma(La~Long*order,
                        log = "xy",
                        data = Data_Animal_de_pd,
                        slope.test = 1/3,
                        robust = T
)
summary(La.SMAlong.order)
plot(La.SMAlong.order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = La.SMAlong.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
La_SMAlong_diet <- sma(La~Long*Diet, 
                       log = "xy",
                       data = Data_Animal_de_pd,
                       slope.test = 1/3, 
                       robust = T
)
summary(La_SMAlong_diet)
plot(La_SMAlong_diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = La_SMAlong_diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#
#### Rmax and LONG ####
Rmax.SMAlong <- sma(Rmax~Long, 
                    log = "xy",
                    data = Data_Animal_de_pd, 
                    slope.test = 1/3, 
                    robust = T
)
summary(Rmax.SMAlong)
plot(Rmax.SMAlong)
# Stance
Rmax.SMAlong.stance <- sma(Rmax~Long*Stance, 
                           log = "xy",
                           data = Data_Animal_de_pd,
                           slope.test = 1/3, 
                           robust = T
)
summary(Rmax.SMAlong.stance)
plot(Rmax.SMAlong.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Rmax.SMAlong.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Rmax.SMAlong.Morphotype <- sma(Rmax~Long*Morphotype, 
                               log = "xy",
                               data = Data_Animal_de_pd, 
                               slope.test = 1/3, 
                               #robust = T
)
summary(Rmax.SMAlong.Morphotype)
plot(Rmax.SMAlong.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Rmax.SMAlong.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Rmax.SMAlong.order <- sma(Rmax~Long*order,
                          log = "xy",
                          data = Data_Animal_de_pd,
                          slope.test = 1/3,
                          robust = T
)
summary(Rmax_SMAlong_order)
plot(Rmax_SMAlong_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Rmax_SMAlong_order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Rmax_SMAlong_diet <- sma(Rmax~Long*Diet, log = "xy" ,data = Data_Animal_de_pd, slope.test = 1/3)
summary(Rmax_SMAlong_diet)
plot(Rmax_SMAlong_diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Rmax_SMAlong_diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#
#### Rmin and LONG ####
Rmin.SMAlong <- sma(Rmin~Long, 
                    log = "xy",
                    data = Data_Animal_de_pd,
                    slope.test = 1/3, 
                    robust = T
)
summary(Rmin.SMAlong)
plot(Rmin.SMAlong)
# Stance
Rmin.SMAlong.stance <- sma(Rmin~Long*Stance, 
                           log = "xy",
                           data = Data_Animal_de_pd,
                           slope.test = 1/3, 
                           robust = T
)
summary(Rmin.SMAlong.stance)
plot(Rmin.SMAlong.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Rmin.SMAlong.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Rmin.SMAlong.Morphotype <- sma(Rmin~Long*Morphotype, 
                               log = "xy",
                               data = Data_Animal_de_pd, 
                               slope.test = 1/3, 
                               robust = T
)
summary(Rmin.SMAlong.Morphotype)
plot(Rmin.SMAlong.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Rmin.SMAlong.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Rmin.SMAlong.order <- sma(Rmin~Long*order, 
                          log = "xy",
                          data = Data_Animal_de_pd, 
                          slope.test = 1/3,
                          robust = T
)
summary(Rmin.SMAlong.order)
plot(Rmin_SMAlong_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Rmin.SMAlong.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Rmin.SMAlong.diet <- sma(Rmin~Long*Diet, 
                         log = "xy",
                         data = Data_Animal_de_pd,
                         slope.test = 1/3,
                         robust = T
)
summary(Rmin.SMAlong.diet)
plot(Rmin.SMAlong.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Rmin.SMAlong.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#
#### DeltaR and LONG ####
DeltaR.SMAlong <- sma(DeltaR~Long, 
                      log = "xy",
                      data = Data_Animal_de_pd,
                      slope.test = 1/3, 
                      robust = T
)
summary(DeltaR.SMAlong)
plot(DeltaR.SMAlong)
# Stance
DeltaR.SMAlong.stance <- sma(DeltaR~Long*Stance, 
                             log = "xy",
                             data = Data_Animal_de_pd,
                             slope.test = 1/3, 
                             robust = T
)
summary(DeltaR.SMAlong.stance)
plot(DeltaR.SMAlong.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = DeltaR.SMAlong.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
DeltaR.SMAlong.Morphotype <- sma(DeltaR~Long*Morphotype, 
                                 log = "xy",
                                 data = Data_Animal_de_pd, 
                                 slope.test = 1/3, 
                                 robust = T
)
summary(DeltaR.SMAlong.Morphotype)
plot(DeltaR.SMAlong.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = DeltaR.SMAlong.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
DeltaR.SMAlong.order <- sma(DeltaR~Long*order, 
                            log = "xy",
                            data = Data_Animal_de_pd, 
                            slope.test = 1/3,
                            robust = T
)
summary(DeltaR.SMAlong.order)
plot(DeltaR_SMAlong_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = DeltaR.SMAlong.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
DeltaR.SMAlong.diet <- sma(DeltaR~Long*Diet, 
                           log = "xy",
                           data = Data_Animal_de_pd,
                           slope.test = 1/3,
                           robust = T
)
summary(DeltaR.SMAlong.diet)
plot(DeltaR.SMAlong.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = DeltaR.SMAlong.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### Dmed and LONG ####
Dmed.SMAlong <- sma(Dmed~Long, 
                    log = "xy",
                    data = Data_Animal_de_pd, 
                    slope.test = 1/3, 
                    robust = T
)
summary(Dmed.SMAlong)
coef(Dmed.SMAlong)
plot(Dmed.SMAlong)
# Stance
Dmed.SMAlong.stance <- sma(Dmed~Long*Stance, 
                           log = "xy",
                           data = Data_Animal_de_pd,
                           slope.test = 1/3, 
                           robust = T
)
summary(Dmed.SMAlong.stance)
coef(Dmed.SMAlong.stance)
plot(Dmed.SMAlong.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Dmed.SMAlong.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Dmed.SMAlong.Morphotype <- sma(Dmed~Long*Morphotype, 
                               log = "xy",
                               data = Data_Animal_de_pd, 
                               slope.test = 1/3, 
                               robust = T
)
summary(Dmed.SMAlong.Morphotype)
plot(Dmed.SMAlong.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Dmed.SMAlong.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Dmed.SMAlong.order <- sma(Dmed~Long*order, 
                          log = "xy",
                          data = Data_Animal_de_pd,
                          slope.test = 1/3,
                          robust = T
)
summary(Dmed.SMAlong.order)
plot(Dmed.SMAlong.order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Dmed.SMAlong.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Dmed_SMAlong_diet <- sma(Dmed~Long*Diet, 
                         log = "xy",
                         data = Data_Animal_de_pd,
                         slope.test = 1/3,
                         robust = T
)
summary(Dmed_SMAlong_diet)
plot(Dmed_SMAlong_diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Dmed_SMAlong_diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#
#### angeqv and LONG ####
angeqv.SMAlong <- sma(angeqv~Long, 
                      log = "xy",
                      data = Data_Animal_de_pd,
                      slope.test = 0.1, 
                      #robust = T
)
summary(angeqv.SMAlong)
plot(angeqv.SMAlong)
# Stance
angeqv.SMAlong.stance <- sma(angeqv~Long*Stance, 
                             log = "xy",
                             data = Data_Animal_de_pd,
                             slope.test = 1/3, 
                             robust = T
)
summary(angeqv.SMAlong.stance)
plot(angeqv.SMAlong.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = angeqv.SMAlong.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
angeqv.SMAlong.Morphotype <- sma(angeqv~Long*Morphotype, 
                                 log = "xy",
                                 data = Data_Animal_de_pd, 
                                 slope.test = 1/3, 
                                 robust = T
)
summary(angeqv.SMAlong.Morphotype)
plot(angeqv.SMAlong.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = angeqv.SMAlong.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
angeqv.SMAlong.order <- sma(angeqv~Long*order, 
                            log = "xy",
                            data = Data_Animal_de_pd, 
                            slope.test = 1/3,
                            robust = T
)
summary(angeqv.SMAlong.order)
plot(angeqv_SMAlong_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = angeqv.SMAlong.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
angeqv.SMAlong.diet <- sma(angeqv~Long*Diet, 
                           log = "xy",
                           data = Data_Animal_de_pd,
                           slope.test = 1/3,
                           robust = T
)
summary(angeqv.SMAlong.diet)
plot(angeqv.SMAlong.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = angeqv.SMAlong.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
#### Pry_Area and LONG ####
Pry_Area.SMAlong <- sma(Pry_Area~Long, 
                        log = "xy",
                        data = Data_Animal_de_pd,
                        slope.test = 1/3, 
                        robust = T
)
summary(Pry_Area.SMAlong)
plot(Pry_Area.SMAlong)
# Stance
Pry_Area.SMAlong.stance <- sma(Pry_Area~Long*Stance, 
                               log = "xy",
                               data = Data_Animal_de_pd,
                               slope.test = 1/3, 
                               robust = T
)
summary(Pry_Area.SMAlong.stance)
plot(Pry_Area.SMAlong.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Pry_Area.SMAlong.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Pry_Area.SMAlong.Morphotype <- sma(Pry_Area~Long*Morphotype, 
                                   log = "xy",
                                   data = Data_Animal_de_pd, 
                                   slope.test = 1/3, 
                                   robust = T
)
summary(Pry_Area.SMAlong.Morphotype)
plot(Pry_Area.SMAlong.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Pry_Area.SMAlong.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Pry_Area.SMAlong.order <- sma(Pry_Area~Long*order, 
                              log = "xy",
                              data = Data_Animal_de_pd, 
                              slope.test = 1/3,
                              robust = T
)
summary(Pry_Area.SMAlong.order)
plot(Pry_Area_SMAlong_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Pry_Area.SMAlong.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Pry_Area.SMAlong.diet <- sma(Pry_Area~Long*Diet, 
                             log = "xy",
                             data = Data_Animal_de_pd,
                             slope.test = 1/3,
                             robust = T
)
summary(Pry_Area.SMAlong.diet)
plot(Pry_Area.SMAlong.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Pry_Area.SMAlong.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### Sf and LONG ####
Sf.SMAlong <- sma(Sf~Long, 
                  log = "xy",
                  data = Data_Animal_de_pd,
                  slope.test = 1/3, 
                  robust = T
)
summary(Sf.SMAlong)
plot(Sf.SMAlong)
# Stance
Sf.SMAlong.stance <- sma(Sf~Long*Stance, 
                         log = "xy",
                         data = Data_Animal_de_pd,
                         slope.test = 1/3, 
                         robust = T
)
summary(Sf.SMAlong.stance)
plot(Sf.SMAlong.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Sf.SMAlong.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Sf.SMAlong.Morphotype <- sma(Sf~Long*Morphotype, 
                             log = "xy",
                             data = Data_Animal_de_pd, 
                             slope.test = 1/3, 
                             robust = T
)
summary(Sf.SMAlong.Morphotype)
plot(Sf.SMAlong.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Sf.SMAlong.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Sf.SMAlong.order <- sma(Sf~Long*order, 
                        log = "xy",
                        data = Data_Animal_de_pd, 
                        slope.test = 1/3,
                        robust = T
)
summary(Sf.SMAlong.order)
plot(Sf_SMAlong_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Sf.SMAlong.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Sf.SMAlong.diet <- sma(Sf~Long*Diet, 
                       log = "xy",
                       data = Data_Animal_de_pd,
                       slope.test = 1/3,
                       robust = T
)
summary(Sf.SMAlong.diet)
plot(Sf.SMAlong.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Sf.SMAlong.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### Mscl_mass_kg and LONG ####
Mscl_mass_kg.SMAlong <- sma(Mscl_mass_kg~Long, 
                            log = "xy",
                            data = Data_Animal_de_pd,
                            slope.test = 1/3, 
                            robust = T
)
summary(Mscl_mass_kg.SMAlong)
plot(Mscl_mass_kg.SMAlong)
# Stance
Mscl_mass_kg.SMAlong.stance <- sma(Mscl_mass_kg~Long*Stance, 
                                   log = "xy",
                                   data = Data_Animal_de_pd,
                                   slope.test = 1/3, 
                                   robust = T
)
summary(Mscl_mass_kg.SMAlong.stance)
plot(Mscl_mass_kg.SMAlong.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Mscl_mass_kg.SMAlong.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Mscl_mass_kg.SMAlong.Morphotype <- sma(Mscl_mass_kg~Long*Morphotype, 
                                       log = "xy",
                                       data = Data_Animal_de_pd, 
                                       slope.test = 1/3, 
                                       robust = T
)
summary(Mscl_mass_kg.SMAlong.Morphotype)
plot(Mscl_mass_kg.SMAlong.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Mscl_mass_kg.SMAlong.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Mscl_mass_kg.SMAlong.order <- sma(Mscl_mass_kg~Long*order, 
                                  log = "xy",
                                  data = Data_Animal_de_pd, 
                                  slope.test = 1/3,
                                  robust = T
)
summary(Mscl_mass_kg.SMAlong.order)
plot(Mscl_mass_kg_SMAlong_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Mscl_mass_kg.SMAlong.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Mscl_mass_kg.SMAlong.diet <- sma(Mscl_mass_kg~Long*Diet, 
                                 log = "xy",
                                 data = Data_Animal_de_pd,
                                 slope.test = 1/3,
                                 robust = T
)
summary(Mscl_mass_kg.SMAlong.diet)
plot(Mscl_mass_kg.SMAlong.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Mscl_mass_kg.SMAlong.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### Mscl_length_m and LONG ####
Mscl_length_m.SMAlong <- sma(Mscl_length_m~Long, 
                             log = "xy",
                             data = Data_Animal_de_pd,
                             slope.test = 1/3, 
                             robust = T
)
summary(Mscl_length_m.SMAlong)
plot(Mscl_length_m.SMAlong)
# Stance
Mscl_length_m.SMAlong.stance <- sma(Mscl_length_m~Long*Stance, 
                                    log = "xy",
                                    data = Data_Animal_de_pd,
                                    slope.test = 1/3, 
                                    robust = T
)
summary(Mscl_length_m.SMAlong.stance)
plot(Mscl_length_m.SMAlong.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Mscl_length_m.SMAlong.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Mscl_length_m.SMAlong.Morphotype <- sma(Mscl_length_m~Long*Morphotype, 
                                        log = "xy",
                                        data = Data_Animal_de_pd, 
                                        slope.test = 1/3, 
                                        robust = T
)
summary(Mscl_length_m.SMAlong.Morphotype)
plot(Mscl_length_m.SMAlong.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Mscl_length_m.SMAlong.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Mscl_length_m.SMAlong.order <- sma(Mscl_length_m~Long*order, 
                                   log = "xy",
                                   data = Data_Animal_de_pd, 
                                   slope.test = 1/3,
                                   robust = T
)
summary(Mscl_length_m.SMAlong.order)
plot(Mscl_length_m_SMAlong_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Mscl_length_m.SMAlong.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Mscl_length_m.SMAlong.diet <- sma(Mscl_length_m~Long*Diet, 
                                  log = "xy",
                                  data = Data_Animal_de_pd,
                                  slope.test = 1/3,
                                  robust = T
)
summary(Mscl_length_m.SMAlong.diet)
plot(Mscl_length_m.SMAlong.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Mscl_length_m.SMAlong.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### MoI and LONG ####
MoI.SMAlong <- sma(MoI~Long, 
                   log = "xy",
                   data = Data_Animal_de_pd,
                   slope.test = 1/3, 
                   robust = T
)
summary(MoI.SMAlong)
plot(MoI.SMAlong)
# Stance
MoI.SMAlong.stance <- sma(MoI~Long*Stance, 
                          log = "xy",
                          data = Data_Animal_de_pd,
                          slope.test = 1/3, 
                          robust = T
)
summary(MoI.SMAlong.stance)
plot(MoI.SMAlong.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI.SMAlong.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
MoI.SMAlong.Morphotype <- sma(MoI~Long*Morphotype, 
                              log = "xy",
                              data = Data_Animal_de_pd, 
                              slope.test = 1/3, 
                              robust = T
)
summary(MoI.SMAlong.Morphotype)
plot(MoI.SMAlong.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI.SMAlong.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
MoI.SMAlong.order <- sma(MoI~Long*order, 
                         log = "xy",
                         data = Data_Animal_de_pd, 
                         slope.test = 1/3,
                         robust = T
)
summary(MoI.SMAlong.order)
plot(MoI_SMAlong_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI.SMAlong.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
MoI.SMAlong.diet <- sma(MoI~Long*Diet, 
                        log = "xy",
                        data = Data_Animal_de_pd,
                        slope.test = 1/3,
                        robust = T
)
summary(MoI.SMAlong.diet)
plot(MoI.SMAlong.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI.SMAlong.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### Limb_mass and LONG ####
Limb_mass.SMAlong <- sma(Limb_mass~Long, 
                         log = "xy",
                         data = Data_Animal_de_pd,
                         slope.test = 1/3, 
                         robust = T
)
summary(Limb_mass.SMAlong)
plot(Limb_mass.SMAlong)
# Stance
Limb_mass.SMAlong.stance <- sma(Limb_mass~Long*Stance, 
                                log = "xy",
                                data = Data_Animal_de_pd,
                                slope.test = 1/3, 
                                robust = T
)
summary(Limb_mass.SMAlong.stance)
plot(Limb_mass.SMAlong.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Limb_mass.SMAlong.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Limb_mass.SMAlong.Morphotype <- sma(Limb_mass~Long*Morphotype, 
                                    log = "xy",
                                    data = Data_Animal_de_pd, 
                                    slope.test = 1/3, 
                                    robust = T
)
summary(Limb_mass.SMAlong.Morphotype)
plot(Limb_mass.SMAlong.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Limb_mass.SMAlong.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Limb_mass.SMAlong.order <- sma(Limb_mass~Long*order, 
                               log = "xy",
                               data = Data_Animal_de_pd, 
                               slope.test = 1/3,
                               robust = T
)
summary(Limb_mass.SMAlong.order)
plot(Limb_mass_SMAlong_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Limb_mass.SMAlong.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Limb_mass.SMAlong.diet <- sma(Limb_mass~Long*Diet, 
                              log = "xy",
                              data = Data_Animal_de_pd,
                              slope.test = 1/3,
                              robust = T
)
summary(Limb_mass.SMAlong.diet)
plot(Limb_mass.SMAlong.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Limb_mass.SMAlong.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### CoM and LONG ####
CoM.SMAlong <- sma(CoM~Long, 
                   log = "xy",
                   data = Data_Animal_de_pd,
                   slope.test = 1/3, 
                   robust = T
)
summary(CoM.SMAlong)
plot(CoM.SMAlong)
# Stance
CoM.SMAlong.stance <- sma(CoM~Long*Stance, 
                          log = "xy",
                          data = Data_Animal_de_pd,
                          slope.test = 1/3, 
                          robust = T
)
summary(CoM.SMAlong.stance)
plot(CoM.SMAlong.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = CoM.SMAlong.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
CoM.SMAlong.Morphotype <- sma(CoM~Long*Morphotype, 
                              log = "xy",
                              data = Data_Animal_de_pd, 
                              slope.test = 1/3, 
                              robust = T
)
summary(CoM.SMAlong.Morphotype)
plot(CoM.SMAlong.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = CoM.SMAlong.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
CoM.SMAlong.order <- sma(CoM~Long*order, 
                         log = "xy",
                         data = Data_Animal_de_pd, 
                         slope.test = 1/3,
                         robust = T
)
summary(CoM.SMAlong.order)
plot(CoM_SMAlong_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = CoM.SMAlong.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
CoM.SMAlong.diet <- sma(CoM~Long*Diet, 
                        log = "xy",
                        data = Data_Animal_de_pd,
                        slope.test = 1/3,
                        robust = T
)
summary(CoM.SMAlong.diet)
plot(CoM.SMAlong.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = CoM.SMAlong.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### MoA_m and LONG ####
MoA_m.SMAlong <- sma(MoA_m~Long, 
                     log = "xy",
                     data = Data_Animal_de_pd,
                     slope.test = 1/3, 
                     robust = T
)
summary(MoA_m.SMAlong)
plot(MoA_m.SMAlong)
# Stance
MoA_m.SMAlong.stance <- sma(MoA_m~Long*Stance, 
                            log = "xy",
                            data = Data_Animal_de_pd,
                            slope.test = 1/3, 
                            robust = T
)
summary(MoA_m.SMAlong.stance)
plot(MoA_m.SMAlong.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoA_m.SMAlong.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
MoA_m.SMAlong.Morphotype <- sma(MoA_m~Long*Morphotype, 
                                log = "xy",
                                data = Data_Animal_de_pd, 
                                slope.test = 1/3, 
                                robust = T
)
summary(MoA_m.SMAlong.Morphotype)
plot(MoA_m.SMAlong.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoA_m.SMAlong.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
MoA_m.SMAlong.order <- sma(MoA_m~Long*order, 
                           log = "xy",
                           data = Data_Animal_de_pd, 
                           slope.test = 1/3,
                           robust = T
)
summary(MoA_m.SMAlong.order)
plot(MoA_m_SMAlong_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoA_m.SMAlong.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
MoA_m.SMAlong.diet <- sma(MoA_m~Long*Diet, 
                          log = "xy",
                          data = Data_Animal_de_pd,
                          slope.test = 1/3,
                          robust = T
)
summary(MoA_m.SMAlong.diet)
plot(MoA_m.SMAlong.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoA_m.SMAlong.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### P_eq and LONG ####
P_eq.SMAlong <- sma(P_eq~Long, 
                    log = "xy",
                    data = Data_Animal_de_pd,
                    slope.test = 1/3, 
                    robust = T
)
summary(P_eq.SMAlong)
plot(P_eq.SMAlong)
# Stance
P_eq.SMAlong.stance <- sma(P_eq~Long*Stance, 
                           log = "xy",
                           data = Data_Animal_de_pd,
                           slope.test = 1/3, 
                           robust = T
)
summary(P_eq.SMAlong.stance)
plot(P_eq.SMAlong.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = P_eq.SMAlong.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
P_eq.SMAlong.Morphotype <- sma(P_eq~Long*Morphotype, 
                               log = "xy",
                               data = Data_Animal_de_pd, 
                               slope.test = 1/3, 
                               robust = T
)
summary(P_eq.SMAlong.Morphotype)
plot(P_eq.SMAlong.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = P_eq.SMAlong.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
P_eq.SMAlong.order <- sma(P_eq~Long*order, 
                          log = "xy",
                          data = Data_Animal_de_pd, 
                          slope.test = 1/3,
                          robust = T
)
summary(P_eq.SMAlong.order)
plot(P_eq_SMAlong_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = P_eq.SMAlong.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
P_eq.SMAlong.diet <- sma(P_eq~Long*Diet, 
                         log = "xy",
                         data = Data_Animal_de_pd,
                         slope.test = 1/3,
                         robust = T
)
summary(P_eq.SMAlong.diet)
plot(P_eq.SMAlong.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = P_eq.SMAlong.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### MoI_G and LONG ####
MoI_G.SMAlong <- sma(MoI_G~Long, 
                     log = "xy",
                     data = Data_Animal_de_pd,
                     slope.test = 1/3, 
                     robust = T
)
summary(MoI_G.SMAlong)
plot(MoI_G.SMAlong)
# Stance
MoI_G.SMAlong.stance <- sma(MoI_G~Long*Stance, 
                            log = "xy",
                            data = Data_Animal_de_pd,
                            slope.test = 1/3, 
                            robust = T
)
summary(MoI_G.SMAlong.stance)
plot(MoI_G.SMAlong.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_G.SMAlong.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
MoI_G.SMAlong.Morphotype <- sma(MoI_G~Long*Morphotype, 
                                log = "xy",
                                data = Data_Animal_de_pd, 
                                slope.test = 1/3, 
                                robust = T
)
summary(MoI_G.SMAlong.Morphotype)
plot(MoI_G.SMAlong.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_G.SMAlong.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
MoI_G.SMAlong.order <- sma(MoI_G~Long*order, 
                           log = "xy",
                           data = Data_Animal_de_pd, 
                           slope.test = 1/3,
                           robust = T
)
summary(MoI_G.SMAlong.order)
plot(MoI_G_SMAlong_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_G.SMAlong.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
MoI_G.SMAlong.diet <- sma(MoI_G~Long*Diet, 
                          log = "xy",
                          data = Data_Animal_de_pd,
                          slope.test = 1/3,
                          robust = T
)
summary(MoI_G.SMAlong.diet)
plot(MoI_G.SMAlong.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_G.SMAlong.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### MoI_H and LONG ####
MoI_H.SMAlong <- sma(MoI_H~Long, 
                     log = "xy",
                     data = Data_Animal_de_pd,
                     slope.test = 1/3, 
                     robust = T
)
summary(MoI_H.SMAlong)
plot(MoI_H.SMAlong)
# Stance
MoI_H.SMAlong.stance <- sma(MoI_H~Long*Stance, 
                            log = "xy",
                            data = Data_Animal_de_pd,
                            slope.test = 1/3, 
                            robust = T
)
summary(MoI_H.SMAlong.stance)
plot(MoI_H.SMAlong.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_H.SMAlong.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
MoI_H.SMAlong.Morphotype <- sma(MoI_H~Long*Morphotype, 
                                log = "xy",
                                data = Data_Animal_de_pd, 
                                slope.test = 1/3, 
                                robust = T
)
summary(MoI_H.SMAlong.Morphotype)
plot(MoI_H.SMAlong.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_H.SMAlong.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
MoI_H.SMAlong.order <- sma(MoI_H~Long*order, 
                           log = "xy",
                           data = Data_Animal_de_pd, 
                           slope.test = 1/3,
                           robust = T
)
summary(MoI_H.SMAlong.order)
plot(MoI_H_SMAlong_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_H.SMAlong.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
MoI_H.SMAlong.diet <- sma(MoI_H~Long*Diet, 
                          log = "xy",
                          data = Data_Animal_de_pd,
                          slope.test = 1/3,
                          robust = T
)
summary(MoI_H.SMAlong.diet)
plot(MoI_H.SMAlong.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_H.SMAlong.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
#### Alpha and LONG ####
Alpha.SMAlong <- sma(Alpha~Long, 
                     log = "xy",
                     data = Data_Animal_de_pd,
                     slope.test = 1/3, 
                     robust = T
)
summary(Alpha.SMAlong)
plot(Alpha.SMAlong)
# Stance
Alpha.SMAlong.stance <- sma(Alpha~Long*Stance, 
                            log = "xy",
                            data = Data_Animal_de_pd,
                            slope.test = 1/3, 
                            robust = T
)
summary(Alpha.SMAlong.stance)
plot(Alpha.SMAlong.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Alpha.SMAlong.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Alpha.SMAlong.Morphotype <- sma(Alpha~Long*Morphotype, 
                                log = "xy",
                                data = Data_Animal_de_pd, 
                                slope.test = 1/3, 
                                robust = T
)
summary(Alpha.SMAlong.Morphotype)
plot(Alpha.SMAlong.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Alpha.SMAlong.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Alpha.SMAlong.order <- sma(Alpha~Long*order, 
                           log = "xy",
                           data = Data_Animal_de_pd, 
                           slope.test = 1/3,
                           robust = T
)
summary(Alpha.SMAlong.order)
plot(Alpha_SMAlong_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Alpha.SMAlong.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Alpha.SMAlong.diet <- sma(Alpha~Long*Diet, 
                          log = "xy",
                          data = Data_Animal_de_pd,
                          slope.test = 1/3,
                          robust = T
)
summary(Alpha.SMAlong.diet)
plot(Alpha.SMAlong.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Alpha.SMAlong.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
#### Omg and LONG ####
Omg.SMAlong <- sma(Omg~Long, 
                   log = "xy",
                   data = Data_Animal_de_pd,
                   slope.test = 1/3, 
                   robust = T
)
summary(Omg.SMAlong)
plot(Omg.SMAlong)
# Stance
Omg.SMAlong.stance <- sma(Omg~Long*Stance, 
                          log = "xy",
                          data = Data_Animal_de_pd,
                          slope.test = 1/3, 
                          robust = T
)
summary(Omg.SMAlong.stance)
plot(Omg.SMAlong.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Omg.SMAlong.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Omg.SMAlong.Morphotype <- sma(Omg~Long*Morphotype, 
                              log = "xy",
                              data = Data_Animal_de_pd, 
                              slope.test = 1/3, 
                              robust = T
)
summary(Omg.SMAlong.Morphotype)
plot(Omg.SMAlong.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Omg.SMAlong.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Omg.SMAlong.order <- sma(Omg~Long*order, 
                         log = "xy",
                         data = Data_Animal_de_pd, 
                         slope.test = 1/3,
                         robust = T
)
summary(Omg.SMAlong.order)
plot(Omg_SMAlong_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Omg.SMAlong.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Omg.SMAlong.diet <- sma(Omg~Long*Diet, 
                        log = "xy",
                        data = Data_Animal_de_pd,
                        slope.test = 1/3,
                        robust = T
)
summary(Omg.SMAlong.diet)
plot(Omg.SMAlong.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Omg.SMAlong.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
#### Vsup and LONG ####
Vsup.SMAlong <- sma(Vsup~Long, 
                    log = "xy",
                    data = Data_Animal_de_pd,
                    slope.test = 1/3, 
                    robust = T
)
summary(Vsup.SMAlong)
plot(Vsup.SMAlong)
# Stance
Vsup.SMAlong.stance <- sma(Vsup~Long*Stance, 
                           log = "xy",
                           data = Data_Animal_de_pd,
                           slope.test = 1/3, 
                           robust = T
)
summary(Vsup.SMAlong.stance)
plot(Vsup.SMAlong.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Vsup.SMAlong.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Vsup.SMAlong.Morphotype <- sma(Vsup~Long*Morphotype, 
                               log = "xy",
                               data = Data_Animal_de_pd, 
                               slope.test = 1/3, 
                               robust = T
)
summary(Vsup.SMAlong.Morphotype)
plot(Vsup.SMAlong.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Vsup.SMAlong.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Vsup.SMAlong.order <- sma(Vsup~Long*order, 
                          log = "xy",
                          data = Data_Animal_de_pd, 
                          slope.test = 1/3,
                          robust = T
)
summary(Vsup.SMAlong.order)
plot(Vsup_SMAlong_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Vsup.SMAlong.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Vsup.SMAlong.diet <- sma(Vsup~Long*Diet, 
                         log = "xy",
                         data = Data_Animal_de_pd,
                         slope.test = 1/3,
                         robust = T
)
summary(Vsup.SMAlong.diet)
plot(Vsup.SMAlong.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Vsup.SMAlong.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### MoI_HC and LONG ####
MoI_HC.SMAlong <- sma(MoI_HC~Long, 
                      log = "xy",
                      data = Data_Animal_de_pd,
                      slope.test = 1/3, 
                      robust = T
)
summary(MoI_HC.SMAlong)
plot(MoI_HC.SMAlong)
# Stance
MoI_HC.SMAlong.stance <- sma(MoI_HC~Long*Stance, 
                             log = "xy",
                             data = Data_Animal_de_pd,
                             slope.test = 1/3, 
                             robust = T
)
summary(MoI_HC.SMAlong.stance)
plot(MoI_HC.SMAlong.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_HC.SMAlong.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
MoI_HC.SMAlong.Morphotype <- sma(MoI_HC~Long*Morphotype, 
                                 log = "xy",
                                 data = Data_Animal_de_pd, 
                                 slope.test = 1/3, 
                                 robust = T
)
summary(MoI_HC.SMAlong.Morphotype)
plot(MoI_HC.SMAlong.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_HC.SMAlong.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
MoI_HC.SMAlong.order <- sma(MoI_HC~Long*order, 
                            log = "xy",
                            data = Data_Animal_de_pd, 
                            slope.test = 1/3,
                            robust = T
)
summary(MoI_HC.SMAlong.order)
plot(MoI_HC_SMAlong_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_HC.SMAlong.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
MoI_HC.SMAlong.diet <- sma(MoI_HC~Long*Diet, 
                           log = "xy",
                           data = Data_Animal_de_pd,
                           slope.test = 1/3,
                           robust = T
)
summary(MoI_HC.SMAlong.diet)
plot(MoI_HC.SMAlong.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_HC.SMAlong.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### Alpha_C and LONG ####
Alpha_C.SMAlong <- sma(Alpha_C~Long, 
                       log = "xy",
                       data = Data_Animal_de_pd,
                       slope.test = 1/3, 
                       robust = T
)
summary(Alpha_C.SMAlong)
plot(Alpha_C.SMAlong)
# Stance
Alpha_C.SMAlong.stance <- sma(Alpha_C~Long*Stance, 
                              log = "xy",
                              data = Data_Animal_de_pd,
                              slope.test = 1/3, 
                              robust = T
)
summary(Alpha_C.SMAlong.stance)
plot(Alpha_C.SMAlong.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Alpha_C.SMAlong.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Alpha_C.SMAlong.Morphotype <- sma(Alpha_C~Long*Morphotype, 
                                  log = "xy",
                                  data = Data_Animal_de_pd, 
                                  slope.test = 1/3, 
                                  robust = T
)
summary(Alpha_C.SMAlong.Morphotype)
plot(Alpha_C.SMAlong.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Alpha_C.SMAlong.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Alpha_C.SMAlong.order <- sma(Alpha_C~Long*order, 
                             log = "xy",
                             data = Data_Animal_de_pd, 
                             slope.test = 1/3,
                             robust = T
)
summary(Alpha_C.SMAlong.order)
plot(Alpha_C_SMAlong_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Alpha_C.SMAlong.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Alpha_C.SMAlong.diet <- sma(Alpha_C~Long*Diet, 
                            log = "xy",
                            data = Data_Animal_de_pd,
                            slope.test = 1/3,
                            robust = T
)
summary(Alpha_C.SMAlong.diet)
plot(Alpha_C.SMAlong.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Alpha_C.SMAlong.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### OmgC and LONG ####
OmgC.SMAlong <- sma(OmgC~Long, 
                    log = "xy",
                    data = Data_Animal_de_pd,
                    slope.test = 1/3, 
                    robust = T
)
summary(OmgC.SMAlong)
plot(OmgC.SMAlong)
# Stance
OmgC.SMAlong.stance <- sma(OmgC~Long*Stance, 
                           log = "xy",
                           data = Data_Animal_de_pd,
                           slope.test = 1/3, 
                           robust = T
)
summary(OmgC.SMAlong.stance)
plot(OmgC.SMAlong.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = OmgC.SMAlong.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
OmgC.SMAlong.Morphotype <- sma(OmgC~Long*Morphotype, 
                               log = "xy",
                               data = Data_Animal_de_pd, 
                               slope.test = 1/3, 
                               robust = T
)
summary(OmgC.SMAlong.Morphotype)
plot(OmgC.SMAlong.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = OmgC.SMAlong.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
OmgC.SMAlong.order <- sma(OmgC~Long*order, 
                          log = "xy",
                          data = Data_Animal_de_pd, 
                          slope.test = 1/3,
                          robust = T
)
summary(OmgC.SMAlong.order)
plot(OmgC_SMAlong_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = OmgC.SMAlong.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
OmgC.SMAlong.diet <- sma(OmgC~Long*Diet, 
                         log = "xy",
                         data = Data_Animal_de_pd,
                         slope.test = 1/3,
                         robust = T
)
summary(OmgC.SMAlong.diet)
plot(OmgC.SMAlong.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = OmgC.SMAlong.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### VsupC and LONG ####
VsupC.SMAlong <- sma(VsupC~Long, 
                     log = "xy",
                     data = Data_Animal_de_pd,
                     slope.test = 1/3, 
                     robust = T
)
summary(VsupC.SMAlong)
plot(VsupC.SMAlong)
# Stance
VsupC.SMAlong.stance <- sma(VsupC~Long*Stance, 
                            log = "xy",
                            data = Data_Animal_de_pd,
                            slope.test = 1/3, 
                            robust = T
)
summary(VsupC.SMAlong.stance)
plot(VsupC.SMAlong.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = VsupC.SMAlong.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
VsupC.SMAlong.Morphotype <- sma(VsupC~Long*Morphotype, 
                                log = "xy",
                                data = Data_Animal_de_pd, 
                                slope.test = 1/3, 
                                robust = T
)
summary(VsupC.SMAlong.Morphotype)
plot(VsupC.SMAlong.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = VsupC.SMAlong.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
VsupC.SMAlong.order <- sma(VsupC~Long*order, 
                           log = "xy",
                           data = Data_Animal_de_pd, 
                           slope.test = 1/3,
                           robust = T
)
summary(VsupC.SMAlong.order)
plot(VsupC_SMAlong_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = VsupC.SMAlong.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
VsupC.SMAlong.diet <- sma(VsupC~Long*Diet, 
                          log = "xy",
                          data = Data_Animal_de_pd,
                          slope.test = 1/3,
                          robust = T
)
summary(VsupC.SMAlong.diet)
plot(VsupC.SMAlong.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = VsupC.SMAlong.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
#### Maxspeed_ms and LONG ####
Maxspeed_ms.SMAlong <- sma(Maxspeed_ms~Long, 
                           log = "xy",
                           data = Data_Animal_de_pd,
                           slope.test = 1/3, 
                           robust = T
)
summary(Maxspeed_ms.SMAlong)
plot(Maxspeed_ms.SMAlong)
# Stance
Maxspeed_ms.SMAlong.stance <- sma(Maxspeed_ms~Long*Stance, 
                                  log = "xy",
                                  data = Data_Animal_de_pd,
                                  slope.test = 1/3, 
                                  robust = T
)
summary(Maxspeed_ms.SMAlong.stance)
plot(Maxspeed_ms.SMAlong.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Maxspeed_ms.SMAlong.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Maxspeed_ms.SMAlong.Morphotype <- sma(Maxspeed_ms~Long*Morphotype, 
                                      log = "xy",
                                      data = Data_Animal_de_pd, 
                                      slope.test = 1/3, 
                                      robust = T
)
summary(Maxspeed_ms.SMAlong.Morphotype)
plot(Maxspeed_ms.SMAlong.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Maxspeed_ms.SMAlong.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Maxspeed_ms.SMAlong.order <- sma(Maxspeed_ms~Long*order, 
                                 log = "xy",
                                 data = Data_Animal_de_pd, 
                                 slope.test = 1/3,
                                 robust = T
)
summary(Maxspeed_ms.SMAlong.order)
plot(Maxspeed_ms.SMAlong.order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Maxspeed_ms.SMAlong.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Maxspeed_ms.SMAlong.diet <- sma(Maxspeed_ms~Long*Diet, 
                                log = "xy",
                                data = Data_Animal_de_pd,
                                slope.test = 1/3,
                                robust = T
)
summary(Maxspeed_ms.SMAlong.diet)
plot(Maxspeed_ms.SMAlong.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Maxspeed_ms.SMAlong.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)


################################################################################
# ALL and LA #################################
#### Rmax and LA ####
Rmax.SMAla <- sma(Rmax~La, 
                  log = "xy",
                  data = Data_Animal_de_pd, 
                  slope.test = 1/3, 
                  robust = T
)
summary(Rmax.SMAla)
plot(Rmax.SMAla)
# Stance
Rmax.SMAla.stance <- sma(Rmax~La*Stance, 
                         log = "xy",
                         data = Data_Animal_de_pd,
                         slope.test = 1/3, 
                         robust = T
)
summary(Rmax.SMAla.stance)
plot(Rmax.SMAla.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Rmax.SMAla.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Rmax.SMAla.Morphotype <- sma(Rmax~La*Morphotype, 
                             log = "xy",
                             data = Data_Animal_de_pd, 
                             slope.test = 1/3, 
                             #robust = T
)
summary(Rmax.SMAla.Morphotype)
plot(Rmax.SMAla.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Rmax.SMAla.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Rmax.SMAla.order <- sma(Rmax~La*order,
                        log = "xy",
                        data = Data_Animal_de_pd,
                        slope.test = 1/3,
                        robust = T
)
summary(Rmax_SMAla_order)
plot(Rmax_SMAla_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Rmax_SMAla_order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Rmax_SMAla_diet <- sma(Rmax~La*Diet, log = "xy" ,data = Data_Animal_de_pd, slope.test = 1/3)
summary(Rmax_SMAla_diet)
plot(Rmax_SMAla_diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Rmax_SMAla_diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#
#### Rmin and LA ####
Rmin.SMAla <- sma(Rmin~La, 
                  log = "xy",
                  data = Data_Animal_de_pd,
                  slope.test = 1/3, 
                  robust = T
)
summary(Rmin.SMAla)
plot(Rmin.SMAla)
# Stance
Rmin.SMAla.stance <- sma(Rmin~La*Stance, 
                         log = "xy",
                         data = Data_Animal_de_pd,
                         slope.test = 1/3, 
                         robust = T
)
summary(Rmin.SMAla.stance)
plot(Rmin.SMAla.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Rmin.SMAla.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Rmin.SMAla.Morphotype <- sma(Rmin~La*Morphotype, 
                             log = "xy",
                             data = Data_Animal_de_pd, 
                             slope.test = 1/3, 
                             robust = T
)
summary(Rmin.SMAla.Morphotype)
plot(Rmin.SMAla.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Rmin.SMAla.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Rmin.SMAla.order <- sma(Rmin~La*order, 
                        log = "xy",
                        data = Data_Animal_de_pd, 
                        slope.test = 1/3,
                        robust = T
)
summary(Rmin.SMAla.order)
plot(Rmin_SMAla_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Rmin.SMAla.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Rmin.SMAla.diet <- sma(Rmin~La*Diet, 
                       log = "xy",
                       data = Data_Animal_de_pd,
                       slope.test = 1/3,
                       robust = T
)
summary(Rmin.SMAla.diet)
plot(Rmin.SMAla.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Rmin.SMAla.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#
#### DeltaR and LA ####
DeltaR.SMAla <- sma(DeltaR~La, 
                    log = "xy",
                    data = Data_Animal_de_pd,
                    slope.test = 1/3, 
                    robust = T
)
summary(DeltaR.SMAla)
plot(DeltaR.SMAla)
# Stance
DeltaR.SMAla.stance <- sma(DeltaR~La*Stance, 
                           log = "xy",
                           data = Data_Animal_de_pd,
                           slope.test = 1/3, 
                           robust = T
)
summary(DeltaR.SMAla.stance)
plot(DeltaR.SMAla.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = DeltaR.SMAla.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
DeltaR.SMAla.Morphotype <- sma(DeltaR~La*Morphotype, 
                               log = "xy",
                               data = Data_Animal_de_pd, 
                               slope.test = 1/3, 
                               robust = T
)
summary(DeltaR.SMAla.Morphotype)
plot(DeltaR.SMAla.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = DeltaR.SMAla.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
DeltaR.SMAla.order <- sma(DeltaR~La*order, 
                          log = "xy",
                          data = Data_Animal_de_pd, 
                          slope.test = 1/3,
                          robust = T
)
summary(DeltaR.SMAla.order)
plot(DeltaR_SMAla_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = DeltaR.SMAla.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
DeltaR.SMAla.diet <- sma(DeltaR~La*Diet, 
                         log = "xy",
                         data = Data_Animal_de_pd,
                         slope.test = 1/3,
                         robust = T
)
summary(DeltaR.SMAla.diet)
plot(DeltaR.SMAla.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = DeltaR.SMAla.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### Dmed and LA ####
Dmed.SMAla <- sma(Dmed~La, 
                  log = "xy",
                  data = Data_Animal_de_pd, 
                  slope.test = 1/3, 
                  robust = T
)
summary(Dmed.SMAla)
coef(Dmed.SMAla)
plot(Dmed.SMAla)
# Stance
Dmed.SMAla.stance <- sma(Dmed~La*Stance, 
                         log = "xy",
                         data = Data_Animal_de_pd,
                         slope.test = 1/3, 
                         robust = T
)
summary(Dmed.SMAla.stance)
coef(Dmed.SMAla.stance)
plot(Dmed.SMAla.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Dmed.SMAla.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Dmed.SMAla.Morphotype <- sma(Dmed~La*Morphotype, 
                             log = "xy",
                             data = Data_Animal_de_pd, 
                             slope.test = 1/3, 
                             robust = T
)
summary(Dmed.SMAla.Morphotype)
plot(Dmed.SMAla.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Dmed.SMAla.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Dmed.SMAla.order <- sma(Dmed~La*order, 
                        log = "xy",
                        data = Data_Animal_de_pd,
                        slope.test = 1/3,
                        robust = T
)
summary(Dmed.SMAla.order)
plot(Dmed.SMAla.order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Dmed.SMAla.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Dmed_SMAla_diet <- sma(Dmed~La*Diet, 
                       log = "xy",
                       data = Data_Animal_de_pd,
                       slope.test = 1/3,
                       robust = T
)
summary(Dmed_SMAla_diet)
plot(Dmed_SMAla_diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Dmed_SMAla_diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#
#### angeqv and LA ####
angeqv.SMAla <- sma(angeqv~La, 
                    log = "xy",
                    data = Data_Animal_de_pd,
                    slope.test = 0.1, 
                    #robust = T
)
summary(angeqv.SMAla)
plot(angeqv.SMAla)
# Stance
angeqv.SMAla.stance <- sma(angeqv~La*Stance, 
                           log = "xy",
                           data = Data_Animal_de_pd,
                           slope.test = 1/3, 
                           robust = T
)
summary(angeqv.SMAla.stance)
plot(angeqv.SMAla.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = angeqv.SMAla.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
angeqv.SMAla.Morphotype <- sma(angeqv~La*Morphotype, 
                               log = "xy",
                               data = Data_Animal_de_pd, 
                               slope.test = 1/3, 
                               robust = T
)
summary(angeqv.SMAla.Morphotype)
plot(angeqv.SMAla.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = angeqv.SMAla.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
angeqv.SMAla.order <- sma(angeqv~La*order, 
                          log = "xy",
                          data = Data_Animal_de_pd, 
                          slope.test = 1/3,
                          robust = T
)
summary(angeqv.SMAla.order)
plot(angeqv_SMAla_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = angeqv.SMAla.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
angeqv.SMAla.diet <- sma(angeqv~La*Diet, 
                         log = "xy",
                         data = Data_Animal_de_pd,
                         slope.test = 1/3,
                         robust = T
)
summary(angeqv.SMAla.diet)
plot(angeqv.SMAla.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = angeqv.SMAla.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
#### Pry_Area and LA ####
Pry_Area.SMAla <- sma(Pry_Area~La, 
                      log = "xy",
                      data = Data_Animal_de_pd,
                      slope.test = 1/3, 
                      robust = T
)
summary(Pry_Area.SMAla)
plot(Pry_Area.SMAla)
# Stance
Pry_Area.SMAla.stance <- sma(Pry_Area~La*Stance, 
                             log = "xy",
                             data = Data_Animal_de_pd,
                             slope.test = 1/3, 
                             robust = T
)
summary(Pry_Area.SMAla.stance)
plot(Pry_Area.SMAla.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Pry_Area.SMAla.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Pry_Area.SMAla.Morphotype <- sma(Pry_Area~La*Morphotype, 
                                 log = "xy",
                                 data = Data_Animal_de_pd, 
                                 slope.test = 1/3, 
                                 robust = T
)
summary(Pry_Area.SMAla.Morphotype)
plot(Pry_Area.SMAla.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Pry_Area.SMAla.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Pry_Area.SMAla.order <- sma(Pry_Area~La*order, 
                            log = "xy",
                            data = Data_Animal_de_pd, 
                            slope.test = 1/3,
                            robust = T
)
summary(Pry_Area.SMAla.order)
plot(Pry_Area_SMAla_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Pry_Area.SMAla.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Pry_Area.SMAla.diet <- sma(Pry_Area~La*Diet, 
                           log = "xy",
                           data = Data_Animal_de_pd,
                           slope.test = 1/3,
                           robust = T
)
summary(Pry_Area.SMAla.diet)
plot(Pry_Area.SMAla.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Pry_Area.SMAla.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### Sf and LA ####
Sf.SMAla <- sma(Sf~La, 
                log = "xy",
                data = Data_Animal_de_pd,
                slope.test = 1/3, 
                robust = T
)
summary(Sf.SMAla)
plot(Sf.SMAla)
# Stance
Sf.SMAla.stance <- sma(Sf~La*Stance, 
                       log = "xy",
                       data = Data_Animal_de_pd,
                       slope.test = 1/3, 
                       robust = T
)
summary(Sf.SMAla.stance)
plot(Sf.SMAla.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Sf.SMAla.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Sf.SMAla.Morphotype <- sma(Sf~La*Morphotype, 
                           log = "xy",
                           data = Data_Animal_de_pd, 
                           slope.test = 1/3, 
                           robust = T
)
summary(Sf.SMAla.Morphotype)
plot(Sf.SMAla.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Sf.SMAla.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Sf.SMAla.order <- sma(Sf~La*order, 
                      log = "xy",
                      data = Data_Animal_de_pd, 
                      slope.test = 1/3,
                      robust = T
)
summary(Sf.SMAla.order)
plot(Sf_SMAla_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Sf.SMAla.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Sf.SMAla.diet <- sma(Sf~La*Diet, 
                     log = "xy",
                     data = Data_Animal_de_pd,
                     slope.test = 1/3,
                     robust = T
)
summary(Sf.SMAla.diet)
plot(Sf.SMAla.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Sf.SMAla.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### Mscl_mass_kg and LA ####
Mscl_mass_kg.SMAla <- sma(Mscl_mass_kg~La, 
                          log = "xy",
                          data = Data_Animal_de_pd,
                          slope.test = 1/3, 
                          robust = T
)
summary(Mscl_mass_kg.SMAla)
plot(Mscl_mass_kg.SMAla)
# Stance
Mscl_mass_kg.SMAla.stance <- sma(Mscl_mass_kg~La*Stance, 
                                 log = "xy",
                                 data = Data_Animal_de_pd,
                                 slope.test = 1/3, 
                                 robust = T
)
summary(Mscl_mass_kg.SMAla.stance)
plot(Mscl_mass_kg.SMAla.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Mscl_mass_kg.SMAla.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Mscl_mass_kg.SMAla.Morphotype <- sma(Mscl_mass_kg~La*Morphotype, 
                                     log = "xy",
                                     data = Data_Animal_de_pd, 
                                     slope.test = 1/3, 
                                     robust = T
)
summary(Mscl_mass_kg.SMAla.Morphotype)
plot(Mscl_mass_kg.SMAla.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Mscl_mass_kg.SMAla.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Mscl_mass_kg.SMAla.order <- sma(Mscl_mass_kg~La*order, 
                                log = "xy",
                                data = Data_Animal_de_pd, 
                                slope.test = 1/3,
                                robust = T
)
summary(Mscl_mass_kg.SMAla.order)
plot(Mscl_mass_kg_SMAla_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Mscl_mass_kg.SMAla.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Mscl_mass_kg.SMAla.diet <- sma(Mscl_mass_kg~La*Diet, 
                               log = "xy",
                               data = Data_Animal_de_pd,
                               slope.test = 1/3,
                               robust = T
)
summary(Mscl_mass_kg.SMAla.diet)
plot(Mscl_mass_kg.SMAla.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Mscl_mass_kg.SMAla.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### Mscl_length_m and LA ####
Mscl_length_m.SMAla <- sma(Mscl_length_m~La, 
                           log = "xy",
                           data = Data_Animal_de_pd,
                           slope.test = 1/3, 
                           robust = T
)
summary(Mscl_length_m.SMAla)
plot(Mscl_length_m.SMAla)
# Stance
Mscl_length_m.SMAla.stance <- sma(Mscl_length_m~La*Stance, 
                                  log = "xy",
                                  data = Data_Animal_de_pd,
                                  slope.test = 1/3, 
                                  robust = T
)
summary(Mscl_length_m.SMAla.stance)
plot(Mscl_length_m.SMAla.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Mscl_length_m.SMAla.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Mscl_length_m.SMAla.Morphotype <- sma(Mscl_length_m~La*Morphotype, 
                                      log = "xy",
                                      data = Data_Animal_de_pd, 
                                      slope.test = 1/3, 
                                      robust = T
)
summary(Mscl_length_m.SMAla.Morphotype)
plot(Mscl_length_m.SMAla.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Mscl_length_m.SMAla.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Mscl_length_m.SMAla.order <- sma(Mscl_length_m~La*order, 
                                 log = "xy",
                                 data = Data_Animal_de_pd, 
                                 slope.test = 1/3,
                                 robust = T
)
summary(Mscl_length_m.SMAla.order)
plot(Mscl_length_m_SMAla_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Mscl_length_m.SMAla.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Mscl_length_m.SMAla.diet <- sma(Mscl_length_m~La*Diet, 
                                log = "xy",
                                data = Data_Animal_de_pd,
                                slope.test = 1/3,
                                robust = T
)
summary(Mscl_length_m.SMAla.diet)
plot(Mscl_length_m.SMAla.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Mscl_length_m.SMAla.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### MoI and LA ####
MoI.SMAla <- sma(MoI~La, 
                 log = "xy",
                 data = Data_Animal_de_pd,
                 slope.test = 1/3, 
                 robust = T
)
summary(MoI.SMAla)
plot(MoI.SMAla)
# Stance
MoI.SMAla.stance <- sma(MoI~La*Stance, 
                        log = "xy",
                        data = Data_Animal_de_pd,
                        slope.test = 1/3, 
                        robust = T
)
summary(MoI.SMAla.stance)
plot(MoI.SMAla.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI.SMAla.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
MoI.SMAla.Morphotype <- sma(MoI~La*Morphotype, 
                            log = "xy",
                            data = Data_Animal_de_pd, 
                            slope.test = 1/3, 
                            robust = T
)
summary(MoI.SMAla.Morphotype)
plot(MoI.SMAla.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI.SMAla.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
MoI.SMAla.order <- sma(MoI~La*order, 
                       log = "xy",
                       data = Data_Animal_de_pd, 
                       slope.test = 1/3,
                       robust = T
)
summary(MoI.SMAla.order)
plot(MoI_SMAla_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI.SMAla.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
MoI.SMAla.diet <- sma(MoI~La*Diet, 
                      log = "xy",
                      data = Data_Animal_de_pd,
                      slope.test = 1/3,
                      robust = T
)
summary(MoI.SMAla.diet)
plot(MoI.SMAla.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI.SMAla.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### Limb_mass and LA ####
Limb_mass.SMAla <- sma(Limb_mass~La, 
                       log = "xy",
                       data = Data_Animal_de_pd,
                       slope.test = 1/3, 
                       robust = T
)
summary(Limb_mass.SMAla)
plot(Limb_mass.SMAla)
# Stance
Limb_mass.SMAla.stance <- sma(Limb_mass~La*Stance, 
                              log = "xy",
                              data = Data_Animal_de_pd,
                              slope.test = 1/3, 
                              robust = T
)
summary(Limb_mass.SMAla.stance)
plot(Limb_mass.SMAla.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Limb_mass.SMAla.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Limb_mass.SMAla.Morphotype <- sma(Limb_mass~La*Morphotype, 
                                  log = "xy",
                                  data = Data_Animal_de_pd, 
                                  slope.test = 1/3, 
                                  robust = T
)
summary(Limb_mass.SMAla.Morphotype)
plot(Limb_mass.SMAla.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Limb_mass.SMAla.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Limb_mass.SMAla.order <- sma(Limb_mass~La*order, 
                             log = "xy",
                             data = Data_Animal_de_pd, 
                             slope.test = 1/3,
                             robust = T
)
summary(Limb_mass.SMAla.order)
plot(Limb_mass_SMAla_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Limb_mass.SMAla.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Limb_mass.SMAla.diet <- sma(Limb_mass~La*Diet, 
                            log = "xy",
                            data = Data_Animal_de_pd,
                            slope.test = 1/3,
                            robust = T
)
summary(Limb_mass.SMAla.diet)
plot(Limb_mass.SMAla.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Limb_mass.SMAla.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### CoM and LA ####
CoM.SMAla <- sma(CoM~La, 
                 log = "xy",
                 data = Data_Animal_de_pd,
                 slope.test = 1/3, 
                 robust = T
)
summary(CoM.SMAla)
plot(CoM.SMAla)
# Stance
CoM.SMAla.stance <- sma(CoM~La*Stance, 
                        log = "xy",
                        data = Data_Animal_de_pd,
                        slope.test = 1/3, 
                        robust = T
)
summary(CoM.SMAla.stance)
plot(CoM.SMAla.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = CoM.SMAla.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
CoM.SMAla.Morphotype <- sma(CoM~La*Morphotype, 
                            log = "xy",
                            data = Data_Animal_de_pd, 
                            slope.test = 1/3, 
                            robust = T
)
summary(CoM.SMAla.Morphotype)
plot(CoM.SMAla.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = CoM.SMAla.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
CoM.SMAla.order <- sma(CoM~La*order, 
                       log = "xy",
                       data = Data_Animal_de_pd, 
                       slope.test = 1/3,
                       robust = T
)
summary(CoM.SMAla.order)
plot(CoM_SMAla_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = CoM.SMAla.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
CoM.SMAla.diet <- sma(CoM~La*Diet, 
                      log = "xy",
                      data = Data_Animal_de_pd,
                      slope.test = 1/3,
                      robust = T
)
summary(CoM.SMAla.diet)
plot(CoM.SMAla.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = CoM.SMAla.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### MoA_m and LA ####
MoA_m.SMAla <- sma(MoA_m~La, 
                   log = "xy",
                   data = Data_Animal_de_pd,
                   slope.test = 1/3, 
                   robust = T
)
summary(MoA_m.SMAla)
plot(MoA_m.SMAla)
# Stance
MoA_m.SMAla.stance <- sma(MoA_m~La*Stance, 
                          log = "xy",
                          data = Data_Animal_de_pd,
                          slope.test = 1/3, 
                          robust = T
)
summary(MoA_m.SMAla.stance)
plot(MoA_m.SMAla.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoA_m.SMAla.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
MoA_m.SMAla.Morphotype <- sma(MoA_m~La*Morphotype, 
                              log = "xy",
                              data = Data_Animal_de_pd, 
                              slope.test = 1/3, 
                              robust = T
)
summary(MoA_m.SMAla.Morphotype)
plot(MoA_m.SMAla.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoA_m.SMAla.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
MoA_m.SMAla.order <- sma(MoA_m~La*order, 
                         log = "xy",
                         data = Data_Animal_de_pd, 
                         slope.test = 1/3,
                         robust = T
)
summary(MoA_m.SMAla.order)
plot(MoA_m_SMAla_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoA_m.SMAla.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
MoA_m.SMAla.diet <- sma(MoA_m~La*Diet, 
                        log = "xy",
                        data = Data_Animal_de_pd,
                        slope.test = 1/3,
                        robust = T
)
summary(MoA_m.SMAla.diet)
plot(MoA_m.SMAla.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoA_m.SMAla.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### P_eq and LA ####
P_eq.SMAla <- sma(P_eq~La, 
                  log = "xy",
                  data = Data_Animal_de_pd,
                  slope.test = 1/3, 
                  robust = T
)
summary(P_eq.SMAla)
plot(P_eq.SMAla)
# Stance
P_eq.SMAla.stance <- sma(P_eq~La*Stance, 
                         log = "xy",
                         data = Data_Animal_de_pd,
                         slope.test = 1/3, 
                         robust = T
)
summary(P_eq.SMAla.stance)
plot(P_eq.SMAla.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = P_eq.SMAla.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
P_eq.SMAla.Morphotype <- sma(P_eq~La*Morphotype, 
                             log = "xy",
                             data = Data_Animal_de_pd, 
                             slope.test = 1/3, 
                             robust = T
)
summary(P_eq.SMAla.Morphotype)
plot(P_eq.SMAla.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = P_eq.SMAla.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
P_eq.SMAla.order <- sma(P_eq~La*order, 
                        log = "xy",
                        data = Data_Animal_de_pd, 
                        slope.test = 1/3,
                        robust = T
)
summary(P_eq.SMAla.order)
plot(P_eq_SMAla_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = P_eq.SMAla.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
P_eq.SMAla.diet <- sma(P_eq~La*Diet, 
                       log = "xy",
                       data = Data_Animal_de_pd,
                       slope.test = 1/3,
                       robust = T
)
summary(P_eq.SMAla.diet)
plot(P_eq.SMAla.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = P_eq.SMAla.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### MoI_G and LA ####
MoI_G.SMAla <- sma(MoI_G~La, 
                   log = "xy",
                   data = Data_Animal_de_pd,
                   slope.test = 1/3, 
                   robust = T
)
summary(MoI_G.SMAla)
plot(MoI_G.SMAla)
# Stance
MoI_G.SMAla.stance <- sma(MoI_G~La*Stance, 
                          log = "xy",
                          data = Data_Animal_de_pd,
                          slope.test = 1/3, 
                          robust = T
)
summary(MoI_G.SMAla.stance)
plot(MoI_G.SMAla.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_G.SMAla.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
MoI_G.SMAla.Morphotype <- sma(MoI_G~La*Morphotype, 
                              log = "xy",
                              data = Data_Animal_de_pd, 
                              slope.test = 1/3, 
                              robust = T
)
summary(MoI_G.SMAla.Morphotype)
plot(MoI_G.SMAla.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_G.SMAla.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
MoI_G.SMAla.order <- sma(MoI_G~La*order, 
                         log = "xy",
                         data = Data_Animal_de_pd, 
                         slope.test = 1/3,
                         robust = T
)
summary(MoI_G.SMAla.order)
plot(MoI_G_SMAla_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_G.SMAla.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
MoI_G.SMAla.diet <- sma(MoI_G~La*Diet, 
                        log = "xy",
                        data = Data_Animal_de_pd,
                        slope.test = 1/3,
                        robust = T
)
summary(MoI_G.SMAla.diet)
plot(MoI_G.SMAla.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_G.SMAla.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### MoI_H and LA ####
MoI_H.SMAla <- sma(MoI_H~La, 
                   log = "xy",
                   data = Data_Animal_de_pd,
                   slope.test = 1/3, 
                   robust = T
)
summary(MoI_H.SMAla)
plot(MoI_H.SMAla)
# Stance
MoI_H.SMAla.stance <- sma(MoI_H~La*Stance, 
                          log = "xy",
                          data = Data_Animal_de_pd,
                          slope.test = 1/3, 
                          robust = T
)
summary(MoI_H.SMAla.stance)
plot(MoI_H.SMAla.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_H.SMAla.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
MoI_H.SMAla.Morphotype <- sma(MoI_H~La*Morphotype, 
                              log = "xy",
                              data = Data_Animal_de_pd, 
                              slope.test = 1/3, 
                              robust = T
)
summary(MoI_H.SMAla.Morphotype)
plot(MoI_H.SMAla.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_H.SMAla.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
MoI_H.SMAla.order <- sma(MoI_H~La*order, 
                         log = "xy",
                         data = Data_Animal_de_pd, 
                         slope.test = 1/3,
                         robust = T
)
summary(MoI_H.SMAla.order)
plot(MoI_H_SMAla_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_H.SMAla.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
MoI_H.SMAla.diet <- sma(MoI_H~La*Diet, 
                        log = "xy",
                        data = Data_Animal_de_pd,
                        slope.test = 1/3,
                        robust = T
)
summary(MoI_H.SMAla.diet)
plot(MoI_H.SMAla.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_H.SMAla.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
#### Alpha and LA ####
Alpha.SMAla <- sma(Alpha~La, 
                   log = "xy",
                   data = Data_Animal_de_pd,
                   slope.test = 1/3, 
                   robust = T
)
summary(Alpha.SMAla)
plot(Alpha.SMAla)
# Stance
Alpha.SMAla.stance <- sma(Alpha~La*Stance, 
                          log = "xy",
                          data = Data_Animal_de_pd,
                          slope.test = 1/3, 
                          robust = T
)
summary(Alpha.SMAla.stance)
plot(Alpha.SMAla.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Alpha.SMAla.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Alpha.SMAla.Morphotype <- sma(Alpha~La*Morphotype, 
                              log = "xy",
                              data = Data_Animal_de_pd, 
                              slope.test = 1/3, 
                              robust = T
)
summary(Alpha.SMAla.Morphotype)
plot(Alpha.SMAla.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Alpha.SMAla.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Alpha.SMAla.order <- sma(Alpha~La*order, 
                         log = "xy",
                         data = Data_Animal_de_pd, 
                         slope.test = 1/3,
                         robust = T
)
summary(Alpha.SMAla.order)
plot(Alpha_SMAla_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Alpha.SMAla.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Alpha.SMAla.diet <- sma(Alpha~La*Diet, 
                        log = "xy",
                        data = Data_Animal_de_pd,
                        slope.test = 1/3,
                        robust = T
)
summary(Alpha.SMAla.diet)
plot(Alpha.SMAla.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Alpha.SMAla.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
#### Omg and LA ####
Omg.SMAla <- sma(Omg~La, 
                 log = "xy",
                 data = Data_Animal_de_pd,
                 slope.test = 1/3, 
                 robust = T
)
summary(Omg.SMAla)
plot(Omg.SMAla)
# Stance
Omg.SMAla.stance <- sma(Omg~La*Stance, 
                        log = "xy",
                        data = Data_Animal_de_pd,
                        slope.test = 1/3, 
                        robust = T
)
summary(Omg.SMAla.stance)
plot(Omg.SMAla.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Omg.SMAla.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Omg.SMAla.Morphotype <- sma(Omg~La*Morphotype, 
                            log = "xy",
                            data = Data_Animal_de_pd, 
                            slope.test = 1/3, 
                            robust = T
)
summary(Omg.SMAla.Morphotype)
plot(Omg.SMAla.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Omg.SMAla.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Omg.SMAla.order <- sma(Omg~La*order, 
                       log = "xy",
                       data = Data_Animal_de_pd, 
                       slope.test = 1/3,
                       robust = T
)
summary(Omg.SMAla.order)
plot(Omg_SMAla_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Omg.SMAla.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Omg.SMAla.diet <- sma(Omg~La*Diet, 
                      log = "xy",
                      data = Data_Animal_de_pd,
                      slope.test = 1/3,
                      robust = T
)
summary(Omg.SMAla.diet)
plot(Omg.SMAla.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Omg.SMAla.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
#### Vsup and LA ####
Vsup.SMAla <- sma(Vsup~La, 
                  log = "xy",
                  data = Data_Animal_de_pd,
                  slope.test = 1/3, 
                  robust = T
)
summary(Vsup.SMAla)
plot(Vsup.SMAla)
# Stance
Vsup.SMAla.stance <- sma(Vsup~La*Stance, 
                         log = "xy",
                         data = Data_Animal_de_pd,
                         slope.test = 1/3, 
                         robust = T
)
summary(Vsup.SMAla.stance)
plot(Vsup.SMAla.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Vsup.SMAla.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Vsup.SMAla.Morphotype <- sma(Vsup~La*Morphotype, 
                             log = "xy",
                             data = Data_Animal_de_pd, 
                             slope.test = 1/3, 
                             robust = T
)
summary(Vsup.SMAla.Morphotype)
plot(Vsup.SMAla.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Vsup.SMAla.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Vsup.SMAla.order <- sma(Vsup~La*order, 
                        log = "xy",
                        data = Data_Animal_de_pd, 
                        slope.test = 1/3,
                        robust = T
)
summary(Vsup.SMAla.order)
plot(Vsup_SMAla_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Vsup.SMAla.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Vsup.SMAla.diet <- sma(Vsup~La*Diet, 
                       log = "xy",
                       data = Data_Animal_de_pd,
                       slope.test = 1/3,
                       robust = T
)
summary(Vsup.SMAla.diet)
plot(Vsup.SMAla.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Vsup.SMAla.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### MoI_HC and LA ####
MoI_HC.SMAla <- sma(MoI_HC~La, 
                    log = "xy",
                    data = Data_Animal_de_pd,
                    slope.test = 1/3, 
                    robust = T
)
summary(MoI_HC.SMAla)
plot(MoI_HC.SMAla)
# Stance
MoI_HC.SMAla.stance <- sma(MoI_HC~La*Stance, 
                           log = "xy",
                           data = Data_Animal_de_pd,
                           slope.test = 1/3, 
                           robust = T
)
summary(MoI_HC.SMAla.stance)
plot(MoI_HC.SMAla.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_HC.SMAla.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
MoI_HC.SMAla.Morphotype <- sma(MoI_HC~La*Morphotype, 
                               log = "xy",
                               data = Data_Animal_de_pd, 
                               slope.test = 1/3, 
                               robust = T
)
summary(MoI_HC.SMAla.Morphotype)
plot(MoI_HC.SMAla.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_HC.SMAla.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
MoI_HC.SMAla.order <- sma(MoI_HC~La*order, 
                          log = "xy",
                          data = Data_Animal_de_pd, 
                          slope.test = 1/3,
                          robust = T
)
summary(MoI_HC.SMAla.order)
plot(MoI_HC_SMAla_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_HC.SMAla.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
MoI_HC.SMAla.diet <- sma(MoI_HC~La*Diet, 
                         log = "xy",
                         data = Data_Animal_de_pd,
                         slope.test = 1/3,
                         robust = T
)
summary(MoI_HC.SMAla.diet)
plot(MoI_HC.SMAla.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = MoI_HC.SMAla.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### Alpha_C and LA ####
Alpha_C.SMAla <- sma(Alpha_C~La, 
                     log = "xy",
                     data = Data_Animal_de_pd,
                     slope.test = 1/3, 
                     robust = T
)
summary(Alpha_C.SMAla)
plot(Alpha_C.SMAla)
# Stance
Alpha_C.SMAla.stance <- sma(Alpha_C~La*Stance, 
                            log = "xy",
                            data = Data_Animal_de_pd,
                            slope.test = 1/3, 
                            robust = T
)
summary(Alpha_C.SMAla.stance)
plot(Alpha_C.SMAla.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Alpha_C.SMAla.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Alpha_C.SMAla.Morphotype <- sma(Alpha_C~La*Morphotype, 
                                log = "xy",
                                data = Data_Animal_de_pd, 
                                slope.test = 1/3, 
                                robust = T
)
summary(Alpha_C.SMAla.Morphotype)
plot(Alpha_C.SMAla.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Alpha_C.SMAla.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Alpha_C.SMAla.order <- sma(Alpha_C~La*order, 
                           log = "xy",
                           data = Data_Animal_de_pd, 
                           slope.test = 1/3,
                           robust = T
)
summary(Alpha_C.SMAla.order)
plot(Alpha_C_SMAla_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Alpha_C.SMAla.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Alpha_C.SMAla.diet <- sma(Alpha_C~La*Diet, 
                          log = "xy",
                          data = Data_Animal_de_pd,
                          slope.test = 1/3,
                          robust = T
)
summary(Alpha_C.SMAla.diet)
plot(Alpha_C.SMAla.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Alpha_C.SMAla.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### OmgC and LA ####
OmgC.SMAla <- sma(OmgC~La, 
                  log = "xy",
                  data = Data_Animal_de_pd,
                  slope.test = 1/3, 
                  robust = T
)
summary(OmgC.SMAla)
plot(OmgC.SMAla)
# Stance
OmgC.SMAla.stance <- sma(OmgC~La*Stance, 
                         log = "xy",
                         data = Data_Animal_de_pd,
                         slope.test = 1/3, 
                         robust = T
)
summary(OmgC.SMAla.stance)
plot(OmgC.SMAla.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = OmgC.SMAla.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
OmgC.SMAla.Morphotype <- sma(OmgC~La*Morphotype, 
                             log = "xy",
                             data = Data_Animal_de_pd, 
                             slope.test = 1/3, 
                             robust = T
)
summary(OmgC.SMAla.Morphotype)
plot(OmgC.SMAla.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = OmgC.SMAla.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
OmgC.SMAla.order <- sma(OmgC~La*order, 
                        log = "xy",
                        data = Data_Animal_de_pd, 
                        slope.test = 1/3,
                        robust = T
)
summary(OmgC.SMAla.order)
plot(OmgC_SMAla_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = OmgC.SMAla.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
OmgC.SMAla.diet <- sma(OmgC~La*Diet, 
                       log = "xy",
                       data = Data_Animal_de_pd,
                       slope.test = 1/3,
                       robust = T
)
summary(OmgC.SMAla.diet)
plot(OmgC.SMAla.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = OmgC.SMAla.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)

#### VsupC and LA ####
VsupC.SMAla <- sma(VsupC~La, 
                   log = "xy",
                   data = Data_Animal_de_pd,
                   slope.test = 1/3, 
                   robust = T
)
summary(VsupC.SMAla)
plot(VsupC.SMAla)
# Stance
VsupC.SMAla.stance <- sma(VsupC~La*Stance, 
                          log = "xy",
                          data = Data_Animal_de_pd,
                          slope.test = 1/3, 
                          robust = T
)
summary(VsupC.SMAla.stance)
plot(VsupC.SMAla.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = VsupC.SMAla.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
VsupC.SMAla.Morphotype <- sma(VsupC~La*Morphotype, 
                              log = "xy",
                              data = Data_Animal_de_pd, 
                              slope.test = 1/3, 
                              robust = T
)
summary(VsupC.SMAla.Morphotype)
plot(VsupC.SMAla.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = VsupC.SMAla.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
VsupC.SMAla.order <- sma(VsupC~La*order, 
                         log = "xy",
                         data = Data_Animal_de_pd, 
                         slope.test = 1/3,
                         robust = T
)
summary(VsupC.SMAla.order)
plot(VsupC_SMAla_order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = VsupC.SMAla.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
VsupC.SMAla.diet <- sma(VsupC~La*Diet, 
                        log = "xy",
                        data = Data_Animal_de_pd,
                        slope.test = 1/3,
                        robust = T
)
summary(VsupC.SMAla.diet)
plot(VsupC.SMAla.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = VsupC.SMAla.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
#### Maxspeed_ms and LA ####
Maxspeed_ms.SMAla <- sma(Maxspeed_ms~La, 
                         log = "xy",
                         data = Data_Animal_de_pd,
                         slope.test = 1/3, 
                         robust = T
)
summary(Maxspeed_ms.SMAla)
plot(Maxspeed_ms.SMAla)
# Stance
Maxspeed_ms.SMAla.stance <- sma(Maxspeed_ms~La*Stance, 
                                log = "xy",
                                data = Data_Animal_de_pd,
                                slope.test = 1/3, 
                                robust = T
)
summary(Maxspeed_ms.SMAla.stance)
plot(Maxspeed_ms.SMAla.stance, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Maxspeed_ms.SMAla.stance[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Morphotype
Maxspeed_ms.SMAla.Morphotype <- sma(Maxspeed_ms~La*Morphotype, 
                                    log = "xy",
                                    data = Data_Animal_de_pd, 
                                    slope.test = 1/3, 
                                    robust = T
)
summary(Maxspeed_ms.SMAla.Morphotype)
plot(Maxspeed_ms.SMAla.Morphotype, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Maxspeed_ms.SMAla.Morphotype[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Order
Maxspeed_ms.SMAla.order <- sma(Maxspeed_ms~La*order, 
                               log = "xy",
                               data = Data_Animal_de_pd, 
                               slope.test = 1/3,
                               robust = T
)
summary(Maxspeed_ms.SMAla.order)
plot(Maxspeed_ms.SMAla.order, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Maxspeed_ms.SMAla.order[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)
# Diet
Maxspeed_ms.SMAla.diet <- sma(Maxspeed_ms~La*Diet, 
                              log = "xy",
                              data = Data_Animal_de_pd,
                              slope.test = 1/3,
                              robust = T
)
summary(Maxspeed_ms.SMAla.diet)
plot(Maxspeed_ms.SMAla.diet, col = c(4,7,6,9, 10, 11, 20),)
legend(x = "topleft",
       legend = Maxspeed_ms.SMAla.diet[["groupsummary"]][["group"]],
       col = c(4,7,6,9, 10, 11, 20),
       cex = 0.75,
       bty = 'n',
       lty = 1)














################################################################################
# for phylo least squares
Long.PGLS = pgls(log(Long)~log(massAvg), data = CompData2, lambda = 'ML')
summary(Long.PGLS)
plot(Long.PGLS)
Long.GLS = gls(log(Long)~log(massAvg), data = Data_Animal_de_pd)
summary(Long.GLS)




# PCA ######################################
#### PCA analysis to check the relation between variables ESPECIES
# PCA analysis (for incomplete data)
# Get the qualitative information
Stc.Data_Animal_de_pd <- Data_Animal_de_pd[,c(10)] #Stance information
Ord.Data_Animal_de_pd <- Data_Animal_de_pd[,c(27)] # Order information
Diet.Data_Animal_de_pd <- Data_Animal_de_pd[,c(15)] # Diet information

# logarithmic transformation
log.Data_Animal_de_pd <- 
  log(Data_Animal_de_pd[,c(4, 31, 32, 33, 34, 35)]) # speeds 28, 29, 30


# IMPUTATION: to complete the data if needed

# check the plausibility of doing imputation
# estimate the number of components from incomplete data
nb <- estim_ncpPCA(log.Data_Animal_de_pd, scale = TRUE) 
# For the uncertainty of imputing values
mi.AnimalPCA_Comp<-MIPCA(log.Data_Animal_de_pd, ncp=nb$ncp, scale = TRUE) 
plot(mi.AnimalPCA_Comp)

# imputation of the values
# iterativePCA algorithm
AnimalPCA.comp <- imputePCA(log.Data_Animal_de_pd, 
                            ncp = nb$ncp, 
                            scale = TRUE, 
                            maxiter = 10000) 
summary(AnimalPCA.comp)
# transforming it as DB: if needed later
AnimalComplete <- as.data.frame(AnimalPCA.comp$completeObs) 

# PCA analysis through package prcomp
AnimalPCA <- prcomp(AnimalPCA.comp$completeObs, 
                    center = TRUE, 
                    scale. = TRUE)
summary(AnimalPCA)
print(AnimalPCA) # returns the standard deviation of each of the PCs
plot(AnimalPCA, type = "l") # variances (y-axis) associated with the PCs (x-axis)
summary(AnimalPCA) # check cumulative variance and take until is 95% (at least)
AnimalPCA$center # mean of the variables
AnimalPCA$scale # standard deviation
AnimalPCA$rotation #  the principal component loading
#AnimalPCA$x #principal component score vectors
AnimalPCA.sdv <- AnimalPCA$sdev # standard deviation
AnimalPCA.var <- AnimalPCA.sdv^2 # variance

AnimalPCA.g <- ggbiplot(AnimalPCA, 
                        obs.scale = 1, 
                        var.scale = 1, 
                        #choices=c(2,3), # to change which components to plot
                        groups = Stc.Data_Animal_de_pd, ellipse = TRUE, 
                        ellipse.prob = 0.68,
                        circle = TRUE) # PCA graphic
AnimalPCA.g <- AnimalPCA.g + scale_color_discrete(name = '')+ theme_minimal()
print(AnimalPCA.g)

# PCA analysis through package FactoMineR
AnimalPCA.MineR <- PCA(AnimalPCA.comp$completeObs) #PCA
summary(AnimalPCA.MineR)
dimdesc(AnimalPCA.MineR) #correlation of each variable on the first dimension
round(AnimalPCA.MineR$ind$cos2,2)
dim_desc(AnimalPCA.MineR)

# Dimensionless
#### PCA analysis to check the relation between DIMENSIONLESS parameters
# get data from the original data frame
lessDim = data.frame(Data_Animal_de_pd$La, 
                     Data_Animal_de_pd$Long, 
                     Data_Animal_de_pd$Rmin, 
                     Data_Animal_de_pd$Rmax, 
                     Data_Animal_de_pd$Rmed, 
                     Data_Animal_de_pd$massAvg)

colnames(lessDim) <- c('La','Long','Rmin','Rmax', 'Rmed', 'massAvg')

# dimensionless parameters, mass/1kg
lessDim['La_LH'] <- lessDim['La']/lessDim['Long']
lessDim['Rmin_LH'] <- lessDim['Rmin']/lessDim['Long']
lessDim['Rmax_LH'] <- lessDim['Rmax']/lessDim['Long']
lessDim['Rmed_LH'] <- lessDim['Rmed']/lessDim['Long']

lessDim <- na.omit(lessDim)

# logarithmic transformation
log.lessDim <- log(lessDim[,c(6, 7, 8, 9, 10)]) # speeds 28, 29, 30

# PCA analysis through package prcomp
lessDimPCA <- prcomp(log.lessDim, center = TRUE, scale. = TRUE) #PCA
summary(lessDimPCA)
print(lessDimPCA) # returns the standard deviation of each of the PCs
plot(lessDimPCA, type = "l") # plot of the variances (y-axis) associated with the PCs (x-axis)
summary(lessDimPCA) # check cumulative variance and take until is 95% (at least)
lessDimPCA$center # mean of the variables
lessDimPCA$scale # standard deviation
lessDimPCA$rotation #  the principal component loading
#lessDimPCA$x #principal component score vectors
lessDimPCA.sdv <- lessDimPCA$sdev # standard deviation
lessDimPCA.var <- lessDimPCA.sdv^2 # variance

lessDimPCA.g <- ggbiplot(lessDimPCA, 
                         obs.scale = 1, 
                         var.scale = 1, 
                         choices=c(1,2),
                         ellipse.prob = 0.68,
                         circle = TRUE) # PCA graphic
lessDimPCA.g <- lessDimPCA.g + scale_color_discrete(name = '') + theme_minimal()
print(lessDimPCA.g)



######### SCRAPES ###############
# to get info of orders
#ddply(Data_Animal_de_pd, .(order, Stance), nrow)
#Data_Animal_de_pd['order']=='Tubulidentata'

#Pandas df for Bones
#Data_Bone_de_pd = pd$read_pickle(Data_Bone_pd)
#summary(Data_Bone_de_pd)

# read data from de .csv format
# .csv df for especies and bones
#Data_Animal <- read.csv(Data_Animal_csv, header = TRUE, sep = ",", row.names = 1)
#summary(Data_Animal)

# get out-'liars'
Data_Animal_de_pd <- 
  Data_Animal_de_pd[!(Data_Animal_de_pd$NCBI %in% c(30554, 9627, 9802, 71006)),]



# Files definition (change the directories accordingly)

#Data_Bone_csv = "/home/kale/Documents/Allometry/DATA/Bone_out.csv"
#Data_Bone_pd = "/home/kale/Documents/Allometry/DATA/Bone_out.pkl"


# Pandas df for Bones
#Data_Bone_de_pd = pd$read_pickle(Data_Bone_pd)
#summary(Data_Bone_de_pd)

# read data from de .csv format .csv df for bones
#Data_Bone <- read.csv(Data_Bone_csv, header = TRUE, sep = ",", row.names = 1)
#summary(Data_Bone)


#  Morphotype and All ####
plot(Morphotype~massAvg, log = "x", data = Data_Animal_de_pd)

plot(Morphotype~Long, data = Data_Animal_de_pd)
Data_Animal_de_pd['Long_esp'] <- Data_Animal_de_pd['Long']/Data_Animal_de_pd['massAvg']
plot(Morphotype~Long_esp, log = "x",data = Data_Animal_de_pd)

Data_Animal_de_pd['Rmax_Rmin'] <- Data_Animal_de_pd['Rmax']/Data_Animal_de_pd['Rmin']
Data_Animal_de_pd['Rmax_esp'] <- Data_Animal_de_pd['Rmax']/Data_Animal_de_pd['massAvg']
plot(Morphotype~Rmax_Rmin, data = Data_Animal_de_pd)
plot(Morphotype~Rmax_esp, log = "x", data = Data_Animal_de_pd)

plot(Morphotype~Rmin, data = Data_Animal_de_pd)
Data_Animal_de_pd['Rmin_esp'] <- Data_Animal_de_pd['Rmin']/Data_Animal_de_pd['massAvg']
plot(Morphotype~Rmin_esp, log = "x", data = Data_Animal_de_pd)

plot(Morphotype~Rmed, data = Data_Animal_de_pd)
Data_Animal_de_pd['Rmed_esp'] <- Data_Animal_de_pd['Rmed']/Data_Animal_de_pd['massAvg']
plot(Morphotype~Rmed_esp, log = "x", data = Data_Animal_de_pd)

plot(Morphotype~angeqv, log="x", data = Data_Animal_de_pd)
Data_Animal_de_pd['angeqv_esp'] <- Data_Animal_de_pd['angeqv']/Data_Animal_de_pd['massAvg']
plot(Morphotype~angeqv_esp, log = "x", data = Data_Animal_de_pd)

plot(Morphotype~Maxspeed_ms, data = Data_Animal_de_pd)
Data_Animal_de_pd['Maxspeed_ms_esp'] <- Data_Animal_de_pd['Maxspeed_ms']/Data_Animal_de_pd['massAvg']
plot(Morphotype~Mscl_Force_esp, log = "x", data = Data_Animal_de_pd)

plot(Morphotype~speedHirt_ms, data = Data_Animal_de_pd)
Data_Animal_de_pd['speedHirt_ms_esp'] <- Data_Animal_de_pd['speedHirt_ms']/Data_Animal_de_pd['massAvg']
plot(Morphotype~speedHirt_ms_esp, log = "x", data = Data_Animal_de_pd)

plot(Morphotype~Pry_Area_m, data = Data_Animal_de_pd)
Data_Animal_de_pd['Pry_Area_m_esp'] <- Data_Animal_de_pd['Pry_Area_m']/Data_Animal_de_pd['massAvg']
plot(Morphotype~Pry_Area_m_esp, log = "x", data = Data_Animal_de_pd)

plot(Morphotype~Sf, data = Data_Animal_de_pd)
Data_Animal_de_pd['Sf_esp'] <- Data_Animal_de_pd['Sf']/Data_Animal_de_pd['massAvg']
plot(Morphotype~Sf_esp, log = "x", data = Data_Animal_de_pd)

plot(Morphotype~Mscl_Force, log="x", data = Data_Animal_de_pd)
Data_Animal_de_pd['Mscl_Force_esp'] <- Data_Animal_de_pd['Mscl_Force']/Data_Animal_de_pd['massAvg']
plot(Morphotype~Mscl_Force_esp, log = "x", data = Data_Animal_de_pd)

plot(Morphotype~MoI, data = Data_Animal_de_pd)
Data_Animal_de_pd['MoI_esp'] <- Data_Animal_de_pd['MoI']/Data_Animal_de_pd['massAvg']
plot(Morphotype~MoI_esp, log = "x", data = Data_Animal_de_pd)

plot(Morphotype~Limb_mass, data = Data_Animal_de_pd)
Data_Animal_de_pd['Limb_mass_esp'] <- Data_Animal_de_pd['Limb_mass']/Data_Animal_de_pd['massAvg']
plot(Morphotype~Limb_mass_esp, log = "x", data = Data_Animal_de_pd)

plot(Morphotype~CoM, data = Data_Animal_de_pd)
Data_Animal_de_pd['CoM_esp'] <- Data_Animal_de_pd['CoM']/Data_Animal_de_pd['massAvg']
plot(Morphotype~Limb_mass_esp, log = "x", data = Data_Animal_de_pd)

plot(Morphotype~MoI_G, data = Data_Animal_de_pd)
Data_Animal_de_pd['MoI_G_esp'] <- Data_Animal_de_pd['MoI_G']/Data_Animal_de_pd['massAvg']
plot(Morphotype~Limb_mass_esp, log = "x", data = Data_Animal_de_pd)

plot(Morphotype~MoA_m, data = Data_Animal_de_pd)
Data_Animal_de_pd['MoA_m_esp'] <- Data_Animal_de_pd['MoA_m']/Data_Animal_de_pd['massAvg']
plot(Morphotype~MoA_m_esp, log = "x", data = Data_Animal_de_pd)

plot(Morphotype~MoI_H, log = "x", data = Data_Animal_de_pd)
Data_Animal_de_pd['MoI_H_esp'] <- Data_Animal_de_pd['MoI_H']/Data_Animal_de_pd['massAvg']
plot(Morphotype~MoI_H_esp, log = "x", data = Data_Animal_de_pd)

plot(Morphotype~Alpha, log="x", data = Data_Animal_de_pd)
Data_Animal_de_pd['Alpha_esp'] <- Data_Animal_de_pd['Alpha']/Data_Animal_de_pd['massAvg']
plot(Morphotype~Alpha_esp, log = "x", data = Data_Animal_de_pd)

plot(Morphotype~Omg, log = "x", data = Data_Animal_de_pd)
Data_Animal_de_pd['Omg_esp'] <- Data_Animal_de_pd['Omg']/Data_Animal_de_pd['massAvg']
plot(Morphotype~Omg_esp, log = "x", data = Data_Animal_de_pd)



Data_Animal_de_pd['VratioP'] <- Data_Animal_de_pd['Vsup']/Data_Animal_de_pd['P_eq']
plot(VratioP~massAvg, log = "xy", data = Data_Animal_de_pd)
