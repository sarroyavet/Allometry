# Clean everything ####
rm(list = ls())

library(reticulate) # import python scripts to R
# use_python("/usr/bin/python3.6", required = T)
library(ape)
library(caper)
library(dplyr)
library(mosaic)
library(phytools)
library(Polychrome)
library(rlang)
library(smatr) # Standardized Major Axis Regression (SMA ou RMA ou Model II)
library(viridis)

pd <- import('pandas')

  # Files definition (change the directories accordingly)
Data_Animal_csv = "/home/kale/Documents/Allometry/DATA/Esp_out.csv"
Data_Animal_pd =  "/home/kale/Documents/Allometry/DATA/Esp_out.pkl"

################################################################################
##############                      READ DATA                     ##############
 # Pandas df for Species
Data_Animal_de_pd = pd$read_pickle(Data_Animal_pd)
class(Data_Animal_de_pd)
summary(Data_Animal_de_pd)

Data_Animal <- read.csv(Data_Animal_csv, header = TRUE, sep = ",", row.names = 1)
class(Data_Animal)
summary(Data_Animal)

# UNIT TRANSFORMATION
 # from kg to gr
Data_Animal_de_pd['massAvg_gr'] <- Data_Animal_de_pd['massAvg']/1000
 # transform speed from km/h to m/s
Data_Animal_de_pd['Maxspeed_ms'] <- Data_Animal_de_pd['Maxspeed']/3.6
Data_Animal_de_pd['speedHirt_ms'] <- Data_Animal_de_pd['speedHirt']/3.6
Data_Animal_de_pd['spdMissing_ms'] <- Data_Animal_de_pd['spdMissing']/3.6
 # transform measurement from mm to m
Data_Animal_de_pd['Long_m'] <- Data_Animal_de_pd['Long']*10^-3
Data_Animal_de_pd['La_m'] <- Data_Animal_de_pd['La']*10^-3
Data_Animal_de_pd['Rmed_m'] <- Data_Animal_de_pd['Rmed']*10^-3
Data_Animal_de_pd['Rmax_m'] <- Data_Animal_de_pd['Rmax']*10^-3
Data_Animal_de_pd['Rmin_m'] <- Data_Animal_de_pd['Rmin']*10^-3

################################################################################
##############         OTHER CALCULATIONS FROM MEASUREMENT        ##############
 # Difference between radii in m
Data_Animal_de_pd['DeltaR'] <- 
  Data_Animal_de_pd['Rmax'] - Data_Animal_de_pd['Rmin'] # Delta radii mm
Data_Animal_de_pd['DeltaR_m'] <- 
  Data_Animal_de_pd['Rmax_m'] - Data_Animal_de_pd['Rmin_m'] # Delta radii m
 # mean diameter
Data_Animal_de_pd['Dmed'] <- Data_Animal_de_pd['Rmed']*2 # Diameter avg mm
Data_Animal_de_pd['Dmed_m'] <- Data_Animal_de_pd['Rmed_m']*2 # Diameter avg m
 # Angle equivalent
Data_Animal_de_pd['angeqv'] <- Data_Animal_de_pd['DeltaR_m'] /
  Data_Animal_de_pd['La_m']# equivalent tangent
 # Projected area
Data_Animal_de_pd['Pry_Area'] <- 
  Data_Animal_de_pd['La']*Data_Animal_de_pd['Dmed'] # Projected area mm2
Data_Animal_de_pd['Pry_Area_m'] <- 
  Data_Animal_de_pd['La_m']*Data_Animal_de_pd['Dmed_m'] # Projected area m2

################################################################################
##############         CALCULATIONS BASED ON THE LITERATURE       ##############
#### Stride frequency ####
Data_Animal_de_pd['Sf'] <- 
  4.70*Data_Animal_de_pd['massAvg']^-0.162 # Heglund & Taylor 1988(Speed, stride, fq)
#### Muscle mass #### 
Data_Animal_de_pd['Mscl_mass_kg'] <- 
  (6.2 * Data_Animal_de_pd['massAvg']^1.11)/1000
#### Muscle characteristic length #### 
Data_Animal_de_pd['Mscl_length_m'] <-
  (18.7 * Data_Animal_de_pd['massAvg']^0.33)/1000
#### Muscle cross sectional area #### 
Data_Animal_de_pd['Mscl_A_m'] <- 
  Data_Animal_de_pd['Mscl_mass_kg']/(1060*Data_Animal_de_pd['Mscl_length_m'])
#### Muscle force #### 
Data_Animal_de_pd['Mscl_Force'] <- 
  0.25*10^6 * Data_Animal_de_pd['Mscl_A_m']
#### Coatham's moment of inertia of the forearm ####
Data_Animal_de_pd['MoI_H_C'] <- 
  (10^(-4.751156)*Data_Animal_de_pd['massAvg']^(1.775433)) # kg/m2
#### Moment of arm of the triceps ####
Data_Animal_de_pd['MoA_m']<- # Alexander, 1981
  (8.7*(10^-3)*Data_Animal_de_pd['massAvg']^(0.41)) # mm *10^-3 to m
################################################################################
# CALCULATIONS FROM DEDUCED EQUATIONS #####

#### Joint pressure #### 
Data_Animal_de_pd['P_eq'] <- 
  Data_Animal_de_pd['Mscl_Force']/Data_Animal_de_pd['Pry_Area_m'] # kg/m2
summary(Data_Animal_de_pd['P_eq'])
Var_P_eq = var(Data_Animal_de_pd['P_eq'], na.rm = T)
sd_P_eq = Var_P_eq^0.5
sd_P_eq

#### Angular acceleration of the elbow COATHAM ####
Data_Animal_de_pd['Alpha_C'] <- 
  Data_Animal_de_pd['Mscl_Force']*Data_Animal_de_pd['MoA_m'] /
  (Data_Animal_de_pd['MoI_H_C']) 
#### Angular vel from Coathams ####
Data_Animal_de_pd['Omg_C'] <- Data_Animal_de_pd['Alpha_C']*0.254/
  (Data_Animal_de_pd['Sf']) # Angular velocity of the joint for only 0.254 of 
                            # the period (assuming it is when we have the 
                            # max speed)

#### Sliding velocity Coathams ####
Data_Animal_de_pd['Vsup_C'] <- Data_Animal_de_pd['Omg_C']*
  Data_Animal_de_pd['Dmed_m']/2 
################################################################################
###################    STANDARD MAJOR AXIS REGRESSION     #####################
# sma simple between two parameters
getsma <- function(P1, P2, dataf, slptst, Robust){
  P1 <- rlang::enexpr(P1)
  P2 <- rlang::enexpr(P2)
  SMAt <- rlang::inject(sma(!!P1 ~ !!P2, 
                            log = "xy",
                            data = dataf,
                            slope.test = slptst, 
                            robust = Robust
  ))
  return(SMAt)
  }

# sma by groups
getsmaGroup <- function(P1, P2, Group, dataf, slptst, Robust){
  P1 <- rlang::enexpr(P1)
  P2 <- rlang::enexpr(P2)
  Group <- rlang::enexpr(Group)
  SMAt <- rlang::inject(sma(!!P1~!!P2*!!Group, 
                            log = "xy",
                            data = Data_Animal_de_pd, 
                            slope.test = slptst,
                            robust = Robust
                            ))
  return(SMAt)
}

# plot groups
plotgrps <- function(sma_res){
  plot(sma_res, col = c(4, 7, 6, 9, 10, 11, 20),)
  legend(x = "topleft",
         legend = sma_res[["groupsummary"]][["group"]],
         col = c(4, 7, 6, 9, 10, 11, 20),
         cex = 0.75,
         bty = 'n',
         lty = 1)
}

#### getSMAs ####
# list of parameters to evaluate
param <- list('massAvg', 'Long', 'La', 'Dmed', 'Rmax', 'Rmin', 'DeltaR', 
              'DeltaR_m', 'angeqv', 'Pry_Area_m', 'Sf', 'Mscl_mass_kg', 
              'Mscl_length_m', 'Mscl_Force', 'P_eq', 'MoI_H_C', 'Alpha_C',
              'Omg_C', 'Vsup_C')
param <- lapply(param, function(x){rlang::sym(x)})
# combination of the parameters for the regressions
comb <- combn(param, 2, simplify = FALSE)

# list of groups
Groups <- list('Stance', 'Morphotype', 'order', 'family', 'Diet')
Groups <- lapply(Groups, function(x){rlang::sym(x)})

# perform the SMA with the parameters combined (paired)
for (p1p2 in comb) {
  p1p2
  P1 <- p1p2[[1]]
  P2 <- p1p2[[2]]
  trai <- try(assign(paste(P2,'.SMA.',P1, sep = ''), 
                     getsma(!!P2, !!P1, Data_Animal_de_pd, 1/3, TRUE)),
              silent = TRUE)
  if (class(trai) != 'try-error') {
    var <- sym(paste(P2,'.SMA.',P1, sep = ''))
    print(var)
    # rlang::inject(summary(!!var))
  }
  for (g in Groups) {
      trai <- try(assign(paste(P2,'.SMA.',P1,'.',g, sep = ''), 
             getsmaGroup(!!P2, !!P1, !!g, Data_Animal_de_pd, 1/3, TRUE)),
             silent = TRUE)
  }
}

# use the code below to find sma for specific things
# test <- getsma(Alpha_C, massAvg, Data_Animal_de_pd, 1/3, TRUE)
# test
# plot(test)

################################################################################
####################             PHYLO_eso              ########################
#### Get data ####
# get the names of the animals with _ as the space
subs <- function(x){
  names <- list()
  for (an in x) {
    # if (an == "Ceratotherium simum cottoni"){
    #   names <- append(names, 'Ceratotherium simum')
    # }else 
      if (an == "Neogale vison") {
      names <- append(names, 'Neovison vison')
    }else if (an == "Equus asinus") {
      names <- append(names, 'Equus africanus')
    }else{
      names <- append(names, an)
    }
  }
  names <- sub(" ", "_", names)
  return(names)
}

Data_Animal_de_pd['Sci_name2'] <- subs(Data_Animal_de_pd$Sci_name)

# read the tree data
treesfile = '/home/kale/Documents/Allometry/DATA/Phylo/output.nex'

Trees = read.nexus(file = treesfile)
summary(Trees)

#### function for several parameters -> pgls ####
getpgls <- function(P1, P2, dataf){
  P1 <- rlang::enexpr(P1)
  P2 <- rlang::enexpr(P2)
  PGLSt <- rlang::inject(pgls(log10(!!P1)~log10(!!P2), data = dataf, lambda = 'ML'))
  return(PGLSt)
}

#### parameters for the regressions ####
# it should be done group by group or it will take to much time
parampgls <- list('Long', 'La', 'Dmed', 'Rmax', 'Rmin', 'DeltaR',
                  'angeqv', 'Pry_Area_m', 'P_eq', 'Vsup_C')
parampgls <- list('Vsup_C')
# parampgls <- list('Long', 'La' )
# parampgls <- list('Rmax', 'Rmin')
# parampgls <- list('DeltaR','Dmed')
# parampgls <- list('angeqv', 'Pry_Area_m')
parampgls <- list('P_eq')
parampgls <- list('Vsup_C')
parampgls <- list('La', 'Dmed', 'angeqv')
parampgls <- list('Long', 'Rmax')
parampgls <- list('Pry_Area_m', 'DeltaR', 'Rmin')

parampgls <- lapply(parampgls, function(x){rlang::sym(x)})

#### Calculate the coefficients for all the trees ####
P1 <- rlang::sym('massAvg')
for (p1p2 in parampgls) {
  int = sym(paste('Int',p1p2, sep = ''))
  slp = sym(paste('slp',p1p2, sep = ''))
  lmb = sym(paste('lmb',p1p2, sep = ''))
  assign(paste('Int',p1p2, sep = ''), list())
  assign(paste('slp',p1p2, sep = ''), list())
  assign(paste('lmb',p1p2, sep = ''), list())
  for (i in 1:length(Trees)) {
    phy <- Trees[[i]]
    # check if all names exists in the tree
    # name.check(phy, data = Data_Animal_de_pd, # check we have data for the tree
             # data.names = Data_Animal_de_pd$Sci_name2)
    # combine phylo with the dataset and ensure structure and ordering
    CompData = comparative.data(phy, 
                                data = Data_Animal_de_pd[,c("massAvg", 
                                                            'Long',
                                                            'Long_m',
                                                            'La',
                                                            'La_m',
                                                            'Rmax',
                                                            'Rmax_m',
                                                            'Rmin',
                                                            'Rmin_m',
                                                            'DeltaR',
                                                            'Dmed',
                                                            'angeqv',
                                                            'Pry_Area_m',
                                                            'P_eq',
                                                            'Vsup_C',
                                                            'Sci_name2')],
                                names.col = 'Sci_name2', 
                                na.omit = TRUE)

    P2 <- p1p2
    trai <- try(assign(paste(P2,'.PGLS.',P1, sep = ''),
                       getpgls(!!P2, !!P1, CompData)), 
                silent = TRUE)
    if (class(trai) != 'try-error') {
      print(p1p2)
      print(i)
      var <- sym(paste(P2,'.PGLS.',P1, sep = ''))
      temp <- rlang::inject(!!var)
      rlang::inject(!!int <- append(!!int,temp[["model"]][["coef"]][["(Intercept)"]]))
      rlang::inject(!!slp <- append(!!slp,temp[["model"]][["coef"]][["log10(massAvg)"]]))
      rlang::inject(!!lmb <- append(!!lmb,temp[["param.CI"]][["lambda"]][["opt"]][["lambda"]]))
    } else {
      print(p1p2)
      print(i)
      print('no')
      var <- sym(paste(P2,'.PGLS.',P1, sep = ''))
      temp <- rlang::inject(!!var)
      rlang::inject(!!int <- append(!!int,NA))
      rlang::inject(!!slp <- append(!!slp,NA))
      rlang::inject(!!lmb <- append(!!lmb,NA))
    }
  }
}

# ##### SAVE ####
# parampgls <- list('Long', 'La', 'Dmed', 'Rmax', 'Rmin', 'DeltaR',
#                   'angeqv', 'Pry_Area_m', 'P_eq', 'Vsup_C')
# # save info, just in case
# for (p1p2 in parampgls) {
#   file1 = paste('/home/kale/Documents/Allometry/DATA/pgls/','Int',p1p2,'.rds', sep = '')
#   file2 = paste('/home/kale/Documents/Allometry/DATA/pgls/','slp',p1p2,'.rds', sep = '')
#   file3 = paste('/home/kale/Documents/Allometry/DATA/pgls/','lmb',p1p2,'.rds', sep = '')
#   int = sym(paste('Int',p1p2, sep = ''))
#   slp = sym(paste('slp',p1p2, sep = ''))
#   lmb = sym(paste('lmb',p1p2, sep = ''))
#   rlang::inject(saveRDS(!!int, file = file1))
#   rlang::inject(saveRDS(!!slp, file = file2))
#   rlang::inject(saveRDS(!!lmb, file = file3))
# }

##### LOAD ####
# parampgls <- list('Long', 'La', 'Dmed', 'Rmax', 'Rmin', 'DeltaR',
#                   'angeqv', 'Pry_Area_m', 'P_eq', 'Vsup_C')
# for (p1p2 in parampgls) {
#   file1 = paste('/home/kale/Documents/Allometry/DATA/pgls/','Int',p1p2,'.rds', sep = '')
#   file2 = paste('/home/kale/Documents/Allometry/DATA/pgls/','slp',p1p2,'.rds', sep = '')
#   file3 = paste('/home/kale/Documents/Allometry/DATA/pgls/','lmb',p1p2,'.rds', sep = '')
#   assign(paste('Int',p1p2, sep = ''), readRDS(file1))
#   assign(paste('slp',p1p2, sep = ''), readRDS(file1))
#   assign(paste('lmb',p1p2, sep = ''), readRDS(file1))
# }

#### Draw the histograms of the coefficients and lambda ####
parampgls <- list('Long', 'La', 'Rmax', 'Rmin', 'DeltaR','Dmed',
                  'angeqv', 'Pry_Area_m', 'P_eq', 'Vsup_C')

for (p1p2 in parampgls) {
  W = 8
  H = 4
  int = sym(paste('Int',p1p2, sep = ''))
  slp = sym(paste('slp',p1p2, sep = ''))
  lmb = sym(paste('lmb',p1p2, sep = ''))
  ffile = paste('/home/kale/Documents/PaperAllometry/images/Rplots/histo/',p1p2,sep = '')
  #
  tmpint <- rlang::inject(na.omit(unlist(!!int)))
  tmpslp <- rlang::inject(na.omit(unlist(!!slp)))
  tmplmb <- rlang::inject(na.omit(unlist(!!lmb)))
  svg(filename = paste(ffile,'_int.svg',sep = ''),
      width = W,
      height = H, 
      family = 'serif')
  histint <- hist(tmpint, 
                  breaks = 20, 
                  freq = TRUE,  
                  col = 'white', 
                  border = rgb(0.27, 0, 0.33, 1),
                  density = 10,
                  xlab = paste('int_', p1p2, sep = ''),
                  main = paste('int_', p1p2, sep = ''))
  rlang::inject(abline(v = mean(unlist(!!int),
                                na.rm = TRUE), 
                              col = rgb(0.05, 0.33, 0, 1), 
                              lwd = 2, 
                              lty = 2
                              ))
  dev.off() 
  #
  svg(filename = paste(ffile,'_slp.svg',sep = ''),
      width = W,
      height = H, 
      family = 'serif')
  histint <- hist(tmpslp, 
                  breaks = 20, 
                  freq = TRUE,  
                  col = 'white', 
                  border = rgb(0.27, 0, 0.33, 1),
                  density = 10,
                  xlab = paste('slp_', p1p2, sep = ''),
                  main = paste('slp_', p1p2, sep = ''))
  rlang::inject(abline(v = mean(unlist(!!slp), 
                                na.rm = TRUE), 
                       col = rgb(0.05, 0.33, 0, 1), 
                       lwd = 2, 
                       lty = 2
  ))
  dev.off()
  #
  svg(filename = paste(ffile,'_lmb.svg',sep = ''),
      width = W,
      height = H, 
      family = 'serif')
  histint <- hist(tmplmb, 
                  breaks = 20, 
                  freq = TRUE,  
                  col = 'white', 
                  border = rgb(0.27, 0, 0.33, 1),
                  density = 10,
                  xlab = paste('lmb_', p1p2, sep = ''),
                  main = paste('lmb_', p1p2, sep = ''))
  rlang::inject(abline(v = mean(unlist(!!lmb), 
                                na.rm = TRUE), 
                       col = rgb(0.05, 0.33, 0, 1), 
                       lwd = 2, 
                       lty = 2
  ))
  dev.off()
}

#### Get the representative tree ####
parampgls <- list('Long', 'La', 'Dmed', 'Rmax', 'Rmin',
                  'angeqv', 'Pry_Area_m', 'P_eq', 'Vsup_C')
aveg <- list() #list of the distances
avtree <- list() #list of averages trees
for (p1p2 in parampgls) {
  int = sym(paste('Int',p1p2, sep = ''))
  slp = sym(paste('slp',p1p2, sep = ''))
  lmb = sym(paste('lmb',p1p2, sep = ''))
  dist <- rlang::inject((unlist(!!int) - mean(unlist(!!int), na.rm = TRUE))^2 +
                        (unlist(!!slp) - mean(unlist(!!slp), na.rm = TRUE))^2 +
                        (unlist(!!lmb) - mean(unlist(!!lmb), na.rm = TRUE))^2)
  closest <- which.min(dist)
  aveg <- append(aveg, dist[closest])
  avtree <- append(avtree, closest)
  print(p1p2)
  print(closest)
  print(dist[closest])
}
MeanTree <- avtree[which.min(unlist(aveg))]

#### List of orders and families ####
GRorders = data.frame()
GRorders <- unique(Data_Animal_de_pd[c('order')])
colsord <- setNames(kelly.colors(length(unique(GRorders$order))),
                    unique(GRorders$order))

GRfamilia = data.frame()
GRfamilia <- unique(Data_Animal_de_pd[c('family')])
colsfamy <- setNames(palette.colors(length(unique(GRfamilia$family)),
                                    palette = 'polychrome'),
                     unique(GRfamilia$family))

#### Get the pgls for the mean tree ####
phy <- Trees[[MeanTree[[1]]]]
parampgls <- list('massAvg', 'Long', 'La', 'Dmed','Rmax', 'Rmin', 'DeltaR',
                  'angeqv', 'Pry_Area_m', 'P_eq', 'Vsup_C')
parampgls <- lapply(parampgls, function(x){rlang::sym(x)})
# combination of the parameters for the regressions
comb <- combn(parampgls, 2, simplify = FALSE)
for (p1p2 in comb) {
  p1p2
  P1 <- p1p2[[1]]
  P2 <- p1p2[[2]]
  # combine phylo with the dataset and ensure structure and ordering
  CompData = comparative.data(phy, 
                              data = Data_Animal_de_pd[,c("massAvg", 
                                                          'Long',
                                                          'Long_m',
                                                          'La',
                                                          'La_m',
                                                          'Rmax',
                                                          'Rmax_m',
                                                          'Rmin',
                                                          'Rmin_m',
                                                          'DeltaR',
                                                          'Dmed',
                                                          'angeqv',
                                                          'Pry_Area_m',
                                                          'P_eq',
                                                          'Vsup_C',
                                                          'order',
                                                          'Sci_name2')],
                              names.col = 'Sci_name2', 
                              na.omit = TRUE)
  # P2 <- rlang::sym(p1p2)
  trai <- try(assign(paste(P2,'.PGLSav.',P1, sep = ''),
                   getpgls(!!P2, !!P1, CompData)), 
            silent = TRUE)
  print(p1p2)
  print(class(trai))
}

# test, run one only
#test2 <- getpgls(Dmed, Long, CompData)

#### PGLS by order ####
for (ord in unlist(GRorders[[1]])) {
  pps <- list('La', 'Dmed','angeqv', 'P_eq', 'Vsup_C')
  ord_data = sym(paste(ord[[1]],'_data', sep = ''))
  rlang::inject(!!ord_data <-
                  Data_Animal_de_pd[Data_Animal_de_pd$order == ord[[1]],])
  comp_ord_data = sym(paste(ord[[1]],'_CompData', sep = ''))
  assign(paste(ord[[1]],'_CompData', sep = ''),
         comparative.data(phy, 
                          data = rlang::inject(!!ord_data)[,c("massAvg", 
                                                              'Long',
                                                              'Long_m',
                                                              'La',
                                                              'La_m',
                                                              'Rmax',
                                                              'Rmax_m',
                                                              'Rmin',
                                                              'Rmin_m',
                                                              'DeltaR',
                                                              'Dmed',
                                                              'angeqv',
                                                              'Pry_Area_m',
                                                              'P_eq',
                                                              'Vsup_C',
                                                              'order',
                                                              'Sci_name2')],
                          names.col = 'Sci_name2', 
                          na.omit = TRUE))
  P1 <- sym('massAvg')
  for (pp in pps) {
    P2 <- sym(pp)
    if (length(unlist(rlang::inject(!!ord_data)['order'])) > 2) {
      print('siiii')
      trai <- assign(paste(ord[[1]],'.PGLSav.',pp, sep = ''),
                     rlang::inject(getpgls(!!P2, !!P1, !!comp_ord_data)))
    }
  }
}

#### Plot the PGLS ####
# Plot function for minor axis in logarithmic scale
minor.ticks.axis <- function(ax,n,t.ratio=0.5, mt, mn, mx,...){
  # major ticks
  iis <- sapply(mt, function(i)log10(i))
  print(iis)
  labels <- sapply(iis,function(i)
    as.expression(bquote(10^ .(i)))
  )
  axis(ax,at = mt,labels = labels,...)
  
  # minor ticks
  for (i in 1:(length(mt) - 1)) {
    minors <- seq(mt[i],mt[i + 1], length = n + 1)
    minors <- minors[-c(1,n + 1)]
    axis(ax, at = minors, tcl = par("tcl")*t.ratio, labels = FALSE)
  }
}

#### Plot DMED ####
plot(CompData$data$massAvg, CompData$data$Dmed, log = 'xy',  frame = FALSE)
atx <- c(0.01, 0.1, 1, 10, 100, 1000, 10000)
aty <- c(0.1, 1, 10, 100)
xlimt <- c(0.01, 10000)
ylimt <- c(0.1, 150)
plot(CompData$data$massAvg, 
     CompData$data$Dmed, 
     xlab = 'mass', 
     ylab = 'Dmed', 
     log  = 'xy', 
     asp  = 1, 
     xaxt = "n", 
     yaxt = "n",
     col  = rgb(0.01, 0.01, 0.4, 1),
     xlim = xlimt,
     ylim = ylimt,
     pch  = 22,
     lwd  = 1.5,
     frame = T,
)
text(CompData$data$massAvg, CompData$data$Dmed, 
     labels = row.names(CompData$data), 
     cex = 0.2)
minor.ticks.axis(1,9, mt = atx, mn = min(atx), mx = max(atx))
minor.ticks.axis(2,9, mt = aty, mn = min(aty), mx = max(aty))
abline(Dmed.PGLSav.massAvg)

for (ord in unlist(GRorders[[1]])) {
  print(ord)
  comp_ord_data = sym(paste(ord[[1]],'_CompData', sep = ''))
  ord_data = sym(paste(ord[[1]],'_data', sep = ''))
  pgls_ord = sym(paste(ord[[1]],'.PGLSav.','Dmed', sep = ''))
  plot(inject(!!comp_ord_data)$data$massAvg,
       inject(!!comp_ord_data)$data$Dmed,
       xlab = 'mass',
       ylab = 'Dmed',
       log  = 'xy',
       asp  = 1,
       xaxt = "n",
       yaxt = "n",
       col  = colsord[ord[[1]]],
       xlim = xlimt,
       ylim = ylimt,
       pch  = 22,
       lwd  = 1.5,
       frame = T,
  )
  # if (length(unlist(rlang::inject(!!ord_data)['order'])) > 2) {
  #   x00 = min(unlist(rlang::inject(!!ord_data)['massAvg']))
  #   x01 = max(unlist(rlang::inject(!!ord_data)['massAvg']))
  #   inty = rlang::inject(!!pgls_ord)[["model"]][["coef"]][["(Intercept)"]]
  #   slpy = rlang::inject(!!pgls_ord)[["model"]][["coef"]][["log10(massAvg)"]]
  #   y00 = 10^inty * x00 ^ slpy  
  #   y01 = 10^inty * x01 ^ slpy  
  #   segments(x0 = x00, y0 = y00, x1 = x01, y1 = y01, col = colsord[ord[[1]]])
  # }
  par(new = T)
}
abline(Dmed.PGLSav.massAvg, col = rgb(0.01, 0.01, 0.4, 1))
minor.ticks.axis(1,9, mt = atx, mn = min(atx), mx = max(atx))
minor.ticks.axis(2,9, mt = aty, mn = min(aty), mx = max(aty))
par(new = F)

#### Plot LA ####
plot(CompData$data$massAvg, CompData$data$La, log = 'xy',  frame = FALSE)
atx <- c(0.01, 0.1, 1, 10, 100, 1000, 10000)
aty <- c(1, 10, 100)
xlimt <- c(0.01, 10000)
ylimt <- c(1, 250)
plot(CompData$data$massAvg, 
     CompData$data$La, 
     xlab = 'mass', 
     ylab = 'La', 
     log  = 'xy', 
     asp  = 1, 
     xaxt = "n", 
     yaxt = "n",
     col  = rgb(0.01, 0.01, 0.4, 1),
     xlim = xlimt,
     ylim = ylimt,
     pch  = 22,
     lwd  = 1.5,
     frame = T,
)
text(CompData$data$massAvg, CompData$data$La, 
     labels = row.names(CompData$data), 
     cex = 0.2)
minor.ticks.axis(1,9, mt = atx, mn = min(atx), mx = max(atx))
minor.ticks.axis(2,9, mt = aty, mn = min(aty), mx = max(aty))
abline(La.PGLSav.massAvg)

for (ord in unlist(GRorders[[1]])) {
  print(ord)
  comp_ord_data = sym(paste(ord[[1]],'_CompData', sep = ''))
  ord_data = sym(paste(ord[[1]],'_data', sep = ''))
  pgls_ord = sym(paste(ord[[1]],'.PGLSav.','La', sep = ''))
  plot(inject(!!comp_ord_data)$data$massAvg,
       inject(!!comp_ord_data)$data$La,
       xlab = 'mass',
       ylab = 'La',
       log  = 'xy',
       asp  = 1,
       xaxt = "n",
       yaxt = "n",
       col  = colsord[ord[[1]]],
       xlim = xlimt,
       ylim = ylimt,
       pch  = 22,
       lwd  = 1.5,
       frame = T,
  )
  # if (length(unlist(rlang::inject(!!ord_data)['order'])) > 2) {
  #   x00 = min(unlist(rlang::inject(!!ord_data)['massAvg']))
  #   x01 = max(unlist(rlang::inject(!!ord_data)['massAvg']))
  #   inty = rlang::inject(!!pgls_ord)[["model"]][["coef"]][["(Intercept)"]]
  #   slpy = rlang::inject(!!pgls_ord)[["model"]][["coef"]][["log10(massAvg)"]]
  #   y00 = 10^inty * x00 ^ slpy  
  #   y01 = 10^inty * x01 ^ slpy  
  #   segments(x0 = x00, y0 = y00, x1 = x01, y1 = y01, col = colsord[ord[[1]]])
  # }
  par(new = T)
}
abline(La.PGLSav.massAvg, col = rgb(0.01, 0.01, 0.4, 1))
minor.ticks.axis(1,9, mt = atx, mn = min(atx), mx = max(atx))
minor.ticks.axis(2,9, mt = aty, mn = min(aty), mx = max(aty))
par(new = F)



#### Plot angeqv ####
plot(CompData$data$massAvg, CompData$data$angeqv, log = 'xy',  frame = FALSE)
atx <- c(0.01, 0.1, 1, 10, 100, 1000, 10000)
aty <- c(0.1, 1)
xlimt <- c(0.01, 10000)
ylimt <- c(0.08, 1.2)
plot(CompData$data$massAvg, 
     CompData$data$angeqv, 
     xlab = 'mass', 
     ylab = 'angeqv', 
     log  = 'xy', 
     asp  = 1, 
     xaxt = "n", 
     yaxt = "n",
     col  = rgb(0.01, 0.01, 0.4, 1),
     xlim = xlimt,
     ylim = ylimt,
     pch  = 22,
     lwd  = 1.5,
     frame = T,
)
text(CompData$data$massAvg, CompData$data$angeqv, 
     labels = row.names(CompData$data), 
     cex = 0.2)
minor.ticks.axis(1,9, mt = atx, mn = min(atx), mx = max(atx))
minor.ticks.axis(2,9, mt = aty, mn = min(aty), mx = max(aty))
abline(angeqv.PGLSav.massAvg)

for (ord in unlist(GRorders[[1]])) {
  print(ord)
  comp_ord_data = sym(paste(ord[[1]],'_CompData', sep = ''))
  ord_data = sym(paste(ord[[1]],'_data', sep = ''))
  pgls_ord = sym(paste(ord[[1]],'.PGLSav.','angeqv', sep = ''))
  plot(inject(!!comp_ord_data)$data$massAvg,
       inject(!!comp_ord_data)$data$angeqv,
       xlab = 'mass',
       ylab = 'angeqv',
       log  = 'xy',
       asp  = 1,
       xaxt = "n",
       yaxt = "n",
       col  = colsord[ord[[1]]],
       xlim = xlimt,
       ylim = ylimt,
       pch  = 22,
       lwd  = 1.5,
       frame = T,
  )
  # if (length(unlist(rlang::inject(!!ord_data)['order'])) > 2) {
  #   x00 = min(unlist(rlang::inject(!!ord_data)['massAvg']))
  #   x01 = max(unlist(rlang::inject(!!ord_data)['massAvg']))
  #   inty = rlang::inject(!!pgls_ord)[["model"]][["coef"]][["(Intercept)"]]
  #   slpy = rlang::inject(!!pgls_ord)[["model"]][["coef"]][["log10(massAvg)"]]
  #   y00 = 10^inty * x00 ^ slpy  
  #   y01 = 10^inty * x01 ^ slpy  
  #   segments(x0 = x00, y0 = y00, x1 = x01, y1 = y01, col = colsord[ord[[1]]])
  # }
  par(new = T)
}
abline(angeqv.PGLSav.massAvg, col = rgb(0.01, 0.01, 0.4, 1))
minor.ticks.axis(1,9, mt = atx, mn = min(atx), mx = max(atx))
minor.ticks.axis(2,9, mt = aty, mn = min(aty), mx = max(aty))
par(new = F)

#### plot P_eq ####
plot(CompData$data$massAvg, CompData$data$P_eq, log = 'xy',  frame = FALSE)
atx <- c(0.01, 0.1, 1, 10, 100, 1000, 10000)
aty <- c(1e6, 1e7)
xlimt <- c(0.01, 10000)
ylimt <- c(4e5, 2e7)
plot(CompData$data$massAvg, 
     CompData$data$P_eq, 
     xlab = 'mass', 
     ylab = 'P_eq', 
     log  = 'xy', 
     asp  = 1, 
     xaxt = "n", 
     yaxt = "n",
     col  = rgb(0.01, 0.01, 0.4, 1),
     xlim = xlimt,
     ylim = ylimt,
     pch  = 22,
     lwd  = 1.5,
     frame = T,
)
text(CompData$data$massAvg, CompData$data$P_eq, 
     labels = row.names(CompData$data), 
     cex = 0.2)
minor.ticks.axis(1,9, mt = atx, mn = min(atx), mx = max(atx))
minor.ticks.axis(2,9, mt = aty, mn = min(aty), mx = max(aty))
abline(P_eq.PGLSav.massAvg)

for (ord in unlist(GRorders[[1]])) {
  print(ord)
  comp_ord_data = sym(paste(ord[[1]],'_CompData', sep = ''))
  ord_data = sym(paste(ord[[1]],'_data', sep = ''))
  pgls_ord = sym(paste(ord[[1]],'.PGLSav.','P_eq', sep = ''))
  plot(inject(!!comp_ord_data)$data$massAvg,
       inject(!!comp_ord_data)$data$P_eq,
       xlab = 'mass',
       ylab = 'P_eq',
       log  = 'xy',
       asp  = 1,
       xaxt = "n",
       yaxt = "n",
       col  = colsord[ord[[1]]],
       xlim = xlimt,
       ylim = ylimt,
       pch  = 22,
       lwd  = 1.5,
       frame = T,
  )
  # if (length(unlist(rlang::inject(!!ord_data)['order'])) > 2) {
  #   x00 = min(unlist(rlang::inject(!!ord_data)['massAvg']))
  #   x01 = max(unlist(rlang::inject(!!ord_data)['massAvg']))
  #   inty = rlang::inject(!!pgls_ord)[["model"]][["coef"]][["(Intercept)"]]
  #   slpy = rlang::inject(!!pgls_ord)[["model"]][["coef"]][["log10(massAvg)"]]
  #   y00 = 10^inty * x00 ^ slpy  
  #   y01 = 10^inty * x01 ^ slpy  
  #   segments(x0 = x00, y0 = y00, x1 = x01, y1 = y01, col = colsord[ord[[1]]])
  # }
  par(new = T)
}
abline(P_eq.PGLSav.massAvg, col = rgb(0.01, 0.01, 0.4, 1))
minor.ticks.axis(1,9, mt = atx, mn = min(atx), mx = max(atx))
minor.ticks.axis(2,9, mt = aty, mn = min(aty), mx = max(aty))
par(new = F)

#### plot Vsup_C ####
plot(CompData$data$massAvg, CompData$data$Vsup_C, log = 'xy',  frame = FALSE)
atx <- c(0.01, 0.1, 1, 10, 100, 1000, 10000)
aty <- c(1, 10)
xlimt <- c(0.01, 10000)
ylimt <- c(0.5, 18)
plot(CompData$data$massAvg, 
     CompData$data$Vsup_C, 
     xlab = 'mass', 
     ylab = 'Vsup_C', 
     log  = 'xy', 
     asp  = 1, 
     xaxt = "n", 
     yaxt = "n",
     col  = rgb(0.01, 0.01, 0.4, 1),
     xlim = xlimt,
     ylim = ylimt,
     pch  = 22,
     lwd  = 1.5,
     frame = T,
)
text(CompData$data$massAvg, CompData$data$Vsup_C, 
     labels = row.names(CompData$data), 
     cex = 0.2)
minor.ticks.axis(1,9, mt = atx, mn = min(atx), mx = max(atx))
minor.ticks.axis(2,9, mt = aty, mn = min(aty), mx = max(aty))
abline(Vsup_C.PGLSav.massAvg)

for (ord in unlist(GRorders[[1]])) {
  print(ord)
  comp_ord_data = sym(paste(ord[[1]],'_CompData', sep = ''))
  ord_data = sym(paste(ord[[1]],'_data', sep = ''))
  pgls_ord = sym(paste(ord[[1]],'.PGLSav.','Vsup_C', sep = ''))
  plot(inject(!!comp_ord_data)$data$massAvg,
       inject(!!comp_ord_data)$data$Vsup_C,
       xlab = 'mass',
       ylab = 'Vsup_C',
       log  = 'xy',
       asp  = 1,
       xaxt = "n",
       yaxt = "n",
       col  = colsord[ord[[1]]],
       xlim = xlimt,
       ylim = ylimt,
       pch  = 22,
       lwd  = 1.5,
       frame = T,
  )
  # if (length(unlist(rlang::inject(!!ord_data)['order'])) > 2) {
  #   x00 = min(unlist(rlang::inject(!!ord_data)['massAvg']))
  #   x01 = max(unlist(rlang::inject(!!ord_data)['massAvg']))
  #   inty = rlang::inject(!!pgls_ord)[["model"]][["coef"]][["(Intercept)"]]
  #   slpy = rlang::inject(!!pgls_ord)[["model"]][["coef"]][["log10(massAvg)"]]
  #   y00 = 10^inty * x00 ^ slpy  
  #   y01 = 10^inty * x01 ^ slpy  
  #   segments(x0 = x00, y0 = y00, x1 = x01, y1 = y01, col = colsord[ord[[1]]])
  # }
  par(new = T)
}
abline(Vsup_C.PGLSav.massAvg, col = rgb(0.01, 0.01, 0.4, 1))
minor.ticks.axis(1,9, mt = atx, mn = min(atx), mx = max(atx))
minor.ticks.axis(2,9, mt = aty, mn = min(aty), mx = max(aty))
par(new = F)

################################################################################
###############             TREE VISUALIZATION              ####################
#### Initial visualization of the tree ####
plot(phy, type = "fan", cex = .5)
plotTree(phy,ftype = "i",fsize = 0.5,lwd = 1)

#### PLOT TREE DEFINITIVE ####
# functions from Revell's blog 
# http://blog.phytools.org/2022/07/a-simple-function-to-plot-radial.html
fanCladogram <- function(tree,use.edge.length=FALSE, ...) {
  if (use.edge.length == FALSE) {
    if (hasArg(power)) power <- list(...)$power
    else power <- 0.6
    tree <- compute.brlen(tree,power = power)
  }
  args <- list(...)
  args$power <- NULL
  args$tree <- tree
  args$type <- "fan"
  args$color <- "transparent"
  do.call(plotTree,args)
  obj <- get("last_plot.phylo",envir = .PlotPhyloEnv)
  X <- cbind(obj$xx[1:Ntip(tree)],
           obj$yy[1:Ntip(tree)])
  rownames(X) <- tree$tip.label
  A <- cbind(obj$xx[1:tree$Nnode + Ntip(tree)],
           obj$yy[1:tree$Nnode + Ntip(tree)])
  rownames(A) <- 1:tree$Nnode + Ntip(tree)
  pp.args <- list(...)
  pp.args$power <- NULL
  pp.args$tree <- tree
  pp.args$X <- X
  pp.args$A <- A
  pp.args$ftype <- "off"
  pp.args$node.size <- c(0,0)
  pp.args$xlim <- obj$x.lim
  pp.args$ylim <- obj$y.lim
  pp.args$xaxt <- "n"
  pp.args$type <- NULL
  pp.args$add <- TRUE
  do.call(phylomorphospace,pp.args)
}
# second part from his answer
# simplified to add bars
fanCladoBar <- function(tree,...) {
  obj <- get("last_plot.phylo",envir = .PlotPhyloEnv)
  X <- cbind(obj$xx[1:Ntip(tree)],
           obj$yy[1:Ntip(tree)])
  rownames(X) <- tree$tip.label
  A <- cbind(obj$xx[1:tree$Nnode + Ntip(tree)],
           obj$yy[1:tree$Nnode + Ntip(tree)])
  rownames(A) <- 1:tree$Nnode + Ntip(tree)
  pp.args <- list(...)
  pp.args$power <- NULL
  pp.args$tree <- tree
  pp.args$X <- X
  pp.args$A <- A
  pp.args$ftype <- "off"
  pp.args$node.size <- c(0,0)
  pp.args$xlim <- obj$x.lim
  pp.args$ylim <- obj$y.lim
  pp.args$xaxt <- "n"
  pp.args$type <- NULL
  pp.args$add <- TRUE
  do.call(phylomorphospace,pp.args)
}

# function to find order a family of each taxa
findspc <- function(valG, 
                    tips, 
                    dataF = Data_Animal_de_pd, 
                    spcname = 'Sci_name2', 
                    val){
  listgr = list()
  for (or in valG) {
    orl = sym(paste(or,'l',sep = ''))
    rlang::inject(!!orl <- list())
    for (spc in tips) {
      valG_spc <- dataF[dataF[spcname] == spc, val]
      if (valG_spc == or) {
        rlang::inject(!!orl <- append(!!orl,spc))
      }
    }
    rlang::inject(listgr[[or]] <- !!orl)
  }
  listspc = list()
  for (spc in tips) {
    valG_spc <- dataF[dataF[spcname] == spc, val]
    listspc <- append(listspc, valG_spc)
  }
  return(list('Gr' = listgr, 'spc' = listspc))
}

# Our tree
ffile = paste('/home/kale/Documents/PaperAllometry/images/Rplots/tree/')
H = 5
W = H*2
svg(filename = paste(ffile,'tree_bars.svg',sep = ''),
    width = W,
    height = H,
    pointsize = 8,
    family = 'serif')

phy.dat <- CompData$phy
phy.dat <- ladderize(phy.dat) # for better visualization
phy.dat <- compute.brlen(phy.dat, power = 0.65) # get rid of branch lengths 

linewidth = 3 # tree width
barwidth = 0.02 # bar width
scaleb = 0.05 # scale of the bars

# data for the heat map
Dmed.dat <- log10(CompData$data$Dmed)
names(Dmed.dat) <- CompData$phy$tip.label
# heat map (plot off)
MapDmed <- contMap(phy.dat,
                   Dmed.dat,
                   plot = F,
                   res = 600,
                   ftype = "i",
                   outline = T)
n <- length(MapDmed$cols)
MapDmed <- setMap(MapDmed,viridis(n))

# data for the arcs
tipis <- MapDmed$tree$tip.label

GRorderslists <- findspc(valG = as.list(GRorders[['order']]),
                         tips = tipis,
                         val = 'order')

GRofamilialists <- findspc(valG = as.list(GRfamilia[['family']]),
                           tips = tipis,
                           val = 'family')

# to color the bars according to faimly
# barcol = list()
# for (spec in 1:length(tipis)) {
#   barcol <- append(barcol,colsfamy[GRofamilialists$spc[[spec]]])
# }
# to color the bars according to order
barcol = list()
for (spec in 1:length(tipis)) {
  barcol <- append(barcol,colsord[GRorderslists$spc[[spec]]])
}
# assign name to colors
barcol <- as.character(barcol)
barcol <- setNames(barcol,tipis)

# bar sizes according to La
La.dat <- log10(CompData$data$La)
names(La.dat) <- CompData$phy$tip.label

# plot the bars, tree color transparent
plotTree.wBars(MapDmed$tree,
               La.dat,
               scale = scaleb,
               type = "fan",
               color = "transparent",
               col = barcol, 
               part = 0.5, 
               outline = T, 
               border = 'white', 
               width = barwidth)

# plot arcs
# para orden
for (i in 1:length(GRorderslists$Gr)) {
  name = names(GRorderslists$Gr)[i]
  tips = as.character(GRorderslists$Gr[[i]])
  if (length(GRorderslists$Gr[[i]]) > 1) {
    nodes = findMRCA(MapDmed$tree,tips)}
  else{
    nodes =  which(MapDmed$tree$tip.label == tips)
  }
  arc.cladelabels(tree = MapDmed$tree,
                  text = name,
                  node = nodes,
                  ln.offset = 1.16,
                  lab.offset = 1.18,
                  mark.node = FALSE,
                  # orientation = 'horizontal'
                  cex = 0.2,
                  col = colsord[[name]],
                  lwd = 2,
                  stretch = 1
  )
}
# # para familia
# for(i in 1:length(GRofamilialists$Gr)){
#   name = names(GRofamilialists$Gr)[i]
#   tips = as.character(GRofamilialists$Gr[[i]])
#   if (length(GRofamilialists$Gr[[i]])>1){
#     nodes = findMRCA(MapDmed$tree,tips)}
#   else{
#     nodes =  which(MapDmed$tree$tip.label==tips)
#   }
#   if (length(nodes)>0){
#     arc.cladelabels(tree = MapDmed$tree,
#                     text = name,
#                     node = nodes,
#                     ln.offset=1.14,
#                     lab.offset=1.16,
#                     mark.node=FALSE,
#                     orientation = 'horizontal',
#                     cex = 0.5,
#     )
#   }
# }

# plot the tree
fanCladoBar(MapDmed$tree, 
            lwd = linewidth,
            colors = MapDmed$cols,
            part = 0.5,
            fsize = 0.8,
            add = TRUE)

# add the color bar of dmed
add.color.bar(0.75,
              cols = MapDmed$cols,
              lims = MapDmed$lims,
              digits = 2,
              x = par()$usr[2] - 2.5,
              y = par()$usr[3] + 0.5,
              prompt = FALSE,
              title = "log10(Dmed)",
              subtitle = "", 
              cex = 0.8, 
              outline = FALSE,
              lwd = linewidth
              )

# based on Revell's blog: to plot the scales of the bar(min, mean, max)
scale <- scaleb
w <- barwidth
# min
xx <- 0.95*par()$usr[2]
yy <- 0.95*par()$usr[4]
# colorb = barcoln[GRofamilialists$spc[[which.min(La.dat)]]]
# colorb = barcoln[GRorderslists$spc[[which.min(La.dat)]]]
colorb = 'grey'
polygon(c(xx,xx - scale*min(La.dat), xx - scale*min(La.dat), xx),
        c(yy - w/2, yy - w/2, yy + w/2, yy + w/2),
        col = colorb, 
        border = F)
text(xx - scale*min(La.dat),yy,
     paste('log10(La) = ',
           round(min(La.dat),1)),
     pos = 2,
     cex = 0.8)
# mean
xx <- 0.95*par()$usr[2]
yy <- 0.9*par()$usr[4]
newladat <- abs(mean(La.dat) - La.dat)
# colorb = barcoln[GRofamilialists$spc[[which.min(newladat)]]]
# colorb = barcoln[GRorderslists$spc[[which.min(La.dat)]]]
colorb = 'grey'
polygon(c(xx,xx - scale*mean(La.dat), xx - scale*mean(La.dat),xx),
        c(yy - w/2, yy - w/2, yy + w/2, yy + w/2),
        col = colorb, 
        border = F)
text(xx - scale*mean(La.dat),yy,
     paste('log10(La) = ',
           round(mean(La.dat),1)),
     pos = 2,
     cex = 0.8)
# max
xx <- 0.95*par()$usr[2]
yy <- 0.85*par()$usr[4]
# colorb = barcoln[GRofamilialists$spc[[which.max(La.dat)]]]
# colorb = barcoln[GRorderslists$spc[[which.max(La.dat)]]]
colorb = 'grey'
polygon(c(xx, xx - scale*max(La.dat), xx - scale*max(La.dat), xx),
        c(yy - w/2, yy - w/2, yy + w/2, yy + w/2),
        col = colorb, 
        border = F)
text(xx - scale*max(La.dat),yy,
     paste('log10(La) = ',
           round(max(La.dat),1)),
     pos = 2,
     cex = 0.8)

# family colors
# legend("topleft",
#        names(barcoln),
#        pch = 22,
#        pt.bg = barcoln,
#        pt.cex = 0.8,
#        cex = 0.8,
#        bty = "n", 
#        x.intersp = 0.2, 
#        ncol = 2, 
#        y.intersp = 0.9)

# order colors
legend("bottomleft",
       names(colsarc),
       pch = 22,
       pt.bg = colsarc,
       pt.cex = 0.8,
       cex = 0.8,
       bty = "o", 
       x.intersp = 0.2, 
       y.intersp = 0.01,
       text.width	= 0.25,
       horiz = T)

dev.off()

# to check the names of each tip
H = 4
W = H*2
svg(filename = paste(ffile,'tree_names.svg',sep = ''),
    width = W,
    height = H,
    pointsize = 8,
    family = 'serif')
plot.new()
fanCladogram(phy.dat,lwd = 6,part = 0.5,fsize = 0.8)
fanCladogram(MapDmed$tree,lwd = 4,colors = MapDmed$cols,part = 0.5,
                      fsize = 0.8,add = TRUE)
dev.off()
