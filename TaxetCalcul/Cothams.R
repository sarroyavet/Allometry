library(smatr) # Standardized Major Axis Regression (SMA ou RMA ou Model II)
library(ggplot2)
library(scales)
library(Hmisc)

# data from Coatham paper
Data2_csv = "/home/kale/OwnCloud/2022/Allom/Iyytemp2sldwrks.csv"

Data2_2 <- read.csv(Data2_csv, header = TRUE, sep = ",", row.names = 1)
summary(Data2_2)

# Plot function for minor axis in logarithmic scale
minor.ticks.axis <- function(ax,n,t.ratio=0.5,mt,mn,mx,...){
  # major ticks
  iis <- sapply(mt, function(i)log10(i))
  print(iis)
  labels <- sapply(iis,function(i)
    as.expression(bquote(10^ .(i)))
  )
  axis(ax,at=mt,labels=labels,...)
  
  # minor ticks
  for (i in 1:(length(mt)-1)){
    minors <- seq(mt[i],mt[i+1], length = n+1)
    minors <- minors[-c(1,n+1)]
    axis(ax,at=minors,tcl=par("tcl")*t.ratio,labels=FALSE)
  }
}


# SMA to calculate the regression of the percentage of the moment of inertia
Per_SMA <- sma(percentage~Body.mass, 
               log = "xy",
               data = Data2_2,
               slope.test = 0.08,
               robust = T
               )
summary(Per_SMA)
plot(Per_SMA)

plot(Per_SMA, 
     xlab= 'log(mass)', 
     ylab = 'log(per) ', 
     log = 'xy', 
     asp = 1, 
     xaxt="n",
     yaxt="n",
     col = rgb(0.01, 0.01, 0.4, 1),
     xlim=c(0.1,1000), 
     ylim=c(0.1, 1.0),
     pch=22,
     lwd=1.5,
)

#text(Data2_2$Body.mass, Data2_2$percentage, labels=row.names(Data2_2))
atx <- c(0.1, 1, 10, 100, 1000)
aty <- c(0.1, 1.0)
minor.ticks.axis(1,9, mt=atx,mn=min(atx), mx = max(atx))
minor.ticks.axis(2,9, mt=aty,mn=min(aty), mx = max(aty))

# Percentage CI
per_mean <- mean(Data2_2$percentage)
per_n <- length(Data2_2$percentage)
per_sd <- sd(Data2_2$percentage)
per_se <- per_sd/sqrt(per_n)
alpha <- 0.05
degrees.freedom <- per_n - 1
per_t.score <- qt(p=alpha/2, df=degrees.freedom,lower.tail=F)
per_error <- per_t.score * per_se
per_lower.bound <- per_mean - per_error
per_upper.bound <- per_mean + per_error



# SMA to calculate the regression of the moment of inertia of the lower arm
Low_SMA <- sma(lower~Body.mass, 
               log = "xy",
               data = Data2_2, 
               slope.test = 0.08, 
               robust = T
               )
summary(Low_SMA)
plot(Low_SMA)

plot(Low_SMA, 
     xlab= 'log(mass)', 
     ylab = 'log(I) ', 
     log = 'xy', 
     asp = 1, 
     xaxt="n", 
     yaxt="n",
     col = rgb(0.01, 0.01, 0.4, 1),
     xlim=c(0.1, 1000), 
     ylim=c(1e-7, 8),
     pch=22,
     lwd=1.5,
     )
text(Data2_2$Body.mass, Data2_2$lower, labels=row.names(Data2_2))
atx <- c(0.1, 1, 10, 100, 1000)
aty <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1)
minor.ticks.axis(1,9, mt=atx,mn=min(atx), mx = max(atx))
minor.ticks.axis(2,9, mt=aty,mn=min(aty), mx = max(aty))


