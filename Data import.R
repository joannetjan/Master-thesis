rm(list = ls())
## In this file, the important packages are installed and the observed wind gusts at the Hamburg Weather Mast and the covariates from the COSMO-REA6 data
## are analyzed.
# Install necessary packages ------------------------------------------------------------------------------------------------------------------------
library("NLP")
library("caret")
library("plyr")
library("readr")
library("ggplot2")
library("GGally")
library("dplyr")
library("mlbench")
library("zoo")
library("xts")
library("forecast")
library("plotly")
library("reshape2")
library("naniar")
library("ggcorrplot")
library("HelpersMG")
library("ismev")
library("evd")
library("evir")
library("quantmod")
library("scoringRules")
library("ensembleMOS")
library("reliaR")
library("gsl")

# Import data files ------------------------------------------------------------------------------------------------------------------------
CV <- read.delim("C:/Users/sky_t/Documents/Joanne/Master scriptie/Data/C_raw.txt", header = FALSE)
CV[1,] = c("Date", CV[1,-17])
CV <- data.frame(CV)
colnames(CV) <- CV[1,]
CV <- CV[-1,]

Z10 <- read.delim("C:/Users/sky_t/Documents/Joanne/Master scriptie/Data/FB010_M60_200401010000-201412312300.txt", header = FALSE)
Z50 <- read.delim("C:/Users/sky_t/Documents/Joanne/Master scriptie/Data/FB050_M60_200401010000-201412312300.txt", header = FALSE)
Z110 <- read.delim("C:/Users/sky_t/Documents/Joanne/Master scriptie/Data/FB110_M60_200401010000-201412312300.txt", header = FALSE)
Z175 <- read.delim("C:/Users/sky_t/Documents/Joanne/Master scriptie/Data/FB175_M60_200401010000-201412312300.txt", header = FALSE)
Z250 <- read.delim("C:/Users/sky_t/Documents/Joanne/Master scriptie/Data/FB250_M60_200401010000-201412312300.txt", header = FALSE)

T10 <- as.numeric(quantile(Z10[Z10!=99999], 0.5))
T50 <- as.numeric(quantile(Z50[Z50!=99999], 0.5))
T110 <- as.numeric(quantile(Z110[Z110!=99999], 0.5))
T175 <- as.numeric(quantile(Z175[Z175!=99999], 0.5))
T250 <- as.numeric(quantile(Z250[Z250!=99999], 0.5))


CV$Date <- seq(from = as.POSIXct("2004-01-01 00:00"), to = as.POSIXct("2014-12-31 23:00"), by = "hour")
# Wind$Date <- seq(from = as.POSIXct("2004-01-01 00:00"), to = as.POSIXct("2014-12-31 23:00"), by = "hour")

Wind <- cbind(CV$Date, Z10, Z50, Z110, Z175, Z250)
colnames(Wind) <- c("Date", 10, 50, 110, 175, 250)
CWind <- Wind 

# Censoring
CWind$`10`[CWind$`10` < T10] <- T10
CWind$`50`[CWind$`50` < T50] <- T50
CWind$`110`[CWind$`110` < T110] <- T110
CWind$`175`[CWind$`175` < T175] <- T175
CWind$`250`[CWind$`250` < T250] <- T250
threshold <- c(T10, T50, T110, T175, T250)
# Initialise variables ------------------------------------------------------------------------------------------------------------------------
N <- dim(Z10)[1]
L <- dim(CV)[2]-1
K <- 3
Z <- c(10, 50, 110, 175, 250)
NZ <- (Z - min(Z))/(max(Z)-min(Z))

F10 <- N - sum(Z10==99999)
F50 <- N - sum(Z50==99999)
F110 <- N - sum(Z110==99999)
F175 <- N - sum(Z175==99999)
F250 <- N - sum(Z250==99999)

# Data analysis ------------------------------------------------------------------------------------------------------------------------
# CV[,-1] <- sapply(CV[,-1], is.numeric)
for (i in seq(L+1)){
  if (i != 1){
    CV[,i] <- as.numeric(CV[,i])
  }
}

df <- cbind(Z10, as.numeric(CV$VMAX_10M))
df <- df[complete.cases(df),]
df <- df[df$V1 != 99999 ,]
dif <- df$V1-df$`as.numeric(CV$VMAX_10M)`
mean(dif)
sd(dif)
hist(dif, breaks = 150)

W2 <- data.frame(Wind)
z <- 2
subW <- W2[,c(1,z)]
subW <- subW[subW[,2] != 99999,]
subW$month <- format(subW$Date, format="%b")
subW$hour <- format(subW$Date, format="%H")
subW$sort_orderM <- as.POSIXlt(subW$Date)$mon
subW$sort_orderH <- as.POSIXlt(subW$Date)$hour

ggplot(subW, aes(x=subW[,1], y = subW[,2])) + geom_line() + xlab("") + ylab("Observed gust at 250M")
ggplot(subW) + geom_boxplot(aes(x=reorder(month, sort_orderM), y=subW[,2])) +xlab("Month") + ylab("Gust at z = 250M" )
ggplot(subW) + geom_boxplot(aes(x=reorder(hour, sort_orderH), y=subW[,2]))+xlab("Hour") + ylab("Gust at z = 250M")


ggplot(subW, aes(x=subW[,1], y = subW[,2])) + geom_line() + xlab("") + ylab("Observed gust at 250M")
ggcorrplot(cor(CV[,(2:17)]), method = "circle")
hist(dif, breaks = 150)

Wind[Wind == 99999] <- NA
CWind[CWind == 99999] <- NA

View(summary(CV))
View(summary(Wind))
sd(Wind$`250`[!is.na(Wind$`250`)])
CV[,2:17] <- scale(CV[,2:17])
CV <- CV[complete.cases(CV),]


ggcorrplot(cor(CV[,(2:17)]), method = "circle")


