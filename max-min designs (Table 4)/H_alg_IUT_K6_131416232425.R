#####################################################################################
###### H-algorithm for min-max design for general Bipartite graph for K=6 for IUT  ##
###### Pairs = (1,3),(1,4), (1,6), (2,3),(2,4), (2,5)                              ##
#####################################################################################



#install.packages("nloptr")
#install.packages("mvtnorm")


library('nloptr')
library(mvtnorm)
rm(list = ls())

# Function for calculating the mean and covariance matrix of Z-ij's and criterion 'Psi'.
psi<- function(design,theta,n,sigma,q)
{
  
  mu13 <- (n*(theta[1]-theta[3]))/(sigma*sqrt((1/design[1])+(1/design[3])))
  mu14 <- (n*(theta[1]-theta[4]))/(sigma*sqrt((1/design[1])+(1/design[4])))
  mu16 <- (n*(theta[1]-theta[6]))/(sigma*sqrt((1/design[1])+(1/design[6])))
  mu23 <- (n*(theta[2]-theta[3]))/(sigma*sqrt((1/design[2])+(1/design[3])))
  mu24 <- (n*(theta[2]-theta[4]))/(sigma*sqrt((1/design[2])+(1/design[4])))
  mu25 <- (n*(theta[2]-theta[5]))/(sigma*sqrt((1/design[2])+(1/design[5])))
  
  p13<- pnorm(q,mu13)-pnorm(-q,mu13)
  p14<- pnorm(q,mu14)-pnorm(-q,mu14)
  p16<- pnorm(q,mu16)-pnorm(-q,mu16)
  p23<- pnorm(q,mu23)-pnorm(-q,mu23)
  p24<- pnorm(q,mu24)-pnorm(-q,mu24)
  p25<- pnorm(q,mu25)-pnorm(-q,mu25)
  
  
  lw1 <- c(-q,-q)
  up1 <- c(q,q)
  
  rho1314 <- sqrt((design[3]*design[4])/((design[1]+design[3])*(design[1]+design[4])))
  s1314<- rbind(c(1,rho1314),c(rho1314,1))
  mu1314<- c(mu13,mu14)
  p1314<- pmvnorm(mean=mu1314, sigma=s1314, lower=lw1, upper=up1)
  
  rho1316 <- sqrt((design[3]*design[6])/((design[1]+design[3])*(design[1]+design[6])))
  s1316<- rbind(c(1,rho1316),c(rho1316,1))
  mu1316<- c(mu13,mu16)
  p1316<- pmvnorm(mean=mu1316, sigma=s1316, lower=lw1, upper=up1)
  
  rho1323 <- sqrt((design[1]*design[2])/((design[1]+design[3])*(design[2]+design[3])))
  s1323<- rbind(c(1,rho1323),c(rho1323,1))
  mu1323<- c(mu13,mu23)
  p1323<- pmvnorm(mean=mu1323, sigma=s1323, lower=lw1, upper=up1)
  
  rho1324 <- 0
  s1324<- rbind(c(1,rho1324),c(rho1324,1))
  mu1324<- c(mu13,mu24)
  p1324<- pmvnorm(mean=mu1324, sigma=s1324, lower=lw1, upper=up1)
  
  rho1325 <- 0
  s1325<- rbind(c(1,rho1325),c(rho1325,1))
  mu1325<- c(mu13,mu25)
  p1325<-  pmvnorm(mean=mu1325, sigma=s1325, lower=lw1, upper=up1)
  
  rho1416 <- sqrt((design[4]*design[6])/((design[1]+design[4])*(design[1]+design[6])))
  s1416<- rbind(c(1,rho1416),c(rho1416,1))
  mu1416<- c(mu14,mu16)
  p1416<- pmvnorm(mean=mu1416, sigma=s1416, lower=lw1, upper=up1)
  
  rho1423 <- 0
  s1423<- rbind(c(1,rho1423),c(rho1423,1))
  mu1423<- c(mu14,mu23)
  p1423<-  pmvnorm(mean=mu1423, sigma=s1423, lower=lw1, upper=up1)
  
  
  rho1424 <- sqrt((design[1]*design[2])/((design[1]+design[4])*(design[2]+design[4])))
  s1424<- rbind(c(1,rho1424),c(rho1424,1))
  mu1424<- c(mu14,mu24)
  p1424<- pmvnorm(mean=mu1424, sigma=s1424, lower=lw1, upper=up1)
  
  rho1425 <- 0
  s1425<- rbind(c(1,rho1425),c(rho1425,1))
  mu1425<- c(mu14,mu25)
  p1425<-  pmvnorm(mean=mu1425, sigma=s1425, lower=lw1, upper=up1)
  
  rho1623 <- 0
  s1623<- rbind(c(1,rho1623),c(rho1623,1))
  mu1623<- c(mu16,mu23)
  p1623<-  pmvnorm(mean=mu1623, sigma=s1623, lower=lw1, upper=up1)
  
  
  rho1624 <- 0
  s1624<- rbind(c(1,rho1624),c(rho1624,1))
  mu1624<- c(mu16,mu24)
  p1624<- pmvnorm(mean=mu1624, sigma=s1624, lower=lw1, upper=up1)
  
  rho1625 <- 0
  s1625<- rbind(c(1,rho1625),c(rho1625,1))
  mu1625<- c(mu16,mu25)
  p1625<-  pmvnorm(mean=mu1625, sigma=s1625, lower=lw1, upper=up1)
  
  
  rho2324 <- sqrt((design[3]*design[4])/((design[2]+design[3])*(design[2]+design[4])))
  s2324<- rbind(c(1,rho2324),c(rho2324,1))
  mu2324<- c(mu23,mu24)
  p2324<- pmvnorm(mean=mu2324, sigma=s2324, lower=lw1, upper=up1)
  
  
  rho2325 <- sqrt((design[3]*design[5])/((design[2]+design[3])*(design[2]+design[5])))
  s2325<- rbind(c(1,rho2325),c(rho2325,1))
  mu2325<- c(mu23,mu25)
  p2325<- pmvnorm(mean=mu2325, sigma=s2325, lower=lw1, upper=up1)
  
  rho2425 <- sqrt((design[4]*design[5])/((design[2]+design[4])*(design[2]+design[5])))
  s2425<- rbind(c(1,rho2425),c(rho2425,1))
  mu2425<- c(mu24,mu25)
  p2425<- pmvnorm(mean=mu2425, sigma=s2425, lower=lw1, upper=up1)
  
  
  
  
  
  lw2 <- c(-q,-q,-q)
  up2 <- c(q,q,q)
  
  s131416<- rbind(c(1,rho1314,rho1316),c(rho1314,1,rho1416),c(rho1316,rho1416,1))
  mu131416<- c(mu13,mu14,mu16)
  p131416<-  pmvnorm(mean=mu131416, sigma=s131416, lower=lw2, upper=up2)
  
  s131423<- rbind(c(1,rho1314,rho1323),c(rho1314,1,rho1423),c(rho1323,rho1423,1))
  mu131423<- c(mu13,mu14,mu23)
  p131423<-  pmvnorm(mean=mu131423, sigma=s131423, lower=lw2, upper=up2)
  
  s131424<- rbind(c(1,rho1314,rho1324),c(rho1314,1,rho1424),c(rho1324,rho1424,1))
  mu131424<- c(mu13,mu14,mu24)
  p131424<-  pmvnorm(mean=mu131424, sigma=s131424, lower=lw2, upper=up2)
  
  s131425<- rbind(c(1,rho1314,rho1325),c(rho1314,1,rho1425),c(rho1325,rho1425,1))
  mu131425<- c(mu13,mu14,mu25)
  p131425<-  pmvnorm(mean=mu131425, sigma=s131425, lower=lw2, upper=up2)
  
  s131623<- rbind(c(1,rho1316,rho1323),c(rho1316,1,rho1623),c(rho1323,rho1623,1))
  mu131623<- c(mu13,mu16,mu23)
  p131623<-  pmvnorm(mean=mu131623, sigma=s131623, lower=lw2, upper=up2)
  
  s131624<- rbind(c(1,rho1316,rho1324),c(rho1316,1,rho1624),c(rho1324,rho1624,1))
  mu131624<- c(mu13,mu16,mu24)
  p131624<-  pmvnorm(mean=mu131624, sigma=s131624, lower=lw2, upper=up2)
  
  s131625<- rbind(c(1,rho1316,rho1325),c(rho1316,1,rho1625),c(rho1325,rho1625,1))
  mu131625<- c(mu13,mu16,mu25)
  p131625<-  pmvnorm(mean=mu131625, sigma=s131625, lower=lw2, upper=up2)
  
  s132324<- rbind(c(1,rho1323,rho1324),c(rho1323,1,rho2324),c(rho1324,rho2324,1))
  mu132324<- c(mu13,mu23,mu24)
  p132324<-  pmvnorm(mean=mu132324, sigma=s132324, lower=lw2, upper=up2)
  
  s132325<- rbind(c(1,rho1323,rho1325),c(rho1323,1,rho2325),c(rho1325,rho2325,1))
  mu132325<- c(mu13,mu23,mu25)
  p132325<-  pmvnorm(mean=mu132325, sigma=s132325, lower=lw2, upper=up2)
  
  s132425<- rbind(c(1,rho1324,rho1325),c(rho1324,1,rho2425),c(rho1325,rho2425,1))
  mu132425<- c(mu13,mu24,mu25)
  p132425<-  pmvnorm(mean=mu132425, sigma=s132425, lower=lw2, upper=up2)
  
  
  s141623<- rbind(c(1,rho1416,rho1423),c(rho1416,1,rho1623),c(rho1423,rho1623,1))
  mu141623<- c(mu14,mu16,mu23)
  p141623<-  pmvnorm(mean=mu141623, sigma=s141623, lower=lw2, upper=up2)
  
  s141624<- rbind(c(1,rho1416,rho1424),c(rho1416,1,rho1624),c(rho1424,rho1624,1))
  mu141624<- c(mu14,mu16,mu24)
  p141624<-  pmvnorm(mean=mu141624, sigma=s141624, lower=lw2, upper=up2)
  
  s141625<- rbind(c(1,rho1416,rho1425),c(rho1416,1,rho1625),c(rho1425,rho1625,1))
  mu141625<- c(mu14,mu16,mu25)
  p141625<-  pmvnorm(mean=mu141625, sigma=s141625, lower=lw2, upper=up2)
  
  s142324<- rbind(c(1,rho1423,rho1424),c(rho1423,1,rho2324),c(rho1424,rho2324,1))
  mu142324<- c(mu14,mu23,mu24)
  p142324<-  pmvnorm(mean=mu142324, sigma=s142324, lower=lw2, upper=up2)
  
  s142325<- rbind(c(1,rho1423,rho1425),c(rho1423,1,rho2325),c(rho1425,rho2325,1))
  mu142325<- c(mu14,mu23,mu25)
  p142325<-  pmvnorm(mean=mu142325, sigma=s142325, lower=lw2, upper=up2)
  
  s142425<- rbind(c(1,rho1424,rho1425),c(rho1424,1,rho2425),c(rho1425,rho2425,1))
  mu142425<- c(mu14,mu24,mu25)
  p142425<-  pmvnorm(mean=mu142425, sigma=s142425, lower=lw2, upper=up2)
  
  s162324<- rbind(c(1,rho1623,rho1624),c(rho1623,1,rho2324),c(rho1624,rho2324,1))
  mu162324<- c(mu16,mu23,mu24)
  p162324<-  pmvnorm(mean=mu162324, sigma=s162324, lower=lw2, upper=up2)
  
  s162325<- rbind(c(1,rho1623,rho1625),c(rho1623,1,rho2325),c(rho1625,rho2325,1))
  mu162325<- c(mu16,mu23,mu25)
  p162325<-  pmvnorm(mean=mu162325, sigma=s162325, lower=lw2, upper=up2)
  
  s162425<- rbind(c(1,rho1624,rho1625),c(rho1624,1,rho2425),c(rho1625,rho2425,1))
  mu162425<- c(mu16,mu24,mu25)
  p162425<-  pmvnorm(mean=mu162425, sigma=s162425, lower=lw2, upper=up2)
  
  s232425<- rbind(c(1,rho2324,rho2325),c(rho2324,1,rho2425),c(rho2325,rho2425,1))
  mu232425<- c(mu23,mu24,mu25)
  p232425<-  pmvnorm(mean=mu232425, sigma=s232425, lower=lw2, upper=up2)
  
  lw3 <- c(-q,-q,-q,-q)
  up3 <- c(q,q,q,q)
  
  s13141623<- rbind(c(1,rho1314,rho1316,rho1323),c(rho1314,1,rho1416,rho1423),c(rho1316,rho1416,1,rho1623),c(rho1323,rho1423,rho1623,1))
  mu13141623<- c(mu13,mu14,mu16,mu23)
  p13141623<-  pmvnorm(mean=mu13141623, sigma=s13141623, lower=lw3, upper=up3)
  
  s13141624<- rbind(c(1,rho1314,rho1316,rho1324),c(rho1314,1,rho1416,rho1424),c(rho1316,rho1416,1,rho1624),c(rho1324,rho1424,rho1624,1))
  mu13141624<- c(mu13,mu14,mu16,mu24)
  p13141624<-  pmvnorm(mean=mu13141624, sigma=s13141624, lower=lw3, upper=up3)
  
  s13141625<- rbind(c(1,rho1314,rho1316,rho1325),c(rho1314,1,rho1416,rho1425),c(rho1316,rho1416,1,rho1625),c(rho1325,rho1425,rho1625,1))
  mu13141625<- c(mu13,mu14,mu16,mu25)
  p13141625<-  pmvnorm(mean=mu13141625, sigma=s13141625, lower=lw3, upper=up3)
  
  s13142324<- rbind(c(1,rho1314,rho1323,rho1324),c(rho1314,1,rho1423,rho1424),c(rho1323,rho1423,1,rho2324),c(rho1324,rho1424,rho2324,1))
  mu13142324<- c(mu13,mu14,mu23,mu24)
  p13142324<-  pmvnorm(mean=mu13142324, sigma=s13142324, lower=lw3, upper=up3)
  
  s13142325<- rbind(c(1,rho1314,rho1323,rho1325),c(rho1314,1,rho1423,rho1425),c(rho1323,rho1423,1,rho2325),c(rho1325,rho1425,rho2325,1))
  mu13142325<- c(mu13,mu14,mu23,mu25)
  p13142325<-  pmvnorm(mean=mu13142325, sigma=s13142325, lower=lw3, upper=up3)
  
  s13142425<- rbind(c(1,rho1314,rho1324,rho1325),c(rho1314,1,rho1424,rho1425),c(rho1324,rho1424,1,rho2425),c(rho1325,rho1425,rho2425,1))
  mu13142425<- c(mu13,mu14,mu24,mu25)
  p13142425<-  pmvnorm(mean=mu13142425, sigma=s13142425, lower=lw3, upper=up3)
  
  s13162324<- rbind(c(1,rho1316,rho1323,rho1324),c(rho1316,1,rho1623,rho1624),c(rho1323,rho1623,1,rho2324),c(rho1324,rho1624,rho2324,1))
  mu13162324<- c(mu13,mu16,mu23,mu24)
  p13162324<-  pmvnorm(mean=mu13162324, sigma=s13162324, lower=lw3, upper=up3)
  
  s13162325<- rbind(c(1,rho1316,rho1323,rho1325),c(rho1316,1,rho1623,rho1625),c(rho1323,rho1623,1,rho2325),c(rho1325,rho1625,rho2325,1))
  mu13162325<- c(mu13,mu16,mu23,mu25)
  p13162325<-  pmvnorm(mean=mu13162325, sigma=s13162325, lower=lw3, upper=up3)
  
  s13162425<- rbind(c(1,rho1316,rho1324,rho1325),c(rho1316,1,rho1624,rho1625),c(rho1324,rho1624,1,rho2425),c(rho1325,rho1625,rho2425,1))
  mu13162425<- c(mu13,mu16,mu24,mu25)
  p13162425<-  pmvnorm(mean=mu13162425, sigma=s13162425, lower=lw3, upper=up3)
  
  s13232425<- rbind(c(1,rho1323,rho1324,rho1325),c(rho1323,1,rho2324,rho2325),c(rho1324,rho2324,1,rho2425),c(rho1325,rho2325,rho2425,1))
  mu13232425<- c(mu13,mu23,mu24,mu25)
  p13232425<-  pmvnorm(mean=mu13232425, sigma=s13232425, lower=lw3, upper=up3)
  
  s14162324<- rbind(c(1,rho1416,rho1423,rho1424),c(rho1416,1,rho1623,rho1624),c(rho1423,rho1623,1,rho2324),c(rho1424,rho1624,rho2324,1))
  mu14162324<- c(mu14,mu16,mu23,mu24)
  p14162324<-  pmvnorm(mean=mu14162324, sigma=s14162324, lower=lw3, upper=up3)
  
  s14162325<- rbind(c(1,rho1416,rho1423,rho1425),c(rho1416,1,rho1623,rho1625),c(rho1423,rho1623,1,rho2325),c(rho1425,rho1625,rho2325,1))
  mu14162325<- c(mu14,mu16,mu23,mu25)
  p14162325<-  pmvnorm(mean=mu14162325, sigma=s14162325, lower=lw3, upper=up3)
  
  s14162425<- rbind(c(1,rho1416,rho1424,rho1425),c(rho1416,1,rho1624,rho1625),c(rho1424,rho1624,1,rho2425),c(rho1425,rho1625,rho2425,1))
  mu14162425<- c(mu14,mu16,mu24,mu25)
  p14162425<-  pmvnorm(mean=mu14162425, sigma=s14162425, lower=lw3, upper=up3)
  
  s14232425<- rbind(c(1,rho1423,rho1424,rho1425),c(rho1423,1,rho2324,rho2325),c(rho1424,rho2324,1,rho2425),c(rho1425,rho2325,rho2425,1))
  mu14232425<- c(mu14,mu23,mu24,mu25)
  p14232425<-  pmvnorm(mean=mu14232425, sigma=s14232425, lower=lw3, upper=up3)
  
  s16232425<- rbind(c(1,rho1623,rho1624,rho1625),c(rho1623,1,rho2324,rho2325),c(rho1624,rho2324,1,rho2425),c(rho1625,rho2325,rho2425,1))
  mu16232425<- c(mu16,mu23,mu24,mu25)
  p16232425<-  pmvnorm(mean=mu16232425, sigma=s16232425, lower=lw3, upper=up3)
  
  lw3 <- c(-q,-q,-q,-q,-q)
  up3 <- c(q,q,q,q,q)
  
  
  s1314162324<- rbind(c(1,rho1314,rho1316,rho1323,rho1324),c(rho1314,1,rho1416,rho1423,rho1424),c(rho1316,rho1416,1,rho1623,rho1624),c(rho1323,rho1423,rho1623,1,rho2324),c(rho1324,rho1424,rho1624,rho2324,1))
  mu1314162324<- c(mu13,mu14,mu16,mu23,mu24)
  p1314162324<-  pmvnorm(mean=mu1314162324, sigma=s1314162324, lower=lw3, upper=up3)
  
  s1314162325<- rbind(c(1,rho1314,rho1316,rho1323,rho1325),c(rho1314,1,rho1416,rho1423,rho1425),c(rho1316,rho1416,1,rho1623,rho1625),c(rho1323,rho1423,rho1623,1,rho2325),c(rho1325,rho1425,rho1625,rho2325,1))
  mu1314162325<- c(mu13,mu14,mu16,mu23,mu25)
  p11314162325<-  pmvnorm(mean=mu1314162325, sigma=s1314162325, lower=lw3, upper=up3)
  
  s1314162425<- rbind(c(1,rho1314,rho1316,rho1324,rho1325),c(rho1314,1,rho1416,rho1424,rho1425),c(rho1316,rho1416,1,rho1624,rho1625),c(rho1324,rho1424,rho1624,1,rho2425),c(rho1325,rho1425,rho1625,rho2425,1))
  mu1314162425<- c(mu13,mu14,mu16,mu24,mu25)
  p1314162425<-  pmvnorm(mean=mu1314162425, sigma=s1314162425, lower=lw3, upper=up3)
  
  
  s1314232425<- rbind(c(1,rho1314,rho1323,rho1324,rho1325),c(rho1314,1,rho1423,rho1424,rho1425),c(rho1323,rho1423,1,rho2324,rho2325),c(rho1324,rho1424,rho2324,1,rho2425),c(rho1325,rho1425,rho2325,rho2425,1))
  mu1314232425<- c(mu13,mu14,mu23,mu24,mu25)
  p1314232425<-  pmvnorm(mean=mu1314232425, sigma=s1314232425, lower=lw3, upper=up3)
  
  s1316232425<- rbind(c(1,rho1316,rho1323,rho1324,rho1325),c(rho1316,1,rho1623,rho1624,rho1625),c(rho1323,rho1623,1,rho2324,rho2325),c(rho1324,rho1624,rho2324,1,rho2425),c(rho1325,rho1625,rho2325,rho2425,1))
  mu1316232425<- c(mu13,mu16,mu23,mu24,mu25)
  p1316232425<-  pmvnorm(mean=mu1316232425, sigma=s1316232425, lower=lw3, upper=up3)
  
  s1416232425<- rbind(c(1,rho1416,rho1423,rho1424,rho1425),c(rho1416,1,rho1623,rho1624,rho1625),c(rho1423,rho1623,1,rho2324,rho2325),c(rho1424,rho1624,rho2324,1,rho2425),c(rho1425,rho1625,rho2325,rho2425,1))
  mu1416232425<- c(mu14,mu16,mu23,mu24,mu25)
  p1416232425<-  pmvnorm(mean=mu1416232425, sigma=s1416232425, lower=lw3, upper=up3)
  
  
  lw4 <- c(-q,-q,-q,-q,-q,-q)
  up4 <- c(q,q,q,q,q,q)
  
  
  s131416232425<- rbind(c(1,rho1314,rho1316,rho1323,rho1324,rho1325),c(rho1314,1,rho1416,rho1423,rho1424,rho1425),c(rho1316,rho1416,1,rho1623,rho1624,rho1625),c(rho1323,rho1423,rho1623,1,rho2324,rho2325),c(rho1324,rho1424,rho1624,rho2324,1,rho2425),c(rho1325,rho1425,rho1625,rho2325,rho2425,1))
  mu131416232425<- c(mu13,mu14,mu16,mu23,mu24,mu25)
  p131416232425<-  pmvnorm(mean=mu131416232425, sigma=s131416232425, lower=lw4, upper=up4)
  
  
  
  # Power of the test  
  p1<- p13+p14+p16+p23+p24+p25
  p2<- p1314+p1316+p1323+p1324+p1325+p1416+p1423+p1424+p1425+p1623+p1624+p1625+p2324+p2325+p2425
  p3<- p131416+p131423+p131424+p131425+p131623+p131624+p131625+p132324+p132325+p132425+p141623+p141624+p141625+p142324+p142325+p142425+p162324+p162325+p162425+p232425
  p4<- p13141623+p13141624+p13141625+p13142324+p13142325+p13142425+p13162324+p13162325+p13162425+p13232425+p14162324+p14162325+p14162425+p14232425+p16232425
  p5<-  p1314162324+p11314162325+p1314162425+p1314232425+p1316232425+p1416232425
  p6<- p131416232425
  p<- -(1-p1+p2-p3+p4-p5+p6)    # Power function of IUT
  return(p[1])
  
}

# Function for calculating the average power B

B<- function(design,theta_vals,prob_vals,n,sigma,q){
  
  s <- 0
  for (i in 1:length(prob_vals)){
    s <- s + prob_vals[i]* psi(design,theta_vals[i,],n,sigma,q)   # Average power associated to LFCs and probabilities
  }
  ret_val<-s
  return(ret_val[1])
}

# Function for compute optimum on the avarage design

get_optimal_on_the_average_design<- function(theta_vals,prob_vals,n,sigma,q)
{
  des_0 <-c(0.1,0.2,0.1,0.2,0.3,0.05)
  #options = optimoptions('fmincon','Display','none')
  const1 <- function(design,theta_vals,prob_vals,n,sigma,q) {
    return( design[1]+design[2]+design[3]+design[4]+design[5]+design[6]-1)
  }
  lb <- c(0,0,0,0,0,0)
  ub <- c(0.5,0.5,0.5,0.5,0.5,0.5)
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-16 )
  opts <- list( "algorithm"= "NLOPT_GN_ISRES",
                "xtol_rel"= 1.0e-16,
                "maxeval"= 20000,
                "local_opts" = local_opts,
                "print_level" = 0 )
  
  
  
  res <- nloptr ( x0 = des_0,
                  theta_vals = theta_vals,
                  prob_vals = prob_vals,
                  n=n,
                  sigma=sigma,
                  q=q,
                  eval_f = B,
                  lb = lb,
                  ub = ub,
                  eval_g_eq = const1,
                  #eval_g_ineq = eval_g0,
                  opts = opts
  )

  return(res$solution)
  
}

#Step 1

step_1<- function(H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,n,sigma,q){

  psi_min <- 100
  #stopping=true
  
  this_psi_neg<- rep(0, length(prob_vals_l))
  this_psi<- rep(0, length(prob_vals_l))
  
  for (i in 1:length(prob_vals_l)){
    this_psi_neg[i] <-  psi(design_l,theta_vals_l[i,],n,sigma,q)
    this_psi[i] <- -this_psi_neg[i]
  }
  
  
  this_psi<- as.numeric(round(this_psi, 3))
  
  
  this_psi_min<- min(this_psi)
  print(this_psi)
  index <- which(this_psi == this_psi_min)
  
  # Stopping Criterion
  
  if(this_psi_min < psi_min)
  {
    psi_min <- this_psi_min
  }
  
  iteration_num
  iteration_num <- iteration_num+1
  psi_min
  B_pi_l
  #if (psi_min < -B_pi_l)
  print("abs(psi_min+B_pi_l)")
  print(abs(psi_min+B_pi_l))
  if((psi_min < -B_pi_l)&(abs(psi_min+B_pi_l)>0.004)){
    stopping <- 0}
  else{
    stopping<- 1}
  
  
  if (stopping==1){
    print("Found minmax design")
    design_l
    print("Least Favourable Distribution")
    print("probs")
    prob_vals_l
    print("thetas")
    theta_vals_l}
  else{
    #theta_l = theta_argmin_psi
    step_2(index,H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,n,sigma,q)
  }
}

#Step 2

step_2<- function(index,H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,n,sigma,q){
  #delta_l = unit mass to theta_l
  step_3(index,H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,n,sigma,q)
}

#Step 3

step_3<- function(index,H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,n,sigma,q)
{
  design_l1 <- design_l
  smallest_B_pi_t_l1 <- -100
  prob_vals_l1 <- prob_vals_l
  theta_vals_l1 <- theta_vals_l
  
  theta_vals_l
  prob_vals_l
  prob_vals_t_l1<- matrix(0, length(H), length(prob_vals_l))
  delta_l <- rep(0, length(prob_vals_l))
  for(i in 1:length(index)){
    delta_l[index[i]] <- 1/length(index)
  }
  print(delta_l)
  for (t_val in 1:length(H)){
    #new T priors
    prob_vals_t_l1[t_val,] =  ((1-H[t_val])*prob_vals_l)+(H[t_val]*delta_l)
  }
  
  # Assign l+1 values to l.
  
  for(i in 1:length(H)){
    design_t_l1 <- get_optimal_on_the_average_design(theta_vals_l, prob_vals_t_l1[i,],n,sigma,q)
    B_pi_t_l1 <- B(design_t_l1, theta_vals_l,prob_vals_t_l1[i,],n,sigma,q)  
    if (-smallest_B_pi_t_l1 > -B_pi_t_l1){
      smallest_B_pi_t_l1 <- B_pi_t_l1
      design_l1 <- design_t_l1
      prob_vals_l1 <- prob_vals_t_l1[i,]
      theta_vals_l1 <- theta_vals_l
    }
  }
  print("prob_vals_l1")
  print(prob_vals_l1)
  print("design_l1")
  print(design_l1)
  step_4(H,h_grid_space,theta_vals_l1,prob_vals_l1,design_l1,B_pi_l,n,sigma,q)
}

# Step 4

step_4<- function(H,h_grid_space,theta_vals_l1,prob_vals_l1,design_l1,B_pi_l,n,sigma,q){
  # STEP 3
  theta_vals_l <- theta_vals_l1
  prob_vals_l <- prob_vals_l1
  design_l <- design_l1
  B_pi_l1 <- B(design_l, theta_vals_l, prob_vals_l,n,sigma,q)
  print(B_pi_l1-B_pi_l)
  
  # Assigning B(l+1) values to B(l)  
  
  if (-B_pi_l1 < -B_pi_l){
    B_pi_l <- B_pi_l1
    print("B_pi_l")
    print(B_pi_l)
    #"Assigned l1 to l"
    
    # return
    step_1(H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,n,sigma,q)
  }
  else{
    h_grid_space <- h_grid_space/2
    print(h_grid_space)
    H <-seq(0, 1, by=h_grid_space)
    
    this_psi_neg_1<- rep(0, length(prob_vals_l))
    this_psi_1<- rep(0, length(prob_vals_l))
    
    for (i in 1:length(prob_vals_l)){
      this_psi_neg_1[i] <-  psi(design_l,theta_vals_l[i,],n,sigma,q)
      this_psi_1[i] <- -this_psi_neg_1[i]
    }
    
    
    this_psi_1<- as.numeric(round(this_psi_1, 3))
    
    this_psi_min_1<- min(this_psi_1)
    index_1 <- which(this_psi_1 == this_psi_min_1)
    
    step_3(index_1,H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,n,sigma,q)
  }
}

N <- 250                    # Total No. of subjects
n <- sqrt(N)
sigma <- 1                  # Standard deviation
al <- 0.05               # significance level alpha
q <- qnorm(1-al)           # critical value
del<-0.6                   # delta
K<-5
iteration_num<-0
B_pi_l <- 100             # Initialization of B_pi
h_grid_space <- 0.5       # grid length
H<- seq(0, 1, by=h_grid_space)      # grid vector


print( "Initial Prior on Theta : ")

# Theta values (LFC's) are their corresponding prior probabilities



the_1<- c(-del,del,0,0,0,0)
the_2<-  c(del,-del,0,0,0,0)
the_3<-  c(-del,-del,0,0,0,0)
the_4<- c(0,0,-del,del,del,del)
the_5<- c(0,0,del,-del,del,del)
the_6<- c(0,0,del,del,-del,del)
the_7<- c(0,0,del,del,del,-del)
the_8<- c(0,0,-del,-del,del,del)
the_9<- c(0,0,-del,del,-del,del)
the_10<- c(0,0,-del,del,del,-del)
the_11<- c(0,0,-del,-del,-del,-del)



theta_vals_l<- rbind(the_1,the_2,the_3,the_4,the_5,the_6,the_7,the_8,the_9,the_10,the_11)

prob_vals_l <- rep(1/11,11)

# Starting optimum on the average design and corresponding power

print( "Starting Design")
design_l <- get_optimal_on_the_average_design(theta_vals_l,prob_vals_l,n,sigma,q)
B_pi_l <- B(design_l,theta_vals_l, prob_vals_l,n,sigma,q)

step_1(H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,n,sigma,q)
