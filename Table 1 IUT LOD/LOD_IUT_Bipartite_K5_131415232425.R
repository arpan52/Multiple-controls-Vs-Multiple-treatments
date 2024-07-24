

#install.packages("nloptr")
library(nloptr)
library(mvtnorm)

IUT_K5_131415232425 <- function(x,mu,n,sigma,q){
  
  
  mu13 <- (n*(mu[1]-mu[3]))/(sigma*sqrt((1/x[1])+(1/x[3])))
  mu14 <- (n*(mu[1]-mu[4]))/(sigma*sqrt((1/x[1])+(1/x[4])))
  mu15 <- (n*(mu[1]-mu[5]))/(sigma*sqrt((1/x[1])+(1/x[5])))
  mu23 <- (n*(mu[2]-mu[3]))/(sigma*sqrt((1/x[2])+(1/x[3])))
  mu24 <- (n*(mu[2]-mu[4]))/(sigma*sqrt((1/x[2])+(1/x[4])))
  mu25 <- (n*(mu[2]-mu[5]))/(sigma*sqrt((1/x[2])+(1/x[5])))
  
  #print(c(mu13,mu14,mu15,mu23,mu24,mu25))
  
  
  p13<- pnorm(q,mu13)-pnorm(-q,mu13)
  p14<- pnorm(q,mu14)-pnorm(-q,mu14)
  p15<- pnorm(q,mu15)-pnorm(-q,mu15)
  p23<- pnorm(q,mu23)-pnorm(-q,mu23)
  p24<-pnorm(q,mu24)-pnorm(-q,mu24)
  p25<-pnorm(q,mu25)-pnorm(-q,mu25)
  
  
  
  lw1 <- c(-q,-q)
  up1 <- c(q,q)
  
  rho1314 <- sqrt((x[3]*x[4])/((x[1]+x[3])*(x[1]+x[4])))
  s1314<- rbind(c(1,rho1314),c(rho1314,1))
  mu1314<- c(mu13,mu14)
  p1314<-  pmvnorm(mean=mu1314, sigma=s1314, lower=lw1, upper=up1)
  #p1213<- p12*p13
  
  rho1315 <- sqrt((x[3]*x[5])/((x[1]+x[3])*(x[1]+x[5])))
  s1315<- rbind(c(1,rho1315),c(rho1315,1))
  mu1315<- c(mu13,mu15)
  p1315<-  pmvnorm(mean=mu1315, sigma=s1315, lower=lw1, upper=up1)
  
  rho1323 <- sqrt((x[1]*x[2])/((x[1]+x[3])*(x[2]+x[3])))
  s1323<- rbind(c(1,rho1323),c(rho1323,1))
  mu1323<- c(mu13,mu23)
  p1323<- pmvnorm(mean=mu1323, sigma=s1323, lower=lw1, upper=up1)
  # p1223<- p12*p23
  
  rho1324 <- 0
  s1324<- rbind(c(1,rho1324),c(rho1324,1))
  mu1324<- c(mu13,mu24)
  p1324<- pmvnorm(mean=mu1324, sigma=s1324, lower=lw1, upper=up1)
  # p1223<- p12*p23
  
  rho1325 <- 0
  s1325<- rbind(c(1,rho1325),c(rho1325,1))
  mu1325<- c(mu13,mu25)
  p1325<- pmvnorm(mean=mu1325, sigma=s1325, lower=lw1, upper=up1)
  
  rho1415 <- sqrt((x[4]*x[5])/((x[1]+x[4])*(x[1]+x[5])))
  s1415<- rbind(c(1,rho1415),c(rho1415,1))
  mu1415<- c(mu14,mu15)
  p1415<-  pmvnorm(mean=mu1415, sigma=s1415, lower=lw1, upper=up1)
  
  rho1423 <- 0
  s1423<- rbind(c(1,rho1423),c(rho1423,1))
  mu1423<- c(mu14,mu23)
  p1423<- pmvnorm(mean=mu1423, sigma=s1423, lower=lw1, upper=up1)
  #p1323<- p13*p23
  
  rho1424 <-  sqrt((x[1]*x[2])/((x[1]+x[4])*(x[2]+x[4])))
  s1424<- rbind(c(1,rho1424),c(rho1424,1))
  mu1424<- c(mu14,mu24)
  p1424<- pmvnorm(mean=mu1424, sigma=s1424, lower=lw1, upper=up1)
  #p1323<- p13*p23
  
  rho1425 <- 0
  s1425<- rbind(c(1,rho1425),c(rho1425,1))
  mu1425<- c(mu14,mu25)
  p1425<- pmvnorm(mean=mu1425, sigma=s1425, lower=lw1, upper=up1)
  
  rho1523 <- 0
  s1523<- rbind(c(1,rho1523),c(rho1523,1))
  mu1523<- c(mu15,mu23)
  p1523<- pmvnorm(mean=mu1523, sigma=s1523, lower=lw1, upper=up1)
  #p1323<- p13*p23
  
  rho1524 <-  0
  s1524<- rbind(c(1,rho1524),c(rho1524,1))
  mu1524<- c(mu15,mu24)
  p1524<- pmvnorm(mean=mu1524, sigma=s1524, lower=lw1, upper=up1)
  #p1323<- p13*p23
  
  rho1525 <- sqrt((x[1]*x[2])/((x[1]+x[5])*(x[2]+x[5])))
  s1525<- rbind(c(1,rho1525),c(rho1525,1))
  mu1525<- c(mu15,mu25)
  p1525<- pmvnorm(mean=mu1525, sigma=s1525, lower=lw1, upper=up1)
  
  rho2324 <- sqrt((x[3]*x[4])/((x[2]+x[3])*(x[2]+x[4])))
  s2324<- rbind(c(1,rho2324),c(rho2324,1))
  mu2324<- c(mu23,mu24)
  p2324<-  pmvnorm(mean=mu2324, sigma=s2324, lower=lw1, upper=up1)
  #p1213<- p12*p13
  
  rho2325 <- sqrt((x[3]*x[5])/((x[2]+x[3])*(x[2]+x[5])))
  s2325<- rbind(c(1,rho2325),c(rho2325,1))
  mu2325<- c(mu23,mu25)
  p2325<-  pmvnorm(mean=mu2325, sigma=s2325, lower=lw1, upper=up1)
  
  rho2425 <- sqrt((x[4]*x[5])/((x[2]+x[4])*(x[2]+x[5])))
  s2425<- rbind(c(1,rho2425),c(rho2425,1))
  mu2425<- c(mu24,mu25)
  p2425<-  pmvnorm(mean=mu2425, sigma=s2425, lower=lw1, upper=up1)
  
  lw2 <- c(-q,-q,-q)
  up2 <- c(q,q,q)
  
  s131415<- rbind(c(1,rho1314,rho1315),c(rho1314,1,rho1415),c(rho1315,rho1415,1))
  mu131415<- c(mu13,mu14,mu15)
  p131415<-  pmvnorm(mean=mu131415, sigma=s131415, lower=lw2, upper=up2)
  
  s131423<- rbind(c(1,rho1314,rho1323),c(rho1314,1,rho1423),c(rho1323,rho1423,1))
  mu131423<- c(mu13,mu14,mu23)
  p131423<-  pmvnorm(mean=mu131423, sigma=s131423, lower=lw2, upper=up2)
  
  s131424<- rbind(c(1,rho1314,rho1324),c(rho1314,1,rho1424),c(rho1324,rho1424,1))
  mu131424<- c(mu13,mu14,mu24)
  p131424<-  pmvnorm(mean=mu131424, sigma=s131424, lower=lw2, upper=up2)
  
  s131425<- rbind(c(1,rho1314,rho1325),c(rho1314,1,rho1425),c(rho1325,rho1425,1))
  mu131425<- c(mu13,mu14,mu25)
  p131425<-  pmvnorm(mean=mu131425, sigma=s131425, lower=lw2, upper=up2)
  
  s131523<- rbind(c(1,rho1315,rho1323),c(rho1315,1,rho1523),c(rho1323,rho1523,1))
  mu131523<- c(mu13,mu15,mu23)
  p131523<-  pmvnorm(mean=mu131523, sigma=s131523, lower=lw2, upper=up2)
  
  s131524<- rbind(c(1,rho1315,rho1324),c(rho1315,1,rho1524),c(rho1324,rho1524,1))
  mu131524<- c(mu13,mu15,mu24)
  p131524<-  pmvnorm(mean=mu131524, sigma=s131524, lower=lw2, upper=up2)
  
  s131525<- rbind(c(1,rho1315,rho1325),c(rho1315,1,rho1525),c(rho1325,rho1525,1))
  mu131525<- c(mu13,mu15,mu25)
  p131525<-  pmvnorm(mean=mu131525, sigma=s131525, lower=lw2, upper=up2)
  
  s132324<- rbind(c(1,rho1323,rho1324),c(rho1323,1,rho2324),c(rho1324,rho2324,1))
  mu132324<- c(mu13,mu23,mu24)
  p132324<-  pmvnorm(mean=mu132324, sigma=s132324, lower=lw2, upper=up2)
  
  s132325<- rbind(c(1,rho1323,rho1325),c(rho1323,1,rho2325),c(rho1325,rho2325,1))
  mu132325<- c(mu13,mu23,mu25)
  p132325<-  pmvnorm(mean=mu132325, sigma=s132325, lower=lw2, upper=up2)
  
  s132425<- rbind(c(1,rho1324,rho1325),c(rho1324,1,rho2425),c(rho1325,rho2425,1))
  mu132425<- c(mu13,mu24,mu25)
  p132425<-  pmvnorm(mean=mu132425, sigma=s132425, lower=lw2, upper=up2)
  
  s141523<- rbind(c(1,rho1415,rho1423),c(rho1415,1,rho1523),c(rho1423,rho1523,1))
  mu141523<- c(mu14,mu15,mu23)
  p141523<-  pmvnorm(mean=mu141523, sigma=s141523, lower=lw2, upper=up2)
  
  s141524<- rbind(c(1,rho1415,rho1424),c(rho1415,1,rho1524),c(rho1424,rho1524,1))
  mu141524<- c(mu14,mu15,mu24)
  p141524<-  pmvnorm(mean=mu141524, sigma=s141524, lower=lw2, upper=up2)
  
  s141525<- rbind(c(1,rho1415,rho1425),c(rho1415,1,rho1525),c(rho1425,rho1525,1))
  mu141525<- c(mu14,mu15,mu25)
  p141525<-  pmvnorm(mean=mu141525, sigma=s141525, lower=lw2, upper=up2)
  
  s142324<- rbind(c(1,rho1423,rho1424),c(rho1423,1,rho2324),c(rho1424,rho2324,1))
  mu142324<- c(mu14,mu23,mu24)
  p142324<-  pmvnorm(mean=mu142324, sigma=s142324, lower=lw2, upper=up2)
  
  s142325<- rbind(c(1,rho1423,rho1425),c(rho1423,1,rho2325),c(rho1425,rho2325,1))
  mu142325<- c(mu14,mu23,mu25)
  p142325<-  pmvnorm(mean=mu142325, sigma=s142325, lower=lw2, upper=up2)
  
  s142425<- rbind(c(1,rho1424,rho1425),c(rho1424,1,rho2425),c(rho1425,rho2425,1))
  mu142425<- c(mu14,mu24,mu25)
  p142425<-  pmvnorm(mean=mu142425, sigma=s142425, lower=lw2, upper=up2)
  
  s152324<- rbind(c(1,rho1523,rho1524),c(rho1523,1,rho2324),c(rho1524,rho2324,1))
  mu152324<- c(mu15,mu23,mu24)
  p152324<-  pmvnorm(mean=mu152324, sigma=s152324, lower=lw2, upper=up2)
  
  s152325<- rbind(c(1,rho1523,rho1525),c(rho1523,1,rho2325),c(rho1525,rho2325,1))
  mu152325<- c(mu15,mu23,mu25)
  p152325<-  pmvnorm(mean=mu152325, sigma=s152325, lower=lw2, upper=up2)
  
  s152425<- rbind(c(1,rho1524,rho1525),c(rho1524,1,rho2425),c(rho1525,rho2425,1))
  mu152425<- c(mu15,mu24,mu25)
  p152425<-  pmvnorm(mean=mu152425, sigma=s152425, lower=lw2, upper=up2)
  
  s232425<- rbind(c(1,rho2324,rho2325),c(rho2324,1,rho2425),c(rho2325,rho2425,1))
  mu232425<- c(mu23,mu24,mu25)
  p232425<-  pmvnorm(mean=mu232425, sigma=s232425, lower=lw2, upper=up2)
  
  lw3 <- c(-q,-q,-q,-q)
  up3 <- c(q,q,q,q)
  
  s13141523<- rbind(c(1,rho1314,rho1315,rho1323),c(rho1314,1,rho1415,rho1423),c(rho1315,rho1415,1,rho1523),c(rho1323,rho1423,rho1523,1))
  mu13141523<- c(mu13,mu14,mu15,mu23)
  p13141523<-  pmvnorm(mean=mu13141523, sigma=s13141523, lower=lw3, upper=up3)
  
  s13141524<- rbind(c(1,rho1314,rho1315,rho1324),c(rho1314,1,rho1415,rho1424),c(rho1315,rho1415,1,rho1524),c(rho1324,rho1424,rho1524,1))
  mu13141524<- c(mu13,mu14,mu15,mu24)
  p13141524<-  pmvnorm(mean=mu13141524, sigma=s13141524, lower=lw3, upper=up3)
  
  s13141525<- rbind(c(1,rho1314,rho1315,rho1325),c(rho1314,1,rho1415,rho1425),c(rho1315,rho1415,1,rho1525),c(rho1325,rho1425,rho1525,1))
  mu13141525<- c(mu13,mu14,mu15,mu25)
  p13141525<-  pmvnorm(mean=mu13141525, sigma=s13141525, lower=lw3, upper=up3)
  
  s13142324<- rbind(c(1,rho1314,rho1323,rho1324),c(rho1314,1,rho1423,rho1424),c(rho1323,rho1423,1,rho2324),c(rho1324,rho1424,rho2324,1))
  mu13142324<- c(mu13,mu14,mu23,mu24)
  p13142324<-  pmvnorm(mean=mu13142324, sigma=s13142324, lower=lw3, upper=up3)
  
  s13142325<- rbind(c(1,rho1314,rho1323,rho1325),c(rho1314,1,rho1423,rho1425),c(rho1323,rho1423,1,rho2325),c(rho1325,rho1425,rho2325,1))
  mu13142325<- c(mu13,mu14,mu23,mu25)
  p13142325<-  pmvnorm(mean=mu13142325, sigma=s13142325, lower=lw3, upper=up3)
  
  s13142425<- rbind(c(1,rho1314,rho1324,rho1325),c(rho1314,1,rho1424,rho1425),c(rho1324,rho1424,1,rho2425),c(rho1325,rho1425,rho2425,1))
  mu13142425<- c(mu13,mu14,mu24,mu25)
  p13142425<-  pmvnorm(mean=mu13142425, sigma=s13142425, lower=lw3, upper=up3)
  
  s13152324<- rbind(c(1,rho1315,rho1323,rho1324),c(rho1315,1,rho1523,rho1524),c(rho1323,rho1523,1,rho2324),c(rho1324,rho1524,rho2324,1))
  mu13152324<- c(mu13,mu15,mu23,mu24)
  p13152324<-  pmvnorm(mean=mu13152324, sigma=s13152324, lower=lw3, upper=up3)
  
  s13152325<- rbind(c(1,rho1315,rho1323,rho1325),c(rho1315,1,rho1523,rho1525),c(rho1323,rho1523,1,rho2325),c(rho1325,rho1525,rho2325,1))
  mu13152325<- c(mu13,mu15,mu23,mu25)
  p13152325<-  pmvnorm(mean=mu13152325, sigma=s13152325, lower=lw3, upper=up3)
  
  s13152425<- rbind(c(1,rho1315,rho1324,rho1325),c(rho1315,1,rho1524,rho1525),c(rho1324,rho1524,1,rho2425),c(rho1325,rho1525,rho2425,1))
  mu13152425<- c(mu13,mu15,mu24,mu25)
  p13152425<-  pmvnorm(mean=mu13152425, sigma=s13152425, lower=lw3, upper=up3)
  
  s13232425<- rbind(c(1,rho1323,rho1324,rho1325),c(rho1323,1,rho2324,rho2325),c(rho1324,rho2324,1,rho2425),c(rho1325,rho2325,rho2425,1))
  mu13232425<- c(mu13,mu23,mu24,mu25)
  p13232425<-  pmvnorm(mean=mu13232425, sigma=s13232425, lower=lw3, upper=up3)
  
  s14152324<- rbind(c(1,rho1415,rho1423,rho1424),c(rho1415,1,rho1523,rho1524),c(rho1423,rho1523,1,rho2324),c(rho1424,rho1524,rho2324,1))
  mu14152324<- c(mu14,mu15,mu23,mu24)
  p14152324<-  pmvnorm(mean=mu14152324, sigma=s14152324, lower=lw3, upper=up3)
  
  s14152325<- rbind(c(1,rho1415,rho1423,rho1425),c(rho1415,1,rho1523,rho1525),c(rho1423,rho1523,1,rho2325),c(rho1425,rho1525,rho2325,1))
  mu14152325<- c(mu14,mu15,mu23,mu25)
  p14152325<-  pmvnorm(mean=mu14152325, sigma=s14152325, lower=lw3, upper=up3)
  
  s14152425<- rbind(c(1,rho1415,rho1424,rho1425),c(rho1415,1,rho1524,rho1525),c(rho1424,rho1524,1,rho2425),c(rho1425,rho1525,rho2425,1))
  mu14152425<- c(mu14,mu15,mu24,mu25)
  p14152425<-  pmvnorm(mean=mu14152425, sigma=s14152425, lower=lw3, upper=up3)
  
  s14232425<- rbind(c(1,rho1423,rho1424,rho1425),c(rho1423,1,rho2324,rho2325),c(rho1424,rho2324,1,rho2425),c(rho1425,rho2325,rho2425,1))
  mu14232425<- c(mu14,mu23,mu24,mu25)
  p14232425<-  pmvnorm(mean=mu14232425, sigma=s14232425, lower=lw3, upper=up3)
  
  s15232425<- rbind(c(1,rho1523,rho1524,rho1525),c(rho1523,1,rho2324,rho2325),c(rho1524,rho2324,1,rho2425),c(rho1525,rho2325,rho2425,1))
  mu15232425<- c(mu15,mu23,mu24,mu25)
  p15232425<-  pmvnorm(mean=mu15232425, sigma=s15232425, lower=lw3, upper=up3)
  
  lw4 <- c(-q,-q,-q,-q,-q)
  up4 <- c(q,q,q,q,q)
  
  s1314152324<- rbind(c(1,rho1314,rho1315,rho1323,rho1324),c(rho1314,1,rho1415,rho1423,rho1424),c(rho1315,rho1415,1,rho1523,rho1524),c(rho1323,rho1423,rho1523,1,rho2324),c(rho1324,rho1424,rho1524,rho2324,1))
  mu1314152324<- c(mu13,mu14,mu15,mu23,mu24)
  p1314152324<-  pmvnorm(mean=mu1314152324, sigma=s1314152324, lower=lw4, upper=up4)                              
  
  s1314152325<- rbind(c(1,rho1314,rho1315,rho1323,rho1325),c(rho1314,1,rho1415,rho1423,rho1425),c(rho1315,rho1415,1,rho1523,rho1525),c(rho1323,rho1423,rho1523,1,rho2325),c(rho1325,rho1425,rho1525,rho2325,1))
  mu1314152325<- c(mu13,mu14,mu15,mu23,mu25)
  p1314152325<-  pmvnorm(mean=mu1314152325, sigma=s1314152325, lower=lw4, upper=up4)
  
  s1314152425<- rbind(c(1,rho1314,rho1315,rho1324,rho1325),c(rho1314,1,rho1415,rho1424,rho1425),c(rho1315,rho1415,1,rho1524,rho1525),c(rho1324,rho1424,rho1524,1,rho2425),c(rho1325,rho1425,rho1525,rho2425,1))
  mu1314152425<- c(mu13,mu14,mu15,mu24,mu25)
  p1314152425<-  pmvnorm(mean=mu1314152425, sigma=s1314152425, lower=lw4, upper=up4)
  
  s1314232425<- rbind(c(1,rho1314,rho1323,rho1324,rho1325),c(rho1314,1,rho1423,rho1424,rho1425),c(rho1323,rho1423,1,rho2324,rho2325),c(rho1324,rho1424,rho2324,1,rho2425),c(rho1325,rho1425,rho2325,rho2425,1))
  mu1314232425<- c(mu13,mu14,mu23,mu24,mu25)
  p1314232425<-  pmvnorm(mean=mu1314232425, sigma=s1314232425, lower=lw4, upper=up4)
  
  s1315232425<- rbind(c(1,rho1315,rho1323,rho1324,rho1325),c(rho1315,1,rho1523,rho1524,rho1525),c(rho1323,rho1523,1,rho2324,rho2325),c(rho1324,rho1524,rho2324,1,rho2425),c(rho1325,rho1525,rho2325,rho2425,1))
  mu1315232425<- c(mu13,mu15,mu23,mu24,mu25)
  p1315232425<-  pmvnorm(mean=mu1315232425, sigma=s1315232425, lower=lw4, upper=up4)
  
  s1415232425<- rbind(c(1,rho1415,rho1423,rho1424,rho1425),c(rho1415,1,rho1523,rho1524,rho1525),c(rho1423,rho1523,1,rho2324,rho2325),c(rho1424,rho1524,rho2324,1,rho2425),c(rho1425,rho1525,rho2325,rho2425,1))
  mu1415232425<- c(mu14,mu15,mu23,mu24,mu25)
  p1415232425<-  pmvnorm(mean=mu1415232425, sigma=s1415232425, lower=lw4, upper=up4)
  
  lw5 <- c(-q,-q,-q,-q,-q,-q)
  up5 <- c(q,q,q,q,q,q)
  
  s131415232425<- rbind(c(1,rho1314,rho1315,rho1323,rho1324,rho1325),c(rho1314,1,rho1415,rho1423,rho1424,rho1425),c(rho1315,rho1415,1,rho1523,rho1524,rho1525),c(rho1323,rho1423,rho1523,1,rho2324,rho2325),c(rho1324,rho1424,rho1524,rho2324,1,rho2425),c(rho1325,rho1425,rho1525,rho2325,rho2425,1))
  
  mu131415232425<- c(mu13,mu14,mu15,mu23,mu24,mu25)
  p131415232425<-  pmvnorm(mean=mu131415232425, sigma=s131415232425, lower=lw5, upper=up5)
  
  p1<- p13+p14+p15+p23+p24+p25
  p2<- p1314+p1315+p1323+p1324+p1325+p1415+p1423+p1424+p1425+p1523+p1524+p1525+p2324+p2325+p2425
  p3<- p131415+p131423+p131424+p131425+p131523+p131524+p131525+p132324+p132325+p132425+p141523+p141524+p141525+p142324+p142325+p142425+p152324+p152325+p152425+p232425
  p4<-  p13141523+ p13141524+ p13141525+ p13142324+ p13142325+ p13142425+p13152324+p13152325+p13152425+p13232425+p14152324+p14152325+p14152425+p14232425+p15232425
  p5<-  p1314152324+ p1314152325+ p1314152425+p1314232425+p1315232425+p1415232425
  return(-(1-(p1-p2+p3-p4+p5-p131415232425)))
  #return(p3)
}

# mu vector
sigma <- 1                    # Standard deviation
N<- 400                       #Total number of subjects
n<- sqrt(N)                   
al <- 0.05                    # significance level
q <- qnorm(1-(al))            # critical value

mu<-c(0,0.1,0.4,0.8,0.6)      # mean vector

# Lower and upper bounds
lb <- c(0,0,0,0,0)
ub <- c(0.35,0.35,0.35,0.35,0.35)



#initial values
x0 <- c(0.2,0.2,0.2,0.2,0.1)

# Set optimization options.
local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
opts <- list( "algorithm"= "NLOPT_GN_ISRES",
              "xtol_rel"= 1.0e-15,
              "maxeval"= 20000,
              "local_opts" = local_opts,
              "print_level" = 0 )


res <- nloptr ( x0 = x0,
                eval_f = IUT_K5_131415232425,
                mu=mu,
                n=n,
                sigma=sigma,
                q=q,
                lb = lb,
                ub = ub,
                #eval_g_eq = const1,
                opts = opts
)

print(res)

