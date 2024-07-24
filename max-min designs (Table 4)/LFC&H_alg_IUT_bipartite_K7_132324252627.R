####################################################################################
###### H-algorithm for min-max design for Bipartite graph for K=7 for IUT         ##
###### pairs = (1,3), (2,3),(2,4), (2,5), (2,6), (2,7)                            ##
#####################################################################################



#install.packages("nloptr")
#install.packages("mvtnorm")


library('nloptr')
library(mvtnorm)
rm(list = ls())

# Initializations

design<- rep(0,7)

# Function for calculating the mean and covariance matrix of Z-ij's and criterion 'Psi'.
psi<- function(x,theta,n,sigma,q){
  
  design[1]<-x[1]
  design[2]<-x[2]
  design[3]<-x[3]
  design[4]<-x[4]
  design[5]<-x[4]
  design[6]<-x[4]
  design[7]<-x[4]
  
  
  mu13 <- (n*(theta[1]-theta[3]))/(sigma*sqrt((1/design[1])+(1/design[3])))
  mu23 <- (n*(theta[2]-theta[3]))/(sigma*sqrt((1/design[2])+(1/design[3])))
  mu24 <- (n*(theta[2]-theta[4]))/(sigma*sqrt((1/design[2])+(1/design[4])))
  mu25 <- (n*(theta[2]-theta[5]))/(sigma*sqrt((1/design[2])+(1/design[5])))
  mu26 <- (n*(theta[2]-theta[6]))/(sigma*sqrt((1/design[2])+(1/design[6])))
  mu27 <- (n*(theta[2]-theta[7]))/(sigma*sqrt((1/design[2])+(1/design[7])))
  
  # Calculating Bivariate, trivariate CDFs and so on...
  
  p13<- pnorm(q,mu13)-pnorm(-q,mu13)
  p23<- pnorm(q,mu23)-pnorm(-q,mu23)
  p24<- pnorm(q,mu24)-pnorm(-q,mu24)
  p25<- pnorm(q,mu25)-pnorm(-q,mu25)
  p26<- pnorm(q,mu26)-pnorm(-q,mu26)
  p27<- pnorm(q,mu27)-pnorm(-q,mu27)
  
  
  lw1 <- c(-q,-q)
  up1 <- c(q,q)
  
  
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
  
  rho1326 <- 0
  s1326<- rbind(c(1,rho1326),c(rho1326,1))
  mu1326<- c(mu13,mu26)
  p1326<-  pmvnorm(mean=mu1326, sigma=s1326, lower=lw1, upper=up1)
  
  rho1327 <- 0
  s1327<- rbind(c(1,rho1327),c(rho1327,1))
  mu1327<- c(mu13,mu27)
  p1327<-  pmvnorm(mean=mu1327, sigma=s1327, lower=lw1, upper=up1)
  
  
  rho2324 <- sqrt((design[3]*design[4])/((design[2]+design[3])*(design[2]+design[4])))
  s2324<- rbind(c(1,rho2324),c(rho2324,1))
  mu2324<- c(mu23,mu24)
  p2324<- pmvnorm(mean=mu2324, sigma=s2324, lower=lw1, upper=up1)
  
  
  rho2325 <- sqrt((design[3]*design[5])/((design[2]+design[3])*(design[2]+design[5])))
  s2325<- rbind(c(1,rho2325),c(rho2325,1))
  mu2325<- c(mu23,mu25)
  p2325<- pmvnorm(mean=mu2325, sigma=s2325, lower=lw1, upper=up1)
  
  rho2326 <- sqrt((design[3]*design[6])/((design[2]+design[3])*(design[2]+design[6])))
  s2326<- rbind(c(1,rho2326),c(rho2326,1))
  mu2326<- c(mu23,mu26)
  p2326<- pmvnorm(mean=mu2326, sigma=s2326, lower=lw1, upper=up1)
  
  rho2327 <- sqrt((design[3]*design[7])/((design[2]+design[3])*(design[2]+design[7])))
  s2327<- rbind(c(1,rho2327),c(rho2327,1))
  mu2327<- c(mu23,mu27)
  p2327<- pmvnorm(mean=mu2327, sigma=s2327, lower=lw1, upper=up1)
  
  rho2425 <- sqrt((design[4]*design[5])/((design[2]+design[4])*(design[2]+design[5])))
  s2425<- rbind(c(1,rho2425),c(rho2425,1))
  mu2425<- c(mu24,mu25)
  p2425<- pmvnorm(mean=mu2425, sigma=s2425, lower=lw1, upper=up1)
  
  rho2426 <- sqrt((design[4]*design[6])/((design[2]+design[4])*(design[2]+design[6])))
  s2426<- rbind(c(1,rho2426),c(rho2426,1))
  mu2426<- c(mu24,mu26)
  p2426<- pmvnorm(mean=mu2426, sigma=s2426, lower=lw1, upper=up1)
  
  rho2427 <- sqrt((design[4]*design[7])/((design[2]+design[4])*(design[2]+design[7])))
  s2427<- rbind(c(1,rho2427),c(rho2427,1))
  mu2427<- c(mu24,mu27)
  p2427<- pmvnorm(mean=mu2427, sigma=s2427, lower=lw1, upper=up1)
  
  rho2526 <- sqrt((design[5]*design[6])/((design[2]+design[5])*(design[2]+design[6])))
  s2526<- rbind(c(1,rho2526),c(rho2526,1))
  mu2526<- c(mu25,mu26)
  p2526<- pmvnorm(mean=mu2526, sigma=s2526, lower=lw1, upper=up1)
  
  rho2527 <- sqrt((design[5]*design[7])/((design[2]+design[5])*(design[2]+design[7])))
  s2527<- rbind(c(1,rho2527),c(rho2527,1))
  mu2527<- c(mu25,mu27)
  p2527<- pmvnorm(mean=mu2527, sigma=s2527, lower=lw1, upper=up1)
  
  rho2627 <- sqrt((design[6]*design[7])/((design[2]+design[6])*(design[2]+design[7])))
  s2627<- rbind(c(1,rho2627),c(rho2627,1))
  mu2627<- c(mu26,mu27)
  p2627<- pmvnorm(mean=mu2627, sigma=s2627, lower=lw1, upper=up1)
  
  lw2 <- c(-q,-q,-q)
  up2 <- c(q,q,q)
  
  
  s132324<- rbind(c(1,rho1323,rho1324),c(rho1323,1,rho2324),c(rho1324,rho2324,1))
  mu132324<- c(mu13,mu23,mu24)
  p132324<-  pmvnorm(mean=mu132324, sigma=s132324, lower=lw2, upper=up2)
  
  s132325<- rbind(c(1,rho1323,rho1325),c(rho1323,1,rho2325),c(rho1325,rho2325,1))
  mu132325<- c(mu13,mu23,mu25)
  p132325<-  pmvnorm(mean=mu132325, sigma=s132325, lower=lw2, upper=up2)
  
  s132326<- rbind(c(1,rho1323,rho1326),c(rho1323,1,rho2326),c(rho1326,rho2326,1))
  mu132326<- c(mu13,mu23,mu26)
  p132326<-  pmvnorm(mean=mu132326, sigma=s132326, lower=lw2, upper=up2)
  
  s132327<- rbind(c(1,rho1323,rho1327),c(rho1323,1,rho2327),c(rho1327,rho2327,1))
  mu132327<- c(mu13,mu23,mu27)
  p132327<-  pmvnorm(mean=mu132327, sigma=s132327, lower=lw2, upper=up2)
  
  s132425<- rbind(c(1,rho1324,rho1325),c(rho1324,1,rho2425),c(rho1325,rho2425,1))
  mu132425<- c(mu13,mu24,mu25)
  p132425<-  pmvnorm(mean=mu132425, sigma=s132425, lower=lw2, upper=up2)
  
  s132426<- rbind(c(1,rho1324,rho1326),c(rho1324,1,rho2426),c(rho1326,rho2426,1))
  mu132426<- c(mu13,mu24,mu26)
  p132426<-  pmvnorm(mean=mu132426, sigma=s132426, lower=lw2, upper=up2)
  
  s132427<- rbind(c(1,rho1324,rho1327),c(rho1324,1,rho2427),c(rho1327,rho2427,1))
  mu132427<- c(mu13,mu24,mu27)
  p132427<-  pmvnorm(mean=mu132427, sigma=s132427, lower=lw2, upper=up2)
  
  s132526<- rbind(c(1,rho1325,rho1326),c(rho1325,1,rho2526),c(rho1326,rho2526,1))
  mu132526<- c(mu13,mu25,mu26)
  p132526<-  pmvnorm(mean=mu132526, sigma=s132526, lower=lw2, upper=up2)
  
  s132527<- rbind(c(1,rho1325,rho1327),c(rho1325,1,rho2527),c(rho1327,rho2527,1))
  mu132527<- c(mu13,mu25,mu27)
  p132527<-  pmvnorm(mean=mu132527, sigma=s132527, lower=lw2, upper=up2)
  
  s132627<- rbind(c(1,rho1326,rho1327),c(rho1326,1,rho2627),c(rho1327,rho2627,1))
  mu132627<- c(mu13,mu26,mu27)
  p132627<-  pmvnorm(mean=mu132627, sigma=s132627, lower=lw2, upper=up2)
  
  
  s232425<- rbind(c(1,rho2324,rho2325),c(rho2324,1,rho2425),c(rho2325,rho2425,1))
  mu232425<- c(mu23,mu24,mu25)
  p232425<-  pmvnorm(mean=mu232425, sigma=s232425, lower=lw2, upper=up2)
  
  s232426<- rbind(c(1,rho2324,rho2326),c(rho2324,1,rho2426),c(rho2326,rho2426,1))
  mu232426<- c(mu23,mu24,mu26)
  p232426<-  pmvnorm(mean=mu232426, sigma=s232426, lower=lw2, upper=up2)
  
  s232427<- rbind(c(1,rho2324,rho2327),c(rho2324,1,rho2427),c(rho2327,rho2427,1))
  mu232427<- c(mu23,mu24,mu27)
  p232427<-  pmvnorm(mean=mu232427, sigma=s232427, lower=lw2, upper=up2)
  
  s232526<- rbind(c(1,rho2325,rho2326),c(rho2325,1,rho2526),c(rho2326,rho2526,1))
  mu232526<- c(mu23,mu25,mu26)
  p232526<-  pmvnorm(mean=mu232526, sigma=s232526, lower=lw2, upper=up2)
  
  s232527<- rbind(c(1,rho2325,rho2327),c(rho2325,1,rho2527),c(rho2327,rho2527,1))
  mu232527<- c(mu23,mu25,mu27)
  p232527<-  pmvnorm(mean=mu232527, sigma=s232527, lower=lw2, upper=up2)
  
  s232627<- rbind(c(1,rho2326,rho2327),c(rho2326,1,rho2627),c(rho2327,rho2627,1))
  mu232627<- c(mu23,mu26,mu27)
  p232627<-  pmvnorm(mean=mu232627, sigma=s232627, lower=lw2, upper=up2)
  
  s242526<- rbind(c(1,rho2425,rho2426),c(rho2425,1,rho2526),c(rho2426,rho2526,1))
  mu242526<- c(mu24,mu25,mu26)
  p242526<-  pmvnorm(mean=mu242526, sigma=s242526, lower=lw2, upper=up2)
  
  s242527<- rbind(c(1,rho2425,rho2427),c(rho2425,1,rho2527),c(rho2427,rho2527,1))
  mu242527<- c(mu24,mu25,mu27)
  p242527<-  pmvnorm(mean=mu242527, sigma=s242527, lower=lw2, upper=up2)
  
  s242627<- rbind(c(1,rho2426,rho2427),c(rho2426,1,rho2627),c(rho2427,rho2627,1))
  mu242627<- c(mu24,mu26,mu27)
  p242627<-  pmvnorm(mean=mu242627, sigma=s242627, lower=lw2, upper=up2)
  
  s252627<- rbind(c(1,rho2526,rho2527),c(rho2526,1,rho2627),c(rho2527,rho2627,1))
  mu252627<- c(mu25,mu26,mu27)
  p252627<-  pmvnorm(mean=mu252627, sigma=s252627, lower=lw2, upper=up2)
  
  lw3 <- c(-q,-q,-q,-q)
  up3 <- c(q,q,q,q)
  
  
  s13232425<- rbind(c(1,rho1323,rho1324,rho1325),c(rho1323,1,rho2324,rho2325),c(rho1324,rho2324,1,rho2425),c(rho1325,rho2325,rho2425,1))
  mu13232425<- c(mu13,mu23,mu24,mu25)
  p13232425<-  pmvnorm(mean=mu13232425, sigma=s13232425, lower=lw3, upper=up3)
  
  s13232426<- rbind(c(1,rho1323,rho1324,rho1326),c(rho1323,1,rho2324,rho2326),c(rho1324,rho2324,1,rho2426),c(rho1326,rho2326,rho2426,1))
  mu13232426<- c(mu13,mu23,mu24,mu26)
  p13232426<-  pmvnorm(mean=mu13232426, sigma=s13232426, lower=lw3, upper=up3)
  
  s13232427<- rbind(c(1,rho1323,rho1324,rho1327),c(rho1323,1,rho2324,rho2327),c(rho1324,rho2324,1,rho2427),c(rho1327,rho2327,rho2427,1))
  mu13232427<- c(mu13,mu23,mu24,mu27)
  p13232427<-  pmvnorm(mean=mu13232427, sigma=s13232427, lower=lw3, upper=up3)
  
  s13232526<- rbind(c(1,rho1323,rho1325,rho1326),c(rho1323,1,rho2325,rho2326),c(rho1325,rho2325,1,rho2526),c(rho1326,rho2326,rho2526,1))
  mu13232526<- c(mu13,mu23,mu25,mu26)
  p13232526<-  pmvnorm(mean=mu13232526, sigma=s13232526, lower=lw3, upper=up3)
  
  s13232526<- rbind(c(1,rho1323,rho1325,rho1326),c(rho1323,1,rho2325,rho2326),c(rho1325,rho2325,1,rho2526),c(rho1326,rho2326,rho2526,1))
  mu13232526<- c(mu13,mu23,mu25,mu26)
  p13232526<-  pmvnorm(mean=mu13232526, sigma=s13232526, lower=lw3, upper=up3)
  
  s13232527<- rbind(c(1,rho1323,rho1325,rho1327),c(rho1323,1,rho2325,rho2327),c(rho1325,rho2325,1,rho2527),c(rho1327,rho2327,rho2527,1))
  mu13232527<- c(mu13,mu23,mu25,mu27)
  p13232527<-  pmvnorm(mean=mu13232527, sigma=s13232527, lower=lw3, upper=up3)
  
  s13232627<- rbind(c(1,rho1323,rho1326,rho1327),c(rho1323,1,rho2326,rho2327),c(rho1326,rho2326,1,rho2627),c(rho1327,rho2327,rho2627,1))
  mu13232627<- c(mu13,mu23,mu26,mu27)
  p13232627<-  pmvnorm(mean=mu13232627, sigma=s13232627, lower=lw3, upper=up3)
  
  s13242526<- rbind(c(1,rho1324,rho1325,rho1326),c(rho1324,1,rho2425,rho2426),c(rho1325,rho2425,1,rho2526),c(rho1326,rho2426,rho2526,1))
  mu13242526<- c(mu13,mu24,mu25,mu26)
  p13242526<-  pmvnorm(mean=mu13242526, sigma=s13242526, lower=lw3, upper=up3)
  
  s13242527<- rbind(c(1,rho1324,rho1325,rho1327),c(rho1324,1,rho2425,rho2427),c(rho1325,rho2425,1,rho2527),c(rho1327,rho2427,rho2527,1))
  mu13242527<- c(mu13,mu24,mu25,mu27)
  p13242527<-  pmvnorm(mean=mu13242527, sigma=s13242527, lower=lw3, upper=up3)
  
  s13242627<- rbind(c(1,rho1324,rho1326,rho1327),c(rho1324,1,rho2426,rho2427),c(rho1326,rho2426,1,rho2627),c(rho1327,rho2427,rho2627,1))
  mu13242627<- c(mu13,mu24,mu26,mu27)
  p13242627<-  pmvnorm(mean=mu13242627, sigma=s13242627, lower=lw3, upper=up3)
  
  s13252627<- rbind(c(1,rho1325,rho1326,rho1327),c(rho1325,1,rho2526,rho2527),c(rho1326,rho2526,1,rho2627),c(rho1327,rho2527,rho2627,1))
  mu13252627<- c(mu13,mu25,mu26,mu27)
  p13252627<-  pmvnorm(mean=mu13252627, sigma=s13252627, lower=lw3, upper=up3)
  
  s23242526<- rbind(c(1,rho2324,rho2325,rho2326),c(rho2324,1,rho2425,rho2426),c(rho2325,rho2425,1,rho2526),c(rho2326,rho2426,rho2526,1))
  mu23242526<- c(mu23,mu24,mu25,mu26)
  p23242526<-  pmvnorm(mean=mu23242526, sigma=s23242526, lower=lw3, upper=up3)
  
  s23242527<- rbind(c(1,rho2324,rho2325,rho2327),c(rho2324,1,rho2425,rho2427),c(rho2325,rho2425,1,rho2527),c(rho2327,rho2427,rho2527,1))
  mu23242527<- c(mu23,mu24,mu25,mu27)
  p23242527<-  pmvnorm(mean=mu23242527, sigma=s23242527, lower=lw3, upper=up3)
  
  s23242627<- rbind(c(1,rho2324,rho2326,rho2327),c(rho2324,1,rho2426,rho2427),c(rho2326,rho2426,1,rho2627),c(rho2327,rho2427,rho2627,1))
  mu23242627<- c(mu23,mu24,mu26,mu27)
  p23242627<-  pmvnorm(mean=mu23242627, sigma=s23242627, lower=lw3, upper=up3)
  
  s23252627<- rbind(c(1,rho2325,rho2326,rho2327),c(rho2325,1,rho2526,rho2527),c(rho2326,rho2526,1,rho2627),c(rho2327,rho2527,rho2627,1))
  mu23252627<- c(mu23,mu25,mu26,mu27)
  p23252627<-  pmvnorm(mean=mu23252627, sigma=s23252627, lower=lw3, upper=up3)
  
  s24252627<- rbind(c(1,rho2425,rho2426,rho2427),c(rho2425,1,rho2526,rho2527),c(rho2426,rho2526,1,rho2627),c(rho2427,rho2527,rho2627,1))
  mu24252627<- c(mu24,mu25,mu26,mu27)
  p24252627<-  pmvnorm(mean=mu24252627, sigma=s24252627, lower=lw3, upper=up3)
  
  lw3 <- c(-q,-q,-q,-q,-q)
  up3 <- c(q,q,q,q,q)
  
  
  s1323242526<- rbind(c(1,rho1323,rho1324,rho1325,rho1326),c(rho1323,1,rho2324,rho2325,rho2326),c(rho1324,rho2324,1,rho2425,rho2426),c(rho1325,rho2325,rho2425,1,rho2526),c(rho1326,rho2326,rho2426,rho2526,1))
  mu1323242526<- c(mu13,mu23,mu24,mu25,mu26)
  p1323242526<-  pmvnorm(mean=mu1323242526, sigma=s1323242526, lower=lw3, upper=up3)
  
  s1323242527<- rbind(c(1,rho1323,rho1324,rho1325,rho1327),c(rho1323,1,rho2324,rho2325,rho2327),c(rho1324,rho2324,1,rho2425,rho2427),c(rho1325,rho2325,rho2425,1,rho2527),c(rho1327,rho2327,rho2427,rho2527,1))
  mu1323242527<- c(mu13,mu23,mu24,mu25,mu27)
  p1323242527<-  pmvnorm(mean=mu1323242527, sigma=s1323242527, lower=lw3, upper=up3)
  
  s1323242627<- rbind(c(1,rho1323,rho1324,rho1326,rho1327),c(rho1323,1,rho2324,rho2326,rho2327),c(rho1324,rho2324,1,rho2426,rho2427),c(rho1326,rho2326,rho2426,1,rho2627),c(rho1327,rho2327,rho2427,rho2627,1))
  mu1323242627<- c(mu13,mu23,mu24,mu26,mu27)
  p1323242627<-  pmvnorm(mean=mu1323242627, sigma=s1323242627, lower=lw3, upper=up3)
  
  s1323252627<- rbind(c(1,rho1323,rho1325,rho1326,rho1327),c(rho1323,1,rho2325,rho2326,rho2327),c(rho1325,rho2325,1,rho2526,rho2527),c(rho1326,rho2326,rho2526,1,rho2627),c(rho1327,rho2327,rho2527,rho2627,1))
  mu1323252627<- c(mu13,mu23,mu25,mu26,mu27)
  p1323252627<-  pmvnorm(mean=mu1323252627, sigma=s1323252627, lower=lw3, upper=up3)
  
  s1324252627<- rbind(c(1,rho1324,rho1325,rho1326,rho1327),c(rho1324,1,rho2425,rho2426,rho2427),c(rho1325,rho2425,1,rho2526,rho2527),c(rho1326,rho2426,rho2526,1,rho2627),c(rho1327,rho2427,rho2527,rho2627,1))
  mu1324252627<- c(mu13,mu24,mu25,mu26,mu27)
  p1324252627<-  pmvnorm(mean=mu1324252627, sigma=s1324252627, lower=lw3, upper=up3)
  
  s2324252627<- rbind(c(1,rho2324,rho2325,rho2326,rho2327),c(rho2324,1,rho2425,rho2426,rho2427),c(rho2325,rho2425,1,rho2526,rho2527),c(rho2326,rho2426,rho2526,1,rho2627),c(rho2327,rho2427,rho2527,rho2627,1))
  mu2324252627<- c(mu23,mu24,mu25,mu26,mu27)
  p2324252627<-  pmvnorm(mean=mu2324252627, sigma=s2324252627, lower=lw3, upper=up3)
  
  lw4 <- c(-q,-q,-q,-q,-q,-q)
  up4 <- c(q,q,q,q,q,q)
  
  s132324252627<- rbind(c(1,rho1323,rho1324,rho1325,rho1326,rho1327),c(rho1323,1,rho2324,rho2325,rho2326,rho2327),c(rho1324,rho2324,1,rho2425,rho2426,rho2427),c(rho1325,rho2325,rho2425,1,rho2526,rho2527),c(rho1326,rho2326,rho2426,rho2526,1,rho2627),c(rho1327,rho2327,rho2427,rho2527,rho2627,1))
  mu132324252627<- c(mu13,mu23,mu24,mu25,mu26,mu27)
  p132324252627<-  pmvnorm(mean=mu132324252627, sigma=s132324252627, lower=lw4, upper=up4)
  
  # Power of the test  
  p1<- p13+p23+p24+p25+p26+p27
  p2<- p1323+p1324+p1325+p1326+p1327+p2324+p2325+p2326+p2327+p2425+p2426+p2427+p2526+p2527+p2627
  p3<- p132324+p132325+p132326+p132327+p132425+p132426+p132427+p132526+p132527+p132627+p232425+p232426+p232427+p232526+p232527+p232627+p242526+p242527+p242627+p252627
  p4<- p13232425+p13232426+p13232427+p13232526+p13232527+p13232627+p13242526+p13242527+p13242627+p13252627+p23242526+p23242527+p23242627+p23252627+p24252627
  p5<- p1323242526+p1323242527+p1323242627+p1323252627+p1324252627+p2324252627
  p6<- p132324252627
  p<- -(1-p1+p2-p3+p4-p5+p6)      # Power of IUT
  return(p[1])
}
# Function for calculating the average power B

B<- function(x,theta_vals,prob_vals,n,sigma,q){
  
  s <- 0
  for (i in 1:length(prob_vals)){
    s <- s + prob_vals[i]* psi(x,theta_vals[i,],n,sigma,q)     # Average power associated to LFCs and probabilities
  }
  ret_val<-s
  return(ret_val[1])
}

# Function for calculating optimum on the avarage design

get_optimal_on_the_average_design<- function(theta_vals,prob_vals,n,sigma,q)
{
  des_0 <-c(0.15,0.2,0.12,0.1)
  #options = optimoptions('fmincon','Display','none')
  const1 <- function(x,theta_vals,prob_vals,n,sigma,q) {
    return( x[1]+x[2]+x[3]+4*x[4]-1)
  }
  const2 <- function(x,theta_vals,prob_vals,n,sigma,q) {
    return( c(x[1]-x[2],x[3]-x[2],x[4]-x[2],x[4]-x[3]))
  }
  
  lb <- c(0,0,0,0)
  ub <- c(0.5,0.5,0.5,0.5)
  
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
                  eval_g_ineq = const2,
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
  index <- which(this_psi == this_psi_min)   # LFC with minimum power
  
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

N <- 62                    # Total No. of subjects
n <- sqrt(N)
sigma <- 5.81              # Standard deviation
al <- 0.05               # significance level alpha
q <- qnorm(1-al)           # critical value
del<- 6.6                  # delta
iteration_num<-0
B_pi_l <- 100             # Initialization of B_pi
h_grid_space <- 0.5       # grid length
H<- seq(0, 1, by=h_grid_space)      # grid vector
0

print( "Initial Prior on Theta : ")

# Theta values (LFC's) are their corresponding prior probabilities



the_1<- c(-del,del,0,0,0,0,0)
the_2<-  c(del,-del,0,0,0,0,0)
the_3<-  c(-del,-del,0,0,0,0,0)
the_4<- c(0,0,-del,del,del,del,del)
the_5<- c(0,0,del,-del,del,del,del)
the_6<- c(0,0,del,del,-del,del,del)
the_7<- c(0,0,del,del,del,-del,del)
the_8<- c(0,0,del,del,del,del,-del)
the_9<- c(0,0,-del,-del,del,del,del)
the_10<- c(0,0,-del,del,-del,del,del)
the_11<- c(0,0,-del,del,del,-del,del)
the_12<- c(0,0,-del,del,del,del,-del)
the_13<- c(0,0,-del,-del,-del,-del,-del)



theta_vals_l<- rbind(the_1,the_2,the_3,the_4,the_5,the_6,the_7,the_8,the_9,the_10,the_11,the_12,the_13)

prob_vals_l <- rep(1/13,13)

# Starting optimum on the average design and corresponding power

print( "Starting Design")
design_l <- get_optimal_on_the_average_design(theta_vals_l,prob_vals_l,n,sigma,q)
B_pi_l <- B(design_l,theta_vals_l, prob_vals_l,n,sigma,q)

step_1(H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,n,sigma,q)