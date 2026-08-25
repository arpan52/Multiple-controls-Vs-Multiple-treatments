rm(list = ls())
set.seed(123)

true.mean <- true.mean <- c(
  5.50,6.26,3.38,3.15,5.52,5.77,5.77,5.52,5.77,5.52) 
K<-length(true.mean)
psd<- 2.29

n_m <- c(34,31,14,14,14,15,15,15,15,15)# c(24,12,10,9,9,9,9,9,9,9)   # max-min design
#n_b <- c(33,27,12,12,12,17,17,17,17,17) #c(20,16,7,7,8,11,10,10,10,10)  #rep(17,10)  #   # balance design
n_b <- rep(18,10)

m<-100000
X.m_mean<-rep(0,K)
X.b_mean<-rep(0,K)

Y.m_mean<-rep(0,K)
Y.b_mean<-rep(0,K)


t_stat_m<-rep(0,m)
t_stat_b<-rep(0,m)


ctr_m<-0
ctr_b<-0


for (i in 1:m) {  #simulating from max-min design
  for (j in 1:K){
    X.m <- rnorm(n_m[j],true.mean[1],psd)  # null data from max--min design
    X.m_mean[j] <- mean(X.m)
    
    X.b <- rnorm(n_b[j],true.mean[1],psd)  # null data from balanced design
    X.b_mean[j] <- mean(X.b)
  }
  
  ## Max-min design
  
  Z_13_m = (X.m_mean[1]-X.m_mean[3])/(psd*sqrt((1/n_m[1])+(1/n_m[3])))
  Z_14_m = (X.m_mean[1]-X.m_mean[4])/(psd*sqrt((1/n_m[1])+(1/n_m[4])))
  Z_15_m = (X.m_mean[1]-X.m_mean[5])/(psd*sqrt((1/n_m[1])+(1/n_m[5])))
  Z_16_m = (X.m_mean[1]-X.m_mean[6])/(psd*sqrt((1/n_m[1])+(1/n_m[6])))
  Z_17_m = (X.m_mean[1]-X.m_mean[7])/(psd*sqrt((1/n_m[1])+(1/n_m[7])))
  Z_18_m = (X.m_mean[1]-X.m_mean[8])/(psd*sqrt((1/n_m[1])+(1/n_m[8])))
  Z_19_m = (X.m_mean[1]-X.m_mean[9])/(psd*sqrt((1/n_m[1])+(1/n_m[9])))
  Z_110_m = (X.m_mean[1]-X.m_mean[10])/(psd*sqrt((1/n_m[1])+(1/n_m[10])))
  
  Z_26_m = (X.m_mean[2]-X.m_mean[6])/(psd*sqrt((1/n_m[2])+(1/n_m[6])))
  Z_27_m = (X.m_mean[2]-X.m_mean[7])/(psd*sqrt((1/n_m[2])+(1/n_m[7])))
  Z_28_m = (X.m_mean[2]-X.m_mean[8])/(psd*sqrt((1/n_m[2])+(1/n_m[8])))
  Z_29_m = (X.m_mean[2]-X.m_mean[9])/(psd*sqrt((1/n_m[2])+(1/n_m[9])))
  Z_210_m = (X.m_mean[2]-X.m_mean[10])/(psd*sqrt((1/n_m[2])+(1/n_m[10])))
  
  t_stat_m[i] <- max(abs(c(
    Z_13_m, Z_14_m, Z_15_m, Z_16_m,
    Z_17_m, Z_18_m, Z_19_m, Z_110_m,
    Z_26_m, Z_27_m, Z_28_m, Z_29_m, Z_210_m
  )))
  
  
  ## Balanced design
  
  Z_13_b = (X.b_mean[1]-X.b_mean[3])/(psd*sqrt((1/n_b[1])+(1/n_b[3])))
  Z_14_b = (X.b_mean[1]-X.b_mean[4])/(psd*sqrt((1/n_b[1])+(1/n_b[4])))
  Z_15_b = (X.b_mean[1]-X.b_mean[5])/(psd*sqrt((1/n_b[1])+(1/n_b[5])))
  Z_16_b = (X.b_mean[1]-X.b_mean[6])/(psd*sqrt((1/n_b[1])+(1/n_b[6])))
  Z_17_b = (X.b_mean[1]-X.b_mean[7])/(psd*sqrt((1/n_b[1])+(1/n_b[7])))
  Z_18_b = (X.b_mean[1]-X.b_mean[8])/(psd*sqrt((1/n_b[1])+(1/n_b[8])))
  Z_19_b = (X.b_mean[1]-X.b_mean[9])/(psd*sqrt((1/n_b[1])+(1/n_b[9])))
  Z_110_b = (X.b_mean[1]-X.b_mean[10])/(psd*sqrt((1/n_b[1])+(1/n_b[10])))
  
  Z_26_b = (X.b_mean[2]-X.b_mean[6])/(psd*sqrt((1/n_b[2])+(1/n_b[6])))
  Z_27_b = (X.b_mean[2]-X.b_mean[7])/(psd*sqrt((1/n_b[2])+(1/n_b[7])))
  Z_28_b = (X.b_mean[2]-X.b_mean[8])/(psd*sqrt((1/n_b[2])+(1/n_b[8])))
  Z_29_b = (X.b_mean[2]-X.b_mean[9])/(psd*sqrt((1/n_b[2])+(1/n_b[9])))
  Z_210_b = (X.b_mean[2]-X.b_mean[10])/(psd*sqrt((1/n_b[2])+(1/n_b[10])))
  
  t_stat_b[i] <- max(abs(c(
    Z_13_b, Z_14_b, Z_15_b, Z_16_b,
    Z_17_b, Z_18_b, Z_19_b, Z_110_b,
    Z_26_b, Z_27_b, Z_28_b, Z_29_b, Z_210_b
  )))
}

q_m<- quantile(t_stat_m,probs = 0.998)
q_b<- quantile(t_stat_b,probs = 0.998)


for (i in 1:m) {
  for (j in 1:K){
    Y.m <- rnorm(n_m[j],true.mean[j],psd)  # alternative data from max--min design
    Y.m_mean[j] <- mean(Y.m)
    
    Y.b <- rnorm(n_b[j],true.mean[j],psd)  # alternative data from balanced design
    Y.b_mean[j] <- mean(Y.b)
  }
  
  ## Max-min design
  
  Z_13_m.A = (Y.m_mean[1]-Y.m_mean[3])/(psd*sqrt((1/n_m[1])+(1/n_m[3])))
  Z_14_m.A = (Y.m_mean[1]-Y.m_mean[4])/(psd*sqrt((1/n_m[1])+(1/n_m[4])))
  Z_15_m.A = (Y.m_mean[1]-Y.m_mean[5])/(psd*sqrt((1/n_m[1])+(1/n_m[5])))
  Z_16_m.A = (Y.m_mean[1]-Y.m_mean[6])/(psd*sqrt((1/n_m[1])+(1/n_m[6])))
  Z_17_m.A = (Y.m_mean[1]-Y.m_mean[7])/(psd*sqrt((1/n_m[1])+(1/n_m[7])))
  Z_18_m.A = (Y.m_mean[1]-Y.m_mean[8])/(psd*sqrt((1/n_m[1])+(1/n_m[8])))
  Z_19_m.A = (Y.m_mean[1]-Y.m_mean[9])/(psd*sqrt((1/n_m[1])+(1/n_m[9])))
  Z_110_m.A = (Y.m_mean[1]-Y.m_mean[10])/(psd*sqrt((1/n_m[1])+(1/n_m[10])))
  
  Z_26_m.A = (Y.m_mean[2]-Y.m_mean[6])/(psd*sqrt((1/n_m[2])+(1/n_m[6])))
  Z_27_m.A = (Y.m_mean[2]-Y.m_mean[7])/(psd*sqrt((1/n_m[2])+(1/n_m[7])))
  Z_28_m.A = (Y.m_mean[2]-Y.m_mean[8])/(psd*sqrt((1/n_m[2])+(1/n_m[8])))
  Z_29_m.A = (Y.m_mean[2]-Y.m_mean[9])/(psd*sqrt((1/n_m[2])+(1/n_m[9])))
  Z_210_m.A = (Y.m_mean[2]-Y.m_mean[10])/(psd*sqrt((1/n_m[2])+(1/n_m[10])))
  
  t_stat_m.A <- max(abs(c(
    Z_13_m.A, Z_14_m.A, Z_15_m.A, Z_16_m.A,
    Z_17_m.A, Z_18_m.A, Z_19_m.A, Z_110_m.A,
    Z_26_m.A, Z_27_m.A, Z_28_m.A, Z_29_m.A, Z_210_m.A
  )))
  
  if(t_stat_m.A>q_m){
    ctr_m <- ctr_m+1
  }
  
  ## Balanced design
  
  Z_13_b.A = (Y.b_mean[1]-Y.b_mean[3])/(psd*sqrt((1/n_b[1])+(1/n_b[3])))
  Z_14_b.A = (Y.b_mean[1]-Y.b_mean[4])/(psd*sqrt((1/n_b[1])+(1/n_b[4])))
  Z_15_b.A = (Y.b_mean[1]-Y.b_mean[5])/(psd*sqrt((1/n_b[1])+(1/n_b[5])))
  Z_16_b.A = (Y.b_mean[1]-Y.b_mean[6])/(psd*sqrt((1/n_b[1])+(1/n_b[6])))
  Z_17_b.A = (Y.b_mean[1]-Y.b_mean[7])/(psd*sqrt((1/n_b[1])+(1/n_b[7])))
  Z_18_b.A = (Y.b_mean[1]-Y.b_mean[8])/(psd*sqrt((1/n_b[1])+(1/n_b[8])))
  Z_19_b.A = (Y.b_mean[1]-Y.b_mean[9])/(psd*sqrt((1/n_b[1])+(1/n_b[9])))
  Z_110_b.A = (Y.b_mean[1]-Y.b_mean[10])/(psd*sqrt((1/n_b[1])+(1/n_b[10])))
  
  Z_26_b.A = (Y.b_mean[2]-Y.b_mean[6])/(psd*sqrt((1/n_b[2])+(1/n_b[6])))
  Z_27_b.A = (Y.b_mean[2]-Y.b_mean[7])/(psd*sqrt((1/n_b[2])+(1/n_b[7])))
  Z_28_b.A = (Y.b_mean[2]-Y.b_mean[8])/(psd*sqrt((1/n_b[2])+(1/n_b[8])))
  Z_29_b.A = (Y.b_mean[2]-Y.b_mean[9])/(psd*sqrt((1/n_b[2])+(1/n_b[9])))
  Z_210_b.A = (Y.b_mean[2]-Y.b_mean[10])/(psd*sqrt((1/n_b[2])+(1/n_b[10])))
  
  t_stat_b.A <- max(abs(c(
    Z_13_b.A, Z_14_b.A, Z_15_b.A, Z_16_b.A,
    Z_17_b.A, Z_18_b.A, Z_19_b.A, Z_110_b.A,
    Z_26_b.A, Z_27_b.A, Z_28_b.A, Z_29_b.A, Z_210_b.A
  )))
  if (t_stat_b.A>q_b){
    ctr_b <- ctr_b+1
  }
  
}


print('CV max-min, CV balanced')
print(c(q_m,q_b))

print('power max-min, power balanced')
print(c(ctr_m/m,ctr_b/m))
