rm(list = ls())
set.seed(123)

true.mean <-  c(6.26,5.50,3.15,5.77,5.52,3.38)  
K<-length(true.mean)
psd<- 2.307

n_m<-  c(18, 26, 15,15,15,20)  # max-min design
#n_b<- c(13,30,15,15,15,20)   # C-optimal design
n_b <- c(17,19,20,18,17,18)  # Balanced design


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
  
  Z_23_m = (X.m_mean[2]-X.m_mean[3])/(psd*sqrt((1/n_m[2])+(1/n_m[3])))
  Z_24_m = (X.m_mean[2]-X.m_mean[4])/(psd*sqrt((1/n_m[2])+(1/n_m[4])))
  Z_25_m = (X.m_mean[2]-X.m_mean[5])/(psd*sqrt((1/n_m[2])+(1/n_m[5])))
  Z_26_m = (X.m_mean[2]-X.m_mean[6])/(psd*sqrt((1/n_m[2])+(1/n_m[6])))
  Z_16_m = (X.m_mean[1]-X.m_mean[6])/(psd*sqrt((1/n_m[1])+(1/n_m[6])))
  
  
  t_stat_m[i]<- max(c(abs(Z_23_m),abs(Z_24_m),abs(Z_25_m),abs(Z_26_m),abs(Z_16_m)))                          #test statistics under the null for max-min design
  
  Z_23_b = (X.b_mean[2]-X.b_mean[3])/(psd*sqrt((1/n_b[2])+(1/n_b[3])))
  Z_24_b = (X.b_mean[2]-X.b_mean[4])/(psd*sqrt((1/n_b[2])+(1/n_b[4])))
  Z_25_b = (X.b_mean[2]-X.b_mean[5])/(psd*sqrt((1/n_b[2])+(1/n_b[5])))
  Z_26_b = (X.b_mean[2]-X.b_mean[6])/(psd*sqrt((1/n_b[2])+(1/n_b[6])))
  Z_16_b = (X.b_mean[1]-X.b_mean[6])/(psd*sqrt((1/n_b[1])+(1/n_b[6])))
  
  
  t_stat_b[i]<- max(c(abs(Z_23_b),abs(Z_24_b),abs(Z_25_b),abs(Z_26_b),abs(Z_16_b)))                          #test statistics under the null for balanced design
  
}

q_m<- quantile(t_stat_m,probs = 0.995)
q_b<- quantile(t_stat_b,probs = 0.995)


for (i in 1:m) {
  for (j in 1:K){
    Y.m <- rnorm(n_m[j],true.mean[j],psd)  # alternative data from max--min design
    Y.m_mean[j] <- mean(Y.m)
    
    Y.b <- rnorm(n_b[j],true.mean[j],psd)  # alternative data from balanced design
    Y.b_mean[j] <- mean(Y.b)
  }
  
  
  Z_23_m.A = (Y.m_mean[2]-Y.m_mean[3])/(psd*sqrt((1/n_m[2])+(1/n_m[3])))
  Z_24_m.A = (Y.m_mean[2]-Y.m_mean[4])/(psd*sqrt((1/n_m[2])+(1/n_m[4])))
  Z_25_m.A = (Y.m_mean[2]-Y.m_mean[5])/(psd*sqrt((1/n_m[2])+(1/n_m[5])))
  Z_26_m.A = (Y.m_mean[2]-Y.m_mean[6])/(psd*sqrt((1/n_m[2])+(1/n_m[6])))
  Z_16_m.A = (Y.m_mean[1]-Y.m_mean[6])/(psd*sqrt((1/n_m[1])+(1/n_m[6])))
  
  
  t_stat_m.A<- max(c(abs(Z_23_m.A),abs(Z_24_m.A),abs(Z_25_m.A),abs(Z_26_m.A),abs(Z_16_m.A)))              #test statistics under the altenative for max-min design
  
  if(t_stat_m.A>q_m){
    ctr_m <- ctr_m+1
  }
  
  Z_23_b.A = (Y.b_mean[2]-Y.b_mean[3])/(psd*sqrt((1/n_b[2])+(1/n_b[3])))
  Z_24_b.A = (Y.b_mean[2]-Y.b_mean[4])/(psd*sqrt((1/n_b[2])+(1/n_b[4])))
  Z_25_b.A = (Y.b_mean[2]-Y.b_mean[5])/(psd*sqrt((1/n_b[2])+(1/n_b[5])))
  Z_26_b.A = (Y.b_mean[2]-Y.b_mean[6])/(psd*sqrt((1/n_b[2])+(1/n_b[6])))
  Z_16_b.A = (Y.b_mean[1]-Y.b_mean[6])/(psd*sqrt((1/n_b[1])+(1/n_b[6])))
  
  
  
  t_stat_b.A<- max(c(abs(Z_23_b.A),abs(Z_24_b.A),abs(Z_25_b.A),abs(Z_26_b.A),abs(Z_16_b.A)))                    #test statistics under the altenative for balanced design
  
  if (t_stat_b.A>q_b){
    ctr_b <- ctr_b+1
  }
  
}


print('CV max-min, CV balanced')
print(c(q_m,q_b))

print('power max-min, power balanced')
print(c(ctr_m/m,ctr_b/m))
