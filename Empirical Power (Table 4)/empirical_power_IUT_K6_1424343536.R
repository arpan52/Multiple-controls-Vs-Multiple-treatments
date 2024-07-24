rm(list = ls())
true.mean <- c(52.3,53.8,57.8,62.1,61.8,63.3)                     # Sample mean
K<-length(true.mean)                                            # No. of comparison groups      
psd<- 5.16                                                       # Pooled standard deviation

n_m<- c(14,14,22,22,14,14)  # max-min design
n_b= c(18,18,18,18,18,18)   # balance design




m<-10000                             #No. of Iterations

# Initializations

X.m_mean<-rep(0,K)
X.b_mean<-rep(0,K)

Y.m_mean<-rep(0,K)
Y.b_mean<-rep(0,K)


t_stat_m<-rep(0,m)
t_stat_b<-rep(0,m)


ctr_m<-0
ctr_b<-0


for (i in 1:m) {  
  for (j in 1:K){
    X.m <- rnorm(n_m[j],true.mean[1],psd)  # null data from max--min design
    X.m_mean[j] <- mean(X.m)
    
    X.b <- rnorm(n_b[j],true.mean[1],psd)  # null data from balanced design
    X.b_mean[j] <- mean(X.b)
  }
  
  Z_14_m = (X.m_mean[1]-X.m_mean[4])/(psd*sqrt((1/n_m[1])+(1/n_m[4])))
  Z_24_m = (X.m_mean[2]-X.m_mean[4])/(psd*sqrt((1/n_m[2])+(1/n_m[4])))
  Z_34_m = (X.m_mean[3]-X.m_mean[4])/(psd*sqrt((1/n_m[3])+(1/n_m[4])))
  Z_35_m = (X.m_mean[3]-X.m_mean[5])/(psd*sqrt((1/n_m[3])+(1/n_m[5])))
  Z_36_m = (X.m_mean[3]-X.m_mean[6])/(psd*sqrt((1/n_m[3])+(1/n_m[6])))
  
  t_stat_m[i]<- min(c(abs(Z_14_m),abs(Z_24_m),abs(Z_34_m),abs(Z_35_m),abs(Z_36_m)))                          #test statistics under the null for max-min design
  
  Z_14_b = (X.b_mean[1]-X.b_mean[4])/(psd*sqrt((1/n_b[1])+(1/n_b[4])))
  Z_24_b = (X.b_mean[2]-X.b_mean[4])/(psd*sqrt((1/n_b[2])+(1/n_b[4])))
  Z_34_b = (X.b_mean[3]-X.b_mean[4])/(psd*sqrt((1/n_b[3])+(1/n_b[4])))
  Z_35_b = (X.b_mean[3]-X.b_mean[5])/(psd*sqrt((1/n_b[3])+(1/n_b[5])))
  Z_36_b = (X.b_mean[3]-X.b_mean[6])/(psd*sqrt((1/n_b[3])+(1/n_b[6])))
  
  
  t_stat_b[i]<- min(c(abs(Z_14_b),abs(Z_24_b),abs(Z_34_b),abs(Z_35_b),abs(Z_35_b)))                          #test statistics under the null for balanced design
  
}

q_m<- quantile(t_stat_b,probs = 0.95)                               # Critical value for Max--min design
q_b<- quantile(t_stat_b,probs = 0.95)                               # Critical value for Balanced design


for (i in 1:m) { 
  for (j in 1:K){
    Y.m <- rnorm(n_m[j],true.mean[j],psd)  # alternative data from max--min design
    Y.m_mean[j] <- mean(Y.m)
    
    Y.b <- rnorm(n_b[j],true.mean[j],psd)  # alternative data from balanced design
    Y.b_mean[j] <- mean(Y.b)
  }
  
  Z_14_m.A = (Y.m_mean[1]-Y.m_mean[4])/(psd*sqrt((1/n_m[1])+(1/n_m[4])))
  Z_24_m.A = (Y.m_mean[2]-Y.m_mean[4])/(psd*sqrt((1/n_m[2])+(1/n_m[4])))
  Z_34_m.A = (Y.m_mean[3]-Y.m_mean[4])/(psd*sqrt((1/n_m[3])+(1/n_m[4])))
  Z_35_m.A = (Y.m_mean[3]-Y.m_mean[5])/(psd*sqrt((1/n_m[3])+(1/n_m[5])))
  Z_36_m.A = (Y.m_mean[3]-Y.m_mean[6])/(psd*sqrt((1/n_m[3])+(1/n_m[6])))
  
  
  t_stat_m.A<- min(c(abs(Z_14_m.A),abs(Z_24_m.A),abs(Z_34_m.A),abs(Z_35_m.A),abs(Z_36_m.A)))              #test statistics under the altenative for max-min design
  
  if(t_stat_m.A>q_m){
    ctr_m <- ctr_m+1
  }
  
  Z_14_b.A = (Y.b_mean[1]-Y.b_mean[4])/(psd*sqrt((1/n_b[1])+(1/n_b[4])))
  Z_24_b.A = (Y.b_mean[2]-Y.b_mean[4])/(psd*sqrt((1/n_b[2])+(1/n_b[4])))
  Z_34_b.A = (Y.b_mean[3]-Y.b_mean[4])/(psd*sqrt((1/n_b[3])+(1/n_b[4])))
  Z_35_b.A = (Y.b_mean[3]-Y.b_mean[5])/(psd*sqrt((1/n_b[3])+(1/n_b[5])))
  Z_36_b.A = (Y.b_mean[3]-Y.b_mean[6])/(psd*sqrt((1/n_b[3])+(1/n_b[6])))
  
  
  t_stat_b.A<- min(c(abs(Z_14_b.A),abs(Z_24_b.A),abs(Z_34_b.A),abs(Z_35_b.A),abs(Z_36_b.A)))                    #test statistics under the altenative for balanced design
  
  if (t_stat_b.A>q_b){
    ctr_b <- ctr_b+1
  }
  
}


print('CV max-min, CV balanced')
print(c(q_m,q_b))

print('power max-min, power balanced')
print(c(ctr_m/m,ctr_b/m))


