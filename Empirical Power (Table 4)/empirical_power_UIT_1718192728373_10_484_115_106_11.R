rm(list = ls())
true.mean <- c(52.7,53.8,53.8,52.3,56.1,52.8,57.8,55.2,56.6,54.8,54.5)     # Sample mean
K<-length(true.mean)                                                      # No. of comparison groups 
psd<- 5.62                                                                  # Pooled standard deviation

n_m<- c(12,8,7,7,7,7,14,14,7,11,11)  # max-min design
n_b= c(10,10,10,10,10,10,9,9,9,9,9)   # balance design

m<-100000                                                 #No. of Iterations

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
  
  Z_17_m = (X.m_mean[1]-X.m_mean[7])/(psd*sqrt((1/n_m[1])+(1/n_m[7])))
  Z_18_m = (X.m_mean[1]-X.m_mean[8])/(psd*sqrt((1/n_m[1])+(1/n_m[8])))
  Z_19_m = (X.m_mean[1]-X.m_mean[9])/(psd*sqrt((1/n_m[1])+(1/n_m[9])))
  Z_27_m = (X.m_mean[2]-X.m_mean[7])/(psd*sqrt((1/n_m[2])+(1/n_m[7])))
  Z_28_m = (X.m_mean[2]-X.m_mean[8])/(psd*sqrt((1/n_m[2])+(1/n_m[8])))
  Z_37_m = (X.m_mean[3]-X.m_mean[7])/(psd*sqrt((1/n_m[3])+(1/n_m[7])))
  Z_310_m = (X.m_mean[3]-X.m_mean[10])/(psd*sqrt((1/n_m[3])+(1/n_m[10])))
  Z_48_m = (X.m_mean[4]-X.m_mean[8])/(psd*sqrt((1/n_m[4])+(1/n_m[8])))
  Z_411_m = (X.m_mean[4]-X.m_mean[11])/(psd*sqrt((1/n_m[4])+(1/n_m[11])))
  Z_510_m = (X.m_mean[5]-X.m_mean[10])/(psd*sqrt((1/n_m[5])+(1/n_m[10])))
  Z_611_m = (X.m_mean[6]-X.m_mean[11])/(psd*sqrt((1/n_m[6])+(1/n_m[11])))
  
  t_stat_m[i]<- max(c(abs(Z_17_m),abs(Z_18_m),abs(Z_19_m),abs(Z_27_m),abs(Z_28_m),abs(Z_37_m),abs(Z_310_m),abs(Z_48_m),abs(Z_411_m),abs(Z_510_m),abs(Z_611_m)))                          #test statistics under the null for max-min design
  
  Z_17_b = (X.b_mean[1]-X.b_mean[7])/(psd*sqrt((1/n_b[1])+(1/n_b[7])))
  Z_18_b = (X.b_mean[1]-X.b_mean[8])/(psd*sqrt((1/n_b[1])+(1/n_b[8])))
  Z_19_b = (X.b_mean[1]-X.b_mean[9])/(psd*sqrt((1/n_b[1])+(1/n_b[9])))
  Z_27_b = (X.b_mean[2]-X.b_mean[7])/(psd*sqrt((1/n_b[2])+(1/n_b[7])))
  Z_28_b = (X.b_mean[2]-X.b_mean[8])/(psd*sqrt((1/n_b[2])+(1/n_b[8])))
  Z_37_b = (X.b_mean[3]-X.b_mean[7])/(psd*sqrt((1/n_b[3])+(1/n_b[7])))
  Z_310_b = (X.b_mean[3]-X.b_mean[10])/(psd*sqrt((1/n_b[3])+(1/n_b[10])))
  Z_48_b = (X.b_mean[4]-X.b_mean[8])/(psd*sqrt((1/n_b[4])+(1/n_b[8])))
  Z_411_b = (X.b_mean[4]-X.b_mean[11])/(psd*sqrt((1/n_b[4])+(1/n_b[11])))
  Z_510_b = (X.b_mean[5]-X.b_mean[10])/(psd*sqrt((1/n_b[5])+(1/n_b[10])))
  Z_611_b = (X.b_mean[6]-X.b_mean[11])/(psd*sqrt((1/n_b[6])+(1/n_b[11])))
  
  t_stat_b[i]<- max(c(abs(Z_17_b),abs(Z_18_b),abs(Z_19_b),abs(Z_27_b),abs(Z_28_b),abs(Z_37_b),abs(Z_310_b),abs(Z_48_b),abs(Z_411_b),abs(Z_510_b),abs(Z_611_b)))                          #test statistics under the null for balanced design
  
}

q_m<- quantile(t_stat_b,probs = 0.95)                                       # Critical value for Max--min design
q_b<- quantile(t_stat_b,probs = 0.95)                                        # Critical value for Max--min design


for (i in 1:m) { 
  for (j in 1:K){
    Y.m <- rnorm(n_m[j],true.mean[j],psd)  # alternative data from max--min design
    Y.m_mean[j] <- mean(Y.m)
    
    Y.b <- rnorm(n_b[j],true.mean[j],psd)  # alternative data from balanced design
    Y.b_mean[j] <- mean(Y.b)
  }
  
  Z_17_m.A = (Y.m_mean[1]-Y.m_mean[7])/(psd*sqrt((1/n_m[1])+(1/n_m[7])))
  Z_18_m.A = (Y.m_mean[1]-Y.m_mean[8])/(psd*sqrt((1/n_m[1])+(1/n_m[8])))
  Z_19_m.A = (Y.m_mean[1]-Y.m_mean[9])/(psd*sqrt((1/n_m[1])+(1/n_m[9])))
  Z_27_m.A = (Y.m_mean[2]-Y.m_mean[7])/(psd*sqrt((1/n_m[2])+(1/n_m[7])))
  Z_28_m.A = (Y.m_mean[2]-Y.m_mean[8])/(psd*sqrt((1/n_m[2])+(1/n_m[8])))
  Z_37_m.A = (Y.m_mean[3]-Y.m_mean[7])/(psd*sqrt((1/n_m[3])+(1/n_m[7])))
  Z_310_m.A = (Y.m_mean[3]-Y.m_mean[10])/(psd*sqrt((1/n_m[3])+(1/n_m[10])))
  Z_48_m.A = (Y.m_mean[4]-Y.m_mean[8])/(psd*sqrt((1/n_m[4])+(1/n_m[8])))
  Z_411_m.A = (Y.m_mean[4]-Y.m_mean[11])/(psd*sqrt((1/n_m[4])+(1/n_m[11])))
  Z_510_m.A = (Y.m_mean[5]-Y.m_mean[10])/(psd*sqrt((1/n_m[5])+(1/n_m[10])))
  Z_611_m.A = (Y.m_mean[6]-Y.m_mean[11])/(psd*sqrt((1/n_m[6])+(1/n_m[11])))
  
  t_stat_m.A<-  max(c(abs(Z_17_m.A),abs(Z_18_m.A),abs(Z_19_m.A),abs(Z_27_m.A),abs(Z_28_m.A),abs(Z_37_m.A),abs(Z_310_m.A),abs(Z_48_m.A),abs(Z_411_m.A),abs(Z_510_m.A),abs(Z_611_m.A)))      #test statistics under the altenative for max-min design
  
  if(t_stat_m.A>q_m){
    ctr_m <- ctr_m+1
  }
  
  Z_17_b.A = (Y.b_mean[1]-Y.b_mean[7])/(psd*sqrt((1/n_b[1])+(1/n_b[7])))
  Z_18_b.A = (Y.b_mean[1]-Y.b_mean[8])/(psd*sqrt((1/n_b[1])+(1/n_b[8])))
  Z_19_b.A = (Y.b_mean[1]-Y.b_mean[9])/(psd*sqrt((1/n_b[1])+(1/n_b[9])))
  Z_27_b.A = (Y.b_mean[2]-Y.b_mean[7])/(psd*sqrt((1/n_b[2])+(1/n_b[7])))
  Z_28_b.A = (Y.b_mean[2]-Y.b_mean[8])/(psd*sqrt((1/n_b[2])+(1/n_b[8])))
  Z_37_b.A = (Y.b_mean[3]-Y.b_mean[7])/(psd*sqrt((1/n_b[3])+(1/n_b[7])))
  Z_310_b.A = (Y.b_mean[3]-Y.b_mean[10])/(psd*sqrt((1/n_b[3])+(1/n_b[10])))
  Z_48_b.A = (Y.b_mean[4]-Y.b_mean[8])/(psd*sqrt((1/n_b[4])+(1/n_b[8])))
  Z_411_b.A = (Y.b_mean[4]-Y.b_mean[11])/(psd*sqrt((1/n_b[4])+(1/n_b[11])))
  Z_510_b.A = (Y.b_mean[5]-Y.b_mean[10])/(psd*sqrt((1/n_b[5])+(1/n_b[10])))
  Z_611_b.A = (Y.b_mean[6]-Y.b_mean[11])/(psd*sqrt((1/n_b[6])+(1/n_b[11])))
  
  t_stat_b.A<-  max(c(abs(Z_17_b.A),abs(Z_18_b.A),abs(Z_19_b.A),abs(Z_27_b.A),abs(Z_28_b.A),abs(Z_37_b.A),abs(Z_310_b.A),abs(Z_48_b.A),abs(Z_411_b.A),abs(Z_510_b.A),abs(Z_611_b.A)))      #test statistics under the altenative for max-min design
  
  if (t_stat_b.A>q_b){
    ctr_b <- ctr_b+1
  }
  
}


print('CV max-min, CV balanced')
print(c(q_m,q_b))

print('power max-min, power balanced')
print(c(ctr_m/m,ctr_b/m))


