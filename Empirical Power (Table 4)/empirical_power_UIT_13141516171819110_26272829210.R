rm(list = ls())
true.mean <- c(52.3,54.5,53.8,53.8,52.7,52.8,55.2,56.6,54.8,57.8)  # Sample mean
K<-length(true.mean)                                               # No. of comparison groups
psd<- 5.62                                                         # Pooled standard deviation

n_m<- c(21,10,8,8,8,8,8,8,8,8)  # max-min design
n_b= c(10,9,10,10,10,10,9,9,9,9)   # balanced design

m<-100000                            #No. of Iterations

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

  t_stat_m[i]<- max(c(abs(Z_13_m),abs(Z_14_m),abs(Z_15_m),abs(Z_16_m),abs(Z_17_m),abs(Z_18_m),abs(Z_19_m),abs(Z_110_m),abs(Z_26_m),abs(Z_27_m),abs(Z_28_m),abs(Z_29_m),abs(Z_210_m)))                          #test statistics under the null for max-min design
  
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
  
  t_stat_b[i]<- max(c(abs(Z_13_b),abs(Z_14_b),abs(Z_15_b),abs(Z_16_b),abs(Z_17_b),abs(Z_18_b),abs(Z_19_b),abs(Z_110_b),abs(Z_26_b),abs(Z_27_b),abs(Z_28_b),abs(Z_29_b),abs(Z_210_b)))                          #test statistics under the null for balanced design

}

q_m<- quantile(t_stat_b,probs = 0.95)      # Critical value for Max--min design
q_b<- quantile(t_stat_b,probs = 0.95)      # Critical value for Balanced design


for (i in 1:m) { 
  for (j in 1:K){
    Y.m <- rnorm(n_m[j],true.mean[j],psd)  # alternative data from max--min design
    Y.m_mean[j] <- mean(Y.m)
    
    Y.b <- rnorm(n_b[j],true.mean[j],psd)  # alternative data from balanced design
    Y.b_mean[j] <- mean(Y.b)
  }

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
  
  t_stat_m.A<- max(c(abs(Z_13_m.A),abs(Z_14_m.A),abs(Z_15_m.A),abs(Z_16_m.A),abs(Z_17_m.A),abs(Z_18_m.A),abs(Z_19_m.A),abs(Z_110_m.A),abs(Z_26_m.A),abs(Z_27_m.A),abs(Z_28_m.A),abs(Z_29_m.A),abs(Z_210_m.A)))              #test statistics under the altenative for max-min design
  
  if(t_stat_m.A>q_m){
    ctr_m <- ctr_m+1
  }
  
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
  
  t_stat_b.A<- max(c(abs(Z_13_b.A),abs(Z_14_b.A),abs(Z_15_b.A),abs(Z_16_b.A),abs(Z_17_b.A),abs(Z_18_b.A),abs(Z_19_b.A),abs(Z_110_b.A),abs(Z_26_b.A),abs(Z_27_b.A),abs(Z_28_b.A),abs(Z_29_b.A),abs(Z_210_b.A)))                    #test statistics under the altenative for balanced design
  
  if (t_stat_b.A>q_b){
    ctr_b <- ctr_b+1
  }
  
}


print('CV max-min, CV balanced')
print(c(q_m,q_b))

print('power max-min, power balanced')
print(c(ctr_m/m,ctr_b/m))


