rm(list = ls())
true.mean <- c(54.5,55.2,65.0,62.1,63.3,64.9,61.8)                      # Sample mean
K<-length(true.mean)                                                    # No. of comparison groups   
psd<- 5.80                                                            # Pooled standard deviation

n_m<- c(7,17,8,7,7,8,8)  # max-min design
n_b= c(9,9,5,9,10,10,10)   # balance design


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
  
  Z_13_m = (X.m_mean[1]-X.m_mean[3])/(psd*sqrt((1/n_m[1])+(1/n_m[3])))
  Z_23_m = (X.m_mean[2]-X.m_mean[3])/(psd*sqrt((1/n_m[2])+(1/n_m[3])))
  Z_24_m = (X.m_mean[2]-X.m_mean[4])/(psd*sqrt((1/n_m[2])+(1/n_m[4])))
  Z_25_m = (X.m_mean[2]-X.m_mean[5])/(psd*sqrt((1/n_m[2])+(1/n_m[5])))
  Z_26_m = (X.m_mean[2]-X.m_mean[6])/(psd*sqrt((1/n_m[2])+(1/n_m[6])))
  Z_27_m = (X.m_mean[2]-X.m_mean[7])/(psd*sqrt((1/n_m[2])+(1/n_m[7])))

  
  t_stat_m[i]<- min(c(abs(Z_13_m),abs(Z_23_m),abs(Z_24_m),abs(Z_25_m),abs(Z_26_m),abs(Z_27_m)))                          #test statistics under the null for max-min design
  
  Z_13_b = (X.b_mean[1]-X.b_mean[3])/(psd*sqrt((1/n_b[1])+(1/n_b[3])))
  Z_23_b = (X.b_mean[2]-X.b_mean[3])/(psd*sqrt((1/n_b[2])+(1/n_b[3])))
  Z_24_b = (X.b_mean[2]-X.b_mean[4])/(psd*sqrt((1/n_b[2])+(1/n_b[4])))
  Z_25_b = (X.b_mean[2]-X.b_mean[5])/(psd*sqrt((1/n_b[2])+(1/n_b[5])))
  Z_26_b = (X.b_mean[2]-X.b_mean[6])/(psd*sqrt((1/n_b[2])+(1/n_b[6])))
  Z_27_b = (X.b_mean[2]-X.b_mean[7])/(psd*sqrt((1/n_b[2])+(1/n_b[7])))
  
  t_stat_b[i]<- min(c(abs(Z_13_b),abs(Z_23_b),abs(Z_24_b),abs(Z_25_b),abs(Z_26_b),abs(Z_27_b)))                          #test statistics under the null for balanced design
  
}

q_m<- quantile(t_stat_b,probs = 0.95)                                 # Critical value for Max--min design
q_b<- quantile(t_stat_b,probs = 0.95)                                  # Critical value for Max--min design


for (i in 1:m) { 
  for (j in 1:K){
    Y.m <- rnorm(n_m[j],true.mean[j],psd)  # alternative data from max--min design
    Y.m_mean[j] <- mean(Y.m)
    
    Y.b <- rnorm(n_b[j],true.mean[j],psd)  # alternative data from balanced design
    Y.b_mean[j] <- mean(Y.b)
  }
  
  Z_13_m.A = (Y.m_mean[1]-Y.m_mean[3])/(psd*sqrt((1/n_m[1])+(1/n_m[3])))
  Z_23_m.A = (Y.m_mean[2]-Y.m_mean[3])/(psd*sqrt((1/n_m[2])+(1/n_m[3])))
  Z_24_m.A = (Y.m_mean[2]-Y.m_mean[4])/(psd*sqrt((1/n_m[2])+(1/n_m[4])))
  Z_25_m.A = (Y.m_mean[2]-Y.m_mean[5])/(psd*sqrt((1/n_m[2])+(1/n_m[5])))
  Z_26_m.A = (Y.m_mean[2]-Y.m_mean[6])/(psd*sqrt((1/n_m[2])+(1/n_m[6])))
  Z_27_m.A = (Y.m_mean[2]-Y.m_mean[7])/(psd*sqrt((1/n_m[2])+(1/n_m[7])))
 
  
  t_stat_m.A<- min(c(abs(Z_13_m.A),abs(Z_23_m.A),abs(Z_24_m.A),abs(Z_25_m.A),abs(Z_26_m.A),abs(Z_27_m.A)))              #test statistics under the altenative for max-min design
  
  if(t_stat_m.A>q_m){
    ctr_m <- ctr_m+1
  }
  
  Z_13_b.A = (Y.b_mean[1]-Y.b_mean[3])/(psd*sqrt((1/n_b[1])+(1/n_b[3])))
  Z_23_b.A = (Y.b_mean[2]-Y.b_mean[3])/(psd*sqrt((1/n_b[2])+(1/n_b[3])))
  Z_24_b.A = (Y.b_mean[2]-Y.b_mean[4])/(psd*sqrt((1/n_b[2])+(1/n_b[4])))
  Z_25_b.A = (Y.b_mean[2]-Y.b_mean[5])/(psd*sqrt((1/n_b[2])+(1/n_b[5])))
  Z_26_b.A = (Y.b_mean[2]-Y.b_mean[6])/(psd*sqrt((1/n_b[2])+(1/n_b[6])))
  Z_27_b.A = (Y.b_mean[2]-Y.b_mean[7])/(psd*sqrt((1/n_b[2])+(1/n_b[7])))

  t_stat_b.A<- min(c(abs(Z_13_b.A),abs(Z_23_b.A),abs(Z_24_b.A),abs(Z_25_b.A),abs(Z_26_b.A),abs(Z_27_b.A)))                    #test statistics under the altenative for balanced design
  
  if (t_stat_b.A>q_b){
    ctr_b <- ctr_b+1
  }
  
}


print('CV max-min, CV balanced')
print(c(q_m,q_b))

print('power max-min, power balanced')
print(c(ctr_m/m,ctr_b/m))


