#####################################################################################
###### H-algorithm for min-max design for a general bipartite graph for K=6 for IUT##
###### Pairs =  (1,4),(2,4),(3,4),(3,5),(3,6)                                      ##
#####################################################################################


# (1,4),(2,4),(3,4),(3,5),(3,6)
# load the following libraries

library('nloptr')
library(mvtnorm)
rm(list = ls())
# initialization

design<- rep(0,6)
# Function for calculating the mean and covariance matrix of Z-ij's and criterion 'Psi'.
psi<- function(x,theta,n,sigma,q){
  
  design[1]<- (1-2*x)/4
  design[2]<- (1-2*x)/4
  design[3]<- x
  design[4]<- x
  design[5]<- (1-2*x)/4
  design[6]<- (1-2*x)/4
  
  
  mu14 <- (n*(theta[1]-theta[4]))/(sigma*sqrt((1/design[1])+(1/design[4])))
  mu24 <- (n*(theta[2]-theta[4]))/(sigma*sqrt((1/design[2])+(1/design[4])))
  mu34 <- (n*(theta[3]-theta[4]))/(sigma*sqrt((1/design[3])+(1/design[4])))
  mu35 <- (n*(theta[3]-theta[5]))/(sigma*sqrt((1/design[3])+(1/design[5])))
  mu36 <- (n*(theta[3]-theta[6]))/(sigma*sqrt((1/design[3])+(1/design[6])))
  
  # Computing bivariate, trivariate CDFs and so on........
  
  p14<- pnorm(q,mu14)-pnorm(-q,mu14)
  p24<- pnorm(q,mu24)-pnorm(-q,mu24)
  p34<- pnorm(q,mu34)-pnorm(-q,mu34)
  p35<- pnorm(q,mu35)-pnorm(-q,mu35)
  p36<- pnorm(q,mu36)-pnorm(-q,mu36)
  
  
  lw1 <- c(-q,-q)
  up1 <- c(q,q)
  
  rho1424 <- sqrt((design[1]*design[2])/((design[1]+design[4])*(design[2]+design[4])))
  s1424<- rbind(c(1,rho1424),c(rho1424,1))
  mu1424<- c(mu14,mu24)
  p1424<- pmvnorm(mean=mu1424, sigma=s1424, lower=lw1, upper=up1)
  
  rho1434 <- sqrt((design[1]*design[3])/((design[1]+design[4])*(design[3]+design[4])))
  s1434<- rbind(c(1,rho1434),c(rho1434,1))
  mu1434<- c(mu14,mu34)
  p1434<- pmvnorm(mean=mu1434, sigma=s1434, lower=lw1, upper=up1)
  
  rho1435 <- 0
  s1435<- rbind(c(1,rho1435),c(rho1435,1))
  mu1435<- c(mu14,mu35)
  p1435<- pmvnorm(mean=mu1435, sigma=s1435, lower=lw1, upper=up1)
  
  rho1436 <- 0
  s1436<- rbind(c(1,rho1436),c(rho1436,1))
  mu1436<- c(mu14,mu36)
  p1436<- pmvnorm(mean=mu1436, sigma=s1436, lower=lw1, upper=up1)
  
  
  rho2434 <- sqrt((design[2]*design[3])/((design[2]+design[4])*(design[3]+design[4])))
  s2434<- rbind(c(1,rho2434),c(rho2434,1))
  mu2434<- c(mu24,mu34)
  p2434<- pmvnorm(mean=mu2434, sigma=s2434, lower=lw1, upper=up1)
  
  rho2435 <- 0
  s2435<- rbind(c(1,rho2435),c(rho2435,1))
  mu2435<- c(mu24,mu35)
  p2435<- pmvnorm(mean=mu2435, sigma=s2435, lower=lw1, upper=up1)
  
  rho2436 <- 0
  s2436<- rbind(c(1,rho2436),c(rho2436,1))
  mu2436<- c(mu24,mu36)
  p2436<- pmvnorm(mean=mu2436, sigma=s2436, lower=lw1, upper=up1)
  
  rho3435 <- sqrt((design[4]*design[5])/((design[3]+design[4])*(design[3]+design[5])))
  s3435<- rbind(c(1,rho3435),c(rho3435,1))
  mu3435<- c(mu34,mu35)
  p3435<- pmvnorm(mean=mu3435, sigma=s3435, lower=lw1, upper=up1)
  
  rho3436 <- sqrt((design[4]*design[6])/((design[3]+design[4])*(design[3]+design[6])))
  s3436<- rbind(c(1,rho3436),c(rho3436,1))
  mu3436<- c(mu34,mu36)
  p3436<- pmvnorm(mean=mu3436, sigma=s3436, lower=lw1, upper=up1)
  
  rho3536 <- sqrt((design[5]*design[6])/((design[3]+design[5])*(design[3]+design[6])))
  s3536<- rbind(c(1,rho3536),c(rho3536,1))
  mu3536<- c(mu35,mu36)
  p3536<- pmvnorm(mean=mu3536, sigma=s3536, lower=lw1, upper=up1)
  
  
  lw2 <- c(-q,-q,-q)
  up2 <- c(q,q,q)
  
  s142434<- rbind(c(1,rho1424,rho1434),c(rho1424,1,rho2434),c(rho1434,rho2434,1))
  mu142434<- c(mu14,mu24,mu34)
  p142434<-  pmvnorm(mean=mu142434, sigma=s142434, lower=lw2, upper=up2)
  
  s142435<- rbind(c(1,rho1424,rho1435),c(rho1424,1,rho2435),c(rho1435,rho2435,1))
  mu142435<- c(mu14,mu24,mu35)
  p142435<-  pmvnorm(mean=mu142435, sigma=s142435, lower=lw2, upper=up2)
  
  s142436<- rbind(c(1,rho1424,rho1436),c(rho1424,1,rho2436),c(rho1436,rho2436,1))
  mu142436<- c(mu14,mu24,mu36)
  p142436<-  pmvnorm(mean=mu142436, sigma=s142436, lower=lw2, upper=up2)
  
  s143435<- rbind(c(1,rho1434,rho1435),c(rho1434,1,rho3435),c(rho1435,rho3435,1))
  mu143435<- c(mu14,mu34,mu35)
  p143435<-  pmvnorm(mean=mu143435, sigma=s143435, lower=lw2, upper=up2)
  
  s143436<- rbind(c(1,rho1434,rho1436),c(rho1434,1,rho3436),c(rho1436,rho3436,1))
  mu143436<- c(mu14,mu34,mu36)
  p143436<-  pmvnorm(mean=mu143436, sigma=s143436, lower=lw2, upper=up2)
  
  s143536<- rbind(c(1,rho1435,rho1436),c(rho1435,1,rho3536),c(rho1436,rho3536,1))
  mu143536<- c(mu14,mu35,mu36)
  p143536<-  pmvnorm(mean=mu143536, sigma=s143536, lower=lw2, upper=up2)
  
  s243435<- rbind(c(1,rho2434,rho2435),c(rho2434,1,rho3435),c(rho2435,rho3435,1))
  mu243435<- c(mu24,mu34,mu35)
  p243435<-  pmvnorm(mean=mu243435, sigma=s243435, lower=lw2, upper=up2)
  
  s243436<- rbind(c(1,rho2434,rho2436),c(rho2434,1,rho3436),c(rho2436,rho3436,1))
  mu243436<- c(mu24,mu34,mu36)
  p243436<-  pmvnorm(mean=mu243436, sigma=s243436, lower=lw2, upper=up2)
  
  s243536<- rbind(c(1,rho2435,rho2436),c(rho2435,1,rho3536),c(rho2436,rho3536,1))
  mu243536<- c(mu24,mu35,mu36)
  p243536<-  pmvnorm(mean=mu243536, sigma=s243536, lower=lw2, upper=up2)
  
  s343536<- rbind(c(1,rho3435,rho3436),c(rho3435,1,rho3536),c(rho3436,rho3536,1))
  mu343536<- c(mu34,mu35,mu36)
  p343536<-  pmvnorm(mean=mu343536, sigma=s343536, lower=lw2, upper=up2)
  
  
  lw3 <- c(-q,-q,-q,-q)
  up3 <- c(q,q,q,q)
  
  s14243435<- rbind(c(1,rho1424,rho1434,rho1435),c(rho1424,1,rho2434,rho2435),c(rho1434,rho2434,1,rho3435),c(rho1435,rho2435,rho3435,1))
  mu14243435<- c(mu14,mu24,mu34,mu35)
  p14243435<-  pmvnorm(mean=mu14243435, sigma=s14243435, lower=lw3, upper=up3)
  
  s14243436<- rbind(c(1,rho1424,rho1434,rho1436),c(rho1424,1,rho2434,rho2436),c(rho1434,rho2434,1,rho3436),c(rho1436,rho2436,rho3436,1))
  mu14243436<- c(mu14,mu24,mu34,mu36)
  p14243436<-  pmvnorm(mean=mu14243436, sigma=s14243436, lower=lw3, upper=up3)
  
  s14243536<- rbind(c(1,rho1424,rho1435,rho1436),c(rho1424,1,rho2435,rho2436),c(rho1435,rho2435,1,rho3536),c(rho1436,rho2436,rho3536,1))
  mu14243536<- c(mu14,mu24,mu35,mu36)
  p14243536<-  pmvnorm(mean=mu14243536, sigma=s14243536, lower=lw3, upper=up3)
  
  s14343536<- rbind(c(1,rho1434,rho1435,rho1436),c(rho1434,1,rho3435,rho3436),c(rho1435,rho3435,1,rho3536),c(rho1436,rho3436,rho3536,1))
  mu14343536<- c(mu14,mu34,mu35,mu36)
  p14343536<-  pmvnorm(mean=mu14343536, sigma=s14343536, lower=lw3, upper=up3)
  
  s24343536<- rbind(c(1,rho2434,rho2435,rho2436),c(rho2434,1,rho3435,rho3436),c(rho2435,rho3435,1,rho3536),c(rho2436,rho3436,rho3536,1))
  mu24343536<- c(mu24,mu34,mu35,mu36)
  p24343536<-  pmvnorm(mean=mu24343536, sigma=s24343536, lower=lw3, upper=up3)
  
  lw3 <- c(-q,-q,-q,-q,-q)
  up3 <- c(q,q,q,q,q)
  
  
  s1424343536<- rbind(c(1,rho1424,rho1434,rho1435,rho1436),c(rho1424,1,rho2434,rho2435,rho2436),c(rho1434,rho2434,1,rho3435,rho3436),c(rho1435,rho2435,rho3435,1,rho3536),c(rho1436,rho2436,rho3436,rho3536,1))
  mu1424343536<- c(mu14,mu24,mu34,mu35,mu36)
  p1424343536<-  pmvnorm(mean=mu1424343536, sigma=s1424343536, lower=lw3, upper=up3)
  
  
  
  # Power of the test  
  p1<- p14+p24+p34+p35+p36
  p2<- p1424+p1434+p1435+p1436+p2434+p2435+p2436+p3435+p3436+p3536
  p3<- p142434+p142435+p142436+p143435+p143436+p143536+p243435+p243436+p243536+p343536
  p4<-  p14243435+p14243436+p14243536+p14343536+p24343536
  p5<- p1424343536
  p<- -(1-p1+p2-p3+p4-p5)   # Power of IUT
  return(p[1])
}


# Function for calculating the average power B

B<- function(x,theta_vals,prob_vals,n,sigma,q){
  
  s <- 0
  for (i in 1:length(prob_vals)){
    s <- s + prob_vals[i]* psi(x,theta_vals[i,],n,sigma,q)   # Average power associated to LFCs and Probabilities
  }
  ret_val<-s
  return(ret_val[1])
}

# Function for calculating optimum on the avarage design

get_optimal_on_the_average_design<- function(theta_vals,prob_vals,n,sigma,q)
{
  des_0 <-c(0.2)
  
  
  lb <- c(0)
  ub <- c(0.5)
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
                  #eval_g_eq = const1,
                  #eval_g_ineq = const2,
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
  index <- which(this_psi == this_psi_min)    # Find LFCs with minimum power 
  
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

N <- 58                    # Total No. of subjects
n <- sqrt(N)
sigma <- 5.16             # Standard deviation
al <- 0.05               # significance level alpha
q <- qnorm(1-al)           # critical value
del<- 5.7                  # delta
iteration_num<-0
B_pi_l <- 100             # Initialization of B_pi
h_grid_space <- 0.5       # grid length
H<- seq(0, 1, by=h_grid_space)      # grid vector
0

print( "Initial Prior on Theta : ")

# Theta values (LFC's) are their corresponding prior probabilities



the_1<- c(-del,del,del,0,0,0)
the_2<-  c(del,-del,del,0,0,0)
the_3<-  c(del,del,-del,0,0,0)
the_4<-  c(del,del,del,0,0,0)
the_5<- c(0,0,0,-del,del,del)
the_6<- c(0,0,0,del,-del,del)
the_7<- c(0,0,0,del,del,-del)
the_8<- c(0,0,0,del,del,del)



theta_vals_l<- rbind(the_1,the_2,the_3,the_4,the_5,the_6,the_7,the_8)

prob_vals_l <- rep(1/8,8)

# Starting optimum on the average design and corresponding power

print( "Starting Design")
design_l <- get_optimal_on_the_average_design(theta_vals_l,prob_vals_l,n,sigma,q)
B_pi_l <- B(design_l,theta_vals_l, prob_vals_l,n,sigma,q)

step_1(H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,n,sigma,q)

