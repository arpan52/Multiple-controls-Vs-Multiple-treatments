#################################################################################################
###### H-algorithm for min-max design for a Complete bipartite graph for K1=2 and K2=4 for UIT ##
#################################################################################################



#install.packages("nloptr")
#install.packages("mvtnorm")


library('nloptr')
library(mvtnorm)
library(tictoc)
tic()

# Initializations
rm(list = ls())
m<-rep(0,6)
# Function for calculating the mean and covariance matrix of Z-ij's and criterion 'Psi'.

psi<- function(design,theta,del,n,sigma,q)
{
  
  t1<- design
  t2<- (1-(2*t1))/4
  
  m[1]<-0
  m[2]<-theta[1]
  m[3]<-del
  m[4]<-theta[2]
  m[5]<-theta[2]
  m[6]<-theta[2]
  
  
  
  de13 <- ((m[1]-m[3])*n)/sigma
  de14 <- ((m[1]-m[4])*n)/sigma
  de15 <- ((m[1]-m[5])*n)/sigma
  de16 <- ((m[1]-m[6])*n)/sigma
  de23 <- ((m[2]-m[3])*n)/sigma
  de24 <- ((m[2]-m[4])*n)/sigma
  de25 <- ((m[2]-m[5])*n)/sigma
  de26 <- ((m[2]-m[6])*n)/sigma
  
  
  
  
  a1 <- t2*t2
  b1 <- (t1+t2)*(t1+t2)
  rho1 <- sqrt(a1/b1)
  
  a2 <- t2*t2
  b2 <- (t1+t2)*(t1+t2)
  rho2 <- sqrt(a2/b2)
  
  a3 <- t2*t2
  b3 <- (t1+t2)*(t1+t2)
  rho3 <- sqrt(a3/b3)
  
  a4 <- t1*t1
  b4 <- (t1+t2)*(t1+t2)
  rho4 <- sqrt(a4/b4)
  
  rho5<- 0
  rho6<- 0
  rho7<-0
  
  
  a8 <- t2*t2
  b8 <- (t1+t2)*(t1+t2)
  rho8 <- sqrt(a8/b8)
  
  a9 <- t2*t2
  b9 <- (t1+t2)*(t1+t2)
  rho9 <- sqrt(a9/b9)
  
  rho10<- 0
  
  a11 <- t1*t1
  b11<- (t1+t2)*(t1+t2)
  rho11<- sqrt(a11/b11)
  
  rho12<-0
  rho13<- 0
  
  a14 <- t2*t2
  b14<- (t1+t2)*(t1+t2)
  rho14<- sqrt(a14/b14)
  
  rho15<-0
  rho16<-0
  
  
  a17 <- t1*t1
  b17<- (t1+t2)*(t1+t2)
  rho17<- sqrt(a17/b17)
  
  rho18<-0
  rho19<-0
  rho20<-0
  rho21<-0
  
  a22 <- t1*t1
  b22<- (t1+t2)*(t1+t2)
  rho22<- sqrt(a22/b22)
  
  
  a23 <- t2*t2
  b23<- (t1+t2)*(t1+t2)
  rho23<- sqrt(a23/b23)
  
  a24 <- t2*t2
  b24<- (t1+t2)*(t1+t2)
  rho24<- sqrt(a24/b24)
  
  a25 <- t2*t2
  b25<- (t1+t2)*(t1+t2)
  rho25<- sqrt(a25/b25)
  
  
  a26 <- t2*t2
  b26<- (t1+t2)*(t1+t2)
  rho26<- sqrt(a26/b26)
  
  a27 <- t2*t2
  b27<- (t1+t2)*(t1+t2)
  rho27<- sqrt(a27/b27)
  
  a28 <- t2*t2
  b28<- (t1+t2)*(t1+t2)
  rho28<- sqrt(a28/b28)
  
  
  mu1 <- (sqrt((t1*t2)/(t1 + t2)))*de13
  mu2 <- (sqrt((t1*t2)/(t1 + t2)))*de14
  mu3 <- (sqrt((t1*t2)/(t1 + t2)))*de15
  mu4 <- (sqrt((t1*t2)/(t1 + t2)))*de16
  mu5 <- (sqrt((t1*t2)/(t1 + t2)))*de23
  mu6 <- (sqrt((t1*t2)/(t1 + t2)))*de24
  mu7 <- (sqrt((t1*t2)/(t1 + t2)))*de25
  mu8 <- (sqrt((t1*t2)/(t1 + t2)))*de26
  
  m <- c(mu1,mu2,mu3,mu4,mu5,mu6,mu7,mu8)
  S <-  rbind(c(1,rho1,rho2,rho3,rho4,rho5,rho6,rho7),c(rho1,1,rho8,rho9,rho10,rho11,rho12,rho13),c(rho2,rho8,1,rho14,rho15,rho16,rho17,rho18),c(rho3,rho9,rho14,1,rho19,rho20,rho21,rho22),c(rho4,rho10,rho15,rho19,1,rho23,rho24,rho25),c(rho5,rho11,rho16,rho20,rho23,1,rho26,rho27),c(rho6,rho12,rho17,rho21,rho24,rho26,1,rho28),c(rho7,rho13,rho18,rho22,rho25,rho27,rho28,1))
  
  lw1 <- c(-q,-q,-q,-q,-q,-q,-q,-q)
  up1 <- c(q,q,q,q,q,q,q,q)
  
  return(-(1-pmvnorm(mean=m, sigma=S, lower=lw1, upper=up1)))   #Power of UIT
}
psi1<- function(theta,design,del,n,sigma,q){
  # Function to compute LFCS
  
  
  t1<- design
  t2<- (1-2*t1)/4
  
  m[1]<-0
  m[2]<-theta[1]
  m[3]<-del
  m[4]<-theta[2]
  m[5]<-theta[2]
  m[6]<-theta[2]
  
  return(-psi(design,mu,del,n,sigma,q))
}
# Function for compute the average power B

B<- function(design,theta_vals,prob_vals,del,n,sigma,q){
  
  s <- 0
  for (i in 1:length(prob_vals)){
    s <- s + prob_vals[i]* psi(design,theta_vals[i,],del,n,sigma,q)     # average power associated to LFCs and Probabilities
  }
  ret_val<-s
  return(ret_val[1])
}

# Function for calculating optimum on the avarage design

get_optimal_on_the_average_design<- function(theta_vals,prob_vals,del,n,sigma,q)
{
  des_0 <-c(0.1)
  lb <- c(0)
  ub <- c(0.35)
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
  opts <- list( "algorithm"= "NLOPT_GN_ISRES",
                "xtol_rel"= 1.0e-15,
                "maxeval"= 20000,
                "local_opts" = local_opts,
                "print_level" = 0 )
  
  # Computing the LFC
  
  res <- nloptr ( x0 = des_0,
                  theta_vals = theta_vals,
                  prob_vals = prob_vals,
                  del=del,
                  n=n,
                  sigma=sigma,
                  q=q,
                  eval_f = B,
                  lb = lb,
                  ub = ub,
                  #eval_g_eq = const1,
                  #eval_g_ineq = eval_g0,
                  opts = opts
  )
  return(res$solution)
  
}

#Step 1

step_1<- function(H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,del,n,sigma,q){

  psi_min <- 100
  #stopping=true
  
  # Lower and upper bounds
  lb <- c(0,-del)
  ub <- c(2*del,del)
  
  the_0 <-c(0.0,0.0)
  eval_g0 <- function(mu,design,del,n,sigma,q) {
    
    g1 <- mu[1]-mu[2]-del
    g2 <- mu[2]-mu[1]-del
    return( c(g1, g2) )
  }
  
  # Set optimization options.
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
  opts <- list( "algorithm"= "NLOPT_GN_ISRES",
                "xtol_rel"= 1.0e-15,
                "maxeval"= 20000,
                "local_opts" = local_opts,
                "print_level" = 0 )
  
  
  res1 <- nloptr (x0 = the_0,
                  design = design_l,
                  del=del,
                  n=n,
                  sigma=sigma,
                  q=q,
                  eval_f = psi1,
                  lb = lb,
                  ub = ub,
                  eval_g_ineq = eval_g0,
                  opts = opts
  )
  
  
  theta_l <- res1$solution
  this_psi_min <-  psi1(theta_l,design_l,del,n,sigma,q)
  
  
  # Stopping Criterion
  
  if(this_psi_min < psi_min)
  {
    psi_min <- this_psi_min
  }
  
  iteration_num
  iteration_num <- iteration_num+1
  print('psi_min')
  print(psi_min)
  print('B_pi_l')
  print(B_pi_l)
  
  #if (psi_min < -B_pi_l)
  print("abs(psi_min+B_pi_l)")
  print(abs(psi_min+B_pi_l))
  if(psi_min < -B_pi_l){
    stopping <- 0}
  else{
    stopping<- 1}
  
  
  if (stopping==1){
    print("Found minmax design")
    print(design_l)
    print("Least Favourable Distribution")
    print("probs")
    print(prob_vals_l)
    print("thetas")
    print(theta_vals_l)}
  else{
    step_3(H,h_grid_space,theta_l,theta_vals_l,prob_vals_l,design_l,B_pi_l,del,n,sigma,q)
  }        
}

#Step 2 & Step 3

step_3<- function(H,h_grid_space,theta_l,theta_vals_l,prob_vals_l,design_l,B_pi_l,del,n,sigma,q)
{
  design_l1 <- design_l
  smallest_B_pi_t_l1 <- -100
  theta_vals_l1 <- rbind(theta_vals_l,theta_l)
  prob_vals_l1 <- c(prob_vals_l,0)
  
  delta_l <- rep(0,length(prob_vals_l1))
  delta_l[length(prob_vals_l1)]<-1
  
  print('delta')
  print(delta_l)
  prob_vals_t_l1<- matrix(0, length(H), length(delta_l))
  for (t_val in 1:length(H)){
    #new T priors
    prob_vals_t_l1[t_val,] =  ((1-H[t_val])*prob_vals_l1)+(H[t_val]*delta_l)
  }
  
  # Assign l+1 values to l.
  
  for(i in 1:length(H)){
    design_t_l1 <- get_optimal_on_the_average_design(theta_vals_l1, prob_vals_t_l1[i,],del,n,sigma,q)
    B_pi_t_l1 <- B(design_t_l1, theta_vals_l1,prob_vals_t_l1[i,],del,n,sigma,q)  
    if (-smallest_B_pi_t_l1 > -B_pi_t_l1){
      smallest_B_pi_t_l1 <- B_pi_t_l1
      design_l1 <- design_t_l1
      prob_vals_l1 <- prob_vals_t_l1[i,]
    }
  }
  print("theta_vals_l1")
  print(theta_vals_l1)
  print("prob_vals_l1")
  print(prob_vals_l1)
  print("design_l1")
  print(design_l1)
  step_4(H,h_grid_space,theta_l,theta_vals_l1,prob_vals_l1,design_l1,B_pi_l,del,n,sigma,q)
}

# Step 4

step_4<- function(H,h_grid_space,theta_l,theta_vals_l1,prob_vals_l1,design_l1,B_pi_l,del,n,sigma,q){
  # STEP 3
  B_pi_l1 <- B(design_l1, theta_vals_l1, prob_vals_l1,del,n,sigma,q)
  print('B_pi_l1-B_pi_l')
  print(B_pi_l1-B_pi_l)
  
  # Assigning B(l+1) values to B(l)  
  
  if (-B_pi_l1 < -B_pi_l){
    B_pi_l <- B_pi_l1
    print("B_pi_l")
    print(B_pi_l)
    
    design_l <- design_l1
    print("design_l")
    print(design_l)
    
    #"Assigned l1 to l"
    theta_vals_l <- theta_vals_l1
    prob_vals_l <- prob_vals_l1
    
    # go to step 1
    step_1(H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,del,n,sigma,q)
  }
  else{
    h_grid_space <- h_grid_space/2
    H <-seq(0, 1, by=h_grid_space)
    # go to step 3
    step_3(H,h_grid_space,theta_l,theta_vals_l,prob_vals_l,design_l,B_pi_l,del,n,sigma,q)
  }
}

N <- 140                    # Total No. of subjects
n <- sqrt(N)
sigma <- 3.141836           # standard deviation
al <- 0.05/16               # significance level alpha
q <- qnorm(1-al)           # critical value
del<- 2.7                  # delta
iteration_num<-0
B_pi_l <- 100             # Initialization of B_pi
h_grid_space <- 0.5       # grid length
H<- seq(0, 1, by=h_grid_space)      # grid vector


print( "Initial Prior on Theta : ")

# Theta values (LFC's) are their corresponding prior probabilities


the_1<- c(0,0)
the_2<- c(del/3,del/2)
theta_vals_l<- rbind(the_1,the_2)
prob_vals_l <- c(0,1)

# Starting optimum on the average design and corresponding power

print( "Starting Design")
design_l <- get_optimal_on_the_average_design(theta_vals_l,prob_vals_l,del,n,sigma,q)
B_pi_l <- B(design_l,theta_vals_l, prob_vals_l,del,n,sigma,q)

step_1(H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,del,n,sigma,q)


toc()