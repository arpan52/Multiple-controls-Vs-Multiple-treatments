#################################################################################################
###### H-algorithm for min-max design for a Complete bipartite graph for K1=2 and K2=3 for UIT ##
#################################################################################################

#install.packages("nloptr")
#install.packages("mvtnorm")


library('nloptr')
library(mvtnorm)
library(tictoc)
tic()


rm(list = ls())
# Initializations

m<-rep(0,5)
# Function for calculating the mean and covariance matrix of Z-ij's and criterion 'Psi'.

psi<- function(design,theta,del,n,sigma,q)
{
  
  m[1]<-0
  m[2]<-theta[1]
  m[3]<-del
  m[4]<- theta[2]
  m[5]<- theta[2]
  
  t1<- design
  t2<- (1-(2*t1))/3
  
  de13 <- ((m[1]-m[3])*n)/sigma
  de14 <- ((m[1]-m[4])*n)/sigma
  de15 <- ((m[1]-m[5])*n)/sigma
  de23 <- ((m[2]-m[3])*n)/sigma
  de24 <- ((m[2]-m[4])*n)/sigma
  de25 <- ((m[2]-m[5])*n)/sigma
  
  
  
  a1 <- t2*t2
  b1 <- (t1+t2)*(t1+t2)
  rho1 <- sqrt(a1/b1)
  
  a2 <- t2*t2
  b2 <- (t1+t2)*(t1+t2)
  rho2 <- sqrt(a2/b2)
  
  a3 <- t1*t1
  b3 <- (t1+t2)*(t1+t2)
  rho3 <- sqrt(a3/b3)
  
  rho4<- 0
  rho5<- 0
  
  a6 <- t2*t2
  b6 <- (t1+t2)*(t1+t2)
  rho6 <- sqrt(a6/b6)
  
  rho7<- 0
  
  a8 <- t1*t1
  b8<- (t1+t2)*(t1+t2)
  rho8<- sqrt(a8/b8)
  
  rho9<-0
  rho10<- 0
  rho11<-0
  
  
  a12 <- t1*t1
  b12<- (t1+t2)*(t1+t2)
  rho12<- sqrt(a12/b12)
  
  
  a13 <- t2*t2
  b13<- (t1+t2)*(t1+t2)
  rho13<- sqrt(a13/b13)
  
  a14 <- t2*t2
  b14<- (t1+t2)*(t1+t2)
  rho14<- sqrt(a14/b14)
  
  a15 <- t2*t2
  b15<- (t1+t2)*(t1+t2)
  rho15<- sqrt(a15/b15)
  
  
  mu1 <- (sqrt((t1*t2)/(t1 + t2)))*de13
  mu2 <- (sqrt((t1*t2)/(t1 + t2)))*de14
  mu3 <- (sqrt((t1*t2)/(t1 + t2)))*de15
  mu4 <- (sqrt((t1*t2)/(t1 + t2)))*de23
  mu5 <- (sqrt((t1*t2)/(t1 + t2)))*de24
  mu6 <- (sqrt((t1*t2)/(t1 + t2)))*de25
  
  meann <- c(mu1,mu2,mu3,mu4,mu5,mu6)
  S <-  rbind(c(1,rho1,rho2,rho3,rho4,rho5),c(rho1,1,rho6,rho7,rho8,rho9),c(rho2,rho6,1,rho10,rho11,rho12),c(rho3,rho7,rho10,1,rho13,rho14),c(rho4,rho8,rho11,rho13,1,rho15),c(rho5,rho9,rho12,rho14,rho15,1))
  
  lw1 <- c(-q,-q,-q,-q,-q,-q)
  up1 <- c(q,q,q,q,q,q)
  
  return(-(1-pmvnorm(mean=meann, sigma=S, lower=lw1, upper=up1)))   #Power of UIT
  
  
}
psi1<- function(theta,design,del,n,sigma,q){
  
  # Function to compute LFCS
  
  
  m[1]<-0
  m[2]<-theta[1]
  m[3]<-del
  m[4]<- theta[2]
  m[5]<- theta[2]
  
  return(-psi(design,mu,del,n,sigma,q))
  
}
# Function for calculating the average power B

B<- function(design,theta_vals,prob_vals,del,n,sigma,q){
  
  s <- 0
  for (i in 1:length(prob_vals)){
    s <- s + prob_vals[i]* psi(design,theta_vals[i,],del,n,sigma,q)   # average power associated to LFCs and Probabilities
  }
  ret_val<-s
  return(ret_val[1])
}

# Function for compute optimum on the avarage design

get_optimal_on_the_average_design<- function(theta_vals,prob_vals,del,n,sigma,q)
{
  des_0 <-c(0.3)
  lb <- c(0)
  ub <- c(0.5)
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
  opts <- list( "algorithm"= "NLOPT_GN_ISRES",
                "xtol_rel"= 1.0e-15,
                "maxeval"= 20000,
                "local_opts" = local_opts,
                "print_level" = 0 )
  
  
  
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
  
  # Computing the LFC
  
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

N <- 200                    # Total No. of subjects
n <- sqrt(N)
sigma <- 1                 # Standard deviation
al <- 0.05/12               # significance level alpha
q <- qnorm(1-al)           # critical value
del<-0.7                   # delta
iteration_num<-0
B_pi_l <- 100             # Initialization of B_pi
h_grid_space <- 0.25       # grid length
H<- seq(0, 1, by=h_grid_space)      # grid vector


print( "Initial Prior on Theta : ")

# Theta values (LFC's) are their corresponding prior probabilities


the_1<- c(0,0)
the_2<- c(del/2,del/2)
theta_vals_l<- rbind(the_1,the_2)
prob_vals_l <- c(1/2,1/2)

# Starting optimum on the average design and corresponding power

print( "Starting Design")
design_l <- get_optimal_on_the_average_design(theta_vals_l,prob_vals_l,del,n,sigma,q)
B_pi_l <- B(design_l,theta_vals_l, prob_vals_l,del,n,sigma,q)

step_1(H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,del,n,sigma,q)


toc()