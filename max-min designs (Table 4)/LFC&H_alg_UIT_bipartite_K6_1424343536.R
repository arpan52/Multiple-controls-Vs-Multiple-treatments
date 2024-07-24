#####################################################################################
###### H-algorithm for min-max design for a general bipartite graph for K=6 for UIT##
###### Pairs =  (1,4),(2,4),(3,4),(3,5),(3,6)                                      ##
#####################################################################################

#install.packages("nloptr")
#install.packages("mvtnorm")


library('nloptr')
library(mvtnorm)
library(tictoc)
tic()


rm(list = ls())
# Initializations

x<-rep(0,6)
mu<-rep(0,6)
# Function for calculating the mean and covariance matrix of Z-ij's and criterion 'Psi'.

psi<- function(design,mu,del,n,sigma,q)
{
  
  # # Calculating mean and covariance matrix of Z_ij's
  x[1]<- design[1]
  x[2]<- design[2]
  x[3]<- design[3]
  x[4]<- design[3]
  x[5]<- design[2]
  x[6]<- design[1]
  
  de14 <- ((mu[1]-mu[4])*n)/sigma
  de24 <- ((mu[2]-mu[4])*n)/sigma
  de34 <- ((mu[3]-mu[4])*n)/sigma
  de35 <- ((mu[3]-mu[5])*n)/sigma
  de36 <- ((mu[3]-mu[6])*n)/sigma
  
  
  
  
  
  
  a1 <- x[1]*x[2]
  b1 <- (x[1]+x[4])*(x[2]+x[4])
  rho1 <- sqrt(a1/b1)
  
  
  a2 <- x[1]*x[3]
  b2 <- (x[1]+x[4])*(x[3]+x[4])
  rho2 <- sqrt(a2/b2)
  
  rho3<- 0
  rho4<- 0
  
  a5 <- x[2]*x[3]
  b5 <- (x[2]+x[4])*(x[3]+x[4])
  rho5 <- sqrt(a5/b5)
  
  rho6<- 0
  rho7<- 0
  
  a8 <- x[4]*x[5]
  b8 <- (x[3]+x[4])*(x[3]+x[5])
  rho8 <- sqrt(a8/b8)
  
  a9 <- x[4]*x[6]
  b9 <- (x[3]+x[4])*(x[3]+x[6])
  rho9 <- sqrt(a9/b9)
  
  
  a10 <- x[5]*x[6]
  b10 <- (x[3]+x[5])*(x[3]+x[6])
  rho10 <- sqrt(a10/b10)
  
  
  mu1 <- (sqrt((x[1]*x[4])/(x[1] + x[4])))*de14
  mu2 <- (sqrt((x[2]*x[4])/(x[2] + x[4])))*de24
  mu3 <- (sqrt((x[3]*x[4])/(x[3] + x[4])))*de34
  mu4 <- (sqrt((x[3]*x[5])/(x[3] + x[5])))*de35
  mu5 <- (sqrt((x[3]*x[6])/(x[3] + x[6])))*de36
  
  
  
  
  m <- c(mu1,mu2,mu3,mu4,mu5)
  
  S <-  rbind(c(1,rho1,rho2,rho3,rho4),c(rho1,1,rho5,rho6,rho7), c(rho2,rho5,1,rho8,rho9),c(rho3,rho6,rho8,1,rho10),c(rho4,rho7,rho9,rho10,1))
  
  
  lw1 <- c(-q,-q,-q,-q,-q)
  up1 <- c(q,q,q,q,q)
  
  
  return(-(1-pmvnorm(mean=m, sigma=S, lower=lw1, upper=up1)))   #Power
}  

psi1<- function(theta,design,del,n,sigma,q){
  # Function to compute LFC
  
  # # Calculating mean and covariance matrix of Z_ij's
  
  mu[1]<- theta[1]
  mu[2]<- theta[2]
  mu[3]<- 0
  mu[4]<- del
  mu[5]<- theta[2]
  mu[6]<- theta[1]
  
  return(-psi(design,mu,del,n,sigma,q))
 
}


B<- function(design,theta_vals,prob_vals,del,n,sigma,q){
  
  s <- 0
  for (i in 1:length(prob_vals)){
    s <- s + prob_vals[i]* psi(design,theta_vals[i,],del,n,sigma,q)  # Average power for LFCs and probabilities
  }
  ret_val<-s
  return(ret_val[1])
}

# Function for computing optimum on the avarage design

get_optimal_on_the_average_design<- function(theta_vals,prob_vals,del,n,sigma,q)
{
  des_0 <-c(0.1,0.1,0.25)
  #options = optimoptions('fmincon','Display','none')
  const1 <- function(design,theta_vals,prob_vals,del,n,sigma,q) {
    return( (2*design[1]+2*design[2]+2*design[3]-1))
  }
  
  const2 <- function(design,theta_vals,prob_vals,del,n,sigma,q) {
    return( c(design[1]-design[3],design[2]-design[3]))
  }
  
  
  lb <- c(0,0,0)
  ub <- c(0.35,0.35,0.35)
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-7 )
  opts <- list( "algorithm"= "NLOPT_GN_ISRES",
                "xtol_rel"= 1.0e-7,
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
                  eval_g_eq = const1,
                  eval_g_ineq = const2,
                  opts = opts
  )
  return(res$solution)
  
}

#Step 1

step_1<- function(H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,del,n,sigma,q){

  psi_min <- 100
  #stopping=true
  
  theta_vals_l<-matrix(0,length(prob_vals_l),6)
  
  
  
  # Lower and upper bounds
  theta_0 <-c(0,0)
  lb <- c(0,0)
  ub <- c(del,del)
  
  # Set optimization options.
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-7 )
  opts <- list( "algorithm"= "NLOPT_GN_ISRES",
                "xtol_rel"= 1.0e-7,
                "maxeval"= 20000,
                "local_opts" = local_opts,
                "print_level" = 0 )
  
  
  res1 <- nloptr (x0 = theta_0,
                  design = design_l,
                  del=del,
                  n=n,
                  sigma=sigma,
                  q=q,
                  eval_f = psi1,
                  lb = lb,
                  ub = ub,
                  opts = opts
  )
  
  
  th1 <- res1$solution
  print('th1')
  print(th1)
  
  theta_vals_l[1,]<- c(th1[1],th1[2],0,del,th1[2],th1[1])
  theta_vals_l[2,] <- c(del,0,0,0,0,0)
  theta_vals_l[3,] <- c(0,del,0,0,0,0)
  theta_vals_l[4,] <- c(0,0,0,0,del,0)
  theta_vals_l[5,] <- c(0,0,0,0,0,del)
  
  print("theta_vals_l")
  print(theta_vals_l)    # Values of LFCs
  
  
  this_psi<- rep(0,length(prob_vals_l))
  for (i in 1:length(prob_vals_l)) {
    this_psi[i] <-  -psi(design_l,theta_vals_l[i,],del,n,sigma,q)
  }
  
  this_psi<- as.numeric(round(this_psi, 2))
  this_psi_min<- min(this_psi)
  
  print("this_psi")
  print(this_psi)
  index <- which(this_psi == this_psi_min)
  
  
  
  
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
  if((psi_min < -B_pi_l)&(abs(psi_min+B_pi_l)>0.004)){
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
    step_3(index,H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,del,n,sigma,q)
  }        
}

#Step 2 & Step 3

step_3<- function(index,H,h_grid_space,theta_vals_l1,prob_vals_l,design_l,B_pi_l,del,n,sigma,q)
{
  design_l1 <- design_l
  smallest_B_pi_t_l1 <- -100
  
  
  prob_vals_t_l1<- matrix(0, length(H), length(prob_vals_l))
  delta_l <- rep(0,length(prob_vals_l))
  for(i in 1:length(prob_vals_l)){
    delta_l[index[i]] <- 1/length(index)
  }
  print(delta_l)
  for (t_val in 1:length(H)){
    #new T priors
    prob_vals_t_l1[t_val,] =  ((1-H[t_val])*prob_vals_l)+(H[t_val]*delta_l)
  }
  # Assign l+1 values to l.
  
  for(i in 1:length(H)){
    design_t_l1 <- get_optimal_on_the_average_design(theta_vals_l1, prob_vals_t_l1[i,],del,n,sigma,q)
    B_pi_t_l1 <- B(design_t_l1,theta_vals_l1,prob_vals_t_l1[i,],del,n,sigma,q)  
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
  step_4(H,h_grid_space,theta_vals_l1,prob_vals_l1,design_l1,B_pi_l,del,n,sigma,q)
}

# Step 4

step_4<- function(H,h_grid_space,theta_vals_l1,prob_vals_l1,design_l1,B_pi_l,del,n,sigma,q){
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
    
    B_pi_l <- B_pi_l1
    print("B_pi_l")
    print(B_pi_l)
    
    design_l <- design_l1
    print("design_l")
    print(design_l)
    
    #"Assigned l1 to l"
    prob_vals_l <- prob_vals_l1
    
    #stopping=true
    
    theta_vals_l<-matrix(0,length(prob_vals_l),6)
    
    
    
    # Lower and upper bounds
    theta_0 <-c(0,0)
    lb <- c(0,0)
    ub <- c(del,del)
    
    # Set optimization options.
    local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-7 )
    opts <- list( "algorithm"= "NLOPT_GN_ISRES",
                  "xtol_rel"= 1.0e-7,
                  "maxeval"= 20000,
                  "local_opts" = local_opts,
                  "print_level" = 0 )
    
    
    res1 <- nloptr (x0 = theta_0,
                    design = design_l,
                    del=del,
                    n=n,
                    sigma=sigma,
                    q=q,
                    eval_f = psi1,
                    lb = lb,
                    ub = ub,
                    opts = opts
    )
    
    
    th1 <- res1$solution
    print('th1')
    print(th1)
    
    theta_vals_l[1,]<- c(th1[1],th1[2],0,del,th1[2],th1[1])
    theta_vals_l[2,] <- c(del,0,0,0,0,0)
    theta_vals_l[3,] <- c(0,del,0,0,0,0)
    theta_vals_l[4,] <- c(0,0,0,0,del,0)
    theta_vals_l[5,] <- c(0,0,0,0,0,del)
    
    print("theta_vals_l")
    print(theta_vals_l)    # The values of LFCs
    
    
    this_psi<- rep(0, length(prob_vals_l))
    for (i in 1:length(prob_vals_l)) {
      this_psi[i] <-  -psi(design_l,theta_vals_l[i,],del,n,sigma,q)
    }
    
    this_psi<- as.numeric(round(this_psi, 2))
    this_psi_min<- min(this_psi)
    
    print("this_psi")
    print(this_psi)
    index <- which(this_psi == this_psi_min)
    
    step_3(index,H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,del,n,sigma,q)
  }
}

N <- 57                    # Total No. of subjects
n <- sqrt(N)
sigma <-5.79                # Standard deviation
al <- 0.05/10               # significance level alpha
q <- qnorm(1-al)           # critical value
del<-4.3                   # delta
iteration_num<-0
B_pi_l <- 100             # Initialization of B_pi
h_grid_space <- 0.25       # grid length
H<- seq(0, 1, by=h_grid_space)      # grid vector


print( "Initial Prior on Theta : ")

# Theta values (LFC's) are their corresponding prior probabilities

the_1<-  c(0,0,0,del,0,0)
the_2<-  c(del,0,0,0,0,0)
the_3<- c(0,del,0,0,0,0)
the_4<- c(0,0,0,0,del,0)
the_5<-  c(0,0,0,0,0,del)

theta_vals_l<- rbind(the_1,the_2,the_3,the_4,the_5)



prob_vals_l <- c(1/5,1/5,1/5,1/5,1/5)

# Starting optimum on the average design and corresponding power

print( "Starting Design")
design_l <- get_optimal_on_the_average_design(theta_vals_l,prob_vals_l,del,n,sigma,q)
B_pi_l <- B(design_l,theta_vals_l, prob_vals_l,del,n,sigma,q)

step_1(H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,del,n,sigma,q)


toc()
