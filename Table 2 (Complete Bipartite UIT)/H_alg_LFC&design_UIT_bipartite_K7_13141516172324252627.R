#################################################################################################
###### H-algorithm for min-max design for a Complete bipartite graph for K1=2 and K2=5 for UIT ##
#################################################################################################

#install.packages("nloptr")
#install.packages("mvtnorm")


library('nloptr')
library(mvtnorm)
library(tictoc)
tic()

# Initializations
rm(list = ls())
mu<-rep(0,7)
# Function for calculating the mean and covariance matrix of Z-ij's and criterion 'Psi'.

psi<- function(design,theta,del,n,sigma,q)
{
  
  
  mu[1]<-0
  mu[2]<-theta[1]
  mu[3]<-del
  mu[4]<-theta[2]
  mu[5]<-theta[2]
  mu[6]<-theta[2]
  mu[7]<-theta[2]
  
  
  
  t1<- design
  t2<- (1-2*design)/5
  
  
  
  de13 <- ((mu[1]-mu[3])*n)/sigma
  de14 <- ((mu[1]-mu[4])*n)/sigma
  de15 <- ((mu[1]-mu[5])*n)/sigma
  de16 <- ((mu[1]-mu[6])*n)/sigma
  de17 <- ((mu[1]-mu[7])*n)/sigma
  de23 <- ((mu[2]-mu[3])*n)/sigma
  de24 <- ((mu[2]-mu[4])*n)/sigma
  de25 <- ((mu[2]-mu[5])*n)/sigma
  de26 <- ((mu[2]-mu[6])*n)/sigma
  de27 <- ((mu[2]-mu[7])*n)/sigma
  
  a1 <- t2*t2
  b1 <- (t1+t2)*(t1+t2)
  rho1 <- sqrt(a1/b1)
  
  a2 <- t2*t2
  b2 <- (t1+t2)*(t1+t2)
  rho2 <- sqrt(a2/b2)
  
  a3 <- t2*t2
  b3 <- (t1+t2)*(t1+t2)
  rho3 <- sqrt(a3/b3)
  
  a4 <- t2*t2
  b4 <- (t1+t2)*(t1+t2)
  rho4 <- sqrt(a4/b4)
  
  a5 <- t1*t1
  b5 <- (t1+t2)*(t1+t2)
  rho5 <- sqrt(a5/b5)
  
  rho6<- 0
  rho7<- 0
  rho8<- 0
  rho9<- 0
  
  a10 <- t2*t2
  b10 <- (t1+t2)*(t1+t2)
  rho10 <- sqrt(a10/b10)
  
  a11 <- t2*t2
  b11 <- (t1+t2)*(t1+t2)
  rho11 <- sqrt(a11/b11)
  
  a12 <- t2*t2
  b12 <- (t1+t2)*(t1+t2)
  rho12 <- sqrt(a12/b12)
  
  rho13<-0
  
  
  a14 <- t1*t1
  b14 <- (t1+t2)*(t1+t2)
  rho14 <- sqrt(a14/b14)
  
  rho15<- 0
  rho16<- 0
  rho17<- 0
  
  
  a18 <- t2*t2
  b18 <- (t1+t2)*(t1+t2)
  rho18 <- sqrt(a18/b18)
  
  a19 <- t2*t2
  b19 <- (t1+t2)*(t1+t2)
  rho19 <- sqrt(a19/b19)
  
  rho20<-0
  rho21<-0
  
  a22 <- t1*t1
  b22 <- (t1+t2)*(t1+t2)
  rho22 <- sqrt(a22/b22)
  
  rho23<-0
  rho24<-0
  
  a25 <- t2*t2
  b25 <- (t1+t2)*(t1+t2)
  rho25 <- sqrt(a25/b25)
  
  rho26<-0
  rho27<-0
  rho28<-0
  
  a29 <- t1*t1
  b29 <- (t1+t2)*(t1+t2)
  rho29 <- sqrt(a29/b29)
  
  rho30<-0
  
  rho31<-0
  rho32<-0
  rho33<-0
  rho34<-0
  
  a35 <- t1*t1
  b35 <- (t1+t2)*(t1+t2)
  rho35 <- sqrt(a35/b35)
  
  
  a36 <- t2*t2
  b36 <- (t1+t2)*(t1+t2)
  rho36 <- sqrt(a36/b36)
  
  a37 <- t2*t2
  b37 <- (t1+t2)*(t1+t2)
  rho37 <- sqrt(a37/b37)
  
  a38 <- t2*t2
  b38 <- (t1+t2)*(t1+t2)
  rho38 <- sqrt(a38/b38)
  
  a39 <- t2*t2
  b39 <- (t1+t2)*(t1+t2)
  rho39 <- sqrt(a39/b39)
  
  a40 <- t2*t2
  b40 <- (t1+t2)*(t1+t2)
  rho40 <- sqrt(a40/b40)
  
  a41 <- t2*t2
  b41 <- (t1+t2)*(t1+t2)
  rho41 <- sqrt(a41/b41)
  
  a42 <- t2*t2
  b42 <- (t1+t2)*(t1+t2)
  rho42 <- sqrt(a42/b42)
  
  a43 <- t2*t2
  b43 <- (t1+t2)*(t1+t2)
  rho43 <- sqrt(a43/b43)
  
  a44 <- t2*t2
  b44 <- (t1+t2)*(t1+t2)
  rho44 <- sqrt(a44/b44)
  
  a45 <- t2*t2
  b45 <- (t1+t2)*(t1+t2)
  rho45 <- sqrt(a45/b45)
  
  
  
  mu1 <- (sqrt((t1*t2)/(t1 + t2)))*de13
  mu2 <- (sqrt((t1*t2)/(t1 + t2)))*de14
  mu3 <- (sqrt((t1*t2)/(t1 + t2)))*de15
  mu4 <- (sqrt((t1*t2)/(t1 + t2)))*de16
  mu5 <- (sqrt((t1*t2)/(t1 + t2)))*de17
  mu6 <- (sqrt((t1*t2)/(t1 + t2)))*de23
  mu7 <- (sqrt((t1*t2)/(t1 + t2)))*de24
  mu8 <- (sqrt((t1*t2)/(t1 + t2)))*de25
  mu9 <- (sqrt((t1*t2)/(t1 + t2)))*de26
  mu10 <- (sqrt((t1*t2)/(t1 + t2)))*de27
  
  m <- c(mu1,mu2,mu3,mu4,mu5,mu6,mu7,mu8,mu9,mu10)
  S <-  rbind(c(1,rho1,rho2,rho3,rho4,rho5,rho6,rho7,rho8,rho9),c(rho1,1,rho10,rho11,rho12,rho13,rho14,rho15,rho16,rho17),c(rho2,rho10,1,rho18,rho19,rho20,rho21,rho22,rho23,rho24),c(rho3,rho11,rho18,1,rho25,rho26,rho27,rho28,rho29,rho30),c(rho4,rho12,rho19,rho25,1,rho31,rho32,rho33,rho34,rho35),c(rho5,rho13,rho20,rho26,rho31,1,rho36,rho37,rho38,rho39),c(rho6,rho14,rho21,rho27,rho32,rho36,1,rho40,rho41,rho42),c(rho7,rho15,rho22,rho28,rho33,rho37,rho40,1,rho43,rho44),c(rho8,rho16,rho23,rho29,rho34,rho38,rho41,rho43,1,rho45),c(rho9,rho17,rho24,rho30,rho35,rho39,rho42,rho44,rho45,1))
  
  lw1 <- c(-q,-q,-q,-q,-q,-q,-q,-q,-q,-q)
  up1 <- c(q,q,q,q,q,q,q,q,q,q)
  
  return(-(1-pmvnorm(mean=m, sigma=S, lower=lw1, upper=up1)))   #Power of UIT
  
}


psi1<- function(theta,design,del,n,sigma,q)
{
  # Function to compute LFCS
  
  mu[1]<-0
  mu[2]<-theta[1]
  mu[3]<-del
  mu[4]<-theta[2]
  mu[5]<-theta[2]
  mu[6]<-theta[2]
  mu[7]<-theta[2]
  
  return(-psi(design,mu,del,n,sigma,q))
  
}


# Function for compute the average power B

B<- function(design,theta_vals,prob_vals,del,n,sigma,q){
  
  s <- 0
  for (i in 1:length(prob_vals)){
    s <- s + prob_vals[i]* psi(design,theta_vals[i,],del,n,sigma,q)            # average power associated to LFCs and Probabilities
  }
  ret_val<-s
  return(ret_val[1])
}

# Function for calculating optimum on the avarage design

get_optimal_on_the_average_design<- function(theta_vals,prob_vals,del,n,sigma,q)
{
  des_0 <-c(0.1)

  lb <- c(0)
  ub <- c(0.4)
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
  print(' theta_l')
  print( theta_l)
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
sigma <- 1                  # standard deviation
al <- 0.05/20               # significance level alpha
q <- qnorm(1-al)           # critical value
del<-1                   # delta
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
