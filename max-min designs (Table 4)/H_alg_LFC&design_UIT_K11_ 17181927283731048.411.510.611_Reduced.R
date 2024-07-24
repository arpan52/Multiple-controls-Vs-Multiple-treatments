######################################################################################
###### H-algorithm for min-max design for a general bipartite graph for K=11 for UIT##
###### Pairs = (1,7),(1,8),(1,9),(2,7),(2,8),(3,7),(3,10),(4,8),(4,11),(5,10),(6,11)##
#####################################################################################

#install.packages("nloptr")
#install.packages("mvtnorm")


library('nloptr')
library(mvtnorm)
library(tictoc)
tic()


rm(list = ls())

# Initializations
mu<-rep(0,11)
x<-rep(0,11)

# Function for calculating the mean and covariance matrix of Z-ij's and criterion 'Psi'.

psi<- function(design,mu,del,n,sigma,q)
{
  
  # # Calculating mean and covariance matrix of Z_ij's
  x[1]<- design[1]
  x[2]<- design[2]
  x[3]<- design[3]
  x[4]<- design[3]
  x[5]<- design[4]
  x[6]<- design[4]
  x[7]<- design[5]
  x[8]<- design[5]
  x[9]<- design[6]
  x[10]<- design[7]
  x[11]<- design[7]
  
  de17 <- ((mu[1]-mu[7])*n)/sigma
  de18 <- ((mu[1]-mu[8])*n)/sigma
  de19 <- ((mu[1]-mu[9])*n)/sigma
  de27 <- ((mu[2]-mu[7])*n)/sigma
  de28 <- ((mu[2]-mu[8])*n)/sigma
  de37 <- ((mu[3]-mu[7])*n)/sigma
  de310 <- ((mu[3]-mu[10])*n)/sigma
  de48 <- ((mu[4]-mu[8])*n)/sigma
  de411 <- ((mu[4]-mu[11])*n)/sigma
  de510 <- ((mu[5]-mu[10])*n)/sigma
  de611 <- ((mu[6]-mu[11])*n)/sigma
  
  
  
  
  
  a1 <- x[7]*x[8]
  b1 <- (x[1]+x[7])*(x[1]+x[8])
  rho1 <- sqrt(a1/b1)
  
  a2 <- x[7]*x[9]
  b2 <- (x[1]+x[7])*(x[1]+x[9])
  rho2 <- sqrt(a2/b2)
  
  a3 <- x[1]*x[2]
  b3 <- (x[1]+x[7])*(x[2]+x[7])
  rho3 <- sqrt(a3/b3)
  
  rho4 <- 0
  
  a5 <- x[1]*x[3]
  b5 <- (x[1]+x[7])*(x[3]+x[7])
  rho5 <- sqrt(a5/b5)
  
  rho6 <- 0
  
  rho7 <- 0
  
  rho8<-0
  rho9<-0
  rho10<-0
  
  
  
  a11 <- x[8]*x[9]
  b11 <- (x[1]+x[8])*(x[1]+x[9])
  rho11 <- sqrt(a11/b11)
  
  
  rho12 <- 0
  
  a13 <- x[1]*x[2]
  b13 <- (x[1]+x[8])*(x[2]+x[8])
  rho13 <- sqrt(a13/b13)
  
  
  rho14 <- 0
  
  rho15 <- 0
  
  a16 <- x[1]*x[4]
  b16 <- (x[1]+x[8])*(x[4]+x[8])
  rho16 <- sqrt(a16/b16)
  
  rho17<-0
  rho18<-0
  rho19<-0
  
  rho20 <- 0
  rho21 <- 0
  rho22 <- 0
  rho23 <- 0
  rho24 <- 0
  rho25 <- 0
  rho26 <- 0
  rho27 <- 0
  
  a28 <- x[7]*x[8]
  b28 <- (x[2]+x[7])*(x[2]+x[8])
  rho28 <- sqrt(a28/b28)
  
  a29 <- x[2]*x[3]
  b29 <- (x[2]+x[7])*(x[3]+x[7])
  rho29 <- sqrt(a29/b29)
  
  rho30 <- 0
  rho31 <- 0
  rho32 <- 0
  rho33 <- 0
  rho34 <- 0
  
  
  rho35 <- 0
  rho36 <- 0
  
  a37 <- x[2]*x[4]
  b37 <- (x[2]+x[8])*(x[4]+x[8])
  rho37 <- sqrt(a37/b37)
  
  rho38 <- 0
  rho39 <- 0
  rho40 <- 0
  
  a41 <- x[7]*x[10]
  b41 <- (x[3]+x[7])*(x[3]+x[10])
  rho41 <- sqrt(a41/b41)
  
  rho42 <- 0
  rho43 <- 0
  rho44 <- 0
  rho45 <- 0
  
  rho46 <-0
  rho47 <- 0
  
  a48 <- x[3]*x[5]
  b48 <- (x[3]+x[10])*(x[5]+x[10])
  rho48 <- sqrt(a48/b48)
  
  rho49 <- 0
  
  
  a50 <- x[8]*x[11]
  b50 <- (x[4]+x[8])*(x[4]+x[11])
  rho50 <- sqrt(a50/b50)
  
  rho51 <-0
  rho52 <- 0
  
  rho53 <- 0
  
  a54 <- x[4]*x[6]
  b54 <- (x[4]+x[11])*(x[6]+x[11])
  rho54 <- sqrt(a54/b54)
  
  rho55 <- 0
  
  mu1 <- (sqrt((x[1]*x[7])/(x[1] + x[7])))*de17
  mu2 <- (sqrt((x[1]*x[8])/(x[1] + x[8])))*de18
  mu3 <- (sqrt((x[1]*x[9])/(x[1] + x[9])))*de19
  mu4 <- (sqrt((x[2]*x[7])/(x[2] + x[7])))*de27
  mu5 <- (sqrt((x[2]*x[8])/(x[2] + x[8])))*de28
  mu6 <- (sqrt((x[3]*x[7])/(x[3] + x[7])))*de37
  mu7 <- (sqrt((x[3]*x[10])/(x[3] + x[10])))*de310
  mu8 <- (sqrt((x[4]*x[8])/(x[4] + x[8])))*de48
  mu9 <- (sqrt((x[4]*x[11])/(x[4] + x[11])))*de411
  mu10 <- (sqrt((x[5]*x[10])/(x[5] + x[10])))*de510
  mu11 <- (sqrt((x[6]*x[11])/(x[6] + x[11])))*de611
  
  
  
  m <- c(mu1,mu2,mu3,mu4,mu5,mu6,mu7,mu8,mu9,mu10,mu11)
  
  S <-  rbind(c(1,rho1,rho2,rho3,rho4,rho5,rho6,rho7,rho8,rho9,rho10),c(rho1,1,rho11,rho12,rho13,rho14,rho15,rho16,rho17,rho18,rho19),c(rho2,rho11,1,rho20,rho21,rho22,rho23,rho24,rho25,rho26,rho27),c(rho3,rho12,rho20,1,rho28,rho29,rho30,rho31,rho32,rho33,rho34),c(rho4,rho13,rho21,rho28,1,rho35,rho36,rho37,rho38,rho39,rho40),c(rho5,rho14,rho22,rho29,rho35,1,rho41,rho42,rho43,rho44,rho45),c(rho6,rho15,rho23,rho30,rho36,rho41,1,rho46,rho47,rho48,rho49),c(rho7,rho16,rho24,rho31,rho37,rho42,rho46,1,rho50,rho51,rho52),c(rho8,rho17,rho25,rho32,rho38,rho43,rho47,rho50,1,rho53,rho54),c(rho9,rho18,rho26,rho33,rho39,rho44,rho48,rho51,rho53,1,rho55),c(rho10,rho19,rho27,rho34,rho40,rho45,rho49,rho52,rho54,rho55,1))
  
  
  lw1 <- c(-q,-q,-q,-q,-q,-q,-q,-q,-q,-q,-q)
  up1 <- c(q,q,q,q,q,q,q,q,q,q,q)
  
  
  return(-(1-pmvnorm(mean=m, sigma=S, lower=lw1, upper=up1)))   #Power of UIT
}  

psi5<- function(theta,design,del,n,sigma,q){
  # LFC associated to delta difference at (3,7)
  
  # Calculating mean and covariance matrix of Z_ij's
  
  mu[1]<- theta[1]
  mu[2]<- theta[2]
  mu[3]<- 0
  mu[4]<- theta[3]
  mu[5]<- theta[4]
  mu[6]<- theta[4]
  mu[7]<- del
  mu[8]<- theta[5]
  mu[9]<- theta[1]
  mu[10]<- theta[6]
  mu[11]<- theta[6]
  return(-psi(design,mu,del,n,sigma,q))
  
}

psi6<- function(theta,design,del,n,sigma,q){
  # LFC associated to delta difference at (3,10)
  
  
  # Calculating mean and covariance matrix of Z_ij's
  
  mu[1]<- theta[1]
  mu[2]<- theta[2]
  mu[3]<- 0
  mu[4]<- theta[3]
  mu[5]<- theta[4]
  mu[6]<- theta[4]
  mu[7]<- theta[5]
  mu[8]<- theta[5]
  mu[9]<- theta[1]
  mu[10]<- del
  mu[11]<- theta[6]
  
  return(-psi(design,mu,del,n,sigma,q))
  
}

psi9<- function(theta,design,del,n,sigma,q){
  
  # LFC associated to delta difference at (5,10)
  
  
  # Calculating mean and covariance matrix of Z_ij's
  
  mu[1]<- theta[1]
  mu[2]<- theta[2]
  mu[3]<- theta[3]
  mu[4]<- theta[4]
  mu[5]<- 0
  mu[6]<- theta[5]
  mu[7]<- theta[6]
  mu[8]<- theta[6]
  mu[9]<- theta[1]
  mu[10]<- del
  mu[11]<- theta[7]
  
  return(-psi(design,mu,del,n,sigma,q))
  
}

psi10<- function(theta,design,del,n,sigma,q){
  
  # LFC associated to delta difference at (6,11)
  
  
  # Calculating mean and covariance matrix of Z_ij's
  
  mu[1]<- theta[1]
  mu[2]<- theta[2]
  mu[3]<- theta[3]
  mu[4]<- theta[4]
  mu[5]<- theta[5]
  mu[6]<- 0
  mu[7]<- theta[6]
  mu[8]<- theta[6]
  mu[9]<- theta[1]
  mu[10]<- theta[7]
  mu[11]<- del
  
  return(-psi(design,mu,del,n,sigma,q))
  
}


B<- function(design,theta_vals,prob_vals,del,n,sigma,q){
  
  s <- 0
  for (i in 1:length(prob_vals)){
    s <- s + prob_vals[i]* psi(design,theta_vals[i,],del,n,sigma,q)
  }
  ret_val<-s
  return(ret_val[1])
}

# Function for computing optimum on the avarage design

get_optimal_on_the_average_design<- function(theta_vals,prob_vals,del,n,sigma,q)
{
  des_0 <- c(0.21436541,0.13570508,0.01503625,0.13115687,0.08633439,0.05772650,0.06357303)
  #options = optimoptions('fmincon','Display','none')
  const1 <- function(design,theta_vals,prob_vals,del,n,sigma,q) {
    return( (design[1]+design[2]+2*design[3]+2*design[4]+2*design[5]+design[6]+2*design[7]-1))
  }
  
  const2 <- function(design,theta_vals,prob_vals,del,n,sigma,q) {
    return( c(design[2]-design[1],design[3]-design[1],design[4]-design[1],design[4]-design[2],design[4]-design[3],design[6]-design[5],design[7]-design[5],design[6]-design[7]))
  }
  
  
  lb <- c(0,0,0,0,0,0,0)
  ub <- c(0.25,0.25,0.25,0.25,0.25,0.25,0.25)
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
  
  this_psi<- rep(0, 5)
  for (i in 1:5) {
    this_psi[i] <-  -psi(design_l,theta_vals_l[i,],del,n,sigma,q)
  }
  
  this_psi<- as.numeric(round(this_psi, 2))
  this_psi_min<- min(this_psi)
  
  print("this_psi")
  print(this_psi)
  index <- which(this_psi == this_psi_min)    # Find LFC with minimum power
  
  
  
  
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
  
  
  prob_vals_t_l1<- matrix(0, length(H), 5)
  delta_l <- rep(0, 5)
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
    
    theta_vals_l <- theta_vals_l1
    
    print(theta_vals_l)  
    
    
    this_psi<- rep(0, 5)
    for (i in 1:5) {
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

N <- 105                    # Total No. of subjects
n <- sqrt(N)
sigma <-5.35                # Standard deviation
al <- 0.05/22               # significance level alpha
q <- qnorm(1-al)           # critical value
del<-5.1                   # delta
iteration_num<-0
B_pi_l <- 100             # Initialization of B_pi
h_grid_space <- 0.25       # grid length
H<- seq(0, 1, by=h_grid_space)      # grid vector


print( "Initial Prior on Theta : ")

# Theta values (LFC's) are their corresponding prior probabilities

the_1<- c(0,0,0,0,0,0,0,0,del,0,0)
the_2<- c(4.54,4.60,0,4.40,2.66,2.66,5.1,4.52,4.54,2.76,2.76)
the_3<- c(0.08,0.04,0,0.15,2.68,2.68,0.06,0.06,0.08,5.1,1.24)
the_4<- c(4.34,4.34,4.43,4.22,0,4.35,4.35,4.35,4.34,5.1,4.28)
the_5<- c(4.36,4.39,4.43,4.41,4.38,0,4.35,4.35,4.36,4.54,5.1)

theta_vals_l<- rbind(the_1,the_2,the_3,the_4,the_5)



prob_vals_l <- c(0.1,0.1,0.2,0.3,0.3)
# Starting optimum on the average design and corresponding power

print( "Starting Design")
design_l <- get_optimal_on_the_average_design(theta_vals_l,prob_vals_l,del,n,sigma,q)
B_pi_l <- B(design_l,theta_vals_l, prob_vals_l,del,n,sigma,q)

theta_vals_l[1,]<- c(0,0,0,0,0,0,0,0,del,0,0)

# Computing the LFC values

# Lower and upper bounds
theta_0 <-c(0,0,0,0,0,0)

lb <- c(0,0,0,0,0,0)
ub <- c(del,del,del,del,del,del)

# Set optimization options.
local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-7 )
opts <- list( "algorithm"= "NLOPT_GN_ISRES",
              "xtol_rel"= 1.0e-7,
              "maxeval"= 20000,
              "local_opts" = local_opts,
              "print_level" = 0 )


res5 <- nloptr (x0 = theta_0,
                design = design_l,
                del=del,
                n=n,
                sigma=sigma,
                q=q,
                eval_f = psi5,
                lb = lb,
                ub = ub,
                opts = opts
)


th5 <- res5$solution
print('th5')
print(th5)

theta_vals_l[2,]<- c(th5[1],th5[2],0,th5[3],th5[4],th5[4],del,th5[5],th5[1],th5[6],th5[6])

# Lower and upper bounds
theta_0 <-c(0,0,0,0,0,0)
lb <- c(0,0,0,0,0,0)
ub <- c(del,del,del,del,del,del)

# Set optimization options.
local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-7 )
opts <- list( "algorithm"= "NLOPT_GN_ISRES",
              "xtol_rel"= 1.0e-7,
              "maxeval"= 20000,
              "local_opts" = local_opts,
              "print_level" = 0 )


res6 <- nloptr (x0 = theta_0,
                design = design_l,
                del=del,
                n=n,
                sigma=sigma,
                q=q,
                eval_f = psi6,
                lb = lb,
                ub = ub,
                opts = opts
)


th6 <- res6$solution
print('th6')
print(th6)

theta_vals_l[3,]<- c(th6[1],th6[2],0,th6[3],th6[4],th6[4],th6[5],th6[5],th6[1],del,th6[6])

# Lower and upper bounds
theta_0 <-c(0,0,0,0,0,0,0)

lb <- c(0,0,0,0,0,0,0)
ub <- c(del,del,del,del,del,del,del)

# Set optimization options.
local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-7 )
opts <- list( "algorithm"= "NLOPT_GN_ISRES",
              "xtol_rel"= 1.0e-7,
              "maxeval"= 20000,
              "local_opts" = local_opts,
              "print_level" = 0 )


res9 <- nloptr (x0 = theta_0,
                design = design_l,
                del=del,
                n=n,
                sigma=sigma,
                q=q,
                eval_f = psi9,
                lb = lb,
                ub = ub,
                opts = opts
)


th9 <- res9$solution
print('th9')
print(th9)

theta_vals_l[4,]<- c(th9[1],th9[2],th9[3],th9[4],0,th9[5],th9[6],th9[6],th9[1],del,th9[7])

# Lower and upper bounds
theta_0 <-c(0,0,0,0,0,0,0)

lb <- c(0,0,0,0,0,0,0)
ub <- c(del,del,del,del,del,del,del)

# Set optimization options.
local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-7 )
opts <- list( "algorithm"= "NLOPT_GN_ISRES",
              "xtol_rel"= 1.0e-7,
              "maxeval"= 20000,
              "local_opts" = local_opts,
              "print_level" = 0 )


res10<- nloptr (x0 = theta_0,
                design = design_l,
                del=del,
                n=n,
                sigma=sigma,
                q=q,
                eval_f = psi10,
                lb = lb,
                ub = ub,
                opts = opts
)


th10 <- res10$solution
print('th10')
print(th10)

theta_vals_l[5,]<- c(th10[1],th10[2],th10[3],th10[4],th10[5],0,th10[6],th10[6],th10[1],th10[7],del)

print("theta_vals_l")
print(theta_vals_l)

step_1(H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,del,n,sigma,q)


toc()
