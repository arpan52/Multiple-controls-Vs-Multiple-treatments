################################################################################################
###### H-algorithm for min-max design for a general bipartite graph for K=10  for UIT         ##
###### Pairs = (1,3),(1,4),(1,5),(1,6),(1,7),(1,8),(1,9),(1,10),(2,6),(2,7),(2,8),(2,9),(2,10)##
################################################################################################




#install.packages("nloptr")
#install.packages("mvtnorm")


library('nloptr')
library(mvtnorm)
library(tictoc)
tic()


rm(list = ls())
# iniatializations

mu<-rep(0,10)
x<-rep(0,10)

# Function for calculating the mean and covariance matrix of Z-ij's and criterion 'Psi'.

psi<- function(design,mu,del,n,sigma,q)
{
  
  # Calculating mean and covariance matrix of Z_ij's
  
  x[1]<- design[1]
  x[2]<- design[2]
  x[3]<- design[3]
  x[4]<- design[3]
  x[5]<- design[3]
  x[6]<- design[4]
  x[7]<- design[4]
  x[8]<- design[4]
  x[9]<- design[4]
  x[10]<- design[4]
  
  de13 <- ((mu[1]-mu[3])*n)/sigma
  de14 <- ((mu[1]-mu[4])*n)/sigma
  de15 <- ((mu[1]-mu[5])*n)/sigma
  de16 <- ((mu[1]-mu[6])*n)/sigma
  de17 <- ((mu[1]-mu[7])*n)/sigma
  de18 <- ((mu[1]-mu[8])*n)/sigma
  de19 <- ((mu[1]-mu[9])*n)/sigma
  de110 <- ((mu[1]-mu[10])*n)/sigma
  de26 <- ((mu[2]-mu[6])*n)/sigma
  de27 <- ((mu[2]-mu[7])*n)/sigma
  de28 <- ((mu[2]-mu[8])*n)/sigma
  de29 <- ((mu[2]-mu[9])*n)/sigma
  de210 <- ((mu[2]-mu[10])*n)/sigma
  
  
  
  
  a1 <- x[3]*x[4]
  b1 <- (x[1]+x[3])*(x[1]+x[4])
  rho1 <- sqrt(a1/b1)
  
  a2 <- x[3]*x[5]
  b2 <- (x[1]+x[3])*(x[1]+x[5])
  rho2 <- sqrt(a2/b2)
  
  a3 <- x[3]*x[6]
  b3 <- (x[1]+x[3])*(x[1]+x[6])
  rho3 <- sqrt(a3/b3)
  
  a4 <- x[3]*x[7]
  b4 <- (x[1]+x[3])*(x[1]+x[7])
  rho4 <- sqrt(a4/b4)
  
  a5 <- x[3]*x[8]
  b5 <- (x[1]+x[3])*(x[1]+x[8])
  rho5 <- sqrt(a5/b5)
  
  a6 <- x[3]*x[9]
  b6 <- (x[1]+x[3])*(x[1]+x[9])
  rho6 <- sqrt(a6/b6)
  
  a7 <- x[3]*x[10]
  b7 <- (x[1]+x[3])*(x[1]+x[10])
  rho7 <- sqrt(a7/b7)
  
  rho8<-0
  rho9<-0
  rho10<-0
  rho11<-0
  rho12<-0
  
  
  a13 <- x[4]*x[5]
  b13 <- (x[1]+x[4])*(x[1]+x[5])
  rho13 <- sqrt(a13/b13)
  
  a14 <- x[4]*x[6]
  b14 <- (x[1]+x[4])*(x[1]+x[6])
  rho14 <- sqrt(a14/b14)
  
  a15 <- x[4]*x[7]
  b15 <- (x[1]+x[4])*(x[1]+x[7])
  rho15 <- sqrt(a15/b15)
  
  a16 <- x[4]*x[8]
  b16 <- (x[1]+x[4])*(x[1]+x[8])
  rho16 <- sqrt(a16/b16)
  
  a17 <- x[4]*x[9]
  b17 <- (x[1]+x[4])*(x[1]+x[9])
  rho17 <- sqrt(a17/b17)
  
  a18 <- x[4]*x[10]
  b18 <- (x[1]+x[4])*(x[1]+x[10])
  rho18 <- sqrt(a18/b18)
  
  rho19<-0
  rho20<-0
  rho21<-0
  rho22<-0
  rho23<-0
  
  a24 <- x[5]*x[6]
  b24 <- (x[1]+x[5])*(x[1]+x[6])
  rho24 <- sqrt(a24/b24)
  
  a25 <- x[5]*x[7]
  b25 <- (x[1]+x[5])*(x[1]+x[7])
  rho25 <- sqrt(a25/b25)
  
  a26 <- x[5]*x[8]
  b26 <- (x[1]+x[5])*(x[1]+x[8])
  rho26 <- sqrt(a26/b26)
  
  a27 <- x[5]*x[9]
  b27 <- (x[1]+x[5])*(x[1]+x[9])
  rho27 <- sqrt(a27/b27)
  
  a28 <- x[5]*x[10]
  b28 <- (x[1]+x[5])*(x[1]+x[10])
  rho28 <- sqrt(a28/b28)
  
  rho29<-0
  rho30<-0
  rho31<-0
  rho32<-0
  rho33<-0
  
  a34 <- x[6]*x[7]
  b34 <- (x[1]+x[6])*(x[1]+x[7])
  rho34 <- sqrt(a34/b34)
  
  a35 <- x[6]*x[8]
  b35 <- (x[1]+x[6])*(x[1]+x[8])
  rho35 <- sqrt(a35/b35)
  
  a36 <- x[6]*x[9]
  b36 <- (x[1]+x[6])*(x[1]+x[9])
  rho36 <- sqrt(a36/b36)
  
  a37 <- x[6]*x[10]
  b37 <- (x[1]+x[6])*(x[1]+x[10])
  rho37 <- sqrt(a37/b37)
  
  a38 <- x[1]*x[2]
  b38 <- (x[1]+x[6])*(x[2]+x[6])
  rho38 <- sqrt(a38/b38)
  
  rho39<-0
  rho40<-0
  rho41<-0
  rho42<-0
  
  a43 <- x[7]*x[8]
  b43 <- (x[1]+x[7])*(x[1]+x[8])
  rho43 <- sqrt(a43/b43)
  
  a44 <- x[7]*x[9]
  b44 <- (x[1]+x[7])*(x[1]+x[9])
  rho44 <- sqrt(a44/b44)
  
  a45 <- x[7]*x[10]
  b45 <- (x[1]+x[7])*(x[1]+x[10])
  rho45 <- sqrt(a45/b45)
  
  rho46<-0
  
  a47 <- x[1]*x[2]
  b47 <- (x[1]+x[7])*(x[2]+x[7])
  rho47 <- sqrt(a47/b47)
  
  
  
  rho48<-0
  rho49<-0
  rho50<-0
  
  a51 <- x[8]*x[9]
  b51 <- (x[1]+x[8])*(x[1]+x[9])
  rho51 <- sqrt(a51/b51)
  
  a52 <- x[8]*x[10]
  b52 <- (x[1]+x[8])*(x[1]+x[10])
  rho52 <- sqrt(a52/b52)
  
  
  rho53<-0
  rho54<-0
  
  a55 <- x[1]*x[2]
  b55 <- (x[1]+x[8])*(x[2]+x[8])
  rho55 <- sqrt(a55/b55)
  
  rho56<-0
  rho57<-0
  
  
  a58 <- x[9]*x[10]
  b58 <- (x[1]+x[9])*(x[1]+x[10])
  rho58 <- sqrt(a58/b58)
  
  
  rho59<-0
  rho60<-0
  rho61<-0
  
  a62 <- x[1]*x[2]
  b62 <- (x[1]+x[9])*(x[2]+x[9])
  rho62 <- sqrt(a62/b62)
  
  rho63<-0
  
  rho64<-0
  rho65<-0
  rho66<-0
  rho67<-0
  
  a68 <- x[1]*x[2]
  b68 <- (x[1]+x[10])*(x[2]+x[10])
  rho68 <- sqrt(a68/b68)
  
  
  
  a69 <- x[6]*x[7]
  b69 <- (x[2]+x[6])*(x[2]+x[7])
  rho69 <- sqrt(a69/b69)
  
  a70 <- x[6]*x[8]
  b70 <- (x[2]+x[6])*(x[2]+x[8])
  rho70 <- sqrt(a70/b70)
  
  a71 <- x[6]*x[9]
  b71 <- (x[2]+x[6])*(x[2]+x[9])
  rho71 <- sqrt(a71/b71)
  
  a72 <- x[6]*x[10]
  b72 <- (x[2]+x[6])*(x[2]+x[10])
  rho72 <- sqrt(a72/b72)
  
  a73 <- x[7]*x[8]
  b73 <- (x[2]+x[7])*(x[2]+x[8])
  rho73 <- sqrt(a73/b73)
  
  a74 <- x[7]*x[9]
  b74 <- (x[2]+x[7])*(x[2]+x[9])
  rho74 <- sqrt(a74/b74)
  
  a75 <- x[7]*x[10]
  b75 <- (x[2]+x[7])*(x[2]+x[10])
  rho75 <- sqrt(a75/b75)
  
  a76 <- x[8]*x[9]
  b76 <- (x[2]+x[8])*(x[2]+x[9])
  rho76 <- sqrt(a76/b76)
  
  a77 <- x[8]*x[10]
  b77 <- (x[2]+x[8])*(x[2]+x[10])
  rho77 <- sqrt(a77/b77)
  
  a78 <- x[9]*x[10]
  b78 <- (x[2]+x[9])*(x[2]+x[10])
  rho78 <- sqrt(a78/b78)
  
  
  mu1 <- (sqrt((x[1]*x[3])/(x[1] + x[3])))*de13
  mu2 <- (sqrt((x[1]*x[4])/(x[1] + x[4])))*de14
  mu3 <- (sqrt((x[1]*x[5])/(x[1] + x[5])))*de15
  mu4 <- (sqrt((x[1]*x[6])/(x[1] + x[6])))*de16
  mu5 <- (sqrt((x[1]*x[7])/(x[1] + x[7])))*de17
  mu6 <- (sqrt((x[1]*x[8])/(x[1] + x[8])))*de18
  mu7 <- (sqrt((x[1]*x[9])/(x[1] + x[9])))*de19
  mu8 <- (sqrt((x[1]*x[10])/(x[1] + x[10])))*de110
  mu9 <- (sqrt((x[2]*x[6])/(x[2] + x[6])))*de26
  mu10 <- (sqrt((x[2]*x[7])/(x[2] + x[7])))*de27
  mu11 <- (sqrt((x[2]*x[8])/(x[2] + x[8])))*de28
  mu12 <- (sqrt((x[2]*x[9])/(x[2] + x[9])))*de29
  mu13 <- (sqrt((x[2]*x[10])/(x[2] + x[10])))*de210
  
  
  m <- c(mu1,mu2,mu3,mu4,mu5,mu6,mu7,mu8,mu9,mu10,mu11,mu12,mu13)
  
  S <-  rbind(c(1,rho1,rho2,rho3,rho4,rho5,rho6,rho7,rho8,rho9,rho10,rho11,rho12),c(rho1,1,rho13,rho14,rho15,rho16,rho17,rho18,rho19,rho20,rho21,rho22,rho23),c(rho2,rho13,1,rho24,rho25,rho26,rho27,rho28,rho29,rho30,rho31,rho32,rho33),c(rho3,rho14,rho24,1,rho34,rho35,rho36,rho37,rho38,rho39,rho40,rho41,rho42),c(rho4,rho15,rho25,rho34,1,rho43,rho44,rho45,rho46,rho47,rho48,rho49,rho50),c(rho5,rho16,rho26,rho35,rho43,1,rho51,rho52,rho53,rho54,rho55,rho56,rho57),c(rho6,rho17,rho27,rho36,rho44,rho51,1,rho58,rho59,rho60,rho61,rho62,rho63),c(rho7,rho18,rho28,rho37,rho45,rho52,rho58,1,rho64,rho65,rho66,rho67,rho68),c(rho8,rho19,rho29,rho38,rho46,rho53,rho59,rho64,1,rho69,rho70,rho71,rho72),c(rho9,rho20,rho30,rho39,rho47,rho54,rho60,rho65,rho69,1,rho73,rho74,rho75),c(rho10,rho21,rho31,rho40,rho48,rho55,rho61,rho66,rho70,rho73,1,rho76,rho77),c(rho11,rho22,rho32,rho41,rho49,rho56,rho62,rho67,rho71,rho74,rho76,1,rho78),c(rho12,rho23,rho33,rho42,rho50,rho57,rho63,rho68,rho72,rho75,rho77,rho78,1))
  
  
  lw1 <- c(-q,-q,-q,-q,-q,-q,-q,-q,-q,-q,-q,-q,-q)
  up1 <- c(q,q,q,q,q,q,q,q,q,q,q,q,q)
  
  
  return(-(1-pmvnorm(mean=m, sigma=S, lower=lw1, upper=up1)))   #Power function
}  

psi1<- function(theta,design,del,n,sigma,q){   # define the power based on the form of the LFCs
  

  mu[1]<- 0
  mu[2]<- theta[1]
  mu[3]<- theta[2]
  mu[4]<- theta[2]
  mu[5]<- theta[2]
  mu[6]<- del
  mu[7]<- theta[3]
  mu[8]<- theta[3]
  mu[9]<- theta[3]
  mu[10]<- theta[3]
  
  return(-psi(design,mu,del,n,sigma,q))
  
}
# Function for calculating the average power B

B<- function(design,theta_vals,prob_vals,del,n,sigma,q){
  
  s <- 0
  for (i in 1:length(prob_vals)){
    s <- s + prob_vals[i]* psi(design,theta_vals[i,],del,n,sigma,q)      # Average power associated to LFCs and probabilities
  }
  ret_val<-s
  return(ret_val[1])
}

# Function for compute optimum on the avarage design

get_optimal_on_the_average_design<- function(theta_vals,prob_vals,del,n,sigma,q)
{
  des_0 <-c(0.1,0.15,0.1,0.07)
  #options = optimoptions('fmincon','Display','none')
  const1 <- function(design,theta_vals,prob_vals,del,n,sigma,q) {
    return( (design[1]+design[2]+3*design[3]+5*design[4]-1))
  }
  
  const2 <- function(design,theta_vals,prob_vals,del,n,sigma,q) {
    return( c(design[2]-design[1],design[3]-design[2],design[4]-design[2],design[3]-design[4]))
  }
  
  
  lb <- c(0,0,0,0)
  ub <- c(0.25,0.25,0.15,0.15)
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

  # Lower and upper bounds
  theta_0 <-c(0,0,0,0)

  lb <- c(0,0,0,0)
  ub <- c(del,del,del,del)
  
  # Set optimization options.
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-7 )
  opts <- list( "algorithm"= "NLOPT_GN_ISRES",
                "xtol_rel"= 1.0e-7,
                "maxeval"= 20000,
                "local_opts" = local_opts,
                "print_level" = 0 )
  # Computing the LFC
  
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
  
  theta_vals_l<-matrix(0,5,10)
  
  th <- res1$solution
  theta_vals_l[4,] <- c(0,th[1],th[2],th[2],th[2],del,th[3],th[3],th[4],th[4])
  theta_vals_l[1,] <- c(0,0,del,0,0,0,0,0,0,0)
  theta_vals_l[2,] <- c(0,0,0,del,0,0,0,0,0,0)
  theta_vals_l[3,] <- c(0,0,0,0,del,0,0,0,0,0)
  

  theta_vals_l[5,]<-  theta_vals_l[4,]
  theta_vals_l[5,1]<-  theta_vals_l[4,2]
  theta_vals_l[5,2]<-  theta_vals_l[4,1]
  theta_vals_l[5,3]<-  theta_vals_l[4,10]
  theta_vals_l[5,4]<-  theta_vals_l[4,10]
  theta_vals_l[5,5]<-  theta_vals_l[4,10]
  theta_vals_l[5,9]<-  theta_vals_l[4,3]
  theta_vals_l[5,10]<-  theta_vals_l[4,3]
  
  print("theta_vals_l")
  print(theta_vals_l)    # Set of computed LFCs
  
  
  this_psi<- rep(0, 5)
  for (i in 1:5) {
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
    
    # Lower and upper bounds
    theta_0 <-c(0,0,0,0)
    lb <- c(0,0,0,0)
    ub <- c(del,del,del,del)
    
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
    
    theta_vals_l<-matrix(0,5,10)
    
    th <- res1$solution
    theta_vals_l[4,]<- c(0,th[1],th[2],th[2],th[2],del,th[3],th[3],th[4],th[4])
    theta_vals_l[1,] <- c(0,0,del,0,0,0,0,0,0,0)
    theta_vals_l[2,] <- c(0,0,0,del,0,0,0,0,0,0)
    theta_vals_l[3,] <- c(0,0,0,0,del,0,0,0,0,0)
    
    
    theta_vals_l[5,]<-  theta_vals_l[4,]
    theta_vals_l[5,1]<-  theta_vals_l[4,2]
    theta_vals_l[5,2]<-  theta_vals_l[4,1]
    theta_vals_l[5,3]<-  theta_vals_l[4,10]
    theta_vals_l[5,4]<-  theta_vals_l[4,10]
    theta_vals_l[5,5]<-  theta_vals_l[4,10]
    theta_vals_l[5,9]<-  theta_vals_l[4,3]
    theta_vals_l[5,10]<-  theta_vals_l[4,3]
    
    
    print("theta_vals_l")
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

N <- 95                    # Total No. of subjects
n <- sqrt(N)
sigma <-5.32               # Standard deviation
al <- 0.05/26               # significance level alpha
q <- qnorm(1-al)           # critical value
del<-4.3                     # delta
iteration_num<-0
B_pi_l <- 100             # Initialization of B_pi
h_grid_space <- 0.5       # grid length
H<- seq(0, 1, by=h_grid_space)      # grid vector


print( "Initial Prior on Theta : ")

# Theta values (LFC's) are their corresponding prior probabilities


the_1<- c(0,0,del,0,0,0,0,0,0,0)
the_2<- c(0,0,0,del,0,0,0,0,0,0)
the_3<- c(0,0,0,0,del,0,0,0,0,0)
the_4<- c(0,0,0,0,0,del,0,0,0,0)
the_5<- c(0,0,0,0,0,0,del,0,0,0)

theta_vals_l<- rbind(the_1,the_2,the_3,the_4,the_5)


p_v<-c(0.3/3,0.35/5,0.35/5)
prob_vals_l <- c(0.1,0.1,0.1,0.3,0.4)

# Starting optimum on the average design and corresponding power

print( "Starting Design")
design_l <- get_optimal_on_the_average_design(theta_vals_l,prob_vals_l,del,n,sigma,q)
B_pi_l <- B(design_l,theta_vals_l, prob_vals_l,del,n,sigma,q)

step_1(H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,del,n,sigma,q)


toc()