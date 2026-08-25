#####################################################################################
###### H-algorithm for min-max design for general Bipartite graph for K=6 for IUT  ##
###### Pairs = (1,3),(1,4), (1,6), (2,3),(2,4), (2,5)                              ##
#####################################################################################



#install.packages("nloptr")
#install.packages("mvtnorm")


library('nloptr')
library(mvtnorm)
rm(list = ls())
x<-rep(0,6)


# Function for calculating the mean and covariance matrix of Z-ij's and criterion 'Psi'.
psi <- function(design,theta,n,sigma,q){
  
  # ------------------------------------------------------------
  # Calculating sample size proportions
  # ------------------------------------------------------------
  
  x[1] <- design[1]
  x[2] <- design[2]
  x[3] <- design[3]
  x[4] <- design[3]
  x[5] <- design[3]
  x[6] <- design[4]
  
  
  # ------------------------------------------------------------
  # Non-centrality parameters for the five comparisons
  # ------------------------------------------------------------
  
  mu16 <- (n*(theta[1]-theta[6]))/
    (sigma*sqrt((1/x[1])+(1/x[6])))
  
  mu23 <- (n*(theta[2]-theta[3]))/
    (sigma*sqrt((1/x[2])+(1/x[3])))
  
  mu24 <- (n*(theta[2]-theta[4]))/
    (sigma*sqrt((1/x[2])+(1/x[4])))
  
  mu25 <- (n*(theta[2]-theta[5]))/
    (sigma*sqrt((1/x[2])+(1/x[5])))
  
  mu26 <- (n*(theta[2]-theta[6]))/
    (sigma*sqrt((1/x[2])+(1/x[6])))
  
  
  # ------------------------------------------------------------
  # Marginal probabilities
  # ------------------------------------------------------------
  
  p16 <- pnorm(q,mu16)-pnorm(-q,mu16)
  p23 <- pnorm(q,mu23)-pnorm(-q,mu23)
  p24 <- pnorm(q,mu24)-pnorm(-q,mu24)
  p25 <- pnorm(q,mu25)-pnorm(-q,mu25)
  p26 <- pnorm(q,mu26)-pnorm(-q,mu26)
  
  
  # ------------------------------------------------------------
  # Correlations between test statistics
  # ------------------------------------------------------------
  
  # (1,6) with (2,3)
  rho1623 <- 0
  
  # (1,6) with (2,4)
  rho1624 <- 0
  
  # (1,6) with (2,5)
  rho1625 <- 0
  
  # (1,6) with (2,6)
  rho1626 <- sqrt(
    (x[1]*x[2])/
      ((x[1]+x[6])*(x[2]+x[6]))
  )
  
  
  # Comparisons sharing treatment 2
  rho2324 <- sqrt(
    (x[3]*x[4])/
      ((x[2]+x[3])*(x[2]+x[4]))
  )
  
  rho2325 <- sqrt(
    (x[3]*x[5])/
      ((x[2]+x[3])*(x[2]+x[5]))
  )
  
  rho2326 <- sqrt(
    (x[3]*x[6])/
      ((x[2]+x[3])*(x[2]+x[6]))
  )
  
  rho2425 <- sqrt(
    (x[4]*x[5])/
      ((x[2]+x[4])*(x[2]+x[5]))
  )
  
  rho2426 <- sqrt(
    (x[4]*x[6])/
      ((x[2]+x[4])*(x[2]+x[6]))
  )
  
  rho2526 <- sqrt(
    (x[5]*x[6])/
      ((x[2]+x[5])*(x[2]+x[6]))
  )
  
  
  # ------------------------------------------------------------
  # Bivariate probabilities
  # ------------------------------------------------------------
  
  lw2 <- c(-q,-q)
  up2 <- c(q,q)
  
  
  s1623 <- rbind(
    c(1,rho1623),
    c(rho1623,1)
  )
  mu1623 <- c(mu16,mu23)
  p1623 <- pmvnorm(
    mean=mu1623,
    sigma=s1623,
    lower=lw2,
    upper=up2
  )
  
  
  s1624 <- rbind(
    c(1,rho1624),
    c(rho1624,1)
  )
  mu1624 <- c(mu16,mu24)
  p1624 <- pmvnorm(
    mean=mu1624,
    sigma=s1624,
    lower=lw2,
    upper=up2
  )
  
  
  s1625 <- rbind(
    c(1,rho1625),
    c(rho1625,1)
  )
  mu1625 <- c(mu16,mu25)
  p1625 <- pmvnorm(
    mean=mu1625,
    sigma=s1625,
    lower=lw2,
    upper=up2
  )
  
  
  s1626 <- rbind(
    c(1,rho1626),
    c(rho1626,1)
  )
  mu1626 <- c(mu16,mu26)
  p1626 <- pmvnorm(
    mean=mu1626,
    sigma=s1626,
    lower=lw2,
    upper=up2
  )
  
  
  s2324 <- rbind(
    c(1,rho2324),
    c(rho2324,1)
  )
  mu2324 <- c(mu23,mu24)
  p2324 <- pmvnorm(
    mean=mu2324,
    sigma=s2324,
    lower=lw2,
    upper=up2
  )
  
  
  s2325 <- rbind(
    c(1,rho2325),
    c(rho2325,1)
  )
  mu2325 <- c(mu23,mu25)
  p2325 <- pmvnorm(
    mean=mu2325,
    sigma=s2325,
    lower=lw2,
    upper=up2
  )
  
  
  s2326 <- rbind(
    c(1,rho2326),
    c(rho2326,1)
  )
  mu2326 <- c(mu23,mu26)
  p2326 <- pmvnorm(
    mean=mu2326,
    sigma=s2326,
    lower=lw2,
    upper=up2
  )
  
  
  s2425 <- rbind(
    c(1,rho2425),
    c(rho2425,1)
  )
  mu2425 <- c(mu24,mu25)
  p2425 <- pmvnorm(
    mean=mu2425,
    sigma=s2425,
    lower=lw2,
    upper=up2
  )
  
  
  s2426 <- rbind(
    c(1,rho2426),
    c(rho2426,1)
  )
  mu2426 <- c(mu24,mu26)
  p2426 <- pmvnorm(
    mean=mu2426,
    sigma=s2426,
    lower=lw2,
    upper=up2
  )
  
  
  s2526 <- rbind(
    c(1,rho2526),
    c(rho2526,1)
  )
  mu2526 <- c(mu25,mu26)
  p2526 <- pmvnorm(
    mean=mu2526,
    sigma=s2526,
    lower=lw2,
    upper=up2
  )
  
  
  # ------------------------------------------------------------
  # Trivariate probabilities
  # ------------------------------------------------------------
  
  lw3 <- c(-q,-q,-q)
  up3 <- c(q,q,q)
  
  
  s162324 <- rbind(
    c(1,rho1623,rho1624),
    c(rho1623,1,rho2324),
    c(rho1624,rho2324,1)
  )
  mu162324 <- c(mu16,mu23,mu24)
  p162324 <- pmvnorm(
    mean=mu162324,
    sigma=s162324,
    lower=lw3,
    upper=up3
  )
  
  
  s162325 <- rbind(
    c(1,rho1623,rho1625),
    c(rho1623,1,rho2325),
    c(rho1625,rho2325,1)
  )
  mu162325 <- c(mu16,mu23,mu25)
  p162325 <- pmvnorm(
    mean=mu162325,
    sigma=s162325,
    lower=lw3,
    upper=up3
  )
  
  
  s162326 <- rbind(
    c(1,rho1623,rho1626),
    c(rho1623,1,rho2326),
    c(rho1626,rho2326,1)
  )
  mu162326 <- c(mu16,mu23,mu26)
  p162326 <- pmvnorm(
    mean=mu162326,
    sigma=s162326,
    lower=lw3,
    upper=up3
  )
  
  
  s162425 <- rbind(
    c(1,rho1624,rho1625),
    c(rho1624,1,rho2425),
    c(rho1625,rho2425,1)
  )
  mu162425 <- c(mu16,mu24,mu25)
  p162425 <- pmvnorm(
    mean=mu162425,
    sigma=s162425,
    lower=lw3,
    upper=up3
  )
  
  
  s162426 <- rbind(
    c(1,rho1624,rho1626),
    c(rho1624,1,rho2426),
    c(rho1626,rho2426,1)
  )
  mu162426 <- c(mu16,mu24,mu26)
  p162426 <- pmvnorm(
    mean=mu162426,
    sigma=s162426,
    lower=lw3,
    upper=up3
  )
  
  
  s162526 <- rbind(
    c(1,rho1625,rho1626),
    c(rho1625,1,rho2526),
    c(rho1626,rho2526,1)
  )
  mu162526 <- c(mu16,mu25,mu26)
  p162526 <- pmvnorm(
    mean=mu162526,
    sigma=s162526,
    lower=lw3,
    upper=up3
  )
  
  
  s232425 <- rbind(
    c(1,rho2324,rho2325),
    c(rho2324,1,rho2425),
    c(rho2325,rho2425,1)
  )
  mu232425 <- c(mu23,mu24,mu25)
  p232425 <- pmvnorm(
    mean=mu232425,
    sigma=s232425,
    lower=lw3,
    upper=up3
  )
  
  
  s232426 <- rbind(
    c(1,rho2324,rho2326),
    c(rho2324,1,rho2426),
    c(rho2326,rho2426,1)
  )
  mu232426 <- c(mu23,mu24,mu26)
  p232426 <- pmvnorm(
    mean=mu232426,
    sigma=s232426,
    lower=lw3,
    upper=up3
  )
  
  
  s232526 <- rbind(
    c(1,rho2325,rho2326),
    c(rho2325,1,rho2526),
    c(rho2326,rho2526,1)
  )
  mu232526 <- c(mu23,mu25,mu26)
  p232526 <- pmvnorm(
    mean=mu232526,
    sigma=s232526,
    lower=lw3,
    upper=up3
  )
  
  
  s242526 <- rbind(
    c(1,rho2425,rho2426),
    c(rho2425,1,rho2526),
    c(rho2426,rho2526,1)
  )
  mu242526 <- c(mu24,mu25,mu26)
  p242526 <- pmvnorm(
    mean=mu242526,
    sigma=s242526,
    lower=lw3,
    upper=up3
  )
  
  
  # ------------------------------------------------------------
  # Four-variate probabilities
  # ------------------------------------------------------------
  
  lw4 <- c(-q,-q,-q,-q)
  up4 <- c(q,q,q,q)
  
  
  s16232425 <- rbind(
    c(1,rho1623,rho1624,rho1625),
    c(rho1623,1,rho2324,rho2325),
    c(rho1624,rho2324,1,rho2425),
    c(rho1625,rho2325,rho2425,1)
  )
  mu16232425 <- c(mu16,mu23,mu24,mu25)
  p16232425 <- pmvnorm(
    mean=mu16232425,
    sigma=s16232425,
    lower=lw4,
    upper=up4
  )
  
  
  s16232426 <- rbind(
    c(1,rho1623,rho1624,rho1626),
    c(rho1623,1,rho2324,rho2326),
    c(rho1624,rho2324,1,rho2426),
    c(rho1626,rho2326,rho2426,1)
  )
  mu16232426 <- c(mu16,mu23,mu24,mu26)
  p16232426 <- pmvnorm(
    mean=mu16232426,
    sigma=s16232426,
    lower=lw4,
    upper=up4
  )
  
  
  s16232526 <- rbind(
    c(1,rho1623,rho1625,rho1626),
    c(rho1623,1,rho2325,rho2326),
    c(rho1625,rho2325,1,rho2526),
    c(rho1626,rho2326,rho2526,1)
  )
  mu16232526 <- c(mu16,mu23,mu25,mu26)
  p16232526 <- pmvnorm(
    mean=mu16232526,
    sigma=s16232526,
    lower=lw4,
    upper=up4
  )
  
  
  s16242526 <- rbind(
    c(1,rho1624,rho1625,rho1626),
    c(rho1624,1,rho2425,rho2426),
    c(rho1625,rho2425,1,rho2526),
    c(rho1626,rho2426,rho2526,1)
  )
  mu16242526 <- c(mu16,mu24,mu25,mu26)
  p16242526 <- pmvnorm(
    mean=mu16242526,
    sigma=s16242526,
    lower=lw4,
    upper=up4
  )
  
  
  s23242526 <- rbind(
    c(1,rho2324,rho2325,rho2326),
    c(rho2324,1,rho2425,rho2426),
    c(rho2325,rho2425,1,rho2526),
    c(rho2326,rho2426,rho2526,1)
  )
  mu23242526 <- c(mu23,mu24,mu25,mu26)
  p23242526 <- pmvnorm(
    mean=mu23242526,
    sigma=s23242526,
    lower=lw4,
    upper=up4
  )
  
  
  # ------------------------------------------------------------
  # Five-variate probability
  # ------------------------------------------------------------
  
  lw5 <- c(-q,-q,-q,-q,-q)
  up5 <- c(q,q,q,q,q)
  
  
  s1623242526 <- rbind(
    c(1,rho1623,rho1624,rho1625,rho1626),
    c(rho1623,1,rho2324,rho2325,rho2326),
    c(rho1624,rho2324,1,rho2425,rho2426),
    c(rho1625,rho2325,rho2425,1,rho2526),
    c(rho1626,rho2326,rho2426,rho2526,1)
  )
  
  mu1623242526 <- c(
    mu16,mu23,mu24,mu25,mu26
  )
  
  p1623242526 <- pmvnorm(
    mean=mu1623242526,
    sigma=s1623242526,
    lower=lw5,
    upper=up5
  )
  
  
  # ------------------------------------------------------------
  # Power of the test
  # Inclusion-exclusion for 5 comparisons
  # ------------------------------------------------------------
  
  p1 <- p16 + p23 + p24 + p25 + p26
  
  p2 <- p1623 +
    p1624 +
    p1625 +
    p1626 +
    p2324 +
    p2325 +
    p2326 +
    p2425 +
    p2426 +
    p2526
  
  p3 <- p162324 +
    p162325 +
    p162326 +
    p162425 +
    p162426 +
    p162526 +
    p232425 +
    p232426 +
    p232526 +
    p242526
  
  p4 <- p16232425 +
    p16232426 +
    p16232526 +
    p16242526 +
    p23242526
  
  p5 <- p1623242526
  
  p <- -(1 - p1 + p2 - p3 + p4 - p5)
  
  return(p[1])
}
# Function for calculating the average power B

B<- function(design,theta_vals,prob_vals,n,sigma,q){
  
  s <- 0
  for (i in 1:length(prob_vals)){
    s <- s + prob_vals[i]* psi(design,theta_vals[i,],n,sigma,q)   # Average power associated to LFCs and probabilities
  }
  ret_val<-s
  return(ret_val[1])
}

# Function for compute optimum on the avarage design

get_optimal_on_the_average_design<- function(theta_vals,prob_vals,n,sigma,q)
{
  des_0 <-c(0.15,0.28,0.15,0.1)
  #options = optimoptions('fmincon','Display','none')
  const1 <- function(design,theta_vals,prob_vals,n,sigma,q) {
    return( design[1]+design[2]+3*design[3]+design[4]-1)
  }
  
  const2 <- function(design,theta_vals,prob_vals,del,n,sigma,q) {
    return( c(design[1]-design[2],design[3]-design[4]))
  }
  
  
  lb <- c(0,0,0,0)
  ub <- c(0.3,0.35,0.3,0.3)
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
                  eval_g_eq = const1,
                  #eval_g_ineq = eval_g0,
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
  index <- which(this_psi == this_psi_min)
  
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

N <- 109                    # Total No. of subjects
n <- sqrt(N)
sigma <- 1                  # Standard deviation
al <- 0.05               # significance level alpha
q <- qnorm(1-al)           # critical value
del<-0.6                   # delta
K<-5
iteration_num<-0
B_pi_l <- 100             # Initialization of B_pi
h_grid_space <- 0.5       # grid length
H<- seq(0, 1, by=h_grid_space)      # grid vector


print( "Initial Prior on Theta : ")

# Theta values (LFC's) are their corresponding prior probabilities



the_1<- c(-del,del,0,0,0,0)
the_2<-  c(del,-del,0,0,0,0)
the_3<-  c(-del,-del,0,0,0,0)
the_4<- c(0,0,-del,del,del,del)
the_5<- c(0,0,del,-del,del,del)
the_6<- c(0,0,del,del,-del,del)
the_7<- c(0,0,del,del,del,-del)
the_8<- c(0,0,-del,-del,del,del)
the_9<- c(0,0,-del,del,-del,del)
the_10<- c(0,0,-del,del,del,-del)
the_11<- c(0,0,-del,-del,-del,-del)



theta_vals_l<- rbind(the_1,the_2,the_3,the_4,the_5,the_6,the_7,the_8,the_9,the_10,the_11)

prob_vals_l <- rep(1/11,11)

# Starting optimum on the average design and corresponding power

print( "Starting Design")
design_l <- get_optimal_on_the_average_design(theta_vals_l,prob_vals_l,n,sigma,q)
B_pi_l <- B(design_l,theta_vals_l, prob_vals_l,n,sigma,q)

step_1(H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,n,sigma,q)
