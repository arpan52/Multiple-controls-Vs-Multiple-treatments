rm(list = ls())

true.mean <- c(84.2,94.4,88.3,90.3,90.1,92.5)#c(84.2,94.4,92.5,90.3,90.1,88.3) #c(97.8, 96.3, 94, 94.2, 93.7, 91.9)
K <- length(true.mean)                                  # No. of comparison groups
psd <- 9.19                                              # Pooled standard deviation

n_m <- c(112,112,84,84,72,72) #c(103,103,77,77,66,66)                       # Max-min design
n_b <-  c(102,102,84,84,59,59)                       # C-optimal design
#n_b <-  c(83,78,79,82,90,80)                     # Balanced design
#n_b<- c(100,100,79,79,67,67)
m <- 100000                                      # No. of iterations


# ---------------------------------------------------------
# Initializations
# ---------------------------------------------------------

Y.m_mean <- rep(0, K)
Y.b_mean <- rep(0, K)

ctr_m <- 0
ctr_b <- 0


# ---------------------------------------------------------
# Critical value
# ---------------------------------------------------------
# One-sided test at alpha = 0.05
# No null simulation and no multiplicity correction

alpha <- 0.05
q <- qnorm(1 - alpha)

print("Critical value")
print(q)


# ---------------------------------------------------------
# Simulation under the alternative
# ---------------------------------------------------------

for (i in 1:m) {
  
  # -------------------------------------------------------
  # Alternative data: Max-min design
  # -------------------------------------------------------
  
  for (j in 1:K) {
    Y.m <- rnorm(n_m[j], true.mean[j], psd)
    Y.m_mean[j] <- mean(Y.m)
  }
  
  
  # -------------------------------------------------------
  # Test statistics for Max-min design
  # Comparisons: (1,3), (1,4), (1,6),
  #              (2,3), (2,4), (2,5)
  # -------------------------------------------------------
  
  Z_13_m <- (Y.m_mean[1] - Y.m_mean[3]) /
    (psd * sqrt((1/n_m[1]) + (1/n_m[3])))
  
  Z_14_m <- (Y.m_mean[1] - Y.m_mean[4]) /
    (psd * sqrt((1/n_m[1]) + (1/n_m[4])))
  
  Z_16_m <- (Y.m_mean[1] - Y.m_mean[6]) /
    (psd * sqrt((1/n_m[1]) + (1/n_m[6])))
  
  Z_23_m <- (Y.m_mean[2] - Y.m_mean[3]) /
    (psd * sqrt((1/n_m[2]) + (1/n_m[3])))
  
  Z_24_m <- (Y.m_mean[2] - Y.m_mean[4]) /
    (psd * sqrt((1/n_m[2]) + (1/n_m[4])))
  
  Z_25_m <- (Y.m_mean[2] - Y.m_mean[5]) /
    (psd * sqrt((1/n_m[2]) + (1/n_m[5])))
  
  
  # Minimum absolute test statistic
  t_stat_m <- min(c(
    abs(Z_13_m),
    abs(Z_14_m),
    abs(Z_16_m),
    abs(Z_23_m),
    abs(Z_24_m),
    abs(Z_25_m)
  ))
  
  
  # Reject H0 if minimum absolute statistic exceeds q
  if (t_stat_m > q) {
    ctr_m <- ctr_m + 1
  }
  
  
  # -------------------------------------------------------
  # Alternative data: C-optimal design
  # -------------------------------------------------------
  
  for (j in 1:K) {
    Y.b <- rnorm(n_b[j], true.mean[j], psd)
    Y.b_mean[j] <- mean(Y.b)
  }
  
  
  # -------------------------------------------------------
  # Test statistics for C-optimal design
  # Comparisons: (1,3), (1,4), (1,6),
  #              (2,3), (2,4), (2,5)
  # -------------------------------------------------------
  
  Z_13_b <- (Y.b_mean[1] - Y.b_mean[3]) /
    (psd * sqrt((1/n_b[1]) + (1/n_b[3])))
  
  Z_14_b <- (Y.b_mean[1] - Y.b_mean[4]) /
    (psd * sqrt((1/n_b[1]) + (1/n_b[4])))
  
  Z_16_b <- (Y.b_mean[1] - Y.b_mean[6]) /
    (psd * sqrt((1/n_b[1]) + (1/n_b[6])))
  
  Z_23_b <- (Y.b_mean[2] - Y.b_mean[3]) /
    (psd * sqrt((1/n_b[2]) + (1/n_b[3])))
  
  Z_24_b <- (Y.b_mean[2] - Y.b_mean[4]) /
    (psd * sqrt((1/n_b[2]) + (1/n_b[4])))
  
  Z_25_b <- (Y.b_mean[2] - Y.b_mean[5]) /
    (psd * sqrt((1/n_b[2]) + (1/n_b[5])))
  
  
  # Minimum absolute test statistic
  t_stat_b <- min(c(
    abs(Z_13_b),
    abs(Z_14_b),
    abs(Z_16_b),
    abs(Z_23_b),
    abs(Z_24_b),
    abs(Z_25_b)
  ))
  
  
  # Reject H0 if minimum absolute statistic exceeds q
  if (t_stat_b > q) {
    ctr_b <- ctr_b + 1
  }
  
}


# ---------------------------------------------------------
# Results
# ---------------------------------------------------------

print("Critical value")
print(q)

print("Power: max-min design, C-optimal design")
print(c(
  ctr_m / m,
  ctr_b / m
))
