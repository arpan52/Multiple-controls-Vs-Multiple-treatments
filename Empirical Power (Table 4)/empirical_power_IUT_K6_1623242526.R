rm(list = ls())

true.mean <- c(84.2,94.4,92.5,90.3,88.3,90.1) #c(97.8,96.3,94,94.2,93.7,91.9) #c(2.91, 4.12, 6.26, 2.61, 3.26, 5.82)   # Population means
K <- length(true.mean)                                  # No. of comparison groups
psd <-  9.19  #7.34                                               # Pooled standard deviation

n_m <- c(63,135,73,73,73,75)                              # Max-min design
#n_b <- c(67,132,66,67,67,93)                              # C-optimal design
n_b <- c(83,78,80,82,79,90)                    # Balanced design
#n_b <- c(75,121,68,68,68,92)
m <- 100000                                               # No. of iterations


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
  # Comparisons: (1,6), (2,3), (2,4), (2,5), (2,6)
  # -------------------------------------------------------
  
  Z_16_m <- (Y.m_mean[1] - Y.m_mean[6]) /
    (psd * sqrt((1/n_m[1]) + (1/n_m[6])))
  
  Z_23_m <- (Y.m_mean[2] - Y.m_mean[3]) /
    (psd * sqrt((1/n_m[2]) + (1/n_m[3])))
  
  Z_24_m <- (Y.m_mean[2] - Y.m_mean[4]) /
    (psd * sqrt((1/n_m[2]) + (1/n_m[4])))
  
  Z_25_m <- (Y.m_mean[2] - Y.m_mean[5]) /
    (psd * sqrt((1/n_m[2]) + (1/n_m[5])))
  
  Z_26_m <- (Y.m_mean[2] - Y.m_mean[6]) /
    (psd * sqrt((1/n_m[2]) + (1/n_m[6])))
  
  
  # Minimum absolute test statistic
  t_stat_m <- min(c(
    abs(Z_16_m),
    abs(Z_23_m),
    abs(Z_24_m),
    abs(Z_25_m),
    abs(Z_26_m)
  ))
  
  
  # Reject H0 if minimum absolute statistic exceeds q
  if (t_stat_m > q) {
    ctr_m <- ctr_m + 1
  }
  
  
  # -------------------------------------------------------
  # Alternative data: Balanced design
  # -------------------------------------------------------
  
  for (j in 1:K) {
    Y.b <- rnorm(n_b[j], true.mean[j], psd)
    Y.b_mean[j] <- mean(Y.b)
  }
  
  
  # -------------------------------------------------------
  # Test statistics for balanced design
  # Comparisons: (1,6), (2,3), (2,4), (2,5), (2,6)
  # -------------------------------------------------------
  
  Z_16_b <- (Y.b_mean[1] - Y.b_mean[6]) /
    (psd * sqrt((1/n_b[1]) + (1/n_b[6])))
  
  Z_23_b <- (Y.b_mean[2] - Y.b_mean[3]) /
    (psd * sqrt((1/n_b[2]) + (1/n_b[3])))
  
  Z_24_b <- (Y.b_mean[2] - Y.b_mean[4]) /
    (psd * sqrt((1/n_b[2]) + (1/n_b[4])))
  
  Z_25_b <- (Y.b_mean[2] - Y.b_mean[5]) /
    (psd * sqrt((1/n_b[2]) + (1/n_b[5])))
  
  Z_26_b <- (Y.b_mean[2] - Y.b_mean[6]) /
    (psd * sqrt((1/n_b[2]) + (1/n_b[6])))
  
  
  # Minimum absolute test statistic
  t_stat_b <- min(c(
    abs(Z_16_b),
    abs(Z_23_b),
    abs(Z_24_b),
    abs(Z_25_b),
    abs(Z_26_b)
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

print("Power: max-min design, balanced design")
print(c(
  ctr_m / m,
  ctr_b / m
))
