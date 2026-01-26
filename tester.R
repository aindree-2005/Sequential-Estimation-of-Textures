
source("sequential_estimation.R")

set.seed(22)
M <- 30
N <- 30
S <- 3
K <- 3  
generate_noise <- function(M, N, S, sd = 1) {
  e <- array(rnorm(M*N*S, 0, sd), dim = c(M,N,S))
  X <- array(0, dim = c(M,N,S))

  for(m in 1:M){
    for(n in 1:N){
      for(s in 1:S){
        X[m,n,s] <- e[m,n,s]
        if(m > 1) X[m,n,s] <- X[m,n,s] + 0.33*e[m-1,n,s]
        if(n > 1) X[m,n,s] <- X[m,n,s] + 0.33*e[m,n-1,s]
        if(s > 1) X[m,n,s] <- X[m,n,s] + 0.33*e[m,n,s-1]
      }
    }
  }
  X
}

generate_synthetic <- function(M, N, S) {
  A <- c(10, 5)
  B <- c(10, 5)
  lambda <- c(2.0, 1.5)
  mu     <- c(2.0, 1.5)
  nu     <- c(2.0, 1.5)

  y <- array(0, dim = c(M,N,S))

  for(k in 1:2){
    for(m in 1:M){
      for(n in 1:N){
        for(s in 1:S){
          phi <- lambda[k]*m + mu[k]*n + nu[k]*s
          y[m,n,s] <- y[m,n,s] +
            A[k]*cos(phi) + B[k]*sin(phi)
        }
      }
    }
  }

  y + generate_noise(M,N,S, sd = 1)
}

y <- generate_synthetic(M,N,S)

# ---- estimate ----
est <- estimate_seq(M,N,S,y,K)

# ---- print results ----
cat("\nEstimated parameters:\n")
for(k in 1:K){
  cat("\nComponent", k, "\n")
  cat("A      =", est$A[k], "\n")
  cat("B      =", est$B[k], "\n")
  cat("lambda =", est$lambda[k], "\n")
  cat("mu     =", est$mu[k], "\n")
  cat("nu     =", est$nu[k], "\n")
}

# ---- true values (for reference) ----
cat("\nTrue parameters:\n")
cat("Comp 1: A=10, B=10, lambda=2.0, mu=2.0, nu=2.0\n")
cat("Comp 2: A=5,  B=5,  lambda=1.5, mu=1.5, nu=1.5\n")
