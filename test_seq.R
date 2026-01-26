source("sequential_estimation.R")

generate_noise <- function(M, N, S, sd = 1) {

  e <- array(rnorm(M*N*S, 0, sd), dim = c(M,N,S))
  X <- array(0, dim = c(M,N,S))

  for(m in 1:M){
    for(n in 1:N){
      for(s in 1:S){

        X[m,n,s] <- e[m,n,s]

        if(m > 1) X[m,n,s] <- X[m,n,s] + 0.33 * e[m-1,n,s]
        if(n > 1) X[m,n,s] <- X[m,n,s] + 0.33 * e[m,n-1,s]
        if(s > 1) X[m,n,s] <- X[m,n,s] + 0.33 * e[m,n,s-1]
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
monte_carlo <- function(M, N, S, K, R = 100) {

  results <- data.frame()

  for(r in 1:R){

    cat("Run", r, "\n")

    y <- generate_synthetic(M,N,S)

    est <- estimate_seq(M,N,S,y,K)

    for(k in 1:K){
      results <- rbind(
        results,
        data.frame(
          run = r,
          component = k,
          A = est$A[k],
          B = est$B[k],
          lambda = est$lambda[k],
          mu = est$mu[k],
          nu = est$nu[k]
        )
      )
    }
  }

  results
}
true_params <- list(
  A  = c(10, 5),
  B  = c(10, 5),
  lambda = c(2.0, 1.5),
  mu     = c(2.0, 1.5),
  nu     = c(2.0, 1.5)
)

mc_summary <- function(M, N, S, K, R = 1) {

  store <- matrix(NA, nrow = R, ncol = 5*K)

  colnames(store) <- as.vector(
    sapply(1:K, function(k)
      c(paste0("A",k), paste0("B",k),
        paste0("lam",k), paste0("mu",k), paste0("nu",k)))
  )

  for(r in 1:R){
    y <- generate_synthetic(M,N,S)
    est <- estimate_seq(M,N,S,y,K)

    for(k in 1:K){
      store[r, (5*(k-1)+1):(5*k)] <-
        c(est$A[k], est$B[k],
          est$lambda[k], est$mu[k], est$nu[k])
    }
  }

  AE  <- colMeans(store)
  VAR <- apply(store, 2, var) * (R-1)/R

  MSE <- rep(NA, length(AE))
  for(k in 1:K){
    idx <- (5*(k-1)+1):(5*k)
    truth <- c(
      if(k <= length(true_params$A)) true_params$A[k] else 0,
      if(k <= length(true_params$B)) true_params$B[k] else 0,
      if(k <= length(true_params$lambda)) true_params$lambda[k] else 0,
      if(k <= length(true_params$mu)) true_params$mu[k] else 0,
      if(k <= length(true_params$nu)) true_params$nu[k] else 0
    )
    MSE[idx] <- colMeans((store[,idx] -
                          matrix(truth, R, 5, byrow=TRUE))^2)
  }

  # remove freq stats for overfitting case (Table 3)
  if(K > length(true_params$A)){
    for(k in (length(true_params$A)+1):K){
      idx <- (5*(k-1)+3):(5*k)
      MSE[idx] <- NA
      VAR[idx] <- NA
    }
  }

  data.frame(
    parameter = colnames(store),
    AE   = AE,
    MSE  = MSE,
    ASYV = VAR
  )
}

table1 <- rbind(
  cbind(M=10, mc_summary(10,10,3,K=1,R=100)),
  cbind(M=20, mc_summary(20,20,3,K=1,R=100)),
  cbind(M=30, mc_summary(30,30,3,K=1,R=100))
)

write.csv(table1, "table1_k1.csv", row.names = FALSE)

table2 <- rbind(
  cbind(M=10, mc_summary(10,10,3,K=2,R=100)),
  cbind(M=20, mc_summary(20,20,3,K=2,R=100)),
  cbind(M=30, mc_summary(30,30,3,K=2,R=100))
)

write.csv(table2, "table2_k2.csv", row.names = FALSE)
table3 <- rbind(
  cbind(M=10, mc_summary(10,10,3,K=3,R=100)),
  cbind(M=20, mc_summary(20,20,3,K=3,R=100)),
  cbind(M=30, mc_summary(30,30,3,K=3,R=100))
)

write.csv(table3, "table3_k3.csv", row.names = FALSE)
