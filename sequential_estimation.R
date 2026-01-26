estimate_seq <- function(M, N, S, y_arr, K) {

  # indices
  grid <- expand.grid(
    m = 1:M,
    n = 1:N,
    s = 1:S
  )

  y <- as.vector(y_arr)
  residual <- y

  # Fourier frequency grids
  lam_grid <- (1:(M-1)) * pi / M
  mu_grid  <- (1:(N-1)) * pi / N
  nu_grid  <- (1:(S-1)) * pi / S

  # storage
  A_hat   <- numeric(K)
  B_hat   <- numeric(K)
  lam_hat <- numeric(K)
  mu_hat  <- numeric(K)
  nu_hat  <- numeric(K)

  # RSS function for optim
  rss3 <- function(theta, yvec) {
    lam <- theta[1]; mu <- theta[2]; nu <- theta[3]
    phi <- lam*grid$m + mu*grid$n + nu*grid$s
    C <- cos(phi); S <- sin(phi)
    Z <- cbind(C,S)
    coef <- solve(t(Z)%*%Z, t(Z)%*%yvec)
    fit <- coef[1]*C + coef[2]*S
    sum((yvec - fit)^2)
  }

  for(k in 1:K){

    best_rss <- Inf
    p1 <- c(NA, NA, NA)

    # -------- GRID SEARCH --------
    for(lam in lam_grid){
      for(mu in mu_grid){
        for(nu in nu_grid){

          phi <- lam*grid$m + mu*grid$n + nu*grid$s
          C <- cos(phi); S <- sin(phi)
          Z <- cbind(C,S)
          coef <- solve(t(Z)%*%Z, t(Z)%*%residual)
          fit <- coef[1]*C + coef[2]*S
          rss <- sum((residual - fit)^2)

          if(rss < best_rss){
            best_rss <- rss
            p1 <- c(lam, mu, nu)
          }
        }
      }
    }
    opt1 <- optim(
      p1, rss3,
      yvec = residual,
      method = "Nelder-Mead",
      control = list(maxit = 500)
    )

    p1 <- pmin(pmax(opt1$par, 1e-3), pi - 1e-3)
    lam <- p1[1]; mu <- p1[2]; nu <- p1[3]
    phi <- lam*grid$m + mu*grid$n + nu*grid$s
    C <- cos(phi); S <- sin(phi)
    Z <- cbind(C,S)
    coef <- solve(t(Z)%*%Z, t(Z)%*%residual)

    A_hat[k]   <- coef[1]
    B_hat[k]   <- coef[2]
    lam_hat[k] <- lam
    mu_hat[k]  <- mu
    nu_hat[k]  <- nu
    residual <- residual - (coef[1]*C + coef[2]*S)
  }

  return(list(
    A = A_hat,
    B = B_hat,
    lambda = lam_hat,
    mu = mu_hat,
    nu = nu_hat
  ))
}
