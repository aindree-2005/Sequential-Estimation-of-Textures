library(pso)

estimate_seq_PSO <- function(M, N, S, y_arr, K){

  grid <- expand.grid(m=1:M, n=1:N, s=1:S)

  y <- as.vector(y_arr)
  residual <- y

  A_hat   <- numeric(K)
  B_hat   <- numeric(K)
  lam_hat <- numeric(K)
  mu_hat  <- numeric(K)
  nu_hat  <- numeric(K)

  rss3 <- function(theta, yvec){
    lam <- theta[1]; mu <- theta[2]; nu <- theta[3]
    phi <- lam*grid$m + mu*grid$n + nu*grid$s
    C  <- cos(phi)
    Sn <- sin(phi)
    Z  <- cbind(C, Sn)
    coef <- solve(t(Z)%*%Z, t(Z)%*%yvec)
    fit <- Z %*% coef
    sum((yvec - fit)^2)
  }

  for(k in 1:K){

    cat("Estimating component", k, "\n")

    rss_idx <- function(p){
      i <- round(p[1])
      j <- round(p[2])
      l <- round(p[3])

      i <- max(1, min(M-1, i))
      j <- max(1, min(N-1, j))
      l <- max(1, min(S-1, l))

      lam <- i*pi/M
      mu  <- j*pi/N
      nu  <- l*pi/S

      phi <- lam*grid$m + mu*grid$n + nu*grid$s
      C  <- cos(phi)
      Sn <- sin(phi)
      Z  <- cbind(C, Sn)
      coef <- solve(t(Z)%*%Z, t(Z)%*%residual)
      fit <- Z %*% coef
      sum((residual - fit)^2)
    }

    pso_res <- psoptim(
      par = c(M/2, N/2, S/2),
      fn  = rss_idx,
      lower = c(1,1,1),
      upper = c(M-1, N-1, S-1),
      control = list(
        maxit = 80,
        s = 40,
        w = 0.7,
        c.p = 1.5,
        c.g = 1.5,
        trace = FALSE
      )
    )

    idx <- round(pso_res$par)
    p1 <- c(idx[1]*pi/M, idx[2]*pi/N, idx[3]*pi/S)

    opt1 <- optim(
      p1, rss3,
      yvec = residual,
      method = "Nelder-Mead",
      control = list(maxit = 200)
    )

    p1 <- pmin(pmax(opt1$par, 1e-6), pi - 1e-6)

    lam <- p1[1]; mu <- p1[2]; nu <- p1[3]

    phi <- lam*grid$m + mu*grid$n + nu*grid$s
    C  <- cos(phi)
    Sn <- sin(phi)
    Z  <- cbind(C, Sn)

    coef <- solve(t(Z)%*%Z, t(Z)%*%residual)

    A_hat[k]   <- coef[1]
    B_hat[k]   <- coef[2]
    lam_hat[k] <- lam
    mu_hat[k]  <- mu
    nu_hat[k]  <- nu

    residual <- residual - Z %*% coef
  }

  list(A=A_hat, B=B_hat, lambda=lam_hat, mu=mu_hat, nu=nu_hat)
}
