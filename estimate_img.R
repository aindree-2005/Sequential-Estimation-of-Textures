library(png)
source("sequential_estimation.R")

img <- readPNG("weave_comp.png")
img <- img[,,1:3] 
M <- dim(img)[1]
N <- dim(img)[2]
S <- 3

y <- img
y <- y - mean(y)   

cat("Image size:", M, "x", N, "x", S, "\n")

lambda_grid <- (1:(M-1)) * pi / M
mu_grid     <- (1:(N-1)) * pi / N
nu_grid     <- (1:(S-1)) * pi / S

I <- array(0, dim = c(M-1, N-1, S-1))

for(si in 1:(S-1)){
  nu <- nu_grid[si]

  for(i in 1:(M-1)){
    lam <- lambda_grid[i]

    for(j in 1:(N-1)){
      mu <- mu_grid[j]

      re <- 0
      im <- 0

      for(m in 1:M){
        for(n in 1:N){
          for(s in 1:S){
            phase <- m*lam + n*mu + s*nu
            re <- re + y[m,n,s] * cos(phase)
            im <- im - y[m,n,s] * sin(phase)
          }
        }
      }

      I[i,j,si] <- (re^2 + im^2) / (M*N*S)
    }
  }
}

K <- 200

vals <- c(I[,,1], I[,,2])
idx <- expand.grid(
  i = 1:(M-1),
  j = 1:(N-1),
  s = 1:(S-1)
)

ord <- order(vals, decreasing = TRUE)
top_idx <- idx[ord[1:K], ]

lambda_init <- lambda_grid[top_idx$i]
mu_init     <- mu_grid[top_idx$j]
nu_init     <- nu_grid[top_idx$s]

y_init <- array(0, dim = c(M,N,S))
y_vec <- as.vector(y)

for(k in 1:K){

  lam <- lambda_init[k]
  mu  <- mu_init[k]
  nu  <- nu_init[k]

  C <- numeric(M*N*S)
  Sg <- numeric(M*N*S)

  idx2 <- 1
  for(m in 1:M){
    for(n in 1:N){
      for(s in 1:S){
        phi <- m*lam + n*mu + s*nu
        C[idx2]  <- cos(phi)
        Sg[idx2] <- sin(phi)
        idx2 <- idx2 + 1
      }
    }
  }

  Z <- cbind(C, Sg)
  coef <- solve(t(Z) %*% Z, t(Z) %*% y_vec)

  idx2 <- 1
  for(m in 1:M){
    for(n in 1:N){
      for(s in 1:S){
        y_init[m,n,s] <- y_init[m,n,s] +
          coef[1]*C[idx2] + coef[2]*Sg[idx2]
        idx2 <- idx2 + 1
      }
    }
  }
}

set.seed(1)
est <- estimate_seq(M, N, S, y, K)
est2 <- estimate_seq(M, N, S, y, K)
y_final <- array(0, dim = c(M,N,S))

for(k in 1:K){
  lam <- est$lambda[k]
  mu  <- est$mu[k]
  nu  <- est$nu[k]
  A   <- est$A[k]
  B   <- est$B[k]

  for(m in 1:M){
    for(n in 1:N){
      for(s in 1:S){
        phi <- m*lam + n*mu + s*nu
        y_final[m,n,s] <- y_final[m,n,s] +
          A*cos(phi) + B*sin(phi)
      }
    }
  }
}
clip01 <- function(x) {
  x[x < 0] <- 0
  x[x > 1] <- 1
  x
}
png("figure5_weave.png", width=1200, height=400)
par(mfrow=c(1,3), mar=c(1,1,2,1))

plot(as.raster(clip01(y + mean(img))), main="(a) Original Texture")
plot(as.raster(clip01(y_init + mean(img))), main="(b) Initial Estimated Texture")
plot(as.raster(clip01(y_final + mean(img))), main="(c) Final Estimated Texture")

dev.off()
