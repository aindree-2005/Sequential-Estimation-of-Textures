library(png)
library(DEoptim)

source("de_seq.R")

img <- readPNG("weave_comp.png")
img <- img[,,1:3]

M <- dim(img)[1]
N <- dim(img)[2]
S <- 3
K <- 200

y <- img
y <- y - mean(y)

cat("Image size:", M, "x", N, "x", S, "\n")

est <- estimate_seq_DE(M, N, S, y, K)

y_rec <- array(0, dim = c(M,N,S))

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
        y_rec[m,n,s] <- y_rec[m,n,s] +
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

png("weave_reconstruction_de.png", width=900, height=400)
par(mfrow=c(1,2), mar=c(1,1,2,1))

plot(
  as.raster(clip01(y + mean(img))),
  main="(a) Original Texture"
)

plot(
  as.raster(clip01(y_rec + mean(img))),
  main=paste("(b) Reconstructed Texture (K =", K, ")")
)

dev.off()
