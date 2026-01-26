library(png)

img <- readPNG("pattern_cell.png")

M <- dim(img)[1]
N <- dim(img)[2]
S <- 3

y <- img

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

png("figure2_periodogram_nu_pi_over_3.png", width = 800, height = 600)

persp(
  lambda_grid,
  mu_grid,
  I[,,1],
  theta = 30, phi = 25,
  xlab = expression(lambda),
  ylab = expression(mu),
  zlab = expression(I(lambda,mu,nu)),
  main = expression(nu == pi/3),
  ticktype = "detailed"
)

dev.off()

png("figure3_periodogram_nu_2pi_over_3.png", width = 800, height = 600)

persp(
  lambda_grid,
  mu_grid,
  I[,,2],
  theta = 30, phi = 25,
  xlab = expression(lambda),
  ylab = expression(mu),
  zlab = expression(I(lambda,mu,nu)),
  main = expression(nu == 2*pi/3),
  ticktype = "detailed"
)

dev.off()
