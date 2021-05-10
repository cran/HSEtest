# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

HSE <- function (Tab){
  # HSE: Testing homogeneity of stratum effects in stratified binary Data
  # Tab: 2kx2 matrix composed of k 2x2 contingency tables across strata 1 to k
  # Output: test statistic and p-value
  k <- nrow(Tab)/2
  delta=sigtilde=sighat=w=N=rep(NA,k)
  for (i in 1:k){
    Tabi <- Tab[(2*i-1):(2*i),]
    N[i] <- sum(Tabi)
    pi12 <- Tabi[1,2]/N[i]
    pi21 <- Tabi[2,1]/N[i]
    delta[i] <- pi12-pi21
    sigtilde[i] <- (pi12+pi21-(pi12-pi21)^2)/N[i]
    sighat[i] <- (pi12+pi21)/N[i]
    if (sigtilde[i]!=0) w[i] <- 1/sigtilde[i]
  }
  if (any(sigtilde==0)) w <- 1
  A <- cbind (diag(k-1), -1)
  deltahat <- sum(w*delta)/sum(w)
  sighat <- sighat-deltahat^2/N
  negindex <- (sighat <= 0)
  sighat[negindex] <- sigtilde[negindex]
  Sigmahat <- A%*%diag(sighat)%*%t(A)
  Adelta <- A%*%delta
  T <- t(Adelta)%*%solve(Sigmahat)%*%Adelta
  p <- 1-pchisq(T, k-1)
  return(c(T=T, p=p))
}
