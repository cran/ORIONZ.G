eapcondef<-function(X, L, PHI, Ree, PSI){

  X <- as.matrix(X)
  L <- as.matrix(L)
  PHI <- as.matrix(PHI)
  Ree <- as.matrix(Ree)
  PSI <- as.matrix(PSI)

  n <- size(X)[1]
  m <- size(L)[1]

  M <- colMeans(X)
  Sx <- apply(X, 2, sd) #sd for each column

  DIF <- X - matrix(1,n,1) %*% M
  Z <- DIF / (matrix(1,n,1) %*% Sx)

  SIG <- (L %*% PHI %*% transpose(L)) + (PSI %*% Ree %*% PSI)
  SIG <- SIG - diag(diag(SIG)) + diag(m)

  th <- transpose(PHI %*% transpose(L) %*% solve(SIG) %*% transpose(Z))
  reli <- diag(PHI %*% transpose(L) %*% solve(SIG) %*% L %*% PHI)

  OUT<-list('th'=th,'reli'=reli)
  invisible(OUT)

}
