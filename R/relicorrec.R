relicorrec<-function(L, PHI, Ree, PSI){

  L <- as.matrix(L)
  PHI <- as.matrix(PHI)
  Ree <- as.matrix(Ree)
  PSI <- as.matrix(PSI)

  m <- size(L)[1]
  D <- diag(m)

  SIG <- (L %*% PHI %*% transpose(L)) + (PSI %*% Ree %*% PSI)
  SIG <- SIG - diag(diag(SIG)) + diag(m)
  SIG2 <- (L %*% PHI %*% transpose(L)) + (PSI %*% D %*% PSI)
  SIG2 <- SIG2 - diag(diag(SIG2)) + diag(m)

  reliunbis <- diag(PHI %*% transpose(L) %*% solve(SIG) %*% L %*% PHI)

  relibis <- diag(PHI %*% transpose(L) %*% solve(SIG2) %*% L %*% PHI)

  corrterm <- reliunbis / relibis

  OUT  <- list("relibis"=relibis, "reliunbis"=reliunbis, "corrterm"=corrterm)
  return(OUT)
}
