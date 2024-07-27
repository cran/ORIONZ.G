eap_grad_obli<-function(X, LAM, PHI, THRES, sigj,grid,index_nodos){

  # 2024 revision

  n <- size(X)[1]
  m <- size(X)[2]

  if (n>m) {ni <- n}
  else { ni <- m}



  ni   <- length(X)
  r    <- ncol(cbind(LAM)) # Works in case LAM is vector or matrix
  nnod <- length(grid)
  hnod <- nrow(index_nodos)


  d<-grid[2]-grid[1]

  L1 <- matrix(0,nnod^r,1)

  nume_th <- matrix(0,r,1)
  nume_se <- matrix(0,r,1)

  TH <- matrix(0,r,1)
  SE <- matrix(0,r,1)
  RELI <- matrix(0,r,1)

  deno <- 0

  p <- matrix(0,ni,1)

  for (h in 1:hnod){

    for (j in 1:ni){
      zi <- matrix(0,r,1)
      x <- 0

      for (i in 1:r){
        zi[i] <- grid[index_nodos[h,i]]
        x <- x + LAM[j,i] * grid[index_nodos[h,i]]
      }

      if ((X[j] -1) != 0){
        p1 <- pgrad(x, THRES[X[j] - 1, j], sigj[j])
      }

      p2 <- 1

      if (THRES[X[j],j] != 0){
        p2 <- pgrad(x, THRES[X[j],j], sigj[j])
      }

      if ((X[j]-1) == 0){
        p[j] <- p2
      }
      else {
        p[j] <- p2 - p1
      }

    }

    L1[h,1] <- 1

    for (j in 1:ni){
      L1[h,1] <- L1[h,1] * p[j]
    }

    # bivariate density for each grid point
    # theta mean and variance are 0 and 1
    term1 <- ordnormulti(zi, PHI)
    term2 <- term1 %*% d %*% d

    # EAP estimations and PSDs

    pp <- (L1[h,1]) * term2

    deno <- deno + pp

    for (i in 1:r){
      nume_th[i] <- nume_th[i] + pp * grid[index_nodos[h,i],1]
      nume_se[i] <- nume_se[i] + pp * grid[index_nodos[h,i],1] * grid[index_nodos[h,i],1]
    }

  }

  for (i in 1:r){
    TH[i,1] <- nume_th[i] / deno
    SE[i,1] <- sqrt(nume_se[i] / deno - (TH[i,1] * TH[i,1]))
    RELI[i,1] <- 1 - (nume_se[i] / deno - (TH[i,1] * TH[i,1]))
  }

  OUT  <- list("TH"=TH, "SE"=SE, "RELI"=RELI)
  return(OUT)
}
