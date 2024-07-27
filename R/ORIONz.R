ORIONz <- function(X, LAM, PHI, THRES, Ree, PSI, model = "linear", disp = TRUE){

  X <- as.matrix(X)
  LAM <- as.matrix(LAM)
  PHI <- as.matrix(PHI)
  if (model == "graded"){
    THRES <- as.matrix(THRES)
    m <- size(THRES)[2]
    THRES <- rbind(THRES,matrix(0,1,m))
  }
  Ree <- as.matrix(Ree)
  PSI <- as.matrix(PSI)
  sigj <- diag(PSI)

  # 1: Calibration

  if (model == "linear"){
    out_calibration <- eapcondef(X, LAM, PHI, Ree, PSI)
  }
  if (model == "graded"){
    out_calibration <- reap_grad_obli(X, LAM, PHI, THRES, sigj, disp)
  }

  # 2: Correction

  out_correction <- relicorrec(LAM, PHI, Ree, PSI)

  # 3: Applying correction

  r<- size(LAM)[2]
  n <- size(X)[1]
  th_corrected <- matrix(NA, n, r)
  if (model == "graded"){
    reli_corrected <- matrix(NA, n, r)
    se_corrected <- matrix(NA, n, r)
  }


  for (i in 1:r){
    th_corrected[,i] <- out_calibration$th[,i] * out_correction$corrterm[i]
    if (model == "graded"){
      reli_corrected[,i] <- out_calibration$reli[,i] * out_correction$corrterm[i]
      se_corrected[,i] <- sqrt(1-reli_corrected[,i])
    }
  }


  if (model=="linear"){
    OUT <- list("th_corrected"=th_corrected, "marginal_reli_corrected" = out_correction$reliunbis)
  }
  if (model=="graded"){
    OUT <- list("th_corrected" = th_corrected, "reli_corrected" = reli_corrected, "se_corrected" = se_corrected,"marginal_reli_corrected" = colMeans(reli_corrected))


  }
  if (disp==TRUE){return(OUT)}
  else {invisible(OUT)}
}
