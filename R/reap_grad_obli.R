reap_grad_obli<-function(X, LAM, PHI, THRES, sigj, disp){

  X <- as.matrix(X)
  LAM <- as.matrix(LAM)
  PHI <- as.matrix(PHI)
  THRES <- as.matrix(THRES)
  sigj <- as.matrix(sigj)

  m<-size(LAM)[1]
  r<-size(LAM)[2]

  #20 nodes
  grid<-transpose(seq(-4,4,(8/20)))

  #grid <- transpose(seq(-4,4,(8/60)))

  # Number of nodes
  ni<-size(grid)[1]
  index_nodos <- Nseriesk(ni,r)

  #########
  n<-size(X)[1]
  m<-size(X)[2]

  th<-numeric()
  se<-numeric()
  reli<-numeric()
  #X<-as.numeric(X)

  ptm_one <- proc.time()

  for (i in 1:n){

    if (min(X)==0){
      X=X+1
    }

    out <- eap_grad_obli(X[i,],LAM,PHI,THRES,sigj,grid,index_nodos)

    th <- rbind(th,t(out$TH))
    se <- rbind(se,t(out$SE))
    reli <- rbind(reli,t(out$RELI))

    if (disp==TRUE){

     # if (i==1){
        compT <- proc.time() - ptm_one
        compT<-compT[3]
        compT<-compT*(n-i)/i

        secondsInAMinute = 60
        secondsInAnHour = 60 * secondsInAMinute
        secondsInADay = 24 * secondsInAnHour

        days <- floor(compT / secondsInADay)

        hourSeconds <- compT %% secondsInADay
        hours <- floor(hourSeconds / secondsInAnHour)

        minuteSeconds <- hourSeconds %% secondsInAnHour
        minutes <- floor(minuteSeconds / secondsInAMinute)

        remainingSeconds <- minuteSeconds %% secondsInAMinute
        seconds <- ceiling(remainingSeconds)

        if (compT > 3600){
          if (days >= 1){ #Very very rare, but just to be sure
            cat("Computing EAP scores Time remaining: +24 hours                                                       \r")
            flush.console()
          }
          else {
            if (hours == 1){
              cat("Computing EAP scores. Time remaining: ", hours,"hour, ",minutes, "minutes and ",seconds, "seconds  \r")
              flush.console()
            }
            else {
              cat("Computing EAP scores. Time remaining: ", hours,"hours, ",minutes, "minutes and ",seconds, "seconds \r")
              flush.console()
            }
          }
        }
        else{
          if (compT >= 60){
            cat("Computing EAP scores. Time remaining: ", minutes, "minutes and ",seconds,"seconds \r")
            flush.console()
          }
          if (compT < 60) {
            cat("Computing EAP scores. Time remaining ",seconds,"seconds                                                                  \r")
            flush.console()
          }
        }
    }
  }

  if (disp==TRUE){
    cat("\r","                                                                                                  ","\r")
  }

  OUT<-list("th"=th,"se"=se,"reli"=reli)
  return(OUT)

}
