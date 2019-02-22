#' Predict from catch-at-age model
#'
#' 
#' @param stoxdata list with input data
#' @param win indicator variable, =1 if windows platform, 0 otherwise
#' @keywords c-program mcmc
#' @export
## Created 30.10.2017
## Modified
##
run.pred <- function(stoxdata,common,win)
{
  nAges<-1+stoxdata$GlobalParameters$maxage-stoxdata$GlobalParameters$minage
  ncat<-nAges
  if(stoxdata$GlobalParameters$CC)ncat<-ncat*2
  
  pc.age<-make.predict.cell(ncat,stoxdata$AgeLength,stoxdata$Landings$AgeLengthCov)
  pc.lga<-make.predict.cell(1,stoxdata$AgeLength,stoxdata$Landings$AgeLengthCov)
  pc.wgl<-make.predict.cell(1,stoxdata$WeightLength,stoxdata$Landings$WeightLengthCov)
  
  
  common.par<-make.caa.params(stoxdata,common)
  data.catch<-make.data.catch(stoxdata,pc.age,pc.wgl)
  dist.cell<-make.dist.cell(pc.age,pc.lga,pc.wgl,stoxdata$GlobalParameters$CC)
  
  
  # TEMPORARY
  common.par.predict<<-common.par
  data.catch<<-data.catch
  dist.cell<<-dist.cell
  
  #write.bin.C.list(common.par,common$inputfolder,"common_par_predict")
  write.ascii.C.list(common.par,common$inputfolder,"common_par_predict_ascii")
  write.bin.C.list(data.catch,common$inputfolder,"data_catch")
  write.bin.C.list(dist.cell,common$inputfolder,"dist_cell")

  if(win){
    caa.call <- system.file("bin/caa_main_predict.exe",package="Reca")
  } else {
    caa.call <- system.file("bin/caa_main_predict",package="Reca")
  }
  s <- system(paste(shQuote(caa.call)," ",common$inputfolder,sep=""), intern=TRUE)
  filename <- paste0(common$inputfolder,"log.txt")
  if(!file.exists(filename)) stop("Error in prediction")
  fp <- file(filename,"r")
  log <- readLines(fp)
  close(fp)
  if(log[length(log)]=="OK")cat("Success in prediction\n")
  else stop("Error in prediction")
  file.remove(filename)  # To make sure that you don't read an old file
  
}

