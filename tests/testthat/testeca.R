#' Tests running eca formatted data. Simply runs eca.estiamate and eca.fit.
#' @param file file with eca formatted data
#' @param runfiledir direcotry to use for eca output
testEcaFormatted <- function(file, runfiledir){
  require(Reca)
  write(paste("Loading from file:", file), stderr())
  load(file)
  GlobalParameters$resultdir <- runfiledir
  if(!(file.exists(GlobalParameters$resultdir))){
    stop(paste("Directory", runfiledir, "does not exist."))
  }
  ## Estimate model
  fit <- eca.estimate(AgeLength,WeightLength,Landings,GlobalParameters)
  pred <- eca.predict(AgeLength,WeightLength,Landings,GlobalParameters)

  return (list(fit=fit, pred=pred))
}

runTest <- function(testfile, testfiles="./testfiles", tmpdir="./tmp"){
  write("\n", stdout())
  write("######", stdout())
  write(paste("Testing:", testfile), stdout())
  write("######", stdout())
  write("\n", stdout())
  
  testcurrent <- function(){
    if(!(file.exists(tmpdir))){
      stop(paste("Directory", tmpdir, "does not exist."))
    }
    if(!(file.exists(testfiles))){
      stop(paste("Directory", testfiles, "does not exist."))
    }
    ecadir <- gsub(".Rdata", "", testfile)
    if (!(file.exists(file.path(tmpdir, ecadir)))){
      dir.create(file.path(tmpdir, ecadir))
    }
    testEcaFormatted(file.path(testfiles,testfile), file.path(tmpdir, ecadir))
  
    
  }
  handle <- function(e){write(paste("Test", testfile, ": fail"), stdout())}
  
  ret <- testcurrent()


#  tryCatch(
#    {testcurrent()
#    write(paste("Test", testfile, ": sucsess"), stdout())
#    ret <- T
#    }, 
#    error = handle,
#    finally={
#      write(paste("Test", testfile, "finished"), stdout())
#      return(ret)
#      }
#    )
  return(TRUE)
}
  

#' @param tmpdir location where tests will generate output
runAllTests <- function(testfiles="./testfiles", tmpdir="./tmp"){
  successes <- c()
  failures <- c()
  if(!(file.exists(tmpdir))){
    stop(paste("Directory", tmpdir, "does not exist."))
  }
  if(!(file.exists(testfiles))){
    stop(paste("Directory", testfiles, "does not exist."))
  }
  for (testfile in list.files(testfiles)){
    success <- runTest(testfile, testfiles, tmpdir)
    if (success){
      successes <- c(successes, testfile)
    }
    else{
      failures <- c(failures, testfile)
    }
  }
  write("\n", stdout())
  write("######", stdout())
  write("Summary:", stdout())
  write("Tests that run successfully:", stdout())
  for (t in successes){
    write(t, stdout())
  }
  write("\n", stdout())
  write("Tests that failed:", stdout())
  for (t in failures){
    write(t, stdout())
  }
  write("######", stdout())
}
# If not already done, clone https://github.com/Sea2Data/Rstox_utils.
# Change directory to the folder "Work" in the local Rstox_utils (on Holmin's computer this is setwd("~/Code/Github/Rstox_utils/Rstox_utils/Work"))
# setwd(dirname(sys.frame(1)$ofile)) #only works when sourcing
if (!file.exists("testfiles")){stop("Set wd to location of script")}

# Run a preliminary test:
test_that("Preliminary test", {
  expect_equal(runTest("herring_2015_tempfixed_gearrandom_100samples.Rdata"), TRUE)
})

# Run the full test, which takes half an hour:
#runAllTests()
