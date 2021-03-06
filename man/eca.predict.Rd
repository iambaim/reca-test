% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eca.main.functions.R
\name{eca.predict}
\alias{eca.predict}
\title{Predict catch-at-age}
\usage{
eca.predict(AgeLength, WeightLength, Landings, GlobalParameters)
}
\arguments{
\item{AgeLength}{A list with elements connected to the age-length model, see \link{eca.estimate}.}

\item{WeightLength}{A list with elements connected to the weight-length model, see \link{eca.estimate}.}

\item{Landings}{A list with elements connected to landings, see Details}

\item{GlobalParameters}{A list with input parameters that are given by the user, see Details}
}
\value{
A list with elements
\item{TotalCount}{array of the total catch (in numbers) with dimensions #length intervals x #age categories x #samples}
\item{MeanLength}{matrix of the mean length in the different age categories with dimensions #age categories x #samples}
\item{MeanWeight}{matrix of the mean weight in the different age categories with dimensions #age categories x #samples}
\item{AgeCategories}{vector with age categories}
\item{LengthIntervalsLog}{vector with length intervals, on log scale}
If a multiple stocks (coastal cod/Atlantic cod) analysis is run, the function returns the same list, but where the
dimension of age categories in all the elements is 2 times the number of age categories in the model (ncat). The first ncat 
valuels correspond to Atlantic cod while the last ncat values correspond to coastal cod.
}
\description{
This function runs the C prediction program and writes results to binary files.
}
\details{
Landings
 \itemize{
 \item LiveWeightKG  - Catch data given in kilograms for each cell.
 \item AgeLengthCov  - Covariate matrix for age-length model specifying which level of the different covariates the 
 corresponding cell in \emph{LiveWeightKG} has, e.g. constant, season, gearfactor, spatial. 
 It also contains a column \emph{midseason} which is the fraction of the year where the catch is taken.
 \item WeightLengthCov  - Covariate matrix for weight-length model specifying which level of the different covariates 
 the corresponding cell in \emph{LiveWeightKG} has, e.g. constant, season, gearfactor, spatial.
 It also contains a column \emph{midseason} which is the fraction of the year where the catch is taken.
 }
 GlobalParameters
 \itemize{
 \item minage - minimum age category
 \item maxage - maximum age category
 \item CC - TRUE if coastal cod is included in the model, FALSE otherwise
 \item resultdir - directory path where the results are stored. 
 Two directories \strong{cfiles} and \strong{resfiles} are in this folder (created if they don't exist). 
 In \strong{cfiles} the temporary binary input files to the c-program are stored, and in \strong{resfiles} the
 binary output files are stored.
 \item fitfile - name of the output files from estimating the model. 
 The files are given the extensions \strong{.mcmc1} and \strong{.mcmc2} for the age and length-given-age model 
 and for the weight-given-length model, respectively.
 \item predictfile - name of the output file from prediction. The file is given the extension \strong{.pred}.
 \item caa.burnin - number of MCMC samples discarded before running the prediction
 \item seed - random seed value
 }
}
\seealso{
\link{eca.estimate} for estimating the catch-at-age model.
}
\keyword{c-program}
\keyword{mcmc}
