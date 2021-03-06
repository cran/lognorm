#' Estimate the parameters of the lognormal approximation to the sum
#'
#' @describeIn estimateSumLognormalSample
#'   In addition to \code{estimateSumLognormal} take care of missing values
#'   and estimate correlation terms.
#' @param mu numeric vector of center parameters of terms at log scale
#' @param sigma numeric vector of scale parameter of terms at log scale
#' @param resLog time series of model-residuals at log scale 
#' to estimate correlation
#' @param effAcf effective autocorrelation
#' coefficients (may provide precomputed for efficiency or if the sample
#' of \code{resLog} is too small) set to 1 to assume uncorrelated sample
#' @param isGapFilled logical vector whether entry is gap-filled 
#' rather than an original measurement, see details
#' @param na.rm neglect terms with NA values in mu or sigma
#' 
#' @details 
#' If there are no gap-filled values, i.e. \code{all(!isGapFilled)} or
#' \code{!length(isGapFilled)} (the default), distribution parameters
#' are estimated using all the samples. Otherwise, the scale parameter
#' (uncertainty) is first estimated using only the non-gapfilled records.
#' 
#' Also use isGapFilled == TRUE for records, where sigma cannot be trusted. 
#' When setting sigma to missing, this is also affecting the expected value.
#' @details If there are only gap-filled records, 
#' assume uncertainty to be 
#' (before v0.1.5: the largest uncertainty of given gap-filled records.)
#' the mean of the given multiplicative standard deviation
#'
#' @return numeric vector with components \code{mu}, \code{sigma}, 
#' and \code{nEff},
#' i.e. the parameters of the lognormal distribution at log scale
#' and the number of effective observations.
#' @export
estimateSumLognormalSample <- function(
  mu, sigma, resLog  
  , effAcf = computeEffectiveAutoCorr(resLog) 
  , isGapFilled = logical(0) 
  , na.rm = TRUE  ##<< 
){
  # only one term, return the parameters
  if (length(mu) == 1 ) return(
    c(mu = as.vector(mu), sigma = as.vector(sigma), nEff = 1)
  )
  if (length(sigma) == 1) sigma <- rep(sigma, length(mu))
  nEff <- computeEffectiveNumObs(resLog, effAcf = effAcf, na.rm = na.rm)
  if (length(isGapFilled) && any(isGapFilled)) {
    isMeasured <- !isGapFilled
    sigmaSum <- if (!sum(isMeasured)) {
      #max(sigma)
      log(mean(exp(sigma), na.rm = TRUE))
    } else {
      # recursive call with only the measured records
      # set non-measured to NA but keep structure for autocorrelation distances
      nTerm <- length(mu)
      muMeas <- sigmaMeas <- resLogMeas <- rep(NA_real_, nTerm)
      muMeas[isMeasured] <- mu[isMeasured]
      sigmaMeas[isMeasured] <- sigma[isMeasured]
      resLogMeas[isMeasured] <- resLog[isMeasured]
      ans1 <- estimateSumLognormalSample(
        muMeas, sigmaMeas, resLogMeas
        ,effAcf = effAcf, isGapFilled = logical(0)
        ,na.rm = na.rm
      )
      ans1["sigma"]
    }
    p <- estimateSumLognormal(
      # sigma will not be used only checked for is.finite()
      mu, sigma = sigma, sigmaSum = sigmaSum, effAcf = effAcf, na.rm = na.rm)
    return(c(p , nEff = nEff))
  }
  p <- estimateSumLognormal( mu, sigma, effAcf = effAcf, na.rm = na.rm)
  return(c(p , nEff = nEff))
}

#' Estimate the parameters of the lognormal approximation to the sum 
#'
#' @describeIn estimateSumLognormalSample
#'   Before calling \code{estimateSumLognormalSample} estimate
#'   lognormal parameters from value and its uncertainty given
#'   on original scale.
#' @param mean numeric vector of expected values
#' @param sigmaOrig numeric vector of standard deviation at original scale
#' @param ... further arguments passed to \code{estimateSumLognormalSample}
#' @export
estimateSumLognormalSampleExpScale <- function(mean, sigmaOrig, ...){
  parlog <- getParmsLognormForMoments(mean, sigmaOrig = sigmaOrig)
  estimateSumLognormalSample( parlog[,"mu"], parlog[,"sigma"], ...)
}  

#' @describeIn estimateSumLognormalSample
#'  Estimate the parameters of the lognormal approximation to the sum
#' @param corr numeric matrix of correlations between the random variables
#' @param sigmaSum numeric scalar: possibility to specify
#' a precomputed scale parameter instead of computing it.
#' @param corrLength integer 
#' scalar: set correlation length to smaller values
#' to speed up computation by neglecting correlations among terms
#' further apart.
#' Set to zero to omit correlations.
#' @param isStopOnNoTerm if no finite estimate is provided then by
#' default NA is returned for the sum.
#' Set this to TRUE to issue an error instead.
#' @param effAcf numeric vector of effective autocorrelation
#' This overrides arguments \code{corr} and \code{corrLength}
#'
#' @references \code{Lo C (2013) WKB approximation for the sum of two 
#' correlated lognormal random variables.
#' Applied Mathematical Sciences, Hikari, Ltd., 7 , 6355-6367 
#' 10.12988/ams.2013.39511}
#' @export
#' @examples
#'   # distribution of the sum of two lognormally distributed random variables
#'   mu1 = log(110)
#'   mu2 = log(100)
#'   sigma1 = log(1.2)
#'   sigma2 = log(1.6)
#'   (coefSum <- estimateSumLognormal( 
#'   c(mu1,mu2), c(sigma1,sigma2) ))
#'   # repeat with correlation
#'   (coefSumCor <- estimateSumLognormal( 
#'   c(mu1,mu2), c(sigma1,sigma2), effAcf = c(1,0.9) ))
#'   # expected value is equal, but variance with correlated variables is larger
#'   getLognormMoments(coefSum["mu"],coefSum["sigma"])
#'   getLognormMoments(coefSumCor["mu"],coefSumCor["sigma"])
estimateSumLognormal <- function(
  mu, sigma, effAcf = c()                
  , corr = Diagonal(length(mu)) 
  , corrLength = if (inherits(corr, "ddiMatrix")) 0 else nTerm  
  , sigmaSum = numeric(0) 
  , isStopOnNoTerm = FALSE 
  , na.rm = isStopOnNoTerm 
){
  lengthMu <- length(mu)
  if (length(sigma) == 1) sigma <- rep(sigma, lengthMu)
  iFinite <- which( is.finite(mu) & is.finite(sigma))
  if (!isTRUE(na.rm) && length(iFinite) != lengthMu) 
    return(c(mu = NA_real_, sigma = NA_real_))
  muFin <- mu[iFinite]
  sigmaFin <- sigma[iFinite]
  nTerm = length(muFin)
  if (nTerm == 0) if (isTRUE(isStopOnNoTerm)) {
    stop( "need to provide at least one term"
          ,", but length of finite mu and sigma was zero")
  } else {
    # if no finite term was given, return NA
    return(c(mu = NA_real_, sigma = NA_real_))
  }
  if (nTerm == 1) return(structure(c(muFin,sigmaFin), names = c("mu","sigma")))
  if (!missing(effAcf) && (length(effAcf) > 1)) {
    corr <- getCorrMatFromAcf(length(mu), effAcf)
    corrLength <- length(effAcf) - 1L
  } 
  corrFin <- corr[iFinite,iFinite]
  #S = exp(muFin) 
  # S denotes the expected value, not mu in Lo13
  S = exp(muFin + sigmaFin*sigmaFin/2) 
  Ssum = sum(S)
  sigma2Eff <- if (length(sigmaSum)) {
    # if sigma of the sum has been pre-specified
    sigmaSum^2
  } else if (corrLength == 0) {
      sumDiag <- sum(sigmaFin*sigmaFin*S*S)/Ssum^2
  } else {
    jStarts <- pmax(1, (1:nTerm) - corrLength)
      jEnds <- pmin(nTerm, (1:nTerm) + corrLength)
    ansi = sapply(1:nTerm, function(i){
      j <- jStarts[i]:jEnds[i]
      ansj <- corrFin[i,j]*sigmaFin[i]*sigmaFin[j]*S[i]*S[j]  #Ssum outside loop
      # ansj2 <- sapply(jStarts[i]:jEnds[i], function(j){
      #   #corrFin[i,j]*sigmaFin[i]*sigmaFin[j]*S[i]/Ssum*S[j]/Ssum
      #   corrFin[i,j]*sigmaFin[i]*sigmaFin[j]*S[i]*S[j]  #Ssum outside loop
      # })
      # sum(ansj2)
      sum(ansj)
    })
    sum(ansi)/Ssum^2
  }
  muSum = log(Ssum) - sigma2Eff/2
  return(c(mu = as.vector(muSum), sigma = as.vector(sqrt(sigma2Eff))))
}


estimateSumLognormalBenchmark <- function(
  ### Estimate the distribution parameters of the lognormal approximation to the sum
  mu      ##<< numeric vector of center parameters
  , sigma  ##<< numeric vector of variance parameter at log sclae
  , corr = diag(nrow = nTerm) # TODO correlations
){
  ##details<< Implements estimation according to
  ## Messica A(2016) A simple low-computation-intensity model for approximating
  ## the distribution function of a sum of non-identical lognormals for
  ## financial applications. 10.1063/1.4964963
  nTerm = length(mu)
  S = exp(mu)
  Ssum = sum(S)
  # diag
  ansi = sapply(1:nTerm, function(i){
    ansj = sapply(i, function(j){
      #corr[i,j]*sigma[i]*sigma[j]*S[i]/Ssum*S[j]/Ssum
      corr[i,j]*sigma[i]*sigma[j]*S[i]*S[j]
    })
    sum(ansj)
  })
  sumDiag = sum(ansi)
  #
  ansi = sapply(1:nTerm, function(i){
    ansj = sapply(1:nTerm, function(j){
      #corr[i,j]*sigma[i]*sigma[j]*S[i]/Ssum*S[j]/Ssum
      corr[i,j]*sigma[i]*sigma[j]*S[i]*S[j]
    })
    sum(ansj)
  })
  sigma2Eff <- sum(ansi)/Ssum^2
  # for performance reasons take out Ssum of the sum
  #sigma2Eff <- sum(ansi)/(nTerm*Ssum)^2
  muSum = log(Ssum) + sigma2Eff/2
  return(c(mu = muSum, sigma = sqrt(sigma2Eff)))
}



