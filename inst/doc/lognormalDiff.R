## ----eval=FALSE, include=FALSE------------------------------------------------
#  # twDev::genVigs()
#  #rmarkdown::render("lognormalDiff.Rmd","md_document")

## ----setup, include=FALSE-----------------------------------------------------
library(knitr)
opts_chunk$set(out.extra = 'style="display:block; margin: auto"'
    #, fig.align = "center"
    #, fig.width = 4.6, fig.height = 3.2
    , fig.width = 6, fig.height = 3.75 #goldener Schnitt 1.6
    , dev.args = list(pointsize = 10)
    , dev = c('png','pdf')
    )
knit_hooks$set(spar = function(before, options, envir) {
    if (before) {
        par( las = 1 )                   #also y axis labels horizontal
        par(mar = c(2.0,3.3,0,0) + 0.3 )  #margins
        par(tck = 0.02 )                          #axe-tick length inside plots             
        par(mgp = c(1.1,0.2,0) )  #positioning of axis title, axis labels, axis
     }
})
library(lognorm) 
if (!require(ggplot2) || !require(dplyr) || !require(tidyr) || !require(purrr)) {
	print("To generate this vignette, ggplot2, dplyr, tidyr, and purrr are required.")
	knit_exit()
}
themeTw <- ggplot2::theme_bw(base_size = 10) + 
  theme(axis.title = element_text(size = 9))

## -----------------------------------------------------------------------------
# generate nSample values of two lognormal random variables
mu1 = log(110)
mu2 = log(100)
sigma1 = 0.25
sigma2 = 0.15
#(coefDiff <- estimateSumLognormal( c(mu1,mu2), c(sigma1,sigma2) ))
(coefDiff <- estimateDiffLognormal(mu1,mu2, sigma1,sigma2, 0))
(expDiff <- getLognormMoments(coefDiff["mu"], coefDiff["sigma"])[,"mean"] - 
    coefDiff["shift"])

## -----------------------------------------------------------------------------
getLognormMoments(coefDiff["mu"], coefDiff["sigma"], shift = coefDiff["shift"])
getLognormMode(coefDiff["mu"], coefDiff["sigma"], shift = coefDiff["shift"])
getLognormMedian(coefDiff["mu"], coefDiff["sigma"], shift = coefDiff["shift"])

## -----------------------------------------------------------------------------
p <- seq(0,1,length.out = 100)[-c(1,100)]
dsPredY <- data.frame(
  var = "y", 
  q_shifted = qlnorm(p, coefDiff["mu"], coefDiff["sigma"] )
) %>%
mutate( 
  q = q_shifted - coefDiff["shift"],
  d = dlnorm(q_shifted, coefDiff["mu"], coefDiff["sigma"])
)

## ----plotUncorr, echo=FALSE, fig.height=2.04, fig.width=3.27------------------
nSample = 2000

ds <- data.frame(
  a = rlnorm(nSample, mu1, sigma1)
  , b = rlnorm(nSample, mu2, sigma2)
) %>%  mutate(
  y = a - b
)
dsw <- gather(ds, "var", "value", a, b, y)
p1 <- dsw %>% filter(var == "y") %>% 
  ggplot(aes(value, color = var)) + geom_density(linetype = "dotted") +
  geom_vline(xintercept = mean(ds$y), linetype = "dotted")
p1 + geom_line(data = dsPredY, aes(q, d, color = var)) +
  geom_vline(aes(xintercept = expDiff, color = var), 
             data = data.frame(var = "y", q = expDiff)) +
themeTw +
theme(legend.position = c(0.98,0.98), legend.justification = c(1,1)) +
theme(axis.title.x = element_blank())

## -----------------------------------------------------------------------------
  mu1 = log(120)
  mu2 = log(60)
  sigma1 = 0.25
  sigma2 = 0.15
  coefDiff <- estimateDiffLognormal( mu1,mu2,sigma1,sigma2, corr = -0.8 )
  pLo <- plnorm(0 + coefDiff["shift"], coefDiff["mu"], coefDiff["sigma"])
  pSample <- pDiffLognormalSample(mu1,mu2,sigma1,sigma2, corr = -0.8)
  c(pLo = as.numeric(pLo), pSample = pSample)

## ----plotSignifZero, echo=FALSE-----------------------------------------------
p <- seq(0,1,length.out = 100)[-c(1,100)]
dsPredY <- tibble(
  var = "y", 
  q_shifted = qlnorm(p, coefDiff["mu"], coefDiff["sigma"] ),
  q = q_shifted - coefDiff["shift"],
  d = dlnorm(q_shifted, coefDiff["mu"], coefDiff["sigma"])
)
dsPredA <- tibble(
  var = "a", 
  q = qlnorm(p, mu1, sigma1),
  d = dlnorm(q, mu1, sigma1)
)
dsPredB <- tibble(
  var = "b", 
  q = qlnorm(p, mu2, sigma2),
  d = dlnorm(q, mu2, sigma2)
)
bind_rows(select(dsPredY, -q_shifted), dsPredA, dsPredB) %>% 
  #filter(var == "a") %>% 
  ggplot(aes(q, d, color = var)) + geom_line() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
themeTw +
theme(legend.position = c(0.98,0.98), legend.justification = c(1,1)) +
theme(axis.title.x = element_blank())

## -----------------------------------------------------------------------------
if (!requireNamespace("mvtnorm")) {
  warning("Remainder of the vignette required mvtnorm installed.")
  knitr::opts_chunk$set(error = TRUE) 
}
corr = 0.8
(coefDiff <- estimateDiffLognormal(mu1,mu2, sigma1,sigma2, corr = corr))
(expDiff <- getLognormMoments(
  coefDiff["mu"], coefDiff["sigma"], shift = coefDiff["shift"])[,"mean"])

## -----------------------------------------------------------------------------
nSample <- 1e5
sigma_vec = c(sigma1, sigma2)
corrM <- setMatrixOffDiagonals(
  diag(nrow = 2), value = corr, isSymmetric = TRUE)
covM <- diag(sigma_vec) %*% corrM %*% diag(sigma_vec)
xObsN <- exp(mvtnorm::rmvnorm(nSample, mean =  c(mu1, mu2), sigma = covM))
head(xObsN)
y = xObsN[,1] - xObsN[,2]

## ----pdfDiffCorr80, echo=FALSE, fig.height=2.04, fig.width=3.27---------------
p <- seq(0,1,length.out = 100)[-c(1,100)]
dsPredY <- data.frame(
  var = "y", 
  q_shifted = qlnorm(p, coefDiff["mu"], coefDiff["sigma"] )
) %>%
  mutate( 
    q = q_shifted - coefDiff["shift"],
    d = dlnorm(q_shifted, coefDiff["mu"], coefDiff["sigma"])
    )
# density plot of the random draws
ggplot(data.frame(y = y), aes(y, color = "random draws")) + 
  geom_density() +
  # line plot of the lognorm density approximation
  geom_line(data = dsPredY, aes(q, d, color = "computed diff")) +
  # expected value
  geom_vline(xintercept = expDiff, linetype = "dashed") +
themeTw +
theme(legend.position = c(0.98,0.98), legend.justification = c(1,1)) +
theme(axis.title.x = element_blank()) +
theme(legend.title = element_blank())

## -----------------------------------------------------------------------------
# generate nSample values of two lognormal random variables
mu1 = log(110)
mu2 = log(100)
sigma1 = 0.15
sigma2 = 0.25
try(coefDiff <- estimateDiffLognormal(mu1,mu2, sigma1,sigma2, 0))

## -----------------------------------------------------------------------------
# note the switch of positions of mu1 and mu2: mu2 - mu1
#(coefDiff <- estimateDiffLognormal(mu1,mu2, sigma1,sigma2, 0)
(coefDiff <- estimateDiffLognormal(mu2,mu1, sigma2,sigma1, 0))
# note the minus sign in front
(expDiff <- -(getLognormMoments(
  coefDiff["mu"], coefDiff["sigma"], shift = coefDiff["shift"])[,"mean"]))

## -----------------------------------------------------------------------------
p <- seq(0,1,length.out = 100)[-c(1,100)]
dsPredY <- data.frame(
  var = "y", 
  q_shifted_neg = qlnorm(p, coefDiff["mu"], coefDiff["sigma"] )
) %>%
  mutate( 
    # note the minus sign in front
    q = -(q_shifted_neg - coefDiff["shift"]),
    d = dlnorm(q_shifted_neg, coefDiff["mu"], coefDiff["sigma"])
    )

## ----plotDiffLargerVar, echo=FALSE, fig.height=2.04, fig.width=3.27-----------
nSample = 2000
ds <- data.frame(
  a = rlnorm(nSample, mu1, sigma1)
  , b = rlnorm(nSample, mu2, sigma2)
) %>%  mutate(
  y = a - b
)
dsw <- gather(ds, "var", "value", a, b, y)
p1 <- dsw %>% filter(var == "y") %>% 
  ggplot(aes(value, color = var)) + geom_density(linetype = "dotted") +
  geom_vline(xintercept = mean(ds$y), linetype = "dotted")
#

p1 + geom_line(data = dsPredY, aes(q, d, color = var)) +
  geom_vline(aes(xintercept = expDiff, color = var), 
             data = data.frame(var = "y", q = expDiff)) +
themeTw +
theme(legend.position = c(0.98,0.98), legend.justification = c(1,1)) +
theme(axis.title.x = element_blank())

