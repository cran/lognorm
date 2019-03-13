## ----eval=FALSE, include=FALSE-------------------------------------------
#  # twDev::genVigs()
#  rmarkdown::render("aggregateCorrelated.Rmd","md_document")

## ----setup, include=FALSE------------------------------------------------
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
if (!require(ggplot2) || !require(dplyr) || !require(tidyr) || !require(purrr) ||
    !require(mvtnorm)
    ) {
  msg = paste0("To generate this vignette, ggplot2, dplyr, tidyr, purr,"
	            ," and mvtnorm are required.")
  print(msg)
	warning(msg)
  knitr::knit_exit()
}
themeTw <- theme_bw(base_size = 10) + 
  theme(axis.title = element_text(size = 9))

