## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#"
)

## ----install, eval = FALSE----------------------------------------------------
#  devtools:::install_github("gunhanb/MetaStan")

## ----dataset------------------------------------------------------------------
library("MetaStan")
data("dat.Berkey1995", package = "MetaStan")
head(dat.Berkey1995)

## ----forestplot, echo=TRUE, results=TRUE--------------------------------------
library(ggplot2)
# Calculating log odds ratios and variances from data
logodds <- function(x) log((x[1] * (x[4] - x[3]))/((x[2] - x[1]) * x[3]))
stdes   <- function(x) sqrt(1/x[1] + 1/(x[2] - x[1]) + 1/x[3] + 1/(x[4] - x[3]))
r_ind   <- apply(cbind(dat.Berkey1995$r2, dat.Berkey1995$n2, 
                 dat.Berkey1995$r1, dat.Berkey1995$n1), 1, logodds)
se_ind  <- apply(cbind(dat.Berkey1995$r2, dat.Berkey1995$n2, 
                 dat.Berkey1995$r1, dat.Berkey1995$n1), 1, stdes)
lower95_ind <- r_ind + qnorm(.025) * se_ind
upper95_ind <- r_ind + qnorm(.975) * se_ind
# Comparison of the results
trials  <- c("1", "2" ,"3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13")
trials <- ordered(trials, levels = trials)

d <- data.frame(x = trials,
                y = r_ind,
                ylo = lower95_ind,
                yhi = upper95_ind)
forest.plot <- ggplot(d, aes(x = x, y = y, ymin = ylo, ymax = yhi)) +
  geom_pointrange() +
  coord_flip() +
  geom_hline(aes(yintercept=0), lty = 2) +
  xlab("Studies") +
  ggtitle("Forest Plot (BCG vaccines)")

plot(forest.plot)

## ----bnhmFit, results="hide"--------------------------------------------------
 data('dat.Berkey1995', package = "MetaStan")
 ## Fitting a Binomial-Normal Hierarchical model using WIP priors
data('dat.Berkey1995', package = "MetaStan")
## Fitting a Binomial-Normal Hierarchical model using WIP priors
dat_MetaStan <- create_MetaStan_dat(dat = dat.Berkey1995,
                                    armVars = c(responders = "r", sampleSize = "n"))

 meta.BCG.stan  <- meta_stan(data = dat_MetaStan,
                           likelihood = "binomial",
                           mu_prior = c(0, 10),
                           theta_prior = c(0, 100),
                           tau_prior = 0.5,
                           tau_prior_dist = "half-normal")

## ----shinystan, eval = FALSE--------------------------------------------------
#  library("shinystan")
#  ## Firstly convert "stan" object to a "shinystan" object
#  bnhm1.BCG.shinystan = as.shinystan(meta.BCG.stan$fit)
#  launch_shinystan(bnhm1.BCG.shinystan)

## ----print--------------------------------------------------------------------
print(meta.BCG.stan)

