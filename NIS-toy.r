source("penGAMfuns.r")
source("NISfuns.r")

install.packages(c("zoo", "fda", "grplasso","mgcv"))


library(zoo)
library(fda)
library(grplasso)
library(mgcv)

funind = 6      ###simulation setting, see NISfuns.r for details.

fit = simuNIS(funind = funind)
