install.packages("devtools")
devtools::install_github("vdorie/aciccomp/2017", upgrade = "never", force = TRUE)


library(aciccomp2017)


X1 <- dgp_2017(24, 52)
X2 <- dgp_2017(24, 49)
