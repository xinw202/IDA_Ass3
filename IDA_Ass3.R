install.packages("mice")
library(mice)
install.packages("ggplot2")
library(ggplot2)

###Q1a
cc(nhanes)
nrow(cc(nhanes))
per_ic <- (nrow(nhanes)-nrow(cc(nhanes)))/nrow(nhanes)
per_ic

###Q1b
imps<- mice(nhanes, printFlag = FALSE, seed = 1)
imps
# predict bmi from age, hyp, and chl by the normal linear regression model
fits<- with(imps, lm(bmi ~ age + hyp + chl))
class(fits)
# result
ests_1 <- pool(fits)
ests_1

###Q1c
# Repeat the analysis for seed2,3,4,5,6
ests_2<-pool(with(mice(nhanes, printFlag = FALSE, seed = 2),lm(bmi ~ age + hyp + chl)))
ests_2
ests_3<-pool(with(mice(nhanes, printFlag = FALSE, seed = 3),lm(bmi ~ age + hyp + chl)))
ests_3
ests_4<-pool(with(mice(nhanes, printFlag = FALSE, seed = 4),lm(bmi ~ age + hyp + chl)))
ests_4
ests_5<-pool(with(mice(nhanes, printFlag = FALSE, seed = 5),lm(bmi ~ age + hyp + chl)))
ests_5
ests_6<-pool(with(mice(nhanes, printFlag = FALSE, seed = 6),lm(bmi ~ age + hyp + chl)))
ests_6

###Q1d
# Repeat the analysis with M = 100 with the same seeds
ests_1_100<-pool(with(mice(nhanes, printFlag = FALSE, seed = 1, m = 100),lm(bmi ~ age + hyp + chl)))
ests_1_100
ests_2_100<-pool(with(mice(nhanes, printFlag = FALSE, seed = 2, m = 100),lm(bmi ~ age + hyp + chl)))
ests_2_100
ests_3_100<-pool(with(mice(nhanes, printFlag = FALSE, seed = 3, m = 100),lm(bmi ~ age + hyp + chl)))
ests_3_100
ests_4_100<-pool(with(mice(nhanes, printFlag = FALSE, seed = 4, m = 100),lm(bmi ~ age + hyp + chl)))
ests_4_100
ests_5_100<-pool(with(mice(nhanes, printFlag = FALSE, seed = 5, m = 100),lm(bmi ~ age + hyp + chl)))
ests_5_100
ests_6_100<-pool(with(mice(nhanes, printFlag = FALSE, seed = 6, m = 100),lm(bmi ~ age + hyp + chl)))
ests_6_100

###Q2
# stochastic regression imputation
load("dataex2.Rdata")
num_sir<-0
n<-nrow(dataex2)
for (i in 1:n){
  imps_sri <- mice(dataex2[,,i], methods = "norm.nob", m = 20, seed = 1, printFlag = FALSE)
  fits_sri <- with(imps_sri, lm(Y~X))
  ests_sri <- pool(fits_sri)
  summary_sri <- summary(ests_sri, conf.int = TRUE)
  if (summary_sri[2,c(7)]<=3 & summary_sri[2,c(8)]>=3){
    num_sir <- num_sir + 1
  }
}
ecp_sri <- num_sir/n
ecp_sri
# bootstrap based version
num_norm_b <- 0
n <- nrow(dataex2)
for (i in 1:n){
  imps_norm_b <- mice(dataex2[,,i], methods = "norm.boot", m = 20, seed = 1, printFlag = FALSE)
  fits_norm_b <- with(imps_norm_b, lm(Y~X))
  ests_norm_b <- pool(fits_norm_b)
  summary_norm_b <- summary(ests_norm_b, conf.int = TRUE)
  if (summary_norm_b[2,c(7)]<=3 & summary_norm_b[2,c(8)]>=3){
    num_norm_b <- num_norm_b + 1
  }
}
ecp_norm_b <- num_norm_b/n
ecp_norm_b


###Q4a
load('dataex4.Rdata')
imps_4a <- mice(dataex4, printFlag = FALSE, seed = 1, m = 50)
fits_4a <- with(imps_4a, lm(y ~ x1 + x2 + x1*x2))
ests_4a <- pool(fits_4a)
summary_4a <- summary(ests_4a, conf.int = TRUE)
summary_4a[,c(1,2,3,7,8)]

###Q4b
x1 <- dataex4$x1
x2 <- dataex4$x2
dataex4$x1x2 <- x1*x2

imps_4b <- mice(dataex4, maxit = 0)
meth <- imps_4b$method
meth["x1x2"] <- "~I(x1*x2)"

pred <- imps_4b$predictorMatrix
pred[c("x1", "x2"), "x1x2"] <- 0
pred
visSeq <- imps_4b$visitSequence
visSeq
imps_pi <- mice(dataex4, method = meth, predictorMatrix = pred, visitSequence = visSeq,
                printFlag = FALSE, m = 50, seed = 1)
fits_pi <- with(imps_pi, lm(y ~ x1 + x2 + x1*x2))
ests_pi <- pool(fits_pi)
summary_pi <- summary(ests_pi, conf.int = TRUE)
summary_pi[,c(1,2,3,7,8)]

###4c
imps_4c <- mice(dataex4, printFlag = FALSE, m = 50, seed = 1)
fits_4c <- with(imps_4c, lm(y ~ x1 + x2 + x1x2))
ests_4c <- pool(fits_4c)
summary_4c <- summary(ests_4c, conf.int = TRUE)
summary_4c[,c(1,2,3,7,8)]


###Q5
load("NHANES2.Rdata")
dim(NHANES2)
str(NHANES2)
summary(NHANES2)
mdpat_mice <- md.pattern(NHANES2)
mdpat_mice
library(JointAI)
md_pattern(NHANES2, pattern = FALSE, color = c('#34111b', '#e30f41'))

par(mar = c(3,3,2,1), mgp = c(2, 0.6, 0))
plot_all(NHANES2, breaks = 30, ncol = 4)

imp0 <- mice(NHANES2, maxit = 0)
imp0

meth <- imp0$method
meth["hgt"] <- "norm"
meth

post <- imp0$post
post["hgt"] <- "imp[[j]][,i] <- squeeze(imp[[j]][,i], c(0, 2.5))"

imp <- mice(NHANES2, method = meth,
            maxit = 20, m = 20, seed = 1, printFlag = FALSE)
imp$loggedEvents

plot(imp, layout = c(4,4))

densityplot(imp)

densityplot(imp, ~SBP|hypten+gender)
densityplot(imp, ~hgt| gender)


require(devtools)
require(reshape2)
require(RColorBrewer)
require(ggplot2)
source_url("https://gist.githubusercontent.com/NErler/0d00375da460dd33839b98faeee2fdab/raw/c6f537ecf80eddcefd94992ec7926aa57d454536/propplot.R")
propplot(imp)

xyplot(imp, hgt ~ wgt | gender, pch = c(1, 20))

fit<- with(imp, lm(wgt~gender+age+hgt+WC))
summary(fit$analyses[[1]])

comp1 <- complete(imp, 1)
plot(fit$analyses[[1]]$fitted.values, residuals(fit$analyses[[1]]),
     xlab = "Fitted values", ylab = "Residuals")

plot(comp1$wgt ~ comp1$gender, xlab = "gender", ylab = "wgt")
plot(comp1$wgt ~ comp1$age, xlab = "age", ylab = "wgt")
plot(comp1$wgt ~ comp1$hgt, xlab = "Height in metres", ylab = "wgt")
plot(comp1$wgt ~ comp1$WC, xlab = "Waist circumference in cm", ylab = "wgt")

qqnorm(rstandard(fit$analyses[[1]]), xlim = c(-4, 4), ylim = c(-6, 6))
qqline(rstandard(fit$analyses[[1]]), col = 2)

pooled_ests <- pool(fit)
summary(pooled_ests, conf.int = TRUE)

pool.r.squared(pooled_ests, adjusted = TRUE)

fit_no_gender <- with(imp, lm(wgt ~ age + hgt + WC))
D1(fit, fit_no_gender)

fit_no_age<- with(imp, lm(wgt ~ gender + hgt + WC))
D1(fit, fit_no_age)

fit_no_hgt<- with(imp, lm(wgt ~ age+ gender + WC))
D1(fit, fit_no_hgt)

fit_no_WC<- with(imp, lm(wgt ~ age + gender + hgt))
D1(fit, fit_no_WC)

