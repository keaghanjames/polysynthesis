library(ordinalNet)
library(ivreg)
library(nloptr)

setwd("~/Dropbox/Lindell/Polysynthesis/plotting language trees")
load("bestmodel")
#poly.iso <- read.csv("data.csv",header=F)
# poly.iso <- as.character(unlist(poly.iso))
# poly.iso <- poly.iso[which(is.element(poly.iso,rownames(data)))]
# data$poly <- 0
# data[poly.iso,"poly"] <- 1
updated <- read.csv('~/Dropbox/Lindell/Polysynthesis/languoid_data_for_analysis.csv')
which((data$id_ISO_lang == updated$ISO) == F)
#data$poly <- updated$poly_definite
data$poly <- updated$pol_borderline

data$poly <- as.numeric(data$poly)
data$poly

#test if polysynhesis langauge has smaller L1 population size
#not accouting for autocorrelation
fit <- lm(lang_L1.POP_lang_tr~poly,data=data,na.action=na.exclude)
summary(fit)
            # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.03655    0.01259   2.903  0.00371 ** 
# poly        -0.67995    0.05430 -12.521  < 2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.9882 on 6509 degrees of freedom
# Multiple R-squared:  0.02352,	Adjusted R-squared:  0.02337 
# F-statistic: 156.8 on 1 and 6509 DF,  p-value: < 2.2e-16

#accounting for spatial and phylogenetic autocorrelation
X <- data$poly
y <- data$lang_L1.POP_lang_tr
Wspy <- Wsp%*%as.numeric(y)
Wnby <- Wnb%*%as.numeric(y)
Wphyy <- Wphy%*%as.numeric(y)
WspX <- Wsp%*%X
WnbX <- Wnb%*%X
WphyX <- Wphy%*%X
fit_sp <- ivreg(y ~ X + Wspy | X + WnbX )
summary(fit_sp)
# Coefficients:
            # Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.23534    0.05035  -4.674 3.01e-06 ***
# X           -0.02617    0.12729  -0.206    0.837    ##no association between L1speaker size and polysynthesis
# Wspy         0.75194    0.13590   5.533 3.27e-08 *** 

# Diagnostic tests:
                  # df1  df2 statistic p-value    
# Weak instruments    1 6508    51.257   9e-13 ***
# Wu-Hausman          1 6507     1.819   0.177    ##pass the test for endogeneity
# Sargan              0   NA        NA      NA    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.8619 on 6508 degrees of freedom
# Multiple R-Squared: 0.2573,	Adjusted R-squared: 0.2571 
# Wald test: 118.4 on 2 and 6508 DF,  p-value: < 2.2e-16 

fit_spphylo <- ivreg(y ~ X + Wspy + Wphyy | X + WnbX + WphyX)
summary(fit_spphylo)
# Coefficients:
            # Estimate Std. Error t value Pr(>|t|)
# (Intercept) -0.18043    0.14688  -1.228    0.219
# X            0.01529    0.09609   0.159    0.874  ##no association between L1speaker size and polysynthesis
# Wspy         0.50921    0.58746   0.867    0.386
# Wphyy        0.38195    0.78138   0.489    0.625

# Diagnostic tests:
                          # df1  df2 statistic p-value    
# Weak instruments (Wspy)     2 6507    44.615  <2e-16 ***
# Weak instruments (Wphyy)    2 6507    56.844  <2e-16 ***
# Wu-Hausman                  2 6505     0.063   0.939    ##pass the test for endogeneity
# Sargan                      0   NA        NA      NA    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.8085 on 6507 degrees of freedom
# Multiple R-Squared: 0.3466,	Adjusted R-squared: 0.3463 
# Wald test: 100.5 on 3 and 6507 DF,  p-value: < 2.2e-16 

#alternative appraoch
# fitspphylo <- function (y,X,Wsp,Wnb,Wphy,a) {
# 	cal <- function (a,y,X,Wsp,Wnb,Wphy) {
# 		W <- a[2]*(a[1]*Wsp+(1-a[1])*Wnb)+(1-a[2])*Wphy
#         Wy <- W%*%as.numeric(y)
#         WX <- W%*%X
#         fit_spphylo <- ivreg(y ~ X + Wy | X + WX )
#         -sum(log(dnorm(x=resid(fit_spphylo),mean=0,sd=fit_spphylo$sigma)))
#     }
#     cal2 <- function (a,y,X,Wsp,Wnb,Wphy) {
# 		W <- a[2]*(a[1]*Wsp+(1-a[1])*Wnb)+(1-a[2])*Wphy
#         Wy <- W%*%as.numeric(y)
#         WX <- W%*%X
#         fit_spphylo <- ivreg(y ~ X + Wy | X + WX )
#         }
#     p.res <- sbplx(a,cal,y=y,X=X,Wsp=Wsp,Wnb=Wnb,Wphy=Wphy,lower=c(0,0),upper=c(1,1))
#   lm.res <- cal2(a=p.res$par,y=y,X=X,Wsp=Wsp,Wnb=Wnb,Wphy=Wphy)
#   list(p.res=p.res,lm.res=lm.res)
#   
# }
# fit_spphylo <- fitspphylo(y=data$lang_L1.POP_lang_tr,X=data$poly,Wsp=Wsp,Wnb=Wnb,Wphy=Wphy,a=c(0.5,0.5))
# fit_spphylo$p.res
# # $par
# # [1] 0.1344412 0.4679864
# 
# # $value
# # [1] 7525.7
# summary(fit_spphylo$lm.res)
# Coefficients:
            # Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.14591    0.02241  -6.511 8.04e-11 ***
# X           -0.08939    0.07773  -1.150     0.25    ##no association between L1speaker size and polysynthesis
# Wy           0.90384    0.09985   9.052  < 2e-16 ***

# Diagnostic tests:
                  # df1  df2 statistic p-value    
# Weak instruments    1 6508   148.450  <2e-16 ***
# Wu-Hausman          1 6507     0.509   0.475    ##pass the test for endogeneity
# Sargan              0   NA        NA      NA    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.7688 on 6508 degrees of freedom
# Multiple R-Squared: 0.4091,	Adjusted R-squared: 0.4089 
# Wald test: 170.5 on 2 and 6508 DF,  p-value: < 2.2e-16 

#conclusion: polysynthetic languages have smaller L1 speaker size, after accounting for spatial and phylogenetic autocorrelation

#test if polysynthetic langauge has higher endangerment levels
#test for parallele assumption
X <- data$poly
X <- as.matrix(X,length(X),1)
y <- data$EGIDS_tr
res1 <- ordinalNet(X,y,standardize=F,family="cumulative",link="probit",nonparallelTerms=F,lambdaVals=0)
res2 <- ordinalNet(X,y,standardize=F,family="cumulative",link="probit",nonparallelTerms=T,lambdaVals=0)
likratio <- res2$loglik-res1$loglik #79.77414
likratio
p.value <- 1-pchisq(likratio,df=res2$nNonzero-res1$nNonzero) #1.631806e-12
p.value
coefs <- res2$coefs
coefs
#	overall  1-2		2-3			3-4		4-5			5-6a		6a-6b		6b-7	7-8a		8a-8b		8b-9	9-10
#-0.6822596 0.3996079 0.5572903 0.2914457 -0.004821173	0.3082924 -0.36416 -0.2599153 -0.1628272 -0.05543166 0.1510533 0.5843285
#conclusion: polysynthesis tends to decrease endangerment level, but its effect on the transition between endangerment levels differs. Since most of these difference are from level 1-6a, group 1-6a

#group 1-6a
y2 <- y
levels(y2)[1:6] <- "6a"
res1 <- ordinalNet(X,y2,standardize=F,family="cumulative",link="probit",nonparallelTerms=F,lambdaVals=0)
res2 <- ordinalNet(X,y2,standardize=F,family="cumulative",link="probit",nonparallelTerms=T,lambdaVals=0)
likratio <- res2$loglik-res1$loglik #32.11455
p.value <- 1-pchisq(likratio,df=res2$nNonzero-res1$nNonzero) #1.551289e-05
p.value
coefs <- res2$coefs
coefs
#	overall 1-6a to 6b		6b to 7		7 to 8a		8a-8b		8b-9		9-10
#-0.7238629 -0.3225696 -0.2183139 -0.1212226 -0.01382595	0.1926597 0.6259353
#conclusion: after grouping 1-6a, polysynthesis tends to decrease endangerment level, but the effect is weaker for higher endangerment level. 

#test if polysynthetic langauge has higher endangerment levels
autoord <- function (y,X,W) {
        Wy <- W%*%as.numeric(y)
        X2 <- W%*%X
        X2 <- cbind(X,X2)
        res <- lm(Wy~X2)$residuals
        X1 <- cbind(Wy,X,res)
        ordinalNet(X1,y,standardize=F,family="cumulative",link="probit",nonparallelTerms=F,lambdaVals=0)
}
W <- a[2]*(a[1]*Wsp+(1-a[1])*Wnb)+(1-a[2])*Wphy
poly.model <- autoord(y=y2,X=as.matrix(cbind(X,data[,intercept])),W=W)
null.model <- autoord(y=y2,X=as.matrix(data[,intercept]),W=W)

poly.model$coefs

# (Intercept):1 (Intercept):2 (Intercept):3 (Intercept):4 (Intercept):5 (Intercept):6                    X
# [1,]     0.6088611      1.230687      1.607527      1.894531      2.582852      2.971053 -0.2846774 0.1134448
# Western_Africa   Oceania      Europe Southern_Asia South_America Northern_America    Africa
# [1,]       0.791551 0.4343945 -0.08631648     0.4795371    -0.4335171       -0.9291352 0.5383821
# South-Eastern_Asia      Arab        Asia Central_America Australia_and_New_Zealand        res
# [1,]          0.1035645 -0.303795 -0.04552283      0.01344796                  -1.25181 -0.4234462
summary(poly.model)
summary(null.model)

likratio <- poly.model$loglik-null.model$loglik #-0.3641632
likratio

p.value <- 1-pchisq(likratio,df=poly.model$nNonzero-null.model$nNonzero) #1
p.value

#conclusion: polysynthetic langauges tend to have higher endangerment level (coefficient=0.1134448), but this is not statistically significant. 


#I think I need to add inteaction terms for each region
data$poly_Western_Africa <- data$poly*data$Western_Africa
data$poly_Oceania<- data$poly*data$Oceania
data$poly_Europe<- data$poly*data$Europe
data$poly_Southern_Asia<- data$poly*data$Southern_Asia
data$poly_South_America<- data$poly*data$South_America
data$poly_Northern_America<- data$poly*data$Northern_America
data$poly_Africa<- data$poly*data$Africa
data$poly_South_Eastern_Asia<- data$poly*data$`South-Eastern_Asia`
data$poly_Arab<- data$poly*data$Arab
data$poly_Asia<- data$poly*data$Asia
data$poly_Central_America <- data$poly*data$Central_America
data$poly_Australia_and_New_Zealand <- data$poly*data$Australia_and_New_Zealand


ncol(data)


#test if polysynthesis has explanatory power to endangerment levels beyond the best model for endangerment levels
autoord.lasso <- function (y,X,W) {
        Wy <- W%*%as.numeric(y)
        X2 <- W%*%X
        X2 <- cbind(X,X2)
        res <- lm(Wy~X2)$residuals
        X1 <- cbind(Wy,X,res)
        ordinalNetCV(X1,y,nFolds=10,alpha=0.5,standardize=F,family="cumulative",link="probit",nonparallelTerms=F,tuneMethod="cvLoglik",maxiterOut=1000)
}
#X <- data[,c(model,data[,'poly'])] #model is the best model we had for language endangerment

X <- data[,c(model,colnames(data)[2808:2820])] #model is the best model we had for language endangerment
X <- as.matrix(X)
W <- a[2]*(a[1]*Wsp+(1-a[1])*Wnb)+(1-a[2])*Wphy
tmp <- which(!is.na(rowSums(X)))
y.tmp <- y2[tmp]
X.tmp <- X[tmp,]
W <- W[tmp,tmp]
best.poly.model <- autoord.lasso(y=y.tmp,X=X.tmp,W=W)
best.poly.model$fit$coefs[,"poly"]
best.poly.model$fit$coefs[20,]

#0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0   #zero estimated effect of polysynthesis to any region.  
#conclusion: polysynthetic langauges do not have explanatory power to higher endangerment levels beyond the best model for endangerment levels
#saveRDS(best.poly.model, 'best.poly.model.rds')
#saveRDS(best.poly.model, 'best.poly.modelwRegionalInteractions.rds')
saveRDS(best.poly.model, 'best.polybord.modelwRegionalInteractions.rds')

best.poly.model$fit$coefs[20,]
# best.poly.model <- readRDS('best.polybord.model.rds')
# write.csv(best.poly.model$fit$coefs, 'endangerment_w_poly__coefficients.csv')


write.csv(best.poly.model$fit$coefs, 'endangerment_w_poly_borderline_coefficients.csv')
#write.csv(best.model$fit$coefs, 'endangerment_no_poly.csv')


summary(best.poly.model)
best.poly.model$fit$loglik[20]


test_stat <- 2 * (best.poly.model$fit$loglik[20]
 - best.model$fit$loglik[20])
p_value <- pchisq(test_stat, 1, lower.tail = FALSE)
test_stat
p_value

test_stat <- 2 * (best.polybord.model$fit$loglik[20]
                  - best.model$fit$loglik[20])
p_value <- pchisq(test_stat, 1, lower.tail = FALSE)
p_value
