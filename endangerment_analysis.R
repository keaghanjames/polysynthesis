library(ordinalNet)
library(ivreg)
library(nloptr)

setwd("") # you will need to set the appropriate path
load("bestmodel")
updated <- read.csv('language_data_for_analysis.csv')

#use the below lines to switch between the Polysynthetic (poly_definite) and Extended (pol_borderline) lists
#data$poly <- updated$poly_definite
data$poly <- updated$pol_borderline

data$poly <- as.numeric(data$poly)

#we add interaction terms for each of the global regions
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


y <- data$EGIDS_tr
levels(y2)[1:6] <- "6a"

#test if polysynthesis has explanatory power to endangerment levels beyond the best model for endangerment levels by refitting
#Bromham et al. 2022's best model with poly and its interaction terms as additional variables
#first we define the function
autoord.lasso <- function (y,X,W) {
        Wy <- W%*%as.numeric(y)
        X2 <- W%*%X
        X2 <- cbind(X,X2)
        res <- lm(Wy~X2)$residuals
        X1 <- cbind(Wy,X,res)
        ordinalNetCV(X1,y,nFolds=10,alpha=0.5,standardize=F,family="cumulative",link="probit",nonparallelTerms=F,tuneMethod="cvLoglik",maxiterOut=1000)
}

#we define our outcome variables
y <- data$EGIDS_tr
levels(y2)[1:6] <- "6a" #and we also collapse the first six levels into a single level to be consistent with Bromham et al. 2022. 

#we extract our predictors
X <- data[,c(model,colnames(data)[2808:2820])] #model is the best model for language endangerment from Bromham et al. 2022.
#this extracts those variables and the additional interaction variables from the data frame

X <- X[,which(colnames(X) != 'island')] #we drop the variable island as we aren't testing it here

X <- as.matrix(X) #X needs to be a matrix for ordinalNetCB

#here we are calculating the weighted sums of our joint phylo, spatial and contact matrix
#the relative contribution of each is determined by the parameter a which we have inherited
#from Bromham et al. 2022 analysis
W <- a[2]*(a[1]*Wsp+(1-a[1])*Wnb)+(1-a[2])*Wphy
tmp <- which(!is.na(rowSums(X))) #just making sure we don't include any rows containing NAs or oridinalNet will fail
y.tmp <- y2[tmp]
X.tmp <- X[tmp,]
W <- W[tmp,tmp]
best.poly.model <- autoord.lasso(y=y.tmp,X=X.tmp,W=W)
best.poly.model$fit$coefs[,"poly"]
best.poly.model$fit$coefs[20,]
#for lasso regression you only need to interpret the coefficients from the last iteration, i.e. row 20

#probably want to save the model object
saveRDS(best.poly.model, 'best.polybord.modelwRegionalInteractions.rds')


#carry out a loglikelihood test to compare the new model to the best model from Bromham et al. 2022
test_stat <- 2 * (best.poly.model$fit$loglik[20]
 - best.model$fit$loglik[20])
p_value <- pchisq(test_stat, 1, lower.tail = FALSE)
test_stat
p_value


