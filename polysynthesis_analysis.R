#The code defines the function autoglm() that takes six arguments: a, y, X, Wsp, Wphy, and Wnb. 
#It also implements the polysynthesis analysis using the Polysynthetic and Extended list

#   - a is a numeric vector with length 2 that is used to compute the weights for the variables in the model.
#   - y is a numeric vector of the response variable.
#   - X is a numeric matrix of predictors.
#   - Wsp, Wphy, and Wnb are numeric matrices that represent the spacial, phylogenetic and contact distances between languages.

#The optim() function is then used to optimize the parameters a for the autoord() function. 
#The optim() function takes the initial values for a as well as the autoglm() function itself and other optimization parameters as inputs, and returns the optimal values of a. These optimal values are then used as inputs to the autoord() function.


setwd('') #you will need to set the path
#lets load in those distance matrices 
load('Wphy.Rdata')
load('Wsp.Rdata')
load('Wnb.Rdata')


#now we read in our data
data <- read.csv('~/Dropbox/Lindell/Polysynthesis/languoid_data_for_analysis.csv')
head(data)

#define a list of our variables of interest - currently its set to run the Extended List analysis ("pol_borderline")
#if you want to run the Polysynthetic List change this to "poly_definite". Each of these variables is described in 
#Table 1 of the manuscript
variables <- c("pol_borderline", "area", "bordering_language_richness", 'L1_pop',
               "bordering_language_richness_perkm", "bordering_language_evenness",
               "roughness", "altitude_range", 
               'islands11km', "small_family", "isolate_or_mono", 
               "vehicularity")

#subset our data frame by these variables 
data <- data[, variables]
data


#in the original Bromham et al. 2022 paper on which this analysis is based, they divided their observations into a training 
#and test set to carry out model validation tests. We didn't carry out those tests for this analysis but we retained the 
#functionality to split the dataset in the event that someone wanted to do it in the future

#the code below is set up such that all observations are included in the training set (1) and no observations are included 
#in the test set (2) and then all subsequent analyses are only carried out on the test set.

#the reason why we retain this is because dividing the covariance matrices is not straight forward and requires a special 
#function from Bromham et al. 2022. Defined below. 

N <- nrow(data)
#idx <- sample.int(N,size=floor(2*N/3))
idx <- 1:N #for all observations


#idx <- 1:N
#and we split the data using those samples
data1 <- data[idx,]
data2 <- data[-idx,]
#and we need to split our covraince matrices as well

divide_co_matrix <- function(m, idx, include = T, N) {
  if(include == T) {
    M1 <- m[idx,idx]
    diag(M1) <- 0
    tmp <- rowSums(M1)
    M1 <- t(M1)/rowSums(M1)
    M1[tmp==0,] <- 1/(length(idx)-1)
    diag(M1) <- 0
  } 
  if(include == F){
    M1 <- m[-idx,-idx]
    diag(M1) <- 0
    tmp <- rowSums(M1)
    M1 <- t(M1)/rowSums(M1)
    M1[tmp==0,] <- 1/(N-length(idx)-1)
    diag(M1) <- 0
  }
  return(M1)
}

Wphy1 <- divide_co_matrix(m = Wphy, idx = idx, N = N)
Wphy2 <- divide_co_matrix(m = Wphy, idx = idx, N = N, include = F)
Wsp1 <- divide_co_matrix(m = Wsp, idx = idx, N = N)
Wsp2 <- divide_co_matrix(m = Wsp, idx = idx, N = N, include = F)
Wnb1 <- divide_co_matrix(m = Wnb, idx = idx, N = N)
Wnb2 <- divide_co_matrix(m = Wnb, idx = idx, N = N, include = F)


#next we need to decide which variables we want to include in the model 
X <- data1[,variables]
X <- as.matrix(X)
X #contains all the predictors
y <- data1$pol_borderline #while y includes the outcome
#again to run just the polysynthetic list alter y to 
#y <- data1$poly_definite
y
y1 <- as.factor(y)

#we need to specify the formula of our model we will be optimising over, and we include two new terms
#weighted_sums and res which are added to the data by the autoglm function
fml <- paste0('pol_borderline ~ 1 + weighted_sums + res')
fml <- formula(fml)
fml

NAs <- apply(X, 1, function(x) {any(is.na(x))})
tmp <- which(NAs == F) #create a vector for all the rows that are free of NAs

#and here we define the autoglm function
autoglm <- function (a, y, X, Wsp, Wphy, Wnb, fml) {
  tmp <- which(!is.na(rowSums(X)))  # identifies rows in X that do not have missing values
  X.tmp <- X[tmp,-which(colnames(X) == as.character(fml)[2])]  # subsets X to exclude rows with missing values
  y.tmp <- y[tmp]  # subsets y to exclude rows with missing values
  W <- a[2]*(a[1]*Wsp + (1-a[1])*Wnb) + (1-a[2])*Wphy  # calculates weights using values in a, Wsp, Wnb, and Wphy
  W <- W[tmp, tmp]  # subsets W to exclude rows and columns with missing values
  Wy <- W %*% as.numeric(y.tmp)  # calculates a weighted sum of y.tmp using W
  X2 <- W %*% X.tmp  # calculates a weighted sum of X.tmp using W
  X2 <- cbind(X.tmp, X2)  # concatenates X.tmp and X2 as columns in a new matrix
  res <- lm(Wy~X2)$residuals  # fits a linear model of Wy on X2 and extracts the residuals
  X1 <- cbind(y.tmp,Wy, X.tmp, res)  # concatenates Wy, X.tmp, and the residuals as columns in a new matrix
  colnames(X1)[1:2] <- c(as.character(fml)[2], 'weighted_sums')
  out2 <- glm(formula = fml, data = as.data.frame(X1), family=binomial(link="probit"))
  #summary(out2) 
  out2 <- -logLik(out2)[[1]]
}

# The code runs optim() to find the values of a that minimize the negative log-likelihood of the 
# glm model, and assigns these values to a. 

#we need to standardise and scale the continuous predictors
X[,2:9] <- scale(X[,2:9])

#we specify an intercept-only model that includes the weighted sums and residuals to get optimise over a
fml_int <- as.formula('pol_borderline ~ 1 + weighted_sums + res')

res <- try(optim(c(0.5, 0.5), autoglm, method="L-BFGS-B", lower=c(0, 0), upper=c(1, 1), y=as.numeric(y1)-1, X=X, Wsp=Wsp1, Wphy=Wphy1, Wnb=Wnb1, fml = fml_int))

#we run the model optimising over a, a gives us the relative contribution of each of the covariance matrices
a <- res$par
a

#these values of a will then be parsed to the rest of our analysis

#lets write a function that only takes fml as an input and returns the BIC

poly_model_maker <- function(fml = fml, a = b) {
  tmp <- which(!is.na(rowSums(X)))  # identifies rows in X that do not have missing values
  X.tmp <- X[tmp,-which(colnames(X) == as.character(fml)[2])]  # subsets X to exclude rows with missing values
  y.tmp <- y[tmp]  # subsets y to exclude rows with missing values
  W <- a[2]*(a[1]*Wsp + (1-a[1])*Wnb) + (1-a[2])*Wphy  # calculates weights using values in a, Wsp, Wnb, and Wphy
  W <- W[tmp, tmp]  # subsets W to exclude rows and columns with missing values
  Wy <- W %*% as.numeric(y.tmp)  # calculates a weighted sum of y.tmp using W
  X2 <- W %*% X.tmp  # calculates a weighted sum of X.tmp using W
  X2 <- cbind(X.tmp, X2)  # concatenates X.tmp and X2 as columns in a new matrix
  res <- lm(Wy~X2)$residuals  # fits a linear model of Wy on X2 and extracts the residuals
  X1 <- cbind(y.tmp,Wy, X.tmp, res)  # concatenates Wy, X.tmp, and the residuals as columns in a new matrix
  colnames(X1)[1:2] <- c(as.character(fml)[2], 'weighted_sums')
  #out <- ordinalNet(X1, y.tmp, standardize=F, family="cumulative", link="probit", nonparallelTerms=F, lambdaVals=0)  
  out2 <- glm(formula = fml, data = as.data.frame(X1), family=binomial(link="probit"))
  saveRDS(out2, paste0('models_extended_list/', gsub(" ", "", as.character(fml)[3]), '.rds'))
  return(BIC(out2))
}

#critically, this function also saves the model object as an RDS file so make sure you have the file structure set up for it

#we can use the a = b argument in the function to fix the values of a, allowing for parallel analysis
b = a
poly_model_maker(fml = fml)


#we just need to make all possible combinations of the model

all_vars <- c(variables)
all_vars <- all_vars[-1]
all_vars
seq <- 1:length(all_vars)
combs <- lapply(seq, function(x) combn(all_vars, x))

model_combos <- lapply(combs, function(x) as.list(as.data.frame(x)))
length(model_combos)
model_combos <- Reduce(c, model_combos)
names(model_combos) <- NULL


for(i in 1:length(model_combos)) {
  model_combos[[i]] <- paste0('pol_borderline ~ ', paste0(c(model_combos[[i]], 'weighted_sums', 'res'), collapse = ' + ' ))
  model_combos[[i]] <- formula(model_combos[[i]])

}
length(model_combos)


#okay, now we can actually run the analysis 
#its quite slow so I use pbmcapply package to run it in parallel 
#if you're on a windows computer you made need to try something else
require(pbmcapply)
models <- pbmclapply(model_combos,
                     function(x) try(poly_model_maker(x)), mc.cores = detectCores()-1)


#we unlist the BICS
BICs <- unlist(models)
BICs

#then we find the order of our models from best fit (lowest BIC) to worst (highest BIC)
model_order <- order(BICs)
model_order
BICs <- BICs[model_order]

#rearrange the model combos list and inspect the best ones
model_combos <- model_combos[model_order]
head(model_combos)

#and then we can read in any model using model_combos you just have to change the number in the [[]]
#as this is calling the ith entry of model_combos, in this case we are calling the first or best fit model
model<- readRDS(paste0('models_extended_list/', gsub(" ", "", as.character(model_combos[[1]][3]), '.rds'))
                     
summary(model) #to inspect the coefficients

# cor_plot
# setwd("..")
# ggsave('correlogram.pdf', dpi = 300)
