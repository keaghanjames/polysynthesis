#The code defines the function autoord() that takes six arguments: a, y, X, Wsp, Wphy, and Wnb. 
#It also implements the polysynthesis analysis using the Polysynthetic and Extended list

#   - a is a numeric vector with length 2 that is used to compute the weights for the variables in the model.
#   - y is a numeric vector of the response variable.
#   - X is a numeric matrix of predictors.
#   - Wsp, Wphy, and Wnb are numeric matrices that represent the spacial, phylogenetic and contact distances between languages.

#The optim() function is then used to optimize the parameters a for the autoord() function. 
#The optim() function takes the initial values for a as well as the autoord() function itself and other optimization parameters as inputs, and returns the optimal values of a. These optimal values are then used as inputs to the autoord() function.


setwd('') #you will need to set the path
#lets load in those covariance matrices 
load('Wphy.Rdata')
load('Wsp.Rdata')
load('Wnb.Rdata')


#now we read in our data
#data$islands11km
data <- read.csv('~/Dropbox/Lindell/Polysynthesis/languoid_data_for_analysis.csv')
head(data)
data$ln_L1_pop <- log(data$L1_pop+1)
variables <- c("pol_borderline", "area", "bordering_language_richness", 'L1_pop',
               "bordering_language_richness_perkm", "bordering_language_evenness",
               "roughness", "altitude_range", 
               'islands11km', "small_family", "isolate_or_mono", 
               "vehicularity")

data <- data[, variables]
data
#data$family_id <- as.factor(data$family_id)
#data$family_id <- as.integer(data$family_id)

# data$documentation <- factor(data$documentation, ordered = T, levels = c('little or none', 'basic', 'detailed'))
# data$documentation <- as.integer(data$documentation)
# data
as.matrix(data)


#we divide our data into the training and test data 
N <- nrow(data)
idx <- sample.int(N,size=floor(2*N/3))
idx <- 1:N #for all observations

idx #we take a random sample of our observations

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

#v <- variable[1:n.var, 1] #we get the list of variables
#X <- data1[,v[which(v %in% colnames(X) == T)]] #and we select them as above but drop the variables that aren't in data1
X <- data1[,variables]
head(X)
ncol(X)
X <- as.matrix(X)
X
y <- data1$pol_borderline
y
y1 <- as.factor(y)
#levels(y1)[1:6] <- "6a"
y1

#we need to specify the formula of our model we will be optimising over, and we include two new terms
#weighted_sums and res which are added to the data by the autoglm function
fml <- paste0('pol_borderline ~ 1 + weighted_sums + res')
fml <- formula(fml)
fml

# fml2 <- paste0('poly_boderline ~ ', paste(colnames(X)[2:ncol(X)], collapse = '+'), '+ weighted_sums + res')
# fml2 <- formula(fml2)
# 
# fml2

NAs <- apply(X, 1, function(x) {any(is.na(x))})
tmp <- which(NAs == F)
tmp
data$islands11km

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
  #out <- ordinalNet(X1, y.tmp, standardize=F, family="cumulative", link="probit", nonparallelTerms=F, lambdaVals=0)  
  out2 <- glm(formula = fml, data = as.data.frame(X1), family=binomial(link="probit"))
  #summary(out2)
  
  #-out$loglik  # returns the negative log-likelihood of the model fit by ordinalNet()
  out2 <- -logLik(out2)[[1]]
}

# The code runs optim() to find the values of a that minimize the negative log-likelihood of the 
# glm model, and assigns these values to a. 

#we need to standerdise the continous predictors
X[,2:9] <- scale(X[,2:9])

#we specify an intercept only model that includes the weighted sums and residuals to get optimise over a
fml_int <- as.formula('pol_borderline ~ 1 + weighted_sums + res')

res <- try(optim(c(0.5, 0.5), autoglm, method="L-BFGS-B", lower=c(0, 0), upper=c(1, 1), y=as.numeric(y1)-1, X=X, Wsp=Wsp1, Wphy=Wphy1, Wnb=Wnb1, fml = fml_int))
a <- res$par
a

#these values of a will then be parsed to the rest of our analysis

#lets write a function that only takes fml as an input

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

#we can use the a = b argument in the function to fix the values of a, allowing for parralel analysis
b = a
poly_model_maker(fml = fml)



#mods <- list(fml, fml)

#we just need to make all possible combinations of the model

all_vars <- c(variables)
all_vars <- all_vars[-1]
all_vars
seq <- 1:length(all_vars)
combs <- lapply(seq, function(x) combn(all_vars, x))
#apply(combs[[12]], 2, function(x) as.list(as.data.frame(x)))

model_combos <- lapply(combs, function(x) as.list(as.data.frame(x)))
length(model_combos)
model_combos <- Reduce(c, model_combos)
names(model_combos) <- NULL
model_combos

for(i in 1:length(model_combos)) {
  model_combos[[i]] <- paste0('pol_borderline ~ ', paste0(c(model_combos[[i]], 'weighted_sums', 'res'), collapse = ' + ' ))
  model_combos[[i]] <- formula(model_combos[[i]])

}
model_combos


require(pbmcapply)
models <- pbmclapply(model_combos,
                     function(x) try(poly_model_maker(x)), mc.cores = detectCores()-1)


BICs <- unlist(models)
BICs
require(qpcR)

compute_bic_weights <- function(bic_values) {
  delta_bic <- bic_values - min(bic_values)
  weights <- exp(-0.5 * delta_bic) / sum(exp(-0.5 * delta_bic))
  return(list(delta_bic, weights))
}


head(compute_bic_weights(BICs))

aw <- akaike.weights(BICs)
aw
require(dplyr)
model_summaries <- tibble(BIC = BICs, delta = aw$deltaAIC, reLL = aw$rel.LL, weights = aw$weights)

model_matrix <- matrix(ncol = length(variables)-1, nrow = 0)
colnames(model_matrix) <- variables[-1]
model_matrix

for(i in 1:length(model_combos)){
  model_matrix <- rbind(model_matrix, colnames(model_matrix) %in% 
                          Reduce(c, strsplit(as.character(model_combos[[i]])[[3]], ' '))
  )
  print(i)
}
head(model_matrix)

row_sums <- rowSums(model_matrix)

model_matrix <- data.frame(apply(model_matrix, 2,function(x) ifelse(x, "X", " ")))
head(model_matrix)
head(model_summaries)
model_matrix <- as_tibble(model_matrix)
model_summaries <- bind_cols(model_matrix, model_summaries)
model_summaries$model_id <- 1:nrow(model_summaries)
model_summaries$Nparam <- row_sums

model_summaries <- arrange(model_summaries, BIC)
model_summaries$model_id
write.csv(model_summaries, 'model_summaries_extended.csv')


significant_predictors <- function(model) {
  # Get the summary of the model
  summary_model <- summary(model)
  
  # Get the p-values for each predictor
  p_values <- summary_model$coefficients[,4]
  
  # Find the indices of significant predictors
  significant_indices <- which(p_values < 0.05)
  
  # Get the names of significant predictors
  significant_predictors <- row.names(summary_model$coefficients)[significant_indices]
  
  # Return the names of significant predictors
  return(significant_predictors)
}

sig_pred <- list()

for(i in 1:length(BICs)) {
  m <- readRDS(paste0('models_extended_list/', gsub(" ", "", as.character(model_combos[[model_summaries$model_id[i]]])[3]), '.rds'))
  sig_pred[[i]] <- significant_predictors(m)
}

sig_pred

#model_summaries <- model_summaries[which(model_summaries$delta < 2),]
model_summaries



for(i in 1:nrow(model_summaries)) {
  model_summaries[i,which(colnames(model_summaries) %in% sig_pred[[i]])] <-  paste0('X*')
  print(i)
}

write.csv(model_summaries, 'best_fits_extended.csv')




require(tidyr)
require(tibble)
i <- 1

summary_by_model <- tibble(variable = vector(),  Estimate= vector(),
                           Std.Error = vector(), 'Pr(>|z|)' = vector(),
                           low = vector(), high = vector(),
                           sig = vector(), model_order = vector(), nPred = vector(), ID = vector())
for(i in 1:length(BICs)){
  m <- readRDS(paste0('models_extended_list/', gsub(" ", "",as.character(model_combos[[model_summaries$model_id[i]]])[3]), '.rds'))
  coef_tbl <- summary(m)$coefficients %>%
    as.data.frame() %>%
    rownames_to_column(var = "variable") %>%
    as_tibble()
  cofint <- confint(m)
  coef_tbl$low = cofint[,1]
  coef_tbl$high = cofint[,2]
  coef_tbl$sig <- F
  coef_tbl$sig[which(coef_tbl$`Pr(>|z|)` < 0.05)] <- T
  coef_tbl$model_order <- i
  coef_tbl$nPred <- nrow(coef_tbl)
  coef_tbl$ID <- i
  summary_by_model <- rbind(summary_by_model, coef_tbl)
  print(i)
}

summary_by_model$best_fit <- F
summary_by_model$best_fit[which(summary_by_model$model_order == 1)] <- T

summary_by_model$fewest_params <- F
summary_by_model$fewest_params[which(summary_by_model$nPred == min(summary_by_model$nPred))] <- T


summary_by_model$Estimate <- as.character(round(summary_by_model$Estimate, 3))
summary_by_model$low <- as.character(round(summary_by_model$low, 3))
summary_by_model$high <- as.character(round(summary_by_model$high, 3))



for(i in 1:nrow(summary_by_model)){
#  for(i in 1:69){
    
  if(summary_by_model$`Pr(>|z|)`[i] <= 0.05) {
    summary_by_model$Estimate[i] <- paste0(summary_by_model$Estimate[i], '*')
  }
  if(summary_by_model$`Pr(>|z|)`[i] <= 0.01) {
    summary_by_model$Estimate[i] <- paste0(summary_by_model$Estimate[i], '*')
  }
  if(summary_by_model$`Pr(>|z|)`[i] <= 0.001) {
    summary_by_model$Estimate[i] <- paste0(summary_by_model$Estimate[i], '*')
  }
  
  summary_by_model$Estimate[i] <- paste0(summary_by_model$Estimate[i], " (", 
                                         summary_by_model$low[i], "-", 
                                         summary_by_model$high[i], ')')
  
}


summary_by_model_wide <- pivot_wider(summary_by_model,id_cols = ID, names_from = variable, values_from = Estimate)
summary_by_model_wide$BIC <- model_summaries$BIC
summary_by_model_wide$delta <- model_summaries$delta
summary_by_model_wide$Nparam <- model_summaries$Nparam

summary_by_model_wide

write.csv(summary_by_model_wide, 'summary_by_model_extended_wide.csv', row.names = F)

summary_by_model_wide[which(summary_by_model_wide$delta < 6)]



require(ggplot2)
summary_by_model
summary_by_model[-which(summary_by_model$variable == 'weighted_sums' | summary_by_model$variable == 'res'),]

plot_df <- summary_by_model[-which(summary_by_model$variable == 'weighted_sums' | summary_by_model$variable == 'res'),]


gg1 <- ggplot(plot_df[plot_df$fewest_params == T,], aes(x = variable, y = Estimate, col = sig)) +
  geom_point() +
  geom_errorbar(aes(ymin = low, ymax = high, width = .5)) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_x_discrete(labels = c('Intercept', 'Bordering LE', 'Bordering LR', 'Island Endemic', 'Isolate or small family', 'LR')) +
  scale_color_manual(values = c('grey', 'black')) +
  ylim(c(-6.5, 1.5)) +
  ylab('Coefficient ± 95%CI') +
  xlab("")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none')
gg1


gg2 <- ggplot(plot_df[plot_df$model_order == 1,], aes(x = variable, y = Estimate, col = sig)) +
  geom_point() +
  geom_errorbar(aes(ymin = low, ymax = high, width = .5)) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_x_discrete(labels = c('Intercept', 'Altitude range' ,'Bordering LE', 'Bordering LR', 'Island Endemic', 
  'Isolate or small family', 'L1 population', 'LR', 'Roughness')) +
  scale_color_manual(values = c('grey', 'black')) +
  ylim(c(-6.5, 1.5)) +
  ylab('Coefficient ± 95%CI') +
  xlab("")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none')
gg2

require(ggpubr)
coef_plot <- ggarrange(gg1,gg2,ncol = 2, widths = c(.6,.9), labels = c('a', 'b'))
coef_plot

ggsave('coefficient_plot.pdf', device = 'pdf', units = 'in', height = 6, width = 10, dpi = 300)




gg_coef <- ggplot(summary_by_model, aes(x = variable, y = Estimate, col = sig, shape = best_fit, size = best_fit,))+
  geom_point() +
  geom_hline(yintercept = 0, lty = 2) +
  scale_shape_manual(values = c(5,4)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

gg_coef

table_df <- summary_by_model[,c('variable', 'Estimate', 'model_order', 'sig')]
table_df$Estimate <- as.character(round(table_df$Estimate, 3))
table_df$Estimate <- paste0(table_df$Estimate, ifelse(table_df$sig, '*', ''))
table_df <- spread(table_df[,-4], variable, Estimate) 
table_df <- cbind(table_df, round(model_summaries[,c('AIC', 'delta', 'weights', 'Nparam')],3))
table_df
write.csv(table_df, 'coefficients_best_models_iso_PS.csv')


# m1 <- readRDS(paste0('models/', as.character(model_combos[[model_summaries$model_id[1]]])[3], '.rds'))
# summary(m1)
# m2 <- readRDS(paste0('models/', as.character(model_combos[[model_summaries$model_id[2]]])[3], '.rds'))
# summary(m2)
# m3 <- readRDS(paste0('models/', as.character(model_combos[[model_summaries$model_id[3]]])[3], '.rds'))
# summary(m3)
# m4 <- readRDS(paste0('models/', as.character(model_combos[[model_summaries$model_id[4]]])[3], '.rds'))
# summary(m4)
# m5 <- readRDS(paste0('models/', as.character(model_combos[[model_summaries$model_id[5]]])[3], '.rds'))
# summary(m5)
# m6 <- readRDS(paste0('models/', as.character(model_combos[[model_summaries$model_id[6]]])[3], '.rds'))
# summary(m6)
# m9 <- readRDS(paste0('models/', as.character(model_combos[[model_summaries$model_id[9]]])[3], '.rds'))
# summary(m9)
# m7$coefficients
# 
# require(stargazer)
# stargazer(m9, type = 'html', ci = T, df = T, single.row = T,summary.logical = T,
#           star.cutoffs = c(0.05, 0.01, 0.001), out = 'm9.html')
# 
# 
# best_model <- readRDS(paste0('models/', as.character(model_combos[[which.min(AICs)]])[3], '.rds'))
# 
# 
# summary(best_model)
# require(car)
# vif(best_model)
# 
# 
# require(sjPlot)
# tab_model(best_model)
# require(stargazer)
# stargazer(best_model, type = 'html', ci = T, df = T, single.row = T,summary.logical = T,
#           star.cutoffs = c(0.05, 0.01, 0.001), out = 'best_model.html')
# require(texreg)
# screenreg(out2, single.row = T, ci.force = T)
# require(ggplot2)
# require(ggcorrplot)
# 
# #at the very least re run the model but 
# 
# 
# cor_mat <- cor(X, use = "complete.obs")
# p.mat <- cor_pmat(X)
# cor_plot <- ggcorrplot(cor_mat, method = 'circle', hc.order = T,
#            type = 'lower', 
#            ggtheme = ggplot2::theme_gray,
#            outline.col = 'white',
#            colors = c("#6D9EC1", "grey", "#E46726"),
#            p.mat = p.mat, pch = 7)
# cor_plot
# setwd("..")
# ggsave('correlogram.pdf', dpi = 300)
