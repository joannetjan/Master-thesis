# Obtaining standard errors  ------------------------------------------------------------------------------------------------------------------------
S <- diag(diag(optim.LL$hessian))
S <- solve(S) # singular
diag(-solve(-optim.LL$hessian)) #standard errors

sqrt(diag(solve(-h2)))
sqrt(diag(solve(optim.LLS$hessian)))
optim.LLS$value

# Performance measures ------------------------------------------------------------------------------------------------------------------------
store.values.CV <- store.values.CV3
aic.bic.values <- matrix(0, 2, (dim(store.values.CV)[1]))

for (i in 1:dim(store.values.CV)[1]) {
  p <- store.values.CV[i,,1]
  se <- store.values.CV[i,,2]
  ll.value <- -store.values.CV[i,1,3]
  t <- abs(p / se)
  significant.parameters <- t > 1.96
  nsp <- sum(significant.parameters)
  AIC.value <- 2*nsp-2*ll.value
  BIC.value <- nsp*log(o)-2*ll.value
  aic.bic.values[,i] <- c(AIC.value, BIC.value)
  }
aic.bic.values

# CRPS scores, one variable
crps_gev.covariates <- function(y, shape, location, scale) {
  # y <- final_yzt.validate
  # shape <- 0
  # location <- store.values.CV1[i,1:17,1]
  # scale <- store.values.CV1[i,18:34,1]
  # 
  final.location <- matrix(0, nrow = N, ncol = 5)
  final.scale <- matrix(0, nrow = N, ncol = 5)
  
  for (i in seq(5)) {
    final.location[,i] <- ar1[i,,]%*%location
    final.scale[,i] <- ar1[i,,]%*%scale
  }
  
  final.scale <- exp(final.scale)
  
  final.location <- c(final.location)
  final.scale <- c(final.scale)
  y <- c(y)
  
  if (!identical(final.location, 0)) y <- y - final.location
  if (!identical(final.scale, 1)) {
    final.scale[final.scale < 0] <- NaN
    y <- y / final.scale
  } 
  
  if (any(ind <- abs(shape) < 1e-12, na.rm = TRUE)) {
    if (length(y) < length(shape)) y <- rep_len(y, length(shape))
    out <- rep_len(NaN, length(y))
    
    out[ind] <- -y[ind] - digamma(1) - log(2) - 2 *
      if (requireNamespace("gsl", quietly = TRUE)) {
        gsl::expint_Ei(-exp(-y[ind]))
      } else {
        warning(paste("The exponential integral is approximated using the 'integrate' function.",
                      "Consider installing the 'gsl' package to leverage a more accurate implementation.",
                      sep = "\n"))
        sapply(-exp(-y[ind]), function(upper) {
          integrate(function(x) exp(x)/x, -Inf, upper)$value
        })
      }
    out[!ind] <- crps_gev(y[!ind], shape[!ind])
  } else {
    x <- 1 + shape * y
    x[x < 0] <- 0
    x <- x^(-1/shape)
    c1 <- 2 * exp(-x) - 1
    out <- (y + 1/shape) * c1 + gamma(1 - shape) / shape *
      (2 * pgamma(x, 1 - shape) - 2^shape)
  }
  
  out <- matrix(out, nrow = N, ncol = 5)
  final.scale <- matrix(final.scale, nrow = N, ncol = 5 )
  return(final.scale * out)
}

ind.matrix <- final_yzt.validate
ind.matrix[ind.matrix > 0] <- 1

CRPS.base <- crps_gev(final_yzt.validate, 0, optim.base$par[1], exp(optim.base$par[2]))
CRPS.base <- matrix(CRPS.base, ncol = 5, nrow = 96432)
CRPS.standard <- crps_gev.covariates(final_yzt.validate, 0, optim.standard$par[1:17], optim.standard$par[18:34])
test <- ((CRPS.standard-CRPS.base)/CRPS.base)*ind.matrix
mean(test[test != 0])

D <- rep(0, N*length(n)*length(Z))
CRPS.values.CV1 <- array(D, c(length(n), N, length(Z)))
CRPS.scores <- rep(0, length(n))
CRPS.scores.standard <- rep(0, length(n))
# CRPS.standard <- crps_gev.covariates(final_yzt.validate, 0, param.standard[1:17], param.standard[18:34])



for (i in 1:length(n)) {
  CRPS.temp <- crps_gev.covariates(final_yzt.validate, 0, store.values.CV1[i,1:17,1], store.values.CV1[i,18:34,1])
  CRPS.values.CV1[i,,] <- CRPS.temp*ind.matrix
  
  temp.scores <- ((CRPS.temp-CRPS.base)/CRPS.base)*ind.matrix
  CRPS.scores[i] <- mean(temp.scores[temp.scores != 0])
  
  temp.scores.standard <- ((CRPS.temp-CRPS.standard)/CRPS.standard)*ind.matrix
  CRPS.scores.standard[i] <- mean(temp.scores.standard[temp.scores.standard != 0])
}

CRPS.scores
CRPS.scores.standard
which.min(CRPS.scores)
c(min(CRPS.scores), min(CRPS.scores.standard))
rbind(CRPS.scores, CRPS.scores.standard, n)
# mean((crps_gev.covariates(final_yzt, 0, store.values.CV1[5,1:17,1], store.values.CV1[5,18:34,1])-CRPS.base)/CRPS.base)
# mean((crps_gev.covariates(final_yzt, 0, optim.ll$par[1:17], optim.ll$par[18:34])-CRPS.base)/CRPS.base)

# CRPS scores, two variable, elastic net
D2 <- rep(0, length(n1)*N*length(n2)*length(Z))
CRPS.values.CV2 <- array(D, c(length(n1)*length(n2), N, length(Z)))
CRPS.scores2 <- rep(0, length(n1)*length(n2))
CRPS.scores2.standard <- CRPS.scores2

for (i in 1:(length(n1)*length(n2))) {
  CRPS.temp <- crps_gev.covariates(final_yzt.validate, 0, store.values.CV2[i,1:17,1], store.values.CV2[i,18:34,1])
  CRPS.values.CV2[i,,] <- CRPS.temp*ind.matrix
  
  temp.scores <- ((CRPS.temp-CRPS.base)/CRPS.base)*ind.matrix
  CRPS.scores2[i] <- mean(temp.scores[temp.scores != 0])
  
  temp.scores.standard <- ((CRPS.temp-CRPS.standard)/CRPS.standard)*ind.matrix
  CRPS.scores2.standard[i] <- mean(temp.scores.standard[temp.scores.standard != 0])
}

CRPS.scores2
CRPS.scores2.standard
c(min(CRPS.scores2), min(CRPS.scores2.standard))
rbind(CRPS.scores2, CRPS.scores2.standard, rep(n1, each = length(n2)), rep(n2, length(n1)))

# CRPS scores, two variable, adaptive lasso
D3 <- rep(0, length(n)*N*length(n.gamma)*length(Z))
CRPS.values.CV3 <- array(D, c(length(n)*length(n.gamma), N, length(Z)))
CRPS.scores3 <- rep(0, length(n)*length(n.gamma))
CRPS.scores3.standard <- CRPS.scores3

for (i in 1:(length(n)*length(n.gamma))) {
  CRPS.temp <- crps_gev.covariates(final_yzt.validate, 0, store.values.CV3[i,1:17,1], store.values.CV3[i,18:34,1])
  CRPS.values.CV3[i,,] <- CRPS.temp*ind.matrix
  
  temp.scores <- ((CRPS.temp-CRPS.base)/CRPS.base)*ind.matrix
  CRPS.scores3[i] <- mean(temp.scores[temp.scores != 0])
  
  temp.scores.standard <- ((CRPS.temp-CRPS.standard)/CRPS.standard)*ind.matrix
  CRPS.scores3.standard[i] <- mean(temp.scores.standard[temp.scores.standard != 0])
}

CRPS.scores3
CRPS.scores3.standard
c(min(CRPS.scores3), min(CRPS.scores3.standard))
which.min(CRPS.scores3)
rbind(CRPS.scores3, CRPS.scores3.standard, rep(n, each = length(n.gamma)), rep(n.gamma, length(n)))

# CRPS scores, three variable, adaptive elastic net
D4 <- rep(0, length(n1)*length(n2)*length(n.gamma)*length(Z)*N)
CRPS.values.CV4 <- array(D, c(length(n1)*length(n.gamma)*length(n2), N, length(Z)))
CRPS.scores4 <- rep(0, length(n1)*length(n.gamma)*length(n2))
CRPS.scores4.standard <- CRPS.scores4

for (i in 1:(length(n1)*length(n2)*length(n.gamma))) {
  CRPS.temp <- crps_gev.covariates(final_yzt.validate, 0, store.values.CV4[i,1:17,1], store.values.CV4[i,18:34,1])
  CRPS.values.CV4[i,,] <- CRPS.temp*ind.matrix
  
  temp.scores <- ((CRPS.temp-CRPS.base)/CRPS.base)*ind.matrix
  CRPS.scores4[i] <- mean(temp.scores[temp.scores != 0])
  
  temp.scores.standard <- ((CRPS.temp-CRPS.standard)/CRPS.standard)*ind.matrix
  CRPS.scores4.standard[i] <- mean(temp.scores.standard[temp.scores.standard != 0])
}

CRPS.scores4
CRPS.scores4.standard
which.min(CRPS.scores4)
c(min(CRPS.scores4), min(CRPS.scores4.standard))

n.1 <- rep(n1, each = length(n2)*length(n.gamma))
n.2 <- rep(rep(n2, each = length(n.gamma)))
n.3 <- rep(n.gamma, length(n1)*length(n2))
rbind(CRPS.scores4,CRPS.scores4.standard,  n.1 ,n.2 , n.3)

# Strongest wind -----------------------------------------------------------------------------------------------------------------------
crps_gev.covariates.change <- function(y, shape, location, scale, ar.change) {
  # y <- final_yzt.validate
  # shape <- 0
  # location <- store.values.CV1[i,1:17,1]
  # scale <- store.values.CV1[i,18:34,1]
  # 
  final.location <- ar.change%*%location
  final.scale <- ar.change%*%scale
  
  final.scale <- exp(final.scale)
  
  final.location <- c(final.location)
  final.scale <- c(final.scale)
  y <- c(y)
  
  if (!identical(final.location, 0)) y <- y - final.location
  if (!identical(final.scale, 1)) {
    final.scale[final.scale < 0] <- NaN
    y <- y / final.scale
  } 
  
  if (any(ind <- abs(shape) < 1e-12, na.rm = TRUE)) {
    if (length(y) < length(shape)) y <- rep_len(y, length(shape))
    out <- rep_len(NaN, length(y))
    
    out[ind] <- -y[ind] - digamma(1) - log(2) - 2 *
      if (requireNamespace("gsl", quietly = TRUE)) {
        gsl::expint_Ei(-exp(-y[ind]))
      } else {
        warning(paste("The exponential integral is approximated using the 'integrate' function.",
                      "Consider installing the 'gsl' package to leverage a more accurate implementation.",
                      sep = "\n"))
        sapply(-exp(-y[ind]), function(upper) {
          integrate(function(x) exp(x)/x, -Inf, upper)$value
        })
      }
    out[!ind] <- crps_gev(y[!ind], shape[!ind])
  } else {
    x <- 1 + shape * y
    x[x < 0] <- 0
    x <- x^(-1/shape)
    c1 <- 2 * exp(-x) - 1
    out <- (y + 1/shape) * c1 + gamma(1 - shape) / shape *
      (2 * pgamma(x, 1 - shape) - 2^shape)
  }
  return(final.scale * out)
}

percent <- 0.05
CRPS.all <- rep(0, length(Z))
CRPS.strong <- rep(0, length(Z))
CRPS.normal <- rep(0, length(Z))
optimal.param.method <- store.values.CV3[4,,]
# optimal.param.method <- optim.standard$par

for (i in 1:length(Z)) {
  swind <- final_yzt.validate[,i]
  sample.df.validate <- length(swind[swind!= 0])
  
  ar.temp <- ar[i,1:sample.df.validate,]
  
  CRPS.base.all <- crps_gev(swind[swind !=0], 0, optim.base$par[1], exp(optim.base$par[2]))
  CRPS.method.all <- crps_gev.covariates.change(swind[swind !=0], 0, optimal.param.method[1:17], optimal.param.method[18:34], ar.temp)
  CRPS.temp <- (CRPS.method.all-CRPS.base.all)/CRPS.base.all
  CRPS.all[i] <- mean(CRPS.temp)
  
  ordered.swind.index <- order(swind, decreasing = TRUE)
  swind <- swind[ordered.swind.index]
  strong <- swind[1:(floor(sample.df.validate*percent))]
  normal <- swind[(length(strong)+1):sample.df.validate]
  
  ar.temp <- ar1[i,ordered.swind.index[1:length(strong)],]
  
  CRPS.base.strong <- crps_gev(strong, 0, optim.base$par[1], exp(optim.base$par[2]))
  CRPS.method.strong <- crps_gev.covariates.change(strong, 0, optimal.param.method[1:17], optimal.param.method[18:34], ar.temp)
  CRPS.temp <- (CRPS.method.strong-CRPS.base.strong)/CRPS.base.strong
  CRPS.strong[i] <- mean(CRPS.temp)
  
  ar.temp <- ar1[i, (length(strong)+1):sample.df.validate, ]
  
  CRPS.base.normal <- crps_gev(normal, 0, optim.base$par[1], exp(optim.base$par[2]))
  CRPS.method.normal <- crps_gev.covariates.change(normal, 0, optimal.param.method[1:17], optimal.param.method[18:34], ar.temp)
  CRPS.temp <- (CRPS.method.normal-CRPS.base.normal)/CRPS.base.normal
  CRPS.normal[i] <- mean(CRPS.temp)
  
}

round(cbind(CRPS.all, CRPS.strong, CRPS.normal), digits = 4)
