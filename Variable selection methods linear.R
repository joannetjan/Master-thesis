# MLE  ------------------------------------------------------------------------------------------------------------------------
set.seed(23)
n.train <- 0.8

# Setting up the dataset with intercept, including covariates 
D <- rep(0, N*(L+1)*length(Z))
ar <- array(D, c(length(Z), N, (L+1)))
ar1 <- array(D, c(length(Z), N, (L+1)))
final_yzt.train <- matrix(0, nrow = N, ncol = 5)
final_yzt.validate <- matrix(0, nrow = N, ncol = 5)

for (i in seq(5)) {
  tempW <- CWind[CWind[,i+1] != 99999, c(1,i+1)]
  df.all <- na.omit(tempW)
  
  tempdf <- merge(x = df.all, y = CV, by = "Date")
  df.train <- tempdf %>% sample_frac(n.train)
  df.train <- df.train[order(df.train[,1], decreasing = FALSE),]
  df.validate <- tempdf[!tempdf[,1] %in% df.train[,1],]
  
  tempCV <- as.matrix(df.train[,3:dim(tempdf)[2]])
  obs <- dim(tempCV)[1]
  
  tempCV1 <- as.matrix(df.validate[,3:dim(tempdf)[2]])
  obs1 <- dim(tempCV1)[1]
  
  final_yzt.train[1:obs,i] <- df.train[,2]
  final_yzt.validate[1:obs1, i] <- df.validate[,2]

  ar[i, 1:obs, 2:(L+1)] <- tempCV
  ar[i, 1:obs, 1] <- rep(1, obs)
  
  ar1[i, 1:obs1, 2:(L+1)] <- tempCV1
  ar1[i, 1:obs1, 1] <- rep(1, obs1)

}

df.wind.all <- as.vector(as.matrix(final_yzt.validate[final_yzt.validate != 0]))
df.wind.train <- as.vector(as.matrix(final_yzt.train[final_yzt.train != 0]))

likelihood_base <- function(parameters) {
  # parameters <- c(rep(0.1,(L*K)), rep(0.01, L*K))
  m <- parameters[1]
  std <- exp(parameters[2])
  
  x_zt <- (df.wind.train-m)/std
  term1 <- -(log(std)*length(x_zt))
  
  if (sum(is.nan(term1)) > 0 | sum(is.infinite(term1)) > 0 | term1 > 0) {
    term1 <- -1e100
  }
  
  term2 <- -sum(x_zt)
  if (sum(is.nan(term2)) > 0 | sum(is.infinite(term2)) > 0 | term2 > 0) {
    term2 <- -1e100
  }
  e_x_zt <- exp(-x_zt)
  
  if (sum(is.infinite(e_x_zt)) > 0) {
    term3 <- -1e100
  } else {
    term3 <- -sum(e_x_zt)
  }
  ll <- term1+term2+term3
  print(-ll)
  print(c(term1, term2, term3))
  
  return(-ll)
}

likelihood_standard <- function(parameters) {
  # parameters <- c(rep(0.1,(L*K)), rep(0.01, L*K))
  m <- parameters[1:(L+1)]
  std <- parameters[(L+2):length(parameters)]
  
  final_mu_zt <- matrix(0, nrow = N, ncol = 5)
  final_sd_zt <- matrix(0, nrow = N, ncol = 5)
  
  for (i in seq(5)) {
    final_mu_zt[,i] <- ar[i,,]%*%m
    final_sd_zt[,i] <- ar[i,,]%*%std
  }
  
  final_sd_zt <- exp(final_sd_zt)
  x_zt <- (final_yzt.train-final_mu_zt)/final_sd_zt
  term1 <- -sum(log(final_sd_zt))
  
  if (sum(is.nan(term1)) > 0 | sum(is.infinite(term1)) > 0 | term1 > 0) {
    term1 <- -1e100
  }
  
  term2 <- -sum(x_zt)
  if (sum(is.nan(term2)) > 0 | sum(is.infinite(term2)) > 0 | term2 >0) {
    term2 <- -1e100
  }
  e_x_zt <- exp(-x_zt)
  e_x_zt[e_x_zt == 1] <- 0
  
  if (sum(is.infinite(e_x_zt)) > 0) {
    term3 <- -1e100
  } else {
    term3 <- -sum(e_x_zt)
  }
  ll <- term1+term2+term3
  print(-ll)
  print(c(term1, term2, term3))
  if ((-ll )< 0 | sum(is.nan(ll)) > 0){
    print(parameters)
  }
  
  return(-ll)
}

# Variable selection method ------------------------------------------------------------------------------------------------------------------------
# LASSO
likelihood_lasso <- function(parameters) {
  # parameters <- c(rep(0.1,(L*K)), rep(0.01, L*K))
  m <- parameters[1:(L+1)]
  std <- parameters[(L+2):length(parameters)]
  
  final_mu_zt <- matrix(0, nrow = N, ncol = 5)
  final_sd_zt <- matrix(0, nrow = N, ncol = 5)
  
  for (i in seq(5)) {
    final_mu_zt[,i] <- ar[i,,]%*%m
    final_sd_zt[,i] <- ar[i,,]%*%std
  }
  
  final_sd_zt <- exp(final_sd_zt)
  x_zt <- (final_yzt.train-final_mu_zt)/final_sd_zt
  term1 <- -sum(log(final_sd_zt))
  
  if (sum(is.nan(term1)) > 0 | sum(is.infinite(term1)) > 0 | term1 > 0) {
    term1 <- -1e100
  }
  
  term2 <- -sum(x_zt)
  if (sum(is.nan(term2)) > 0 | sum(is.infinite(term2)) > 0 | term2 >0) {
    term2 <- -1e100
  }
  e_x_zt <- exp(-x_zt)
  e_x_zt[e_x_zt == 1] <- 0
  
  if (sum(is.infinite(e_x_zt)) > 0) {
    term3 <- -1e100
  } else {
    term3 <- -sum(e_x_zt)
  }
  
  lmbd <- first.lmbd*o
  term4 <- -lmbd*sum(abs(parameters[c(-1,-18)]))
  ll <- term1+term2+term3+term4
  print(-ll)
  print(c(term1, term2, term3, term4))
  if ((-ll )< 0 | sum(is.nan(ll)) > 0){
    print(parameters)
  }
  
  return(-ll)
}

# Ridge
likelihood_ridge <- function(parameters) {
  # parameters <- c(rep(0.1,(L*K)), rep(0.01, L*K))
  m <- parameters[1:(L+1)]
  std <- parameters[(L+2):length(parameters)]
  
  final_mu_zt <- matrix(0, nrow = N, ncol = 5)
  final_sd_zt <- matrix(0, nrow = N, ncol = 5)
  
  for (i in seq(5)) {
    final_mu_zt[,i] <- ar[i,,]%*%m
    final_sd_zt[,i] <- ar[i,,]%*%std
  }
  
  final_sd_zt <- exp(final_sd_zt)
  x_zt <- (final_yzt.train-final_mu_zt)/final_sd_zt
  term1 <- -sum(log(final_sd_zt))
  
  if (sum(is.nan(term1)) > 0 | sum(is.infinite(term1)) > 0 | term1 > 0) {
    term1 <- -1e100
  }
  
  term2 <- -sum(x_zt)
  if (sum(is.nan(term2)) > 0 | sum(is.infinite(term2)) > 0 | term2 >0) {
    term2 <- -1e100
  }
  e_x_zt <- exp(-x_zt)
  e_x_zt[e_x_zt == 1] <- 0
  
  if (sum(is.infinite(e_x_zt)) > 0) {
    term3 <- -1e100
  } else {
    term3 <- -sum(e_x_zt)
  }
  
  term4 <- -second.lmbd*o*sum((parameters[c(-1,-18)])^2)
  ll <- term1+term2+term3+term4
  print(-ll)
  print(c(term1, term2, term3, term4))
  if ((-ll )< 0 | sum(is.nan(ll)) > 0){
    print(parameters)
  }
  
  return(-ll)
}

#Elastic-net
likelihood_elastic.net<- function(parameters) {
  # parameters <- c(rep(0.1,(L*K)), rep(0.01, L*K))
  m <- parameters[1:(L+1)]
  std <- parameters[(L+2):length(parameters)]
  
  final_mu_zt <- matrix(0, nrow = N, ncol = 5)
  final_sd_zt <- matrix(0, nrow = N, ncol = 5)
  
  for (i in seq(5)) {
    final_mu_zt[,i] <- ar[i,,]%*%m
    final_sd_zt[,i] <- ar[i,,]%*%std
  }
  
  final_sd_zt <- exp(final_sd_zt)
  x_zt <- (final_yzt.train-final_mu_zt)/final_sd_zt
  term1 <- -sum(log(final_sd_zt))
  
  if (sum(is.nan(term1)) > 0 | sum(is.infinite(term1)) > 0 | term1 > 0) {
    term1 <- -1e100
  }
  
  term2 <- -sum(x_zt)
  if (sum(is.nan(term2)) > 0 | sum(is.infinite(term2)) > 0 | term2 >0) {
    term2 <- -1e100
  }
  e_x_zt <- exp(-x_zt)
  e_x_zt[e_x_zt == 1] <- 0
  
  if (sum(is.infinite(e_x_zt)) > 0) {
    term3 <- -1e100
  } else {
    term3 <- -sum(e_x_zt)
  }
  print(first.lmbd, second.lmbd)
  term4 <- -(first.lmbd*o*sum(abs(parameters[c(-1,-18)]))+ second.lmbd*o*sum((parameters[c(-1,-18)])^2))
  print(parameters)
  if (sum(is.nan(term4)) > 0 | sum(is.infinite(term4)) > 0 | term4 > 0) {
    term4 <- -1e100
  }
  
  ll <- term1+term2+term3+term4
  print(-ll)
  print(c(term1, term2, term3, term4))
  if ((-ll )< 0 | sum(is.nan(ll)) > 0){
    print(parameters)
  }
  
  return(-ll)
}

#Adaptive LASSO
likelihood_adaptive.LASSO <- function(parameters) {
  # parameters <- c(rep(0.1,(L*K)), rep(0.01, L*K))
  m <- parameters[1:(L+1)]
  std <- parameters[(L+2):length(parameters)]
  
  final_mu_zt <- matrix(0, nrow = N, ncol = 5)
  final_sd_zt <- matrix(0, nrow = N, ncol = 5)
  
  for (i in seq(5)) {
    final_mu_zt[,i] <- ar[i,,]%*%m
    final_sd_zt[,i] <- ar[i,,]%*%std
  }
  
  final_sd_zt <- exp(final_sd_zt)
  x_zt <- (final_yzt.train-final_mu_zt)/final_sd_zt
  term1 <- -sum(log(final_sd_zt))
  
  if (sum(is.nan(term1)) > 0 | sum(is.infinite(term1)) > 0 | term1 > 0) {
    term1 <- -1e100
  }
  
  term2 <- -sum(x_zt)
  if (sum(is.nan(term2)) > 0 | sum(is.infinite(term2)) > 0 | term2 >0) {
    term2 <- -1e100
  }
  e_x_zt <- exp(-x_zt)
  e_x_zt[e_x_zt == 1] <- 0
  
  if (sum(is.infinite(e_x_zt)) > 0) {
    term3 <- -1e100
  } else {
    term3 <- -sum(e_x_zt)
  }
  
  w <- (abs(ini.par[c(-1,-18)]))^(-third.gamma)
  term4 <- -first.lmbd*o*sum(w*abs(parameters[c(-1,-18)]))
  ll <- term1+term2+term3+term4
  
  print(-ll)
  print(c(term1, term2, term3, term4))
  
  if ((-ll )< 0 | sum(is.nan(ll)) > 0){
    print(parameters)
  }
  
  return(-ll)
}

#Adaptive Elastic-Net
likelihood_adaptive.EN <- function(parameters) {
  # parameters <- c(rep(0.1,(L*K)), rep(0.01, L*K))
  m <- parameters[1:(L+1)]
  std <- parameters[(L+2):length(parameters)]
  
  final_mu_zt <- matrix(0, nrow = N, ncol = 5)
  final_sd_zt <- matrix(0, nrow = N, ncol = 5)
  
  for (i in seq(5)) {
    final_mu_zt[,i] <- ar[i,,]%*%m
    final_sd_zt[,i] <- ar[i,,]%*%std
  }
  
  final_sd_zt <- exp(final_sd_zt)
  x_zt <- (final_yzt.train-final_mu_zt)/final_sd_zt
  term1 <- -sum(log(final_sd_zt))
  
  if (sum(is.nan(term1)) > 0 | sum(is.infinite(term1)) > 0 | term1 > 0) {
    term1 <- -1e100
  }
  
  term2 <- -sum(x_zt)
  if (sum(is.nan(term2)) > 0 | sum(is.infinite(term2)) > 0 | term2 >0) {
    term2 <- -1e100
  }
  e_x_zt <- exp(-x_zt)
  e_x_zt[e_x_zt == 1] <- 0
  
  if (sum(is.infinite(e_x_zt)) > 0) {
    term3 <- -1e100
  } else {
    term3 <- -sum(e_x_zt)
  }
  
  w <- (abs(ini.par.EN[c(-1,-18)]))^(-third.gamma)
  term4 <- -(first.lmbd*o*sum(w*abs(parameters[c(-1,-18)]))+ second.lmbd*o*sum((parameters[c(-1,-18)])^2))
  
  ll <- term1+term2+term3+term4
  
  print(-ll)
  print(c(term1, term2, term3, term4))
  
  if ((-ll ) < 0 | sum(is.nan(ll)) > 0){
    print(parameters)
  }
  
  return(-ll)
}


# first.lmbd <- 0.02
# second.lmbd <- 0.02
# third.gamma <- 1

o <- sum(final_yzt.train > 0)
parameters <- c(rep(0.1,(L+1)*2))
likelihood_standard(parameters)
likelihood_lasso(parameters)
likelihood_ridge(parameters)
likelihood_elastic.net(parameters)
likelihood_adaptive.LASSO(parameters)
likelihood_adaptive.EN(parameters)

start_time <- Sys.time()
optim.standard <- optim(par=c(rep(0.1,(L+1)*2)), fn = likelihood_standard, method = "BFGS", control = list(maxit = 1e6), hessian = TRUE)
end_time <- Sys.time()
end_time - start_time

# Cross-validation ------------------------------------------------------------------------------------------------------------------------
# one variable
begin.value <- 0.085
end.value <- 0.089
step <- 0.001

n <- seq(begin.value, end.value, by = step)
A <- rep(0, length(n)*((L+1)*2)*3)
store.values.CV1 <- array(D, c(length(n), (L+1)*2,3))

# 
likelihood <- likelihood_lasso
start_time <- Sys.time()
for (i in 1:length(n)) {
  first.lmbd <- second.lmbd <- n[i]
  # second.lmbd <- n[i]
  # third.gamma <- n[i]
  optim.CV <- optim(par=param.standard, fn = likelihood, method = "BFGS", control = list(maxit = 1e6), hessian = TRUE)
  store.values.CV1[i,,1] <- optim.CV$par
  store.values.CV1[i,,2] <- sqrt(diag(solve(optim.CV$hessian)))
  store.values.CV1[i,1,3] <- optim.CV$value
  store.values.CV1[i,2,3] <- n[i]
}
end_time <- Sys.time()
end_time - start_time

# store.values.LASSO <- store.values.CV1
store.value.Ridge <- store.values.cv1
x# store.value.EN

# CV Elastic net
begin.value1 <- 0.001
end.value1 <- 0.001
step1 <- 0.001

begin.value2 <- 0.106
end.value2 <- 0.115
step2 <- 0.001

n1 <- seq(begin.value1, end.value1, by = step1)
n2 <- seq(begin.value2, end.value2, by = step2)

# n1 <- c(0.15, 0.2, 0.005, 0.005)
# n2 <- c(0.06, 0.07, 0.104, 0.06)
# n1 <- c(0.005, 0.006, 0.007)
# n2 <- c(0.104, 0.105, 0.106)

A <- rep(0, length(n1)*length(n2)*((L+1)*2)*3)
store.values.CV2 <- array(D, c(length(n1)*length(n2), (L+1)*2,3))

likelihood <- likelihood_elastic.net
count <- 1
start_time <- Sys.time()
for (i in 1:length(n1)) {
  first.lmbd <-  n1[i]
  for (j in 1:length(n2)) {
    second.lmbd <- n2[j]
    optim.CV <- optim(par=param.standard, fn = likelihood, method = "BFGS", control = list(maxit = 1e6), hessian = TRUE)
    store.values.CV2[count,,1] <- optim.CV$par
    store.values.CV2[count,,2] <- sqrt(diag(solve(optim.CV$hessian)))
    store.values.CV2[count,1,3] <- optim.CV$value
    store.values.CV2[count,2,3] <- n1[i]
    store.values.CV2[count,3,3] <- n2[j]
    count <- count + 1
  }
}
end_time <- Sys.time()
end_time - start_time

store.values.EN <- store.values.CV2

# CV adaptive lasso
ini.par <- store.values.LASSO[which.min(CRPS.scores),,1]

begin.value <- 0.069
end.value <- 0.071
step <- 0.002

begin.value.gamma <- 0.087
end.value.gamma <- 0.089
step.gamma <- 0.002

n <- seq(begin.value, end.value, by = step)
n.gamma <- seq(begin.value.gamma, end.value.gamma, by = step.gamma)

A <- rep(0, length(n.gamma)*length(n)*((L+1)*2)*3)
store.values.CV3 <- array(D, c(length(n.gamma)*length(n), (L+1)*2,3))

likelihood <- likelihood_adaptive.LASSO
count <- 1
start_time <- Sys.time()
for (i in 1:length(n)) {
  first.lmbd <-  n[i]
  for (j in 1:length(n.gamma)) {
    tryCatch({
    third.gamma <- n.gamma[j]
    optim.CV <- optim(par=param.standard, fn = likelihood, method = "BFGS", control = list(maxit = 1e6), hessian = TRUE)
    store.values.CV3[count,,1] <- optim.CV$par
    store.values.CV3[count,,2] <- sqrt(diag(solve(optim.CV$hessian)))
    store.values.CV3[count,1,3] <- optim.CV$value
    store.values.CV3[count,2,3] <- n[i]
    store.values.CV3[count,3,3] <- n.gamma[j]
    count <- count + 1
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}
end_time <- Sys.time()
end_time - start_time

store.values.alasso <- store.values.CV3

# CV adaptive elastic net
ini.par.EN <- store.values.EN[which.min(CRPS.scores2),,1]

begin.value1 <- 0.022
end.value1 <- 0.024
step1 <- 0.001

begin.value2 <- 0.067
end.value2 <- 0.067
step2 <- 0.002

begin.value.gamma <- 0.013
end.value.gamma <- 0.013
step.gamma <- 0.001

n1 <- seq(begin.value1, end.value1, by = step1)
n2 <- seq(begin.value2, end.value2, by = step2)
n.gamma <- seq(begin.value.gamma, end.value.gamma, by = step.gamma)

A <- rep(0, length(n.gamma)*length(n1)*length(n1)*((L+1)*2)*3)
store.values.CV4 <- array(D, c(length(n.gamma)*length(n1)*length(n2), (L+1)*2,3))

likelihood <- likelihood_adaptive.EN
count <- 1
start_time <- Sys.time()
for (i in 1:length(n1)) {
  first.lmbd <-  n1[i]
  for (j in 1:lendgth(n2)) {
    second.lmbd <- n2[j]
    for (k in 1:length(n.gamma)) {q
      third.gamma <- n.gamma[k]
      optim.CV <- optim(par=param.standard, fn = likelihood, method = "BFGS", control = list(maxit = 1e6), hessian = TRUE)
      store.values.CV4[count,,1] <- optim.CV$par
      store.values.CV4[count,,2] <- sqrt(diag(solve(optim.CV$hessian)))
      store.values.CV4[count,1,3] <- optim.CV$value
      store.values.CV4[count,2,3] <- n1[i]
      store.values.CV4[count,3,3] <- n.gamma[k]
      store.values.CV4[count,4,3] <- n2[j]
      count <- count + 1
    }
  }

}
end_time <- Sys.time()
end_time - start_time

store.values.adaptive.EN <- store.values.CV4
