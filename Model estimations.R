## A single estimation is done for each variable selection method as well as for the baseline model and the standard model

# MLE  ------------------------------------------------------------------------------------------------------------------------
first.lmbd <- 0.02
second.lmbd <- 0.02
final.third.gamma <- 1
o <- 53865 + 56682 + 54046 + 53322 + 50740
parameters <- c(rep(0.1,(L+1)*2))

likelihood_base(c(0.1, 2))
test <- gumbel(df.wind.all)
likelihood_standard(parameters)
likelihood_lasso(parameters)
likelihood_ridge(parameters)
likelihood_elastic.net(parameters)
likelihood_adaptive.LASSO(parameters)
likelihood_adaptive.EN(parameters)

start_time <- Sys.time()
likelihood_lasso(parameters)
end_time <- Sys.time()
end_time - start_time

final.Ridge <- c(0.12)
final.Lasso <- c(0.0086)
final.EN <- c(0.007, 0.106)
final.adaptive.Lasso <- c(0.072, 0.088)
final.adaptive.EN <- c(0.024, 0.067, 0.013)

start_time <- Sys.time()
optim.ll <- optim(par=parameters, fn = ll, method = "BFGS", control = list(maxit = 1e6), hessian = TRUE)
end_time <- Sys.time()
end_time - start_time

ll <- likelihood_elastic.net
first.lmbd <- final.EN[1]
second.lmbd <- final.EN[2]

ll <- likelihood_adaptive.LASSO
first.lmbd <- final.adaptive.Lasso[1]
third.gamma <- final.adaptive.Lasso[2]

start_time <- Sys.time()
optim.ll.lasso <- optim(par=parameters, fn = ll, method = "BFGS", control = list(maxit = 1e6), hessian = TRUE)
end_time <- Sys.time()
end_time - start_time

optim.base <- optim(par=c(0.1, 2), fn = likelihood_base, method = "BFGS", control = list(maxit = 1e6), hessian = TRUE)
optim.standard <- optim(par=parameters, fn = likelihood_standard, method = "BFGS", control = list(maxit = 1e6), hessian = TRUE)

first.lmbd <- final.Lasso
optim.Lasso <- optim(par=parameters, fn = likelihood_lasso, method = "BFGS", control = list(maxit = 1e6), hessian = TRUE)
optim.Ridge <- optim(par=parameters, fn = likelihood_ridge, method = "BFGS", control = list(maxit = 1e6), hessian = TRUE)
optim.EN <- optim(par=parameters, fn = likelihood_elastic.net, method = "BFGS", control = list(maxit = 1e6), hessian = TRUE)

ini.par.Lasso <- optim.Lasso$par
ini.par.EN <- optim.EN$par
optim.adaptive.LASSO <- optim(par=parameters, fn = likelihood_adaptive.LASSO, method = "BFGS", control = list(maxit = 1e6), hessian = TRUE)
optim.adaptive.EN <- optim(par=parameters, fn = likelihood_adaptive.EN, method = "BFGS", control = list(maxit = 1e6), hessian = TRUE)


