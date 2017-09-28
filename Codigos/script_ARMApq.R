
# Define função do modelo ARMA(p,q) (sem intercepto)
# Argumentos da função:
# params: vetor de parâmetros 
# data: série temporal
# np: ordem AR(p)
# nq: ordem MA(q)
ARMApq <- function(params,data,np,nq) {
  T <- length(data)
  errors <- integer(T)
  m <- 1+max(np,nq)
  for (t in m:T) {
    errors[t] <- data[t]
    for (i in 1:np) {
      errors[t] <- errors[t] - params[i]*data[t-i]
    }
    for (i in 1:nq) {
      errors[t] <- errors[t] - params[np+i]*errors[t-i]
    }
    errors[t] <- errors[t];
  }
  sigma2 <- t(errors)%*%errors/(T-m-1)
  #verossim <- 0.5*(2*log(sqrt(sigma2)) + log(sigma2) + errors^2/sigma2 + log(2*pi))
  verossim <- 0.5*(sum(log(sigma2)) + sum((errors^2)/sigma2)  +  T*log(2*pi));
  return(sum(verossim))
}

# Exemplo com dados simulados através de um processo ARMA(2,2) estacionario SEM constante
T <- 5000
set.seed(1234)
# y_t = 0.5y_t-1 + 0+4y_t-2 + e_t - 0.2e_t-1 + 0.3e_t-2
sim.data <- arima.sim(n = T, list(ar = c(0.5, -0.4), ma = c(-0.2, 0.3)),n.start=1000)

# Estima parâmetros
init.params <- c(0.1,0.1,0.1,0.1)
np <- 2
nq <- 2
resultados <- nlm(ARMApq, init.params, sim.data, np, nq, hessian=TRUE)
print(resultados$estimate)

# Testa com outro otimizador
resultados2 <- optim(par=init.params, ARMApq, data=sim.data, np=2, nq=2, method = "BFGS")
print(resultados2$par)

# Compara com função arima do R
resultados3 <- arima(sim.data,order = c(2,0,2),method="ML",include.mean=FALSE)
print(resultados3$coef)

# Testa função da biblioteca optimx
# install.packages("optimx")
library(optimx)
resultados4 <- optimx(par=init.params, ARMApq, data=sim.data, np=2, nq=2)
print(resultados4)

