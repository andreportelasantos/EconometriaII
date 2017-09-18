# Define função do modelo ARMA(p,q)-GARCH(r,s) (sem intercepto na média)
# Argumentos da função:
# params: vetor de parâmetros 
# data: série temporal
# np: ordem AR
# nq: ordem MA
# nr: ordem ARCH
# ns: ordem GARCH
ARMAGARCH <- function(params,data,np,nq,nr,ns) {
  T <- length(data)
  errors <- integer(T)
  sigma2 <- integer(T)
  m <- 1+max(np,nq,nr,ns)
  sigma2[1:m] <- var(data)
  for (t in m:T) {
    # Recursão AR
    errors[t] <- data[t]
    for (i in 1:np) {
      errors[t] <- errors[t] - params[i]*data[t-i]
    }
    # Recursão MA
    for (i in 1:nq) {
      errors[t] <- errors[t] - params[np+i]*errors[t-i]
    }
    errors[t] <- errors[t];
    # Recursão ARCH
    sigma2[t] <- params[np+nq+1]
    for (i in 1:nr) {
      sigma2[t] <- sigma2[t] + params[np+nq+1+i]*errors[t-i]^2
    }
    # Recursão GARCH
    for (i in 1:ns) {
      sigma2[t] <- sigma2[t] + params[np+nq+nr+1+i]*sigma2[t-i]
    }
    sigma2[t] <- sigma2[t]
  }
  verossim <- 0.5*(2*log(sqrt(sigma2)) + log(sigma2) + errors^2/sigma2 + log(2*pi))
  return(list(LLF=sum(verossim),sigma2=sigma2,residuals=errors/sqrt(sigma2)))
}

# Baixa série de preços das ações da Apple (AAPL)
# install.packages("BatchGetSymbols")
library(BatchGetSymbols)
library(ggplot2)
my.ticker <- c('AAPL')
first.date <- Sys.Date()-1500
last.date <- Sys.Date()
l.out <- BatchGetSymbols(tickers = my.ticker,first.date = first.date,last.date = last.date)
returns <- data.frame(retornos=diff(log(l.out$df.tickers$price.adjusted))*100,datas=l.out$df.tickers$ref.date[2:l.out$df.control$total.obs])

# Estima modelo ARMA(1,1)-GARCH(1,1) 
init.params <- c(0.1,0.1,0.1,0.1,0.5)
# Define restrições
ui <- rbind(c(0,0,1,0,0),c(0,0,0,1,0),c(0,0,0,0,1),c(0,0,0,-1,-1))
ci <- c(0,0,0,-1)
ARMAGARCH.optim <- function(params,data,np,nq,nr,ns) {
  ARMAGARCH(params,data,np,nq,nr,ns)$LLF
}
resultados <- constrOptim(init.params,ARMAGARCH.optim,data=returns[,1],np=1,nq=1,nr=1,ns=1,grad=NULL,ui=ui,ci=ci)
options(scipen=999)
print(resultados$par)

# Retorna variáveis e avalia qualidade do ajuste
resultados.finais <- ARMAGARCH(resultados$par,data=returns[,1],np=1,nq=1,nr=1,ns=1)
df <- data.frame(returns,sigma2=resultados.finais$sigma2,residuals=resultados.finais$residuals)
# Faz gráfico 
require(gridExtra)
require(forecast)
p1 <- ggplot(data = returns, aes(x = datas, y = retornos))
p1 <- p1 + geom_line()
p1 <- p1 + labs(x = 'Dates', y = 'Retornos')

p2 <- ggplot(data = df, aes(x = datas, y = sqrt(sigma2)))
p2 <- p2 + geom_line()
p2 <- p2 + labs(x = 'Dates', y = 'Desvio padrão condicional')

p3 <- ggAcf(df$residuals^2, main="ACF do quadrado dos resíduos padronizados")

grid.arrange(p1, p2, p3, ncol=1)
