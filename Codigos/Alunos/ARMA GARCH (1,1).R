# Define função do modelo ARMA(p,q)-GARCH(r,s) (sem intercepto na mÉdia)
# Argumentos da função:
# params: vetor de parametros 
# data: série temporal
# np: ordem AR
# nq: ordem MA
# nr: ordem ARCH
# ns: ordem GARCH
ARMAGARCH <- function(params,data,np,nq,nr,ns,dens) {
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
  if (dens=="t") {
    v<-params[np+nq+nr+ns+2]
    verossim <- -(T*log(gamma((v+1)/2)/(sqrt((v-2)*pi)*gamma(v/2)))-sum(log(sigma2))/2-sum(((v+1)/2)*log(1+errors^2/(sigma2*(v-2)))))
  } else if (dens=="GED") {
    v <- params[np+nq+nr+ns+2]
    lambda <- (2^(-2/v) * gamma(1/v)/gamma(3/v))^0.5
    verossim <- -T*log(v) + sum(0.5 * abs(errors/(lambda))^v) + T*log(lambda * 2^(1+1/v) * gamma(1/v))
  } else if (dens=="PED") {
    lambda<-params[np+nq+nr+ns+2]
    gam<-params[np+nq+nr+ns+3]
    verossim<- -(T*log(1-gam^2)-T*log(2)-T*log(gamma(1+1/lambda))-T*log(lambda)/lambda-sum((1/lambda)*((abs(data-gam*abs(data)))^lambda)))
  } else {
    verossim <- 0.5*(sum(log(sigma2)) + sum((errors^2)/sigma2)  +  T*log(2*pi));
  }
  return(list(LLF=sum(verossim),sigma2=sigma2,residuals=errors/sqrt(sigma2)))
}


# Baixa série de preços da Apple
# install.packages("BatchGetSymbols")
library(BatchGetSymbols)
library(ggplot2)
my.ticker <- c("AAPL")
first.date <- Sys.Date()-1500
last.date <- Sys.Date()
l.out <- BatchGetSymbols(tickers = my.ticker,first.date = first.date,last.date = last.date)
returns <- data.frame(retornos=diff(log(l.out$df.tickers$price.adjusted))*100,datas=l.out$df.tickers$ref.date[2:l.out$df.control$total.obs])

# Estima modelo ARMA,(1,1)-GARCH(1,1) - NORMAL
init.params <- c(0.1,0.1,0.1,0.1,0.5)
# Define restrições
ui <- rbind(c(0,0,1,0,0),c(0,0,0,1,0),c(0,0,0,0,1),c(0,0,0,-1,-1))
ci <- c(0,0,0,-1)
ARMAGARCH.optim <- function(params,data,np,nq,nr,ns,dens) {
  ARMAGARCH(params,data,np,nq,nr,ns,dens)$LLF
}
resultados <- constrOptim(init.params,ARMAGARCH.optim,data=returns[,1],np=1,nq=1,nr=1,ns=1,grad=NULL,ui=ui,ci=ci,dens="normal")
options(scipen=999)
print(resultados$par)
resultados.finais.normal <- ARMAGARCH(resultados$par,data=returns[,1],np=1,nq=1,nr=1,ns=1,dens="normal")
df.normal <- data.frame(returns,sigma2=resultados.finais.normal$sigma2,residuals=resultados.finais.normal$residuals)
library(fGarch)
garchFit(formula=~arma(1,1)+garch(1,1),data=returns[,1],cond.dist=c("norm"),include.mean=F)

# Estima modelo ARMA(1,1)-GARCH(1,1) -T de Student
init.params <- c(0.1,0.1,0.1,0.1,0.5,10)
# Define restrições
ui <- rbind(c(0,0,1,0,0,0),c(0,0,0,1,0,0),c(0,0,0,0,1,0),c(0,0,0,-1,-1,0),c(0,0,0,0,0,1),c(0,0,0,0,0,-1))
ci <- c(0,0,0,-1,2+.Machine$double.eps,-30)
ARMAGARCH.optim <- function(params,data,np,nq,nr,ns,dens) {
  ARMAGARCH(params,data,np,nq,nr,ns,dens)$LLF
}
resultados <- constrOptim(init.params,ARMAGARCH.optim,data=returns[,1],np=1,nq=1,nr=1,ns=1,grad=NULL,ui=ui,ci=ci,dens="t")
options(scipen=999)
print(resultados$par)
resultados.finais.t <- ARMAGARCH(resultados$par,data=returns[,1],np=1,nq=1,nr=1,ns=1,dens="t")
df.t <- data.frame(returns,sigma2=resultados.finais.t$sigma2,residuals=resultados.finais.t$residuals)

# Estima modelo ARMA(1,1)-GARCH(1,1) -Generalized Exponential Distribution
init.params <- c(0.1,0.1,0.1,0.1,0.1,2)
# Define restrições
ui <- rbind(c(0,0,1,0,0,0),c(0,0,0,1,0,0),c(0,0,0,0,1,0),c(0,0,0,-1,-1,0),c(0,0,0,0,0,1))
ci <- c(0,0,0,-1,1)
ARMAGARCH.optim <- function(params,data,np,nq,nr,ns,dens) {
  ARMAGARCH(params,data,np,nq,nr,ns,dens)$LLF
}
resultados <- constrOptim(init.params,ARMAGARCH.optim,data=returns[,1],np=1,nq=1,nr=1,ns=1,grad=NULL,ui=ui,ci=ci,dens="GED")
options(scipen=999)
print(resultados$par)
resultados.finais.GED <- ARMAGARCH(resultados$par,data=returns[,1],np=1,nq=1,nr=1,ns=1,dens="GED")
df.GED <- data.frame(returns,sigma2=resultados.finais.GED$sigma2,residuals=resultados.finais.GED$residuals)

# Estima modelo ARMA(1,1)-GARCH(1,1) -Pòwer Exponential Distribution
init.params <- c(0.1,0.1,0.1,0.1,0.1,3,0.1)
# Define restrições
ui <- rbind(c(0,0,1,0,0,0,0),c(0,0,0,1,0,0,0),c(0,0,0,0,1,0,0),c(0,0,0,0,0,1,0),c(0,0,0,-1,-1,0,0))
ci <- c(0,0,0,0,-1)
ARMAGARCH.optim <- function(params,data,np,nq,nr,ns,dens) {
  ARMAGARCH(params,data,np,nq,nr,ns,dens)$LLF
}
resultados <- constrOptim(init.params,ARMAGARCH.optim,data=returns[,1],np=1,nq=1,nr=1,ns=1,grad=NULL,ui=ui,ci=ci,dens="PED")
options(scipen=999)
print(resultados$par)
resultados.finais.PED <- ARMAGARCH(resultados$par,data=returns[,1],np=1,nq=1,nr=1,ns=1,dens="PED")
df.PED <- data.frame(returns,sigma2=resultados.finais.PED$sigma2,residuals=resultados.finais.PED$residuals)

# Retorna variáveis e avalia qualidade do ajuste

# Faz gráfico 
require(gridExtra)
require(forecast)
require(ggplot2)
p1.normal <- ggplot(data = returns, aes(x = datas, y = retornos))
p1.normal <- p1.normal + geom_line(linetype="dashed") + geom_line(data = returns-df.normal$residuals, colour = "red")
p1.normal <- p1.normal + labs(x = 'Dates', y = 'Log-Retornos Previstos')
p1.normal <- p1.normal + ggtitle("Standard Normal")

p2.normal <- ggplot(data = df.normal, aes(x = datas, y = sqrt(sigma2)))
p2.normal <- p2.normal + geom_line(colour = "gold4")
p2.normal <- p2.normal + labs(x = 'Dates', y = 'Desvio padrão condicional')

p3.normal <- ggAcf(df.normal$residuals^2, main="ACF do quadrado dos resíduos padronizados")

p1.t <- ggplot(data = returns, aes(x = datas, y = retornos))
p1.t <- p1.t + geom_line(linetype="dashed") + geom_line(data = returns-df.t$residuals, colour = "red")
p1.t <- p1.t + labs(x = 'Dates', y = 'Log-Retornos Previstos')
p1.t <- p1.t + ggtitle("t Student")

p2.t <- ggplot(data = df.t, aes(x = datas, y = sqrt(sigma2)))
p2.t <- p2.t + geom_line(colour = "gold4")
p2.t <- p2.t + labs(x = 'Dates', y = 'Desvio padrão condicional')

p3.t <- ggAcf(df.t$residuals^2, main="ACF do quadrado dos resíduos padronizados")

p1.GED <- ggplot(data = returns, aes(x = datas, y = retornos))
p1.GED <- p1.GED + geom_line(linetype="dashed") + geom_line(data = returns-df.GED$residuals, colour = "red")
p1.GED <- p1.GED + labs(x = 'Dates', y = 'Log-Retornos Previstos')
p1.GED <- p1.GED + ggtitle("Generalized Error Distribution")

p2.GED <- ggplot(data = df.GED, aes(x = datas, y = sqrt(sigma2)))
p2.GED <- p2.GED + geom_line(colour = "gold4")
p2.GED <- p2.GED + labs(x = 'Dates', y = 'Desvio padrão condicional')

p3.GED <- ggAcf(df.GED$residuals^2, main="ACF do quadrado dos resíduos padronizados")

p1.PED <- ggplot(data = returns, aes(x = datas, y = retornos))
p1.PED <- p1.PED + geom_line(linetype="dashed") + geom_line(data = returns-df.PED$residuals, colour = "red")
p1.PED <- p1.PED + labs(x = 'Dates', y = 'Log-Retornos Previstos')
p1.PED <- p1.GED + ggtitle("Power Exponential Distribution")

p2.PED <- ggplot(data = df.PED, aes(x = datas, y = sqrt(sigma2)))
p2.PED <- p2.PED + geom_line(colour = "gold4")
p2.PED <- p2.PED + labs(x = 'Dates', y = 'Desvio padrão condicional',main='Power Exponential Distribution')

p3.PED <- ggAcf(df.PED$residuals^2, main="ACF do quadrado dos resíduos padronizados")

grid.arrange(p1.normal, p1.t, p1.GED, p1.PED,p2.normal, p2.t, p2.GED, p2.PED, p3.normal, p3.t, p3.GED, p3.PED, ncol=4, nrow=3)

RESUMO<-c(sum(df.normal$residuals^2),sum(df.t$residuals^2),sum(df.GED$residuals^2),sum(df.PED$residuals^2))
