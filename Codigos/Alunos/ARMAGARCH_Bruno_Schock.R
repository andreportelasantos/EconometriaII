# Define função do modelo ARMA(p,q)-GARCH(r,s) (sem intercepto na média)
# Argumentos da função:
# params: vetor de parâmetros 
# data: série temporal
# np: ordem AR
# nq: ordem MA
# nr: ordem ARCH
# ns: ordem GARCH
library(BatchGetSymbols)
library(ggplot2)
require(gridExtra)
require(forecast)
library(fGarch)


#ARMAGARCH
{ARMAGARCH <- function(params,data,np,nq,nr,ns, distr) {
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
  if(distr=="STT"){
  v <- params[np+nq+nr+ns+2]
  
    verossimt <- -log(gamma((v+1)/2)) + log(gamma((v/2))) + 0.5*log(pi) + 0.5*log(v-2) + log(sigma2) + ((v+1)/2)*log(1+(errors^2)/(sigma2*(v-2)))
  return(list(LLF=sum(verossimt),sigma2=sigma2, residuals=errors/sqrt(sigma2)))}
  
  else if(distr=="STTpré"){
    v <- 4
    verossimt2 <- -log(gamma((v+1)/2)) + log(gamma((v/2))) + 0.5*log(pi) + 0.5*log(v-2) + log(sigma2) + ((v+1)/2)*log(1+(errors^2)/(sigma2*(v-2)))
    return(list(LLF=sum(verossimt2),sigma2=sigma2, residuals=errors/sqrt(sigma2)))}
  
  else if(distr=="NORMAL"){
    
    verossim <- 0.5*(2*log(sqrt(sigma2)) + log(sigma2) + errors^2/sigma2 + log(2*pi))
    return(list(LLF=sum(verossim),sigma2=sigma2,residuals=errors/sqrt(sigma2)))}
  
  
  else if (distr == "GEDpré"){
    v <- 1.5
    lambda <- (2^(-2/v) * gamma(1/v)/gamma(3/v))^0.5
    verossimgedp <- -log(v) + 0.5 * abs(errors/(lambda))^v + log(lambda * 2^(1+1/v) * gamma(1/v))
    return(list(LLF=sum(verossimgedp),sigma2=sigma2, residuals=errors/sqrt(sigma2)))}

 else (distr == "GED")
v <- params[np+nq+nr+ns+2]
lambda <- (2^(-2/v) * gamma(1/v)/gamma(3/v))^0.5
verossimged <- -log(v) + 0.5*abs(errors/lambda)^v + log(lambda*2^(1+1/2)*gamma(1/v))
return(list(LLF=sum(verossimged),sigma2=sigma2, residuals=errors/sqrt(sigma2)))}
  }

#Teste
my.ticker <- c('AAPL')
first.date <- Sys.Date()-1500
last.date <- Sys.Date()
l.out <- BatchGetSymbols(tickers = my.ticker,first.date = first.date,last.date = last.date)
returns <- data.frame(retornos=diff(log(l.out$df.tickers$price.adjusted))*100,datas=l.out$df.tickers$ref.date[2:l.out$df.control$total.obs])


  
#Normal
{# Estima modelo ARMA(1,1)-GARCH(1,1) para a normal 
init.params <- c(0.1,0.1,0.1,0.1,0.5)
# Define restrições
ui <- rbind(c(0,0,1,0,0),c(0,0,0,1,0),c(0,0,0,0,1),c(0,0,0,-1,-1))
ci <- c(0,0,0,-1)
ARMAGARCH.optim <- function(params,data,np,nq,nr,ns, distr) {
  ARMAGARCH(params,data,np,nq,nr,ns, distr)$LLF
}
resultados <- constrOptim(init.params,ARMAGARCH.optim,data=returns[,1],np=1,nq=1,nr=1,ns=1,grad=NULL,ui=ui,ci=ci,distr="NORMAL")
options(scipen=999)
resnorm <- print(resultados$par)

#Teste
library(BatchGetSymbols)
library(ggplot2)
my.ticker <- c('AAPL')
first.date <- Sys.Date()-1500
last.date <- Sys.Date()
l.out <- BatchGetSymbols(tickers = my.ticker,first.date = first.date,last.date = last.date)
returns <- data.frame(retornos=diff(log(l.out$df.tickers$price.adjusted))*100,datas=l.out$df.tickers$ref.date[2:l.out$df.control$total.obs])

# Retorna variáveis e avalia qualidade do ajuste
resultados.finais <- ARMAGARCH(resultados$par,data=returns[,1],np=1,nq=1,nr=1,ns=1, distr="NORMAL")
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

Normal <- grid.arrange(p1, p2, p3, ncol=1)
grid.arrange(p1, p2, p3, ncol=1)
} 

#T com parâmetros pré-definidos
{# Estima modelo ARMA(1,1)-GARCH(1,1) para a T pré 
  init.params <- c(0.1,0.1,0.1,0.1,0.5)
  # Define restrições
  ui <- rbind(c(0,0,1,0,0),c(0,0,0,1,0),c(0,0,0,0,1),c(0,0,0,-1,-1))
  ci <- c(0,0,0,-1)
  ARMAGARCH.optim <- function(params,data,np,nq,nr,ns, distr) {
    ARMAGARCH(params,data,np,nq,nr,ns, distr)$LLF
  }
  resultados <- constrOptim(init.params,ARMAGARCH.optim,data=returns[,1],np=1,nq=1,nr=1,ns=1,grad=NULL,ui=ui,ci=ci,distr="STTpré")
  options(scipen=999)
  restpré <- print(resultados$par)
  
  #Teste
  library(BatchGetSymbols)
  library(ggplot2)
  my.ticker <- c('AAPL')
  first.date <- Sys.Date()-1500
  last.date <- Sys.Date()
  l.out <- BatchGetSymbols(tickers = my.ticker,first.date = first.date,last.date = last.date)
  returns <- data.frame(retornos=diff(log(l.out$df.tickers$price.adjusted))*100,datas=l.out$df.tickers$ref.date[2:l.out$df.control$total.obs])
  
  # Retorna variáveis e avalia qualidade do ajuste
  resultados.finais <- ARMAGARCH(resultados$par,data=returns[,1],np=1,nq=1,nr=1,ns=1, distr="STTpré")
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
  
STTpré <-  grid.arrange(p1, p2, p3, ncol=1)
grid.arrange(p1, p2, p3, ncol=1)
}

#GED com parâmetros pré-definidos
{# Estima modelo ARMA(1,1)-GARCH(1,1) para a T pré 
  init.params <- c(0.1,0.1,0.1,0.1,0.5)
  # Define restrições
  ui <- rbind(c(0,0,1,0,0),c(0,0,0,1,0),c(0,0,0,0,1),c(0,0,0,-1,-1))
  ci <- c(0,0,0,-1)
  ARMAGARCH.optim <- function(params,data,np,nq,nr,ns, distr) {
    ARMAGARCH(params,data,np,nq,nr,ns, distr)$LLF
  }
  resultados <- constrOptim(init.params,ARMAGARCH.optim,data=returns[,1],np=1,nq=1,nr=1,ns=1,grad=NULL,ui=ui,ci=ci,distr="GEDpré")
  options(scipen=999)
  resgedpré <- print(resultados$par)
  
  #Teste
  library(BatchGetSymbols)
  library(ggplot2)
  my.ticker <- c('AAPL')
  first.date <- Sys.Date()-1500
  last.date <- Sys.Date()
  l.out <- BatchGetSymbols(tickers = my.ticker,first.date = first.date,last.date = last.date)
  returns <- data.frame(retornos=diff(log(l.out$df.tickers$price.adjusted))*100,datas=l.out$df.tickers$ref.date[2:l.out$df.control$total.obs])
  
  # Retorna variáveis e avalia qualidade do ajuste
  resultados.finais <- ARMAGARCH(resultados$par,data=returns[,1],np=1,nq=1,nr=1,ns=1, distr="GEDpré")
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
  
GEDpré <- grid.arrange(p1, p2, p3, ncol=1)
grid.arrange(p1, p2, p3, ncol=1)
}

#T com parâmetros a serem estimados
{# Para a t-student
init.params <- c(0.1,0.1,0.1,0.1,0.5,3)
# Define restrições
ui <- rbind(c(0,0,1,0,0,0),c(0,0,0,1,0,0),c(0,0,0,0,1,0),c(0,0,0,-1,-1,0),c(0,0,0,0,0,1),c(0,0,0,0,0,-1))
ci <- c(0,0,0,-1,2+.Machine$double.eps,-30)
ARMAGARCH.optim <- function(params,data,np,nq,nr,ns,distr) {
  ARMAGARCH(params,data,np,nq,nr,ns,distr)$LLF
}
resultados <- constrOptim(init.params,ARMAGARCH.optim,data=returns[,1],np=1,nq=1,nr=1,ns=1,grad=NULL,ui=ui,ci=ci,distr="STT")
options(scipen=999)
rest <- print(resultados$par)
resultados.finais.t <- ARMAGARCH(resultados$par,data=returns[,1],np=1,nq=1,nr=1,ns=1,distr="STT")
df.t <- data.frame(returns,sigma2=resultados.finais.t$sigma2,residuals=resultados.finais.t$residuals)

#Teste
library(BatchGetSymbols)
library(ggplot2)
my.ticker <- c('AAPL')
first.date <- Sys.Date()-1500
last.date <- Sys.Date()
l.out <- BatchGetSymbols(tickers = my.ticker,first.date = first.date,last.date = last.date)
returns <- data.frame(retornos=diff(log(l.out$df.tickers$price.adjusted))*100,datas=l.out$df.tickers$ref.date[2:l.out$df.control$total.obs])

# Retorna variáveis e avalia qualidade do ajuste
resultados.finais <- ARMAGARCH(resultados$par,data=returns[,1],np=1,nq=1,nr=1,ns=1, distr="STT")
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

STT <- grid.arrange(p1, p2, p3, ncol=1)
grid.arrange(p1, p2, p3, ncol=1)
}

#GED com parâmetros a serem estimados
{# Para a GED
  init.params <- c(0.1,0.1,0.1,0.1,0.5,5)
  # Define restrições
  ui <- rbind(c(0,0,1,0,0,0),c(0,0,0,1,0,0),c(0,0,0,0,1,0),c(0,0,0,-1,-1,0),c(0,0,0,0,0,1),c(0,0,0,0,0,-1))
  ci <- c(0,0,0,-1,2+.Machine$double.eps,-30)
  ARMAGARCH.optim <- function(params,data,np,nq,nr,ns,distr) {
    ARMAGARCH(params,data,np,nq,nr,ns,distr)$LLF
  }
  resultados <- constrOptim(init.params,ARMAGARCH.optim,data=returns[,1],np=1,nq=1,nr=1,ns=1,grad=NULL,ui=ui,ci=ci,distr="GED")
  options(scipen=999)
  resged <- print(resultados$par)
  resultados.finais.t <- ARMAGARCH(resultados$par,data=returns[,1],np=1,nq=1,nr=1,ns=1,distr="GED")
  df.t <- data.frame(returns,sigma2=resultados.finais.t$sigma2,residuals=resultados.finais.t$residuals)
  
  #Teste
  library(BatchGetSymbols)
  library(ggplot2)
  my.ticker <- c('AAPL')
  first.date <- Sys.Date()-1500
  last.date <- Sys.Date()
  l.out <- BatchGetSymbols(tickers = my.ticker,first.date = first.date,last.date = last.date)
  returns <- data.frame(retornos=diff(log(l.out$df.tickers$price.adjusted))*100,datas=l.out$df.tickers$ref.date[2:l.out$df.control$total.obs])
  
  # Retorna variáveis e avalia qualidade do ajuste
  resultados.finais <- ARMAGARCH(resultados$par,data=returns[,1],np=1,nq=1,nr=1,ns=1, distr="GED")
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
  
GED <-  grid.arrange(p1, p2, p3, ncol=1)
}

#Arranjo com a Normal
grid.arrange(Normal, STT, STTpré, GED, GEDpré, ncol=5, nrow=1)  

#Arranjo sem a normal
grid.arrange(STT, STTpré, GED, GEDpré, ncol=4, nrow=1)


#Comparativo do pacote fGarch com a Normal
m1 <- garchFit(~arma(1,1)+garch(1,1), data=returns[,1], trace=F, cond.dist = "norm", description = FALSE)
coef(m1)
resnorm

#Comparativo do pacote fGarch com a t-pré
m2 <- garchFit(~arma(1,1)+garch(1,1), data=returns[,1], trace=F, cond.dist = "std", description = FALSE)
coef(m2)
restpré

#Comparativo do pacote fGarch com a GED-pré
m3 <- garchFit(~arma(1,1)+garch(1,1), data=returns[,1], trace=F, cond.dist = "ged", description = FALSE)
coef(m3)
resgedpré

#Comparativo do pacote fGarch com a T
coef(m2)
rest

#Comparativo do pacote fGarch com a GED
coef(m3)
resged
