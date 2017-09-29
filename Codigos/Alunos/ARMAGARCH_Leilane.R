list.of.packages <- c("BatchGetSymbols","ggplot2","fGarch","gridExtra","forecast","DescTools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(BatchGetSymbols)
library(ggplot2)
library(fGarch)
library(gridExtra)
library(forecast)
library(DescTools)
options(scipen=999)


# Função para estimar parâmetros de um modelo ARMA(1,1)-GARCH(1,1) sem constante na média condicional utilizando máxima verossimilhança condicional

# Argumentos da função:
# params: vetor de parâmetros 
# data: série temporal
# np: ordem AR
# nq: ordem MA
# nr: ordem ARCH
# ns: ordem GARCH
# distr: distribuição, pode ser 'gauss' (default), 't' (t de Student) e 'ged' (Generalized Error Distribution)
# gl: graus de liberdade, 0 para ser estimado com os outros parâmetros

ARMAGARCH <- function(params,data,np,nq,nr,ns,distr="gauss") {
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
  }
  if(distr=='ged') {
    gl <- params[np+nq+nr+ns+2]
    lamb <- (2^(-2/gl) * gamma(1/gl)/gamma(3/gl))^(0.5);
    verossim <- -(log(gl) - 0.5*(abs(errors/lamb)^gl) - log(lamb) - log(2^(1+1/gl)) - lgamma(1/gl))
  } else if(distr=='t') {
    gl <- params[np+nq+nr+ns+2]
    verossim <- -(lgamma((gl+1)/2)-lgamma(gl/2)-0.5*log((gl-2)*pi)) + 0.5*(gl+1)*log(1+(errors^2/((gl-2)*sigma2))) + 0.5*log(sigma2)
  }
#  else if(distr=='t'&& gl!=0) {       # t de Student com gl pré-especificado
#    verossim <- -sum(0.5*(gl+1)*log(1+(errors^2)/((gl-2)*sigma2)) + 0.5*log(sigma2))
#  }
  else {    # se distr == 'gauss' ou qualquer outro valor não programado
    verossim <- 0.5*(log(sigma2) + (errors^2)/sigma2  +  log(2*pi));
  }
  return(list(LLF=sum(verossim),sigma2=sigma2,residuals=errors/sqrt(sigma2)))
}

# Parâmetros iniciais da estimação
init.params <- c(0.1,0.1,0.1,0.1,0.5)
# Define restrições
ui <- rbind(c(0,0,1,0,0),c(0,0,0,1,0),c(0,0,0,0,1),c(0,0,0,-1,-1))
ci <- c(0,0,0,-1)

# Parâmetros iniciais da estimação e restrições para t de Student e GED sem gl pré-especificado
init.params2 <- c(0.1,0.1,0.1,0.1,0.5,4)
ui2 <- rbind(c(0,0,1,0,0,0),c(0,0,0,1,0,0),c(0,0,0,0,1,0),c(0,0,0,-1,-1,0),c(0,0,0,0,0,1))
ci2 <- c(0,0,0,-1,3)
ci3 <- c(0,0,0,-1,0) # valor da restrição para GED


# Helper para maximizar a função verossimilhança condicional
ARMAGARCH.optim <- function(params,data,np,nq,nr,ns,distr) {
  ARMAGARCH(params,data,np,nq,nr,ns,distr)$LLF
}


## Estima parâmetros do modelo ARMA(1,1)-GARCH(1,1) para as séries de retorno das ações da Aaple (AAPL)

# Baixa últimos 1500 dados mais recentes
my.ticker <- c('AAPL')
first.date <- Sys.Date()-1500
last.date <- Sys.Date()
l.out <- BatchGetSymbols(tickers = my.ticker,first.date = first.date,last.date = last.date)
returns <- data.frame(retornos=diff(log(l.out$df.tickers$price.adjusted))*100,datas=l.out$df.tickers$ref.date[2:l.out$df.control$total.obs])


# distr = 'gauss'
# Exibe parâmetros estimados
resultados.gauss <- constrOptim(init.params,ARMAGARCH.optim,data=returns[,1],np=1,nq=1,nr=1,ns=1,distr="gauss",grad=NULL,ui=ui,ci=ci)
print(resultados.gauss$par)

# Retorna variáveis e avalia qualidade do ajuste
resultados.finais.gauss <- ARMAGARCH(resultados.gauss$par,data=returns[,1],np=1,nq=1,nr=1,ns=1,distr="gauss")
df <- data.frame(returns,sigma2=resultados.finais.gauss$sigma2,residuals=resultados.finais.gauss$residuals)


# distr == 't'
# Exibe parâmetros estimados
resultados.t <- constrOptim(init.params2,ARMAGARCH.optim,data=returns[,1],np=1,nq=1,nr=1,ns=1,distr="t",grad=NULL,ui=ui2,ci=ci2)
print(resultados.t$par)

# Retorna variáveis e avalia qualidade do ajuste
resultados.finais.t <- ARMAGARCH(resultados.t$par,data=returns[,1],np=1,nq=1,nr=1,ns=1,distr="t")
df <- data.frame(returns,sigma2=resultados.finais.t$sigma2,residuals=resultados.finais.t$residuals)


# distr == 'ged'
# Exibe parâmetros estimados
resultados.ged <- constrOptim(init.params2,ARMAGARCH.optim,data=returns[,1],np=1,nq=1,nr=1,ns=1,distr="ged",grad=NULL,ui=ui2,ci=ci3)
print(resultados.ged$par)

# Retorna variáveis e avalia qualidade do ajuste
resultados.finais.ged <- ARMAGARCH(resultados.ged$par,data=returns[,1],np=1,nq=1,nr=1,ns=1,distr="ged")
df <- data.frame(returns,sigma2=resultados.finais.ged$sigma2,residuals=resultados.finais.ged$residuals)


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


## Compara com os parâmetros estimados através da biblitoeca fGarch do R

# Estima modelo
resultados.finais.gauss2 <- garchFit(formula = ~ arma(1,1) + garch(1, 1), data = returns[,1], cond.dist = c("QMLE"), include.mean = F)
resultados.finais.t2 <- garchFit(formula = ~ arma(1,1) + garch(1, 1), data = returns[,1], cond.dist = c("std"), include.mean = F)
resultados.finais.ged2 <- garchFit(formula = ~ arma(1,1) + garch(1, 1), data = returns[,1], cond.dist = c("ged"), include.mean = F)

# Compara a verossimilhança obtida com código próprio e com pacote do R
LLF.meu.codigo.gauss <- print(-ARMAGARCH(resultados.gauss$par,data=returns[,1],np=1,nq=1,nr=1,ns=1,distr="gauss")$LLF)
LLF.fGarch.gauss <- print(-ARMAGARCH(c(-0.88144668,0.89127409,0.29517193,0.08502336,0.77621936),data=returns[,1],np=1,nq=1,nr=1,ns=1,distr="gauss")$LLF)

LLF.meu.codigo.t <- print(-ARMAGARCH(resultados.t$par,data=returns[,1],np=1,nq=1,nr=1,ns=1,distr="t")$LLF)
LLF.fGarch.t <- print(-ARMAGARCH(c(-0.9061262,0.9133388,0.1836010,0.1005357,0.8258584,3.8675290),data=returns[,1],np=1,nq=1,nr=1,ns=1,distr="t")$LLF)

LLF.meu.codigo.ged <- print(-ARMAGARCH(resultados.ged$par,data=returns[,1],np=1,nq=1,nr=1,ns=1,distr="ged")$LLF)
LLF.fGarch.ged <- print(-ARMAGARCH(c(-0.4866899,0.4578924,0.2096673,0.0946025,0.8057347,1.0893742),data=returns[,1],np=1,nq=1,nr=1,ns=1,distr="ged")$LLF)