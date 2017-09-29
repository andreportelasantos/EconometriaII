#Clear Workspace
rm(list=ls())

list.of.packages <- c("BatchGetSymbols","ggplot2","fGarch","gridExtra","forecast")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(BatchGetSymbols)
library(ggplot2)
library(fGarch)
library(gridExtra)
library(forecast)

# Função para estimar parâmetros de um modelo ARMA(1,1)-GARCH(1,1) sem constante na média condicional utilizando máxima verossimilhança condicional

# Argumentos da função:
# params: vetor de parâmetros 
# data: série temporal
# np: ordem AR
# nq: ordem MA
# nr: ordem ARCH
# ns: ordem GARCH

#Distribuição Estatística para os EMV

distribuicao=1


# Baixa últimos 1500 dados mais recentes
my.ticker <- c('AAPL')
first.date <- Sys.Date()-1500
last.date <- Sys.Date()
l.out <- BatchGetSymbols(tickers = my.ticker,first.date = first.date,last.date = last.date)
returns <- data.frame(retornos=diff(log(l.out$df.tickers$price.adjusted))*100,datas=l.out$df.tickers$ref.date[2:l.out$df.control$total.obs])



# 1 distribuição normal
# 2 distribuição t-student
# 3 distribuição normal generalizada

ARMAGARCH <- function(params,data,np,nq,nr,ns,nu) {
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
  
  #Distribuição Normal
  if(distribuicao==1){
    verossim <- 0.5*(sum(log(sigma2)) + sum((errors^2)/sigma2)  +  T*log(2*pi));
    return(list(LLF=(verossim),sigma2=sigma2,residuals=errors/sqrt(sigma2)))
    
    #Distribuição t-Student
  } else if (distribuicao == 2) {
    ratio1 <- ((gamma((nu+1)*0.5))/gamma((nu)*0.5))
    verossim <-  log(ratio1) - 0.5*log(pi*(nu-2)) - 0.5*log(sigma2) - 0.5*(nu+1)*(log((1+((errors)^(2))/(sigma2)*(nu-2))))
    return(list(LLF=-sum(verossim),sigma2=sigma2,residuals=errors/sqrt(sigma2)))
    
    #Distribuição normal generalizada
  }
  else if (distribuicao == 3) {
    beta = ((2^(-2/nu) * gamma(1/nu)/gamma(3/nu))^(0.5))
    verossim <- (log(nu)/(beta*(2^(1+1/nu))*gamma(1/nu))) - 0.5*log(sigma2)- 0.5*(abs(errors/(sqrt(sigma2)*beta))^nu);
    return(list(LLF=-sum(verossim),sigma2=sigma2,residuals=errors/sqrt(sigma2)))

  }
  else if (!(distribuicao<1 && distribuicao>3)) {stop ("Distribuição Inválida")}

  # return(list(LLF=-sum(verossim),sigma2=sigma2,residuals=errors/sqrt(sigma2)))

}



if(distribuicao==1){
  init.params <- c(0.1,0.1,0.1,0.1,0.5)
  ui <- rbind( c(0,0,1,0,0),c(0,0,0,1,0),c(0,0,0,0,1),c(0,0,0,-1,-1) )
  ci <- c(0,0,0,-1)
  ARMAGARCH.optim <- function(params,data,np,nq,nr,ns) { ARMAGARCH(params,data,np,nq,nr,ns)$LLF }
  resultados <- constrOptim(init.params,ARMAGARCH.optim,data=returns[,1],np=1,nq=1,nr=1,ns=1,grad=NULL,ui=ui,ci=ci)
  modelo_FIT <-garchFit(formula = ~ arma(1,1) + garch(1, 1), data = returns[,1], cond.dist = c("QMLE"), include.mean = F)

  
  } else if (distribuicao==2){
  init.params    <- c(0.1,0.1,0.1,0.1,0.5,3)
  ui  <- rbind( c(0,0,1,0,0,0),c(0,0,0,1,0,0),c(0,0,0,0,1,0),c(0,0,0,-1,-1,0),c(0,0,0,0,0,1) )
  ci   <- c(0,0,0,-1,2)
  ARMAGARCH.optim <- function(params,data,np,nq,nr,ns,nu) { ARMAGARCH(params,data,np,nq,nr,ns,nu)$LLF }
  resultados <- constrOptim(init.params,ARMAGARCH.optim,data=returns[,1],np=1,nq=1,nr=1,ns=1,nu=3,grad=NULL,ui=ui,ci=ci)
  modelo_FIT <-garchFit(formula = ~ arma(1,1) + garch(1, 1), data = returns[,1], cond.dist = c("std"),  include.mean = F)

  
  
  } else if(distribuicao==3){
    init.params <- c(0.1,0.1,0.1,0.1,0.5,2) 
    ui  <- rbind( c(0,0,1,0,0,0),c(0,0,0,1,0,0),c(0,0,0,0,1,0),c(0,0,0,-1,-1,0),c(0,0,0,0,0,1) )
    ci  <- c(0,0,0,-1,0)
    ARMAGARCH.optim <- function(params,data,np,nq,nr,ns,nu) { ARMAGARCH(params,data,np,nq,nr,ns,nu)$LLF }
    resultados <- constrOptim(init.params,ARMAGARCH.optim,data=returns[,1],np=1,nq=1,nr=1,ns=1,nu=3,grad=NULL,ui=ui,ci=ci)
    modelo_FIT <- garchFit(formula = ~ arma(1,1) + garch(1, 1), data = returns[,1], cond.dist = c("ged"),  include.mean = F)

}

# Compara com os parâmetros estimados através da biblitoeca fGarch do R
print(paste("ar1",resultados$par[1],"ma1",resultados$par[2],"omega",resultados$par[3],"alpha1",resultados$par[4],"beta1",resultados$par[5],"nu",resultados$par[6]))
print(modelo_FIT@fit$coef)

# Compara a verossimilhança obtida com código próprio e com pacote do R
LLF.meu.codigo <- print(-ARMAGARCH(resultados$par,data=returns[,1],np=1,nq=1,nr=1,ns=1, nu=3)$LLF)
LLF.fGarch <- print(-ARMAGARCH(c(0.0206345,-0.0062683,0.2994673,0.0862977,0.7732459),data=returns[,1],np=1,nq=1,nr=1,ns=1, nu=3)$LLF)

# Retorna variáveis e avalia qualidade do ajuste
resultados.finais <- ARMAGARCH(resultados$par,data=returns[,1],np=1,nq=1,nr=1,ns=1,nu=3)
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






