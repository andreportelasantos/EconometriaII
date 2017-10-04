# Argumentos da função:
# params: vetor de parâmetros 
# data: série temporal
# np: ordem AR
# nq: ordem MA
# nr: ordem ARCH
# ns: ordem GARCH
ARMAGARCH <- function(params,data,np,nq,nr,ns, distrib) {
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
  if (distrib == "normal"){
    verossim <- 0.5*(sum(log(sigma2)) + sum((errors^2)/sigma2)  +  (T-m)*log(2*pi));
    return(list(LLF=verossim,sigma2=sigma2,residuals=errors/sqrt(sigma2)))
  }
  if (distrib == "t-esp-gl"){
    v <- 4
    
    verossim <- ((v+1)/2)*log(1+(errors^2)/sigma2*(v-2) + 0.5*log(sigma2))
    return(list(LLF=sum(verossim),sigma2=sigma2, residuals=errors/sqrt(sigma2)))
  }
  if (distrib == "t"){
    v <- params[np+nq+nr+ns+2]
    verossim <- -(log(gamma((v+1)/2)/(sqrt((v-2)*pi)*gamma(v/2)))-log(sigma2)/2-((v+1)/2)*log(1+errors^2/(sigma2*(v-2))))
    return(list(LLF=sum(verossim),sigma2=sigma2,residuals=errors/sqrt(sigma2)))
  }
  if (distrib == "GED"){
    v <- 1.09005185
    lambda <- function(v){
      (2^(-2/v) * gamma(1/v) / gamma(3/v))^0.5
    }
    #verossim <- -log(v) + 0.5 * sum(abs(errors/(lambda(v)))^v) - log(lambda(v) * 2^(1+1/v) * gamma(1/v))
    verossim <- -(T-m)*log(v/lambda(v)) + 0.5*sum(abs(errors/(sqrt(sigma2)*lambda(v)))^v) + (T-m)*(1+1/v)*log(2)+(T-m)*log(gamma(1/v)) + 0.5*sum(log(sigma2))
    
    return(list(LLF=verossim,sigma2=sigma2, residuals=errors/sqrt(sigma2)))
  }
  
}

#************************************
# CARREGANDO OS DADOS
#************************************
# Baixa últimos 1500 dados mais recentes
library(BatchGetSymbols)
my.ticker <- c('AAPL')
first.date <- Sys.Date()-1500
last.date <- Sys.Date()
l.out <- BatchGetSymbols(tickers = my.ticker,first.date = first.date,last.date = last.date)
returns <- data.frame(retornos=diff(log(l.out$df.tickers$price.adjusted))*100,datas=l.out$df.tickers$ref.date[2:l.out$df.control$total.obs])

#*****************************************************************************************************
#Estima parâmetros do modelo ARMA(1,1)-GARCH(1,1) NORMAL
#*****************************************************************************************************
  
# Parâmetros iniciais da estimação
init.params <- c(0.1,0.1,0.1,0.1,0.5)
# Define restrições
ui <- rbind(c(0,0,1,0,0),c(0,0,0,1,0),c(0,0,0,0,1),c(0,0,0,-1,-1))
ci <- c(0,0,0,-1)
# Helper para maximizar a função verossimilhança condicional
ARMAGARCH.optim <- function(params,data,np,nq,nr,ns, distrib) {
  ARMAGARCH(params,data,np,nq,nr,ns, distrib)$LLF
}

# Exibe parâmetros estimados
resultados <- constrOptim(init.params,ARMAGARCH.optim,data=returns[,1],np=1,nq=1,nr=1,ns=1, distrib = "normal",grad=NULL,ui=ui,ci=ci)
print(resultados$par)

# Retorna variáveis e avalia qualidade do ajuste
resultados.finais1 <- ARMAGARCH(resultados$par,data=returns[,1],np=1,nq=1,nr=1,ns=1, distrib = "normal")
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

#*****************************************************************************************************
#Estima parâmetros do modelo ARMA(1,1)-GARCH(1,1) T-STUDENT (GL = 4)
#*****************************************************************************************************

# Parâmetros iniciais da estimação
init.params <- c(0.1,0.1,0.1,0.1,0.5)
# Define restrições
ui <- rbind(c(0,0,1,0,0),c(0,0,0,1,0),c(0,0,0,0,1),c(0,0,0,-1,-1))
ci <- c(0,0,0,-1)
# Helper para maximizar a função verossimilhança condicional
ARMAGARCH.optim <- function(params,data,np,nq,nr,ns, distrib) {
  ARMAGARCH(params,data,np,nq,nr,ns, distrib)$LLF
}

# Exibe parâmetros estimados
resultados <- constrOptim(init.params,ARMAGARCH.optim,data=returns[,1],np=1,nq=1,nr=1,ns=1, distrib = "t-esp-gl",grad=NULL,ui=ui,ci=ci)
print(resultados$par)

# Retorna variáveis e avalia qualidade do ajuste
resultados.finais2 <- ARMAGARCH(resultados$par,data=returns[,1],np=1,nq=1,nr=1,ns=1, distrib = "t-esp-gl")
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

#*****************************************************************************************************
#Estima parâmetros do modelo ARMA(1,1)-GARCH(1,1) T-STUDENT 
#*****************************************************************************************************

# Parâmetros iniciais da estimação
init.params <- c(0.1,0.1,0.1,0.1,0.5,10)
# Define restrições
ui <- rbind(c(0,0,1,0,0,0),c(0,0,0,1,0,0),c(0,0,0,0,1,0),c(0,0,0,-1,-1,0),c(0,0,0,0,0,1),c(0,0,0,0,0,-1))
ci <- c(0,0,0,-1,2+.Machine$double.eps,-30)
# Helper para maximizar a função verossimilhança condicional
ARMAGARCH.optim <- function(params,data,np,nq,nr,ns, distrib) {
  ARMAGARCH(params,data,np,nq,nr,ns, distrib)$LLF
}

# Exibe parâmetros estimados
resultados <- constrOptim(init.params,ARMAGARCH.optim,data=returns[,1],np=1,nq=1,nr=1,ns=1, distrib = "t",grad=NULL,ui=ui,ci=ci)
print(resultados$par)

# Retorna variáveis e avalia qualidade do ajuste
resultados.finais3 <- ARMAGARCH(resultados$par,data=returns[,1],np=1,nq=1,nr=1,ns=1, distrib = "t")
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


#*****************************************************************************************************
#Estima parâmetros do modelo ARMA(1,1)-GARCH(1,1) GED 
#*****************************************************************************************************

# Parâmetros iniciais da estimação
init.params <- c(0.1,0.1,0.1,0.1,0.5)
# Define restrições
ui <- rbind(c(0,0,1,0,0),c(0,0,0,1,0),c(0,0,0,0,1),c(0,0,0,-1,-1))
ci <- c(0,0,0,-1)
# Helper para maximizar a função verossimilhança condicional
ARMAGARCH.optim <- function(params,data,np,nq,nr,ns, distrib) {
  ARMAGARCH(params,data,np,nq,nr,ns, distrib)$LLF
}

# Exibe parâmetros estimados
resultados <- constrOptim(init.params,ARMAGARCH.optim,data=returns[,1],np=1,nq=1,nr=1,ns=1, distrib = "GED",grad=NULL,ui=ui,ci=ci)
print(resultados$par)

# Retorna variáveis e avalia qualidade do ajuste
resultados.finais4 <- ARMAGARCH(resultados$par,data=returns[,1],np=1,nq=1,nr=1,ns=1, distrib = "GED")
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

#_______________________________________________________





#**********************************************************************
#Compara com os parâmetros estimados através da biblitoeca fGarch do R
#**********************************************************************

# Estima modelo
#install.packages("fGarch")
library(fGarch)
resultados.finais2 <- garchFit(formula = ~ arma(1,1) + garch(1, 1), data = returns[,1], cond.dist = c("ged"), include.mean = F)

# Compara a verossimilhança obtida com código próprio e com pacote do R
LLF.meu.codigo <- print(-ARMAGARCH(resultados$par,data=returns[,1],np=1,nq=1,nr=1,ns=1)$LLF)
#LLF.fGarch <- print(-ARMAGARCH(c(0.0206345,-0.0062683,0.2994673,0.0862977,0.7732459),data=returns[,1],np=1,nq=1,nr=1,ns=1)$LLF)

