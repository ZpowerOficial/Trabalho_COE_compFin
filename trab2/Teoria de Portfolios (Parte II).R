# Carregando as bibliotecas necess?rias
library(zoo)
library(tseries)
library(PerformanceAnalytics)
library(tidyverse)
library(rportfolios)

# Definindo o per?odo de an?lise
START = "2013-08-01"
END = "2018-09-01"

# Escolhendo os c?digos das a??es
codigos = c("ABEV3", "BBAS3", "CIEL3")

# Criando uma estrutura de dados zoo para armazenar os pre?os de fechamento
close.zoo = zoo()
# modificando o ?ndice para armazenar datas no formato m?s/ano
index(close.zoo) = as.yearmon(index(close.zoo))

# Obtendo os pre?os de fechamento das a??es
for(codigo in codigos) {
  aux = get.hist.quote(
    instrument=paste(codigo, 'SA', sep = '.'), 
    start=START, end=END, 
    quote="AdjClose", provider="yahoo", compression="m", 
    retclass="zoo", quiet=TRUE)
  index(aux) = as.yearmon(index(aux)) 
  close.zoo = merge(close.zoo, aux)
}

# Nomeando as colunas com os c?digos das a??es
colnames(close.zoo) = codigos

# Transformando pre?os de fechamento em retornos mensais
close.df = data.frame(close.zoo)
ret.df = apply(log(close.df), 2, diff)
head(ret.df)

# Calculando estat?sticas amostrais
# Retornos m?dios
(mu.vec = apply(ret.df, 2, mean))
# Volatilidades
(sig.vec = apply(ret.df, 2, sd))

# Calculando a matriz de covari?ncias
(sigma.mat = cov(ret.df))

# GERANDO PORTF?LIOS ALEAT?RIOS 
NPORT = 1e4
x.mat = matrix(0, nrow = NPORT, ncol = length(codigos))
colnames(x.mat) = codigos
head(x.mat)

# Gerando pesos aleat?rios para os ativos
for (i in 1:NPORT) {
  # gera um n?mero entre 1 e 2, correspondente ao somat?rio das posi??es longas
  x.long = runif(1, 1, 2)
  x.mat[i, ] = random.longshort(
    n = length(codigos), 
    x.t.long = x.long, 
    x.t.short = x.long - 1)
}

# Calculando retorno e volatilidade para cada carteira aleat?ria
mu.p = x.mat %*% mu.vec
head(mu.p)

# Inicializa um vetor para armazenar as volatilidades
sig.p <- rep(0, NPORT)
head(sig.p)

# Loop for para calcular a volatilidade para cada linha da matriz
for (i in 1:NPORT) {
  # Calcula a volatilidade para a i-?sima linha
  sig.p[i] <- sqrt(t(x.mat[i,]) %*% sigma.mat %*% x.mat[i,])
}

# Plotando os portf?lios aleat?rios
plot(sig.p, mu.p, pch = 16, col = 'black',
     xlim = c(0, max(sig.p)),
     xlab = expression(sigma[p]), ylab = expression(mu[p]))

# Marcando os ativos A, B e C no gr?fico
points(sig.vec, mu.vec, pch=16, cex = 1.5, col='red')
text(sig.vec, mu.vec, codigos, pos = 4, cex = 0.7, col = 'red')

# Determinando os pesos no portf?lio de vari?ncia m?nima
# Teoria de Portf?lios (slide 79)
sigma.inv = solve(sigma.mat)
m.vec = rowSums(sigma.inv)/sum(sigma.inv)

# Adicionando o ponto do portf?lio de vari?ncia m?nima no gr?fico
(mu.gmvp = t(m.vec) %*% mu.vec)
(sig.gmvp = sqrt(t(m.vec) %*% sigma.mat %*% m.vec))
points(sig.gmvp, mu.gmvp, pch = 16, col = 'magenta', cex = 2.0)
text(sig.gmvp, mu.gmvp, 'GMVP', pos = 2, cex = 0.7)


############################################################
# DEMO 1
############################################################
#' Definindo uma fun??o para os pesos e as caracter?sticas 
#' de um portf?lio na fronteira, dado o retorno alvo
#' slide 5
portfolio_fronteira <- function(mu_alvo, mu_vec, sigma_mat) {
  # Criando vetores de zeros e uns
  DIM = nrow(sigma_mat)
  zero_vec = rep(0, DIM)
  one_vec = rep(1, DIM)
  
  # Construindo a matriz A
  A = matrix(0, nrow = DIM+2, ncol = DIM+2)
  A[1:DIM, 1:DIM] = 2*sigma_mat
  A[DIM+1, 1:DIM] = mu.vec
  A[DIM+2, 1:DIM] = one_vec
  A[1:DIM, DIM+1] = mu.vec
  A[1:DIM, DIM+2] = one_vec

  # Construindo o vetor b
  b = c(zero_vec, mu_alvo, 1)
  
  # Calculando a inversa da matriz A
  A_inv = solve(A)
  
  # Calculando o vetor z
  z = A_inv %*% b
  
  # Extraindo o vetor x (removendo as duas ?ltimas linhas)
  x_vec = head(z, -2)
  
  # Calculando as caracter?sticas do portf?lio
  mu_p = t(x_vec) %*% mu_vec
  sig_p = sqrt(t(x_vec) %*% sigma_mat %*% x_vec)
  
  # Retornando as caracter?sticas do portf?lio
  return(list(x_vec = x_vec, mu_p = mu_p, sig_p = sig_p))
}

# Calculando o portfolio na fronteira com o mesmo retorno de BBAS3
portfolio_BBAS3 = portfolio_fronteira(mu.vec['BBAS3'], mu.vec, sigma.mat)

# Calculando o portfolio na fronteira com o mesmo retorno de CIEL3
portfolio_CIEL3 = portfolio_fronteira(mu.vec['CIEL3'], mu.vec, sigma.mat)

# Obtendo covari?ncia e correla??o entre os portf?lios (vai ser ?til mais tarde)
sig_xy = t(portfolio_BBAS3$x_vec) %*% sigma.mat %*% portfolio_CIEL3$x_vec
rho_xy = sig_xy / portfolio_BBAS3$sig_p / portfolio_CIEL3$sig_p

# Plotando os pontos no gr?fico
points(
  x = c(portfolio_BBAS3$sig_p, portfolio_CIEL3$sig_p),
  y = c(portfolio_BBAS3$mu_p, portfolio_CIEL3$mu_p),
  pch = 16, col = 'green', cex = 1.5
)

# Adicionando r?tulos aos pontos no gr?fico
text(
  x = c(portfolio_BBAS3$sig_p, portfolio_CIEL3$sig_p),
  y = c(portfolio_BBAS3$mu_p, portfolio_CIEL3$mu_p),
  labels = c('BBAS3 Portfolio', 'CIEL3 Portfolio'), pos = 2, cex = 0.7
)

# Imprimindo a correla??o entre os portf?lios
cat("Correla??o entre BBAS3 Portfolio e CIEL3 Portfolio:", rho_xy, "\n")


############################################################
# DEMO 2
############################################################
# a fronteira da bala de Markowitz como combina??o convexa
# de dois outros portf?lios sobre a fronteira
p.x = portfolio_BBAS3
p.y = portfolio_CIEL3

# Definindo a propor??o alfa para a combina??o convexa
alfa = 0.5

# Calculando o vetor de peso para o portf?lio Z como combina??o convexa de X e Y
z.vec = alfa * p.x$x_vec + (1 - alfa) * p.y$x_vec

# Calculando o retorno esperado do portf?lio Z
mu.p.z = t(z.vec) %*% mu.vec

# Calculando a vari?ncia e o desvio padr?o do portf?lio Z
sig2.p.z = t(z.vec) %*% sigma.mat %*% z.vec
sig.p.z = sqrt(sig2.p.z)

# Tra?ando o gr?fico
plot(0, xlim = c(0, 0.2), ylim=c(0, 0.012), type = 'n',  
     xlab = expression(sigma[p]), ylab = expression(mu[p]))

# Marcando os ativos A, B e C no gr?fico
points(sig.vec, mu.vec, pch=16, cex = 1.5, col='red')
text(sig.vec, mu.vec, codigos, pos = 4, cex = 0.7)

# Marcando o portf?lio de vari?ncia m?nima (GMVP) no gr?fico
points(sig.gmvp, mu.gmvp, pch = 16, col = 'green', cex = 1.5)
text(sig.gmvp, mu.gmvp, 'GMVP', pos = 2, cex = 0.7)

# Marcando os portf?lios X e Y no gr?fico
points(x = c(p.x$sig_p, p.y$sig_p), y = c(p.x$mu_p, p.y$mu_p),
       pch = 16, col = 'green', cex = 1.5)
text(x = c(p.x$sig_p, p.y$sig_p), y = c(p.x$mu_p, p.y$mu_p),
     labels = c('X', 'Y'), pos = 2, cex = 0.7)

# Marcando o portf?lio Z no gr?fico
points(sig.p.z, mu.p.z, pch = 16, col = 'magenta', cex = 1.5)
text(sig.p.z, mu.p.z, 'E1', pos = 2, cex = 0.7)


############################################################
# DEMO 3
############################################################
# Definindo o objetivo de retorno esperado para o portf?lio Z
mu.alvo = 0.012

# Calculando a propor??o alfa para a combina??o convexa que atinge o retorno desejado
alfa = as.numeric((mu.alvo - p.y$mu_p)/(p.x$mu_p - p.y$mu_p))

# Calculando o vetor de peso para o portf?lio Z
z.vec = alfa * x.vec + (1 - alfa) * y.vec

# Calculando as caracter?sticas do portf?lio Z para incluir no gr?fico
(mu.p.z = t(z.vec) %*% mu.vec)
sig.p.z = sqrt(t(z.vec) %*% sigma.mat %*% z.vec)

# Plotando o ponto do portf?lio Z no gr?fico
points(sig.p.z, mu.p.z, pch = 16, col = 'green', cex = 1.5)
text(sig.p.z, mu.p.z, 'E2', pos = 2, cex = 0.7)


############################################################
# DEMO 4
############################################################
# determina a fronteira eficiente como uma combina??o convexa 
# do portfolio de vari?ncia m?nima global (GMVP) e do portf?lio eficiente
# com retorno igual a max(mu.vec)
#
# PASSO 1: determina o GMVP
# uma f?rmula f?cil de guardar...
sigma.inv = solve(sigma.mat)
(m.vec = rowSums(sigma.inv)/sum(sigma.inv))
# determina as caracter?sticas do GMVP
(mu.p.m = t(m.vec) %*% mu.vec)
sig2.p.m = t(m.vec) %*% sigma.mat %*% m.vec
(sig.p.m = sqrt(sig2.p.m))

# PASSO 2: encontre o ativo com o maior retorno esperado
(mu.0 = max(mu.vec))
# usa a fun??o que escrevemos anteriormente para determinar
# os pesos e caracter?sticas do portfolio com retorno igual a mu.0
X = portfolio_fronteira(mu.0, mu.vec, sigma.mat)

# covariancia entre os portfolios
(sig.m.x = t(m.vec) %*% sigma.mat %*% X$x_vec)

#
# PASSO 3
# determina a fronteira eficiente como uma combina??o convexa dos dois portfolios
alfa.vec = seq(from = 2, to = -2, by = -0.1)
(n.pontos = length(alfa.vec))

# aloca espa?o para os par?metros dos portf?lios na fronteira
mu.p.fronteira = NULL
sig.p.fronteira = NULL
for (i in 1:n.pontos) {
  alfa = alfa.vec[i]
  aux = alfa*m.vec + (1 - alfa)*X$x_vec
  mu.p.fronteira[i] = t(aux) %*% mu.vec
  sig.p.fronteira[i] = sqrt(t(aux) %*% sigma.mat %*% aux)
}
#
# PASSO 4
# tra?a o gr?fico
plot(0, xlim=c(0, 0.22), ylim=c(-0.0061, 0.023), type = 'n', 
     xlab = expression(sigma[p]), ylab = expression(mu[p]))
# marca os ativos A, B e C no gr?fico
points(sig.vec, mu.vec, pch=16, cex = 1.5, col='red')
text(sig.vec, mu.vec, codigos, pos = 4, cex = 0.7)
# tra?a a fronteira
points(sig.p.fronteira, mu.p.fronteira, pch = 16, col = 'green')


############################################################
# DEMO 5
############################################################
# calcula o portf?lio tangente
# slide 23
(DIM = length(codigos))

# rendimento m?dio mensal do CDI no per?odo
rf = 0.8896/100

# inverso da matriz sigma
sigma.inv = solve(sigma.mat)

# vetor auxiliar
one.vec = rep(1, DIM)

# calculando o vetor tangente usando a f?rmula no slide 23
(t.vec = sigma.inv %*% (mu.vec - rf * one.vec) / as.numeric(
  t(one.vec) %*% sigma.inv %*% (mu.vec - rf * one.vec))
)

# calcula a m?dia e o desvio padr?o do portf?lio tangente
(mu.p.t = t(t.vec) %*% mu.vec)
sig2.p.t = t(t.vec) %*% sigma.mat %*% t.vec
(sig.p.t = sqrt(sig2.p.t))

# tra?a o portf?lio tangente no gr?fico
points(sig.p.t, mu.p.t, pch = 16, col = 'magenta', cex=1.2)
text(sig.p.t, mu.p.t, 'tangente', pos = 4, cex = 0.7)

#' Tra?a a reta dos investimentos eficientes (
#' reta que passa pelo ativo livre de risco e pelo 
#' portf?lio tangente)
abline(a = rf, b = (mu.p.t - rf)/sig.p.t, col = 'blue')
points(0, rf, pch = 16, col = 'magenta')
text(0, rf, 'rf', pos = 4, cex = 0.7)

#' como o coeficiente angular da reta ? negativo, inverte o coeficiente angular
#' Interpreta??o: Os investimentos eficientes envolvem encurtar o 
#' portf?lio tangente para investir no ativo livre de risco
abline(a = rf, b = (rf - mu.p.t)/sig.p.t, lty=2)

############################################################
# DEMO 6
############################################################
# Usa as f?rmulas no slide 27
# Definindo o risco alvo (SD) desejado
sig.alvo = 0.02

# Calculando o peso do portf?lio tangente
# Nota: ajustando o sinal para compensar o ?ndice de Sharpe negativo
(x.t = -sig.alvo / sig.p.t)

# Calculando o peso do ativo livre de risco
(x.f = 1 - x.t)

# Calculando o retorno esperado do investimento eficiente (slide 27)
(mu.p.e = rf + x.t * (mu.p.t - rf))

# Calculando o risco do portf?lio eficiente
# Nota: ajustando o sinal para compensar o ?ndice de Sharpe negativo
(sig.p.e = -x.t * sig.p.t)

# Tra?ando o ponto no gr?fico
points(sig.p.e, mu.p.e, pch = 16)

# Adicionando linhas de orienta??o ao ponto no gr?fico
segments(sig.p.e, -0.01, sig.p.e, mu.p.e, lty = 3)
axis(side = 1, at = sig.p.e)
segments(-0.01, mu.p.e, sig.p.e, mu.p.e, lty = 3)
axis(side = 2, at = sprintf(mu.p.e, fmt='%.3f'))


############################################################
# DEMO 7
############################################################
# Use as f?rmulas no slide 27
# Encontre o portf?lio eficiente com retorno alvo = 0.02
mu.alvo = 0.02

# Calculando o peso do portf?lio tangente usando as f?rmulas do slide 27
(x.t = (mu.alvo - rf)/(mu.p.t - rf))

# peso do ativo de renda fixa
(x.f = 1 - x.t)

# retorno esperado do portf?lio eficiente
(mu.p.e = rf + x.t*(mu.p.t - rf))

# Calculando o risco do investimento eficiente
# Ajustando o sinal para compensar o ?ndice de Sharpe negativo
(sig.p.e = -x.t * sig.p.t)

# tra?a o ponto no gr?fico
points(sig.p.e, mu.p.e, pch=16)

# Adicionando linhas de orienta??o ao ponto no gr?fico
segments(sig.p.e, -0.01, sig.p.e, mu.p.e, lty = 3)
segments(-0.01, mu.p.e, sig.p.e, mu.p.e, lty = 3)

# Se eu tenho $100,000.00 para investir, quanto vou alocar em cada ativo?
W0 = 100000

# Alocando no ativo de renda fixa
allocation_rf = W0 * x.f
cat("Aloca??o em ativo de renda fixa: $", round(allocation_rf, 2), "\n")

# Alocando em ativos de risco
allocation_risk = W0 * as.numeric(x.t) * t.vec
cat("Aloca??o em ativos de risco:\n")
print(allocation_risk)


############################################################
# DEMO 8
############################################################
# Determine o investimento eficiente que tem o mesmo retorno esperado de BBAS3
# slide 27
(x.t = (mu.vec['BBAS3'] - rf)/(mu.p.t - rf))
(x.f = 1 - x.t)

# caracter?sticas do portf?lio eficiente
mu.p.e = rf + x.t*(mu.p.t - rf)
cat('Retorno do portf?lio eficiente:', mu.p.e, '\n')

# volatilidade do portf?lio eficiente
# ajustando para o ?ndice de Sharpe negativo do portf?lio tangente
sig.p.e = -x.t * sig.p.t
cat('Volatilidade do portf?lio eficiente:', sig.p.e, '\n')

# tra?a o ponto no gr?fico
points(sig.p.e, mu.p.e, pch=16, col='magenta')

# Compare o VaR.05 para BBAS3 e para o portf?lio eficiente 
# VaR a 5% para o Banco do Brasil
W0 = 100000   # capital investido
q.BBAS.05 = mu.vec['BBAS3'] + sig.vec['BBAS3']*qnorm(0.05)
cat('quantil 5% para BBAS3:', q.BBAS.05, '\n')
VaR.BBAS.05 = W0 * q.BBAS.05
cat('Valor em Risco a 5% para BBAS3:', VaR.BBAS.05, '\n')

# VaR a 5% para o portf?lio eficiente
q.e.05 = mu.p.e + sig.p.e * qnorm(0.0125)
cat('quantil 5% para o portf?lio eficiente:', q.e.05, '\n')
VaR.e.05 = W0 * q.e.05
cat('Valor em Risco a 5% para o portf?lio eficiente:', VaR.e.05, '\n')


############################################################
# DEMO 9
############################################################
# Define o diret?rio de trabalho e carrega as fun??es do arquivo "portfolio.r"
getwd()
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
dir()
source("portfolio.r")

# Cria um portfolio com pesos iguais
pesos_iguais_port = getPortfolio(er = mu.vec, cov.mat = sigma.mat, weights = rep(1/3, 3))
plot(pesos_iguais_port, main = 'Portf?lio com Pesos Iguais', col = 'blue')

# Obt?m o GMVP (Global Minimum Variance Portfolio) para um portf?lio
# formado por ABEV3, BBAS3 e CIEL3
gmvp_port = globalMin.portfolio(er = mu.vec, cov.mat = sigma.mat)
plot(gmvp_port, col = 'blue', main = 'Portf?lio de Vari?ncia M?nima Global')

#' Encontra um portf?lio eficiente, isto ?, um portf?lio sobre a  
#' fronteira da bala de Markowitz, dado um retorno alvo (o retorno esperado de BBAS3)
ef_port = efficient.portfolio(er = mu.vec, cov.mat = sigma.mat, target.return = mu.vec['BBAS3'])
plot(ef_port, col = 'blue', main = 'Portf?lio Eficiente com Retorno Alvo')

# Determina o portf?lio tangente
# Taxa livre de risco
# baixando os juros na marra. :))))))))
rf = 0.17e-2
tan_port = tangency.portfolio(er = mu.vec, cov.mat = sigma.mat, risk.free = rf)
plot(tan_port, col = 'blue', main = 'Portf?lio Tangente')

# Determina a fronteira eficiente
ef_fronteira = efficient.frontier(
  er = mu.vec, cov.mat = sigma.mat, nport = 30, alpha.min = -2, alpha.max = 2
)
# caracter?sticas (ER e SD) dos portf?lios na fronteira
ef_fronteira
# pesos dos ativos nos portf?lios na fronteira
ef_fronteira$weights
# somente os retornos dos portf?lios na fronteira
ef_fronteira$er
# somente o desvio padr?o dos portf?lios na fronteira
ef_fronteira$sd
# plota a fronteira eficiente
plot(ef_fronteira, main = 'Fronteira Eficiente')


############################################################
# DEMO 10
############################################################
# ajusta o nome das vari?veis para aqueles usados pelo GPT
efficient_frontier = ef_fronteira
gmvp_portfolio = gmvp_port
tan_portfolio = tan_port
risk_free_rate = rf

# Plota a fronteira eficiente
plot(efficient_frontier, plot.assets = TRUE, col = 'blue', pch = 16)

# Marca o GMVP (Global Minimum Variance Portfolio) no gr?fico
points(gmvp_portfolio$sd, gmvp_portfolio$er, col = 'green', pch = 16, cex = 2)
text(gmvp_portfolio$sd, gmvp_portfolio$er, labels = 'GMVP', pos = 2)

# Marca o portf?lio tangente no gr?fico
points(tan_portfolio$sd, tan_portfolio$er, col = 'red', pch = 16, cex = 2)
text(tan_portfolio$sd, tan_portfolio$er, labels = 'TANGENTE', pos = 2)

# Tra?a a reta dos investimentos eficientes
sharpe_ratio_tan = (tan_portfolio$er - risk_free_rate) / tan_portfolio$sd
abline(a = risk_free_rate, b = sharpe_ratio_tan, col = 'green', lwd = 2)


############################################################
# DEMO 10a
############################################################
# Determina a fronteira eficiente sem vendas a descoberto
library(quadprog) # usada por portfolio.r
efficient_frontier_no_shorts = efficient.frontier(
  er = mu.vec, cov.mat = sigma.mat, shorts = FALSE,
  nport = 30
)

# Adiciona os pontos da fronteira eficiente sem vendas a descoberto ao gr?fico existente
points(
  efficient_frontier_no_shorts$sd, efficient_frontier_no_shorts$er, pch = 16
)

# Carrega a biblioteca 'zoom'
library(zoom)

# Chama a fun??o zoom() para ampliar o gr?fico
zm()


############################################################
# DEMO 11
############################################################
# Determina o GMVP sem vendas a descoberto usando o Excel
# exporta dados para o excel
# copia mu.vec para o clipboard
write.table(
  mu.vec, file = 'clipboard', sep = '\t', 
  row.names = TRUE, col.names = F, dec = ',')
# copia sigma.mat para o clipboard
write.table(sigma.mat, file = 'clipboard', sep = '\t', dec = ',')


############################################################
# DEMO 12
############################################################
# usando as fun??es do Zivot
# GMVP permitindo vendas a descoberto
gmvp_portfolio_shorts_allowed = globalMin.portfolio(er = mu.vec, cov.mat = sigma.mat)
gmvp_portfolio_shorts_allowed$weights
#
# GMVP sem vendas a descoberto
gmvp_portfolio_no_shorts = globalMin.portfolio(er = mu.vec, cov.mat = sigma.mat, shorts = FALSE)
gmvp_portfolio_no_shorts$weights
#
# Usando solve.QP() para calcular o GMVP (Global Minimum Variance Portfolio)
# F?rmulas no Slide 58

# Carregando a biblioteca quadprog
library(quadprog)

# Obtendo o n?mero de ativos
n <- length(codigos)

# Criando vetores ?teis
one.vec <- rep(1, n)
zero.vec <- rep(0, n)

# Criando a matriz A e o vetor b para as restri??es (equa??es e inequa??es)
A <- cbind(one.vec, diag(n))
b <- c(1, zero.vec)

# Resolvendo o problema de otimiza??o quadr?tica
gmvp_no_shorts <- solve.QP(
  Dmat = 2 * sigma.mat, dvec = zero.vec,
  Amat = A, bvec = b, meq = 1
)

# Obtendo os pesos dos ativos no portf?lio
weights_gmvp_no_shorts <- gmvp_no_shorts$solution
names(weights_gmvp_no_shorts) <- codigos
weights_gmvp_no_shorts

# Calculando vari?ncia e volatilidade do portf?lio
variance_gmvp_no_shorts <- gmvp_no_shorts$value
volatility_gmvp_no_shorts <- sqrt(variance_gmvp_no_shorts)

# Calculando retorno esperado
(expected_return_gmvp_no_shorts <- t(weights_gmvp_no_shorts) %*% mu.vec)

# Adicionando o ponto do GMVP no gr?fico existente
points(volatility_gmvp_no_shorts, expected_return_gmvp_no_shorts, pch = 16, col = 'magenta', cex = 2)


############################################################
# DEMO 12a
############################################################
#' GMVP (Global Minimum Variance Portfolio) sem vendas a descoberto, 
#' com no m?ximo 50% alocado em cada ativo
#' usando solve.QP()
#' F?rmulas no Slide 60 sem a restri??o do retorno alvo

# Carregando a biblioteca quadprog
library(quadprog)

# Obtendo o n?mero de ativos
n <- length(codigos)

# Criando vetores ?teis
one.vec <- rep(1, n)
zero.vec <- rep(0, n)

# Criando a matriz A e o vetor b para as restri??es
A <- cbind(one.vec, diag(n), -diag(n))
M <- 0.5
b <- c(1, zero.vec, -M * one.vec)

# Resolvendo o problema de otimiza??o quadr?tica com restri??es
gmvp_no_shorts_max_weight <- solve.QP(
  Dmat = 2 * sigma.mat, dvec = zero.vec, Amat = A, bvec = b, meq = 1
)

# Obtendo os resultados do GMVP com restri??es de peso
weights_gmvp_no_shorts_max_weight <- gmvp_no_shorts_max_weight$solution
variance_gmvp_no_shorts_max_weight <- gmvp_no_shorts_max_weight$value

# Exibindo os resultados
gmvp_no_shorts_max_weight
weights_gmvp_no_shorts_max_weight
variance_gmvp_no_shorts_max_weight

# Adicionando o ponto do GMVP no gr?fico existente
points(sqrt(variance_gmvp_no_shorts_max_weight), t(weights_gmvp_no_shorts_max_weight) %*% mu.vec,
       pch = 16, col = 'red', cex = 2)


############################################################
# DEMO 13
############################################################
# A fronteira eficiente
# Usando as fun??es do Zivot
source("portfolio.r")

# Calcula a fronteira eficiente com vendas a descoberto
efficient_frontier_shorts_allowed <- efficient.frontier(
  er = mu.vec, cov.mat = sigma.mat,
  alpha.max = 1, alpha.min = -1, nport = 20,
  shorts = TRUE
)

# Obt?m os pesos dos ativos na fronteira com vendas a descoberto
weights_efficient_frontier_shorts_allowed <- efficient_frontier_shorts_allowed$weights
weights_efficient_frontier_shorts_allowed

# Calcula a fronteira eficiente sem vendas a descoberto
efficient_frontier_no_shorts <- efficient.frontier(
  er = mu.vec, cov.mat =  sigma.mat,
  alpha.max = 1, alpha.min = -1, nport = 20,
  shorts = FALSE
)

# Obt?m os pesos dos ativos na fronteira sem vendas a descoberto
weights_efficient_frontier_no_shorts <- efficient_frontier_no_shorts$weights
weights_efficient_frontier_no_shorts

# Tra?a os gr?ficos
plot(efficient_frontier_shorts_allowed, pch = 16, col = 'blue')
points(efficient_frontier_no_shorts$sd, efficient_frontier_no_shorts$er, pch = 16)



# Ou, usando diretamente a solve.QP()
# GMVP (Global Minimum Variance Portfolio) sem vendas a descoberto
# Slide 57

# Obtendo o n?mero de ativos
n_ativos <- length(codigos)

# Criando vetores ?teis
zero_vec <- rep(0, n_ativos)
one_vec <- rep(1, n_ativos)
D <- 2 * sigma.mat
d <- zero_vec

# Calcula o GMVP (Global Minimum Variance Portfolio)
A <- cbind(one_vec, diag(n_ativos))
b <- c(1, zero_vec)
gmvp_no_shorts <- solve.QP(Dmat = D, dvec = d, Amat = A, bvec = b, meq = 1)

# Obt?m os pesos dos ativos no GMVP
weights_gmvp_no_shorts <- gmvp_no_shorts$solution

# Determina a fronteira da bala de Markowitz sem vendas a descoberto
mu_alvo_min <- t(weights_gmvp_no_shorts) %*% mu.vec
mu_alvo_max <- max(mu.vec)
mu_vals <- seq(as.numeric(mu_alvo_min), as.numeric(mu_alvo_max), length.out = 20)
n_portfolios <- length(mu_vals)

# Prepara os vetores de sa?da
w_mat <- matrix(0, n_portfolios, n_ativos)
sd_vals <- rep(0, n_portfolios)

# Resolve um problema de converg?ncia do algoritmo
# quando mu_alvo = max(mu.vec)
w_mat[n_portfolios, which.max(mu.vec)] <- 1
sd_vals[n_portfolios] <- sig.vec[which.max(mu.vec)]

# Resolve v?rias vezes o algoritmo de Markowitz
A <- cbind(mu.vec, one_vec, diag(n_ativos))
for (i in 1:(n_portfolios - 1)) {
  mu_alvo <- mu_vals[i]
  b <- c(mu_alvo, 1, zero_vec)
  qp_out <- solve.QP(Dmat = D, dvec = d, Amat = A, bvec = b, meq = 2)
  w_mat[i, ] <- qp_out$solution
  sd_vals[i] <- sqrt(qp_out$value)
}

# Tra?a os gr?ficos
plot(0, xlim = range(efficient_frontier_shorts_allowed$sd, sd_vals),
     ylim = range(efficient_frontier_shorts_allowed$er, mu_vals),
     type = 'n', xlab = expression(sigma[p]), ylab = expression(mu[p]))

points(efficient_frontier_shorts_allowed$sd, efficient_frontier_shorts_allowed$er, pch = 16, col = 'blue')
points(sd_vals, mu_vals, pch = 16)


############################################################
# DEMO 14
############################################################
# Exemplo: Um portf?lio n?o realiz?vel sem vendas a descoberto
# Tentando um retorno alvo de 0.01
# Slide 57

# Obtendo o n?mero de ativos
n_ativos <- length(codigos)

# Criando vetores ?teis
D <- 2 * sigma.mat
d <- zero.vec

# Calculando o portf?lio n?o realiz?vel
A <- cbind(mu.vec, one.vec, diag(n_ativos))
mu_alvo = 0.01
b <- c(mu_alvo, 1, zero.vec)

# Resolvendo o problema de otimiza??o quadr?tica
qp_out_non_feasible <- solve.QP(Dmat = D, dvec = d, Amat = A, bvec = b, meq = 2)


############################################################
# DEMO 14a
############################################################
# Impondo pesos m?ximos aos ativos, por exemplo, 50%
# Calcula o GMVP (Global Minimum Variance Portfolio)
# Slide 60 sem a primeira coluna da matriz A

# Definindo o limite m?ximo para os pesos dos ativos
M <- 0.5

# Obtendo o n?mero de ativos
n_ativos <- length(mu.vec)

# Criando vetores ?teis
D <- 2 * sigma.mat
d <- zero.vec

# Calculando o GMVP com o limite de participa??o dos ativos
A <- cbind(one.vec, diag(n_ativos), -diag(n_ativos))
b <- c(1, zero.vec, -M * one.vec)
gmvp <- solve.QP(Dmat = D, dvec = d, Amat = A, bvec = b, meq = 1)

# Obtendo os pesos dos ativos no GMVP
(weights_gmvp <- gmvp$solution)

# Calculando os pontos sobre a fronteira
mu_alvo_min <- t(weights_gmvp) %*% mu.vec
# Encontrando empiricamente um alvo m?ximo realiz?vel
mu_alvo_max <- 0.0077
mu_vals <- seq(as.numeric(mu_alvo_min), as.numeric(mu_alvo_max), length.out = 20)

# Alocando espa?o para as sa?das
sig_vals <- rep(0, 20)
w_mat <- matrix(0, 20, n_ativos)

# La?o para encontrar os pontos sobre a fronteira (Slide 60)
A <- cbind(mu.vec, one.vec, diag(n_ativos), -diag(n_ativos))
for (i in 1:20) {
  b <- c(mu_vals[i], 1, zero.vec, -M * one.vec)
  qp_out <- solve.QP(Dmat = D, dvec = d, Amat = A, bvec = b, meq = 2)
  w_mat[i, ] <- qp_out$solution
  sig_vals[i] <- sqrt(qp_out$value)
}

# Plota o gr?fico
points(sig_vals, mu_vals, pch = 16, col = 'red')

# Zoom no gr?fico
library(zoom)
zm()


############################################################
# DEMO 15
############################################################
# Elipses de confian?a 95% para mu e sigma
library(ellipse)
n.obs = nrow(ret.df)

# SE te?rico do retorno esperado
se.mu.vec = sig.vec/sqrt(n.obs)
rbind(mu.vec, se.mu.vec)

# SE te?rico da volatilidade
se.sig.vec = sig.vec/sqrt(2*n.obs)
rbind(sig.vec, se.sig.vec)

# determina e plota as elipses
# ABEV3
ellipse.ABEV = ellipse::ellipse(
  x = 0, # correla??o entre mu e sigma
  level=0.95, # intervalo de confian?a
  scale = c(se.sig.vec['ABEV3'], se.mu.vec['ABEV3']), # desvios padr?o de sigma e mu
  centre = c(sig.vec['ABEV3'], mu.vec['ABEV3'])) # o centro da elipse
# BBAS3
ellipse.BBAS = ellipse::ellipse(
  x = 0, 
  scale = c(se.sig.vec['BBAS3'], se.mu.vec['BBAS3']),
  centre = c(sig.vec['BBAS3'], mu.vec['BBAS3']))
# CIEL3
ellipse.CIEL = ellipse::ellipse(
  x = 0, 
  scale = c(se.sig.vec['CIEL3'], se.mu.vec['CIEL3']),
  centre = c(sig.vec['CIEL3'], mu.vec['CIEL3']))
# tra?a os gr?ficos
plot(sig.vec, mu.vec, pch=16, 
     xlim = range(ellipse.ABEV[,'x'], ellipse.BBAS[,'x'], ellipse.CIEL[,'x']),
     ylim = range(ellipse.ABEV[,'y'], ellipse.BBAS[,'y'], ellipse.CIEL[,'y']),
     ylab=expression(mu[p]), xlab=expression(sigma[p]))
text(sig.vec, mu.vec, labels = codigos, pos = 4)
points(ellipse.ABEV, type = 'l')
points(ellipse.BBAS, type = 'l')
points(ellipse.CIEL, type = 'l')

############################################################
# DEMO 16
############################################################
# Bootstrapping o GMVP
gmvp = globalMin.portfolio(er = mu.vec, cov.mat = sigma.mat)
gmvp
plot(gmvp, col = 'slateblue1')
#
# aloca espa?o para os resultados
n.boot = 1e4
mu.gmvp.boot = rep(0, n.boot)
sd.gmvp.boot = rep(0, n.boot)
w.gmvp.boot = matrix(0, n.boot, n.ativos)
colnames(w.gmvp.boot) = codigos
head(w.gmvp.boot)
#
# bootstrap
for (i in 1:n.boot) {
  boot.idx = sample(n.obs, replace=TRUE)
  ret.df.boot = ret.df[boot.idx, ] 
  mu.vec.boot = apply(ret.df.boot, 2, mean)
  sigma.mat.boot = cov(ret.df.boot) 
  gmvp = globalMin.portfolio(mu.vec.boot, sigma.mat.boot)
  mu.gmvp.boot[i] = gmvp$er
  sd.gmvp.boot[i] = gmvp$sd
  w.gmvp.boot[i, ] = gmvp$weights
}
# plota os resultados
# gerando o slide da apresenta??o
#
# valor esperado do retorno do portf?lio
layout(mat = matrix(1:6, 2, 3, byrow = T))
hist(mu.gmvp.boot, freq = F, col = 'slateblue1', main = 'er', xlab = 'Valor')
points(density(x = mu.gmvp.boot), type = 'l')

# desvio padr?o do portf?lio
hist(sd.gmvp.boot, freq = F, col = 'slateblue1', main = 'sd', xlab = 'Valor')
points(density(x = sd.gmvp.boot), type = 'l')

# peso de ABEV3 no GMVP
hist(w.gmvp.boot[, 'ABEV3'], freq = F, col = 'slateblue1', main = 'ABEV3', xlab = 'Valor')
points(density(x = w.gmvp.boot[, 'ABEV3']), type = 'l')

# peso de BBAS3 no GMVP
hist(w.gmvp.boot[, 'BBAS3'], freq = F, col = 'slateblue1', main = 'BBAS', xlab = 'Valor')
points(density(x = w.gmvp.boot[, 'BBAS3']), type = 'l')

# peso de CIEL3 no GMVP
hist(w.gmvp.boot[, 'CIEL3'], freq = F, col = 'slateblue1', main = 'CIEL', xlab = 'Valor')
points(density(x = w.gmvp.boot[, 'CIEL3']), type = 'l')

layout(1)


############################################################
# DEMO 16a
############################################################
# summary statistics
out = matrix(0, 5, 3)
rownames(out) = c('er', 'sd', 'ABEV3', 'BBAS3', 'CIEL3')
colnames(out) = c('Observado', 'M?dia', 'SE')
out['er',] = c(gmvp.port$er, mean(mu.gmvp.boot), sd(mu.gmvp.boot))
out['sd',] = c(gmvp.port$sd, mean(sd.gmvp.boot), sd(sd.gmvp.boot))
out['ABEV3',] = c(gmvp.port$weights['ABEV3'], mean(w.gmvp.boot[, 'ABEV3']), sd(w.gmvp.boot[, 'ABEV3']))
out['BBAS3',] = c(gmvp.port$weights['BBAS3'], mean(w.gmvp.boot[, 'BBAS3']), sd(w.gmvp.boot[, 'BBAS3']))
out['CIEL3',] = c(gmvp.port$weights['CIEL3'], mean(w.gmvp.boot[, 'CIEL3']), sd(w.gmvp.boot[, 'CIEL3']))
out


############################################################
# DEMO 16b
############################################################
# graficamente: a incerteza do GMVP
# para n?o cortar o label. Vide https://stackoverflow.com/questions/8100765/y-axis-label-falling-outside-graphics-window
mar.default = c(5, 4, 4, 2) + 0.1
par(mar=mar.default+c(0, 0.5, 0, 0))
plot(0, 0,  
     ylim=range(mu.gmvp.boot, mu.vec), 
     xlim=range(sd.gmvp.boot, sig.vec), 
     ylab=expression(mu[p]),
     xlab=expression(sigma[p]), 
     pch=16, col="blue", cex.lab=1.5, type='n')  
# plot bootstrapped GMVP
points(sd.gmvp.boot, mu.gmvp.boot, col="green", pch=16)
points(gmvp.port$sd, gmvp.port$er, pch=16, cex=1.5, col="magenta")
points(sig.vec, mu.vec, pch=16, cex=1.5)
text(sig.vec, mu.vec, labels=codigos, pos=4)
text(gmvp.port$sd, gmvp.port$er, labels="GMVP", pos=2)


############################################################
# DEMO 17
############################################################
# Bootstrapping o portf?lio eficiente com retorno alvo = 0.01
ef.01 = efficient.portfolio(er = mu.vec, cov.mat = sigma.mat,
                            target.return = 0.01)
ef.01
# bootstrap
n.boot = 1e4
# aloca espa?o para os resultados
sd.ef.boot = rep(0, n.boot)
w.ef.boot = matrix(0, n.boot, n.ativos)
colnames(w.ef.boot) = codigos
# bootstrap
for (i in 1:n.boot) {
  boot.idx = sample(n.obs, replace=TRUE)
  ret.boot = ret.df[boot.idx, ] 
  mu.boot = apply(ret.boot, 2, mean)
  cov.boot = cov(ret.boot) 
  ef.boot = efficient.portfolio(mu.boot, cov.boot, target.return = 0.01)
  sd.ef.boot[i] = ef.boot$sd
  w.ef.boot[i, ] = ef.boot$weights
}

# plota os resultados
# gerando o slide da apresenta??o
#
layout(mat = matrix(1:4, 2, 2, byrow = T))
hist(sd.ef.boot, freq = F, col = 'slateblue1', main = 'sd', xlab = 'Valor')
hist(w.ef.boot[, 'ABEV3'], freq = F, col = 'slateblue1', main = 'ABEV3', xlab = 'Valor')
hist(w.ef.boot[, 'BBAS3'], freq = F, col = 'slateblue1', main = 'BBAS3', xlab = 'Valor')
hist(w.ef.boot[, 'CIEL3'], freq = F, col = 'slateblue1', main = 'CIEL3', xlab = 'Valor')
layout(1)
#
# voc? n?o deveria remover os "outliers", mas...
layout(mat = matrix(1:4, 2, 2, byrow = T))
x = sd.ef.boot
hist(x[!x %in% boxplot.stats(x)$out], freq = F, col = 'slateblue1', main = 'sd', xlab = 'Valor')
x = w.ef.boot[, 'ABEV3']
hist(x[!x %in% boxplot.stats(x)$out], freq = F, col = 'slateblue1', main = 'ABEV3', xlab = 'Valor')
x = w.ef.boot[, 'BBAS3']
hist(x[!x %in% boxplot.stats(x)$out], freq = F, col = 'slateblue1', main = 'BRFS3', xlab = 'Valor')
x = w.ef.boot[, 'CIEL3']
hist(x[!x %in% boxplot.stats(x)$out], freq = F, col = 'slateblue1', main = 'CIEL3', xlab = 'Valor')
layout(1)


############################################################
# DEMO 17a
############################################################
# summary statistics
out = matrix(0, 4, 3)
rownames(out) = c('sd', 'ABEV3', 'BBAS3', 'CIEL3')
colnames(out) = c('Observado', 'M?dia', 'SE')
out['sd',] = c(ef.01$sd, mean(sd.ef.boot), sd(sd.ef.boot))
out['ABEV3',] = c(ef.01$weights['ABEV3'], mean(w.ef.boot[, 'ABEV3']), sd(w.ef.boot[, 'ABEV3']))
out['BBAS3',] = c(ef.01$weights['BBAS3'], mean(w.ef.boot[, 'BBAS3']), sd(w.ef.boot[, 'BBAS3']))
out['CIEL3',] = c(ef.01$weights['CIEL3'], mean(w.ef.boot[, 'CIEL3']), sd(w.ef.boot[, 'CIEL3']))
out

# muitos outliers...
# talvez fosse melhor fazer...
out = matrix(0, 4, 3)
rownames(out) = c('sd', codigos)
colnames(out) = c('Observado', 'M?dia', 'SE')
x = sd.ef.boot
not.outliers = !x %in% boxplot.stats(x)$out
out['sd',] = c(ef.01$sd, mean(x[not.outliers]), sd(x[not.outliers]))
for (codigo in codigos) {
  x = w.ef.boot[, codigo]
  not.outliers = !x %in% boxplot.stats(x)$out
  out[codigo,] = c(ef.01$weights[codigo], mean(x[not.outliers]), sd(x[not.outliers]))
}  
out


############################################################
# DEMO 18
############################################################
# bootstrap efficient frontier
# fronteira eficiente a partir das amostras
ef = efficient.frontier(mu.vec, sigma.mat, nport = 20, alpha.max = 1, alpha.min = -2)
#
# tra?a o gr?fico
plot(ef$sd, ef$er, type="b", pch='o', col="blue",
     ylim=c(-0.007982799, 0.131248878), xlim=c(0.03060674, 0.47784044),
     ylab=expression(mu[p]), xlab=expression(sigma[p]))
#
# bootstrap efficient frontier
set.seed(4)
for(i in 1:50) {
  idx = sample(x = n.obs, replace = T)
  ret.boot = ret.df[idx, ]
  er.boot = apply(ret.boot, 2, mean)
  cov.boot = cov(ret.boot)
  ef.boot = efficient.frontier(er = er.boot, cov.mat = cov.boot, nport = 20, alpha.max = 1, alpha.min = -2)
  points(ef.boot$sd, ef.boot$er, type="l", col=i)
}


############################################################
# DEMO 19
############################################################
# rolling 36-month global minimum variance portfolios
#
rollGMVP = function(ret.roll) {
  mu.roll = apply(ret.roll, 2, mean)
  sigma.roll = cov(ret.roll)
  gmvp = globalMin.portfolio(er = mu.roll, cov.mat = sigma.roll)
  ans = c(gmvp$er, gmvp$sd, gmvp$weights)
  return(ans)
}

roll.gmvp = rollapply(
  data = ret.df,   # s?rie de observa??es
  width=36,        # tamanho da janela
  by.column=FALSE, # if TRUE, FUN is applied to each column separately
  FUN=rollGMVP)    # a fun??o a ser aplicada
colnames(roll.gmvp) = c('er', 'sd', codigos)
head(roll.gmvp)
#
# plota os pesos dos ativos no GMVP
plot(as.zoo(roll.gmvp[, codigos]), main="", 
     plot.type="single", col=1:3, lwd=3, ylab="pesos")    
legend(x="topleft", legend=codigos, lty=rep(1, 3), col=1:3, lwd=3)
#
# plota m?dias e sds "rolantes" do GMVP
plot(as.zoo(roll.gmvp[, c('er', 'sd')]), 
     ylim=c(-0.007456047, 0.055248210),
     main = '', plot.type="single", 
     ylab="%", 
     col=c("black","blue"), 
     lwd=3)
legend(x="topleft", legend=c("m?dia rolante", "sd rolante"), lty=rep(1, 2),
       col=c("black", "blue"), lwd=3)
# tra?a o intervalo de confian?a do desvio padr?o.
# usa os dados do slide 73 para tra?ar os intervalos de confian?a
abline(h=0.043736534+2*0.005755838, col = 'blue')
abline(h=0.043736534-2*0.005755838, col = 'blue')
# tra?a o intervalo de confian?a do retorno esperado
# usa os dados do slide 73 para tra?ar os intervalos de confian?a
abline(h=0.004481293+2*0.005968670, col = 'black')
abline(h=0.004481293-2*0.005968670, col = 'black')


############################################################
# DEMO 20
############################################################
# rolando o portf?lio eficiente com retorno alvo = 0.01
rollFunction = function(ret.roll) {
  mu.roll = apply(ret.roll, 2, mean)
  sigma.roll = cov(ret.roll)
  eff = efficient.portfolio(er = mu.roll, cov.mat = sigma.roll, target.return = 0.01)
  ans = c(eff$er, eff$sd, eff$weights)
  return(ans)
}

roll.eff = rollapply(
  data = ret.df, 
  width=36, 
  by.column=FALSE, 
  FUN=rollFunction)

colnames(roll.eff) = c('er', 'sd', codigos)
head(roll.eff)
#
# plota os pesos "rolantes" no portf?lio eficiente
plot(as.zoo(roll.eff[, codigos]), main="", 
     plot.type="single", col=1:3, lwd=3, ylab="pesos")    
abline(h=0)
legend(x="bottomleft", legend=codigos, lty=rep(1, 3), col=1:3, lwd=3)
#
# plota as m?dias e desvios padr?es "rolantes" do portf?lio eficiente
plot(as.zoo(roll.eff[, c('er', 'sd')]), plot.type="single", ylab="%", main="", 
     ylim = c(0, 0.22), col=c("black","blue"), lwd=3)
legend(x="topright", legend=c("m?dia rolante", "sd rolante"), lty=rep(1, 2),
       col=c("black", "blue"), lwd=3)
# tra?a o intervalo de confian?a do desvio padr?o "rolante"
# usa os dados do slide 78 para tra?ar os intervalos de confian?a
abline(h=0.062171112+2*0.02299613, col = 'blue', lty = 3)
abline(h=0.062171112-2*0.02299613, col = 'blue', lty = 3)


