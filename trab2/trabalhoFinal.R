# Integrantes do Grupo
# Christopher Henry Pinheiro Gonzaga DRE: 122047223
# Victor Wohlers Cardoso DRE: 119157174


# Carregar bibliotecas
library(readxl)
library(quadprog)
library(rportfolios)
source("portfolio.r")
getwd()
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# (i) Retirando informacoes da planilha
dados <- read_excel("dados.xlsx", sheet = "close.df")
dados_numericos <- dados[sapply(dados, is.numeric)]
close.df <- data.frame(dados_numericos)
dados_codigos <- colnames(close.df)

# (ii, iii, iv, v) Calculando retorno cc, esperanca e volatilidade media e covariancias
(ret.df = apply(log(close.df), 2, diff))
(mu.vec = apply(ret.df, 2, mean))
(sig.vec = apply(ret.df, 2, sd))
(sigma.mat = cov(ret.df))

# (vi Opcional) Criando grafico risco x retorno de carteiras aleatórias
POPU = 1e4
x.mat = matrix(0, nrow = POPU, ncol = length(dados_codigos))

for (i in 1:POPU) {
  # gera um número entre 1 e 2, correspondente ao somatório das posições longas
  x.long = runif(1, 1, 2)
  x.mat[i, ] = random.longshort(
    n = length(dados_codigos), 
    x.t.long = x.long, 
    x.t.short = x.long - 1)
}

mu.p = x.mat %*% mu.vec
sig.p = rep(0,POPU)

# Loop for para calcular a volatilidade para cada linha da matriz
for (i in 1:POPU) {
  # Calcula a volatilidade para a i-?sima linha
  sig.p[i] <- sqrt(t(x.mat[i,]) %*% sigma.mat %*% x.mat[i,])
}

plot(sig.p, mu.p, pch = 16, col = "lightgrey",
     xlim = c(0, max(sig.p)), ylim = c(-0.0035, 0.008),
     xlab = expression(sigma[p]), ylab = expression(mu[p]))

# Marcando os 16 ativos no gráifico
points(sig.vec, mu.vec, pch=16, cex = 0.8, col='black')
text(sig.vec, mu.vec, dados_codigos, pos = 4, cex = 0.65, col = 'black')

# (vii) Determinando GMVP pela formula do slide
sigma.inv = solve(sigma.mat)
gmvp.formula = rowSums(sigma.inv)/sum(sigma.inv)

# (viii) Calculando GMVP pelo portfolio.r
gmvp.zivot = globalMin.portfolio(er = mu.vec ,cov.mat = sigma.mat, shorts = TRUE)

# (x Opcional) Calcular parametros e acrescentar no grafico
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
alfa.vec = seq(from = 4, to = -5.5, by = -0.1)
n.pontos = length(alfa.vec)

mu.0 = max(mu.vec)
X = portfolio_fronteira(mu.0, mu.vec, sigma.mat)

mu.p.fron = NULL
sig.p.fron = NULL
for (i in 1:n.pontos) {
  alfa = alfa.vec[i]
  aux = alfa*gmvp.formula + (1 - alfa)*X$x_vec
  mu.p.fron[i] = t(aux) %*% mu.vec
  sig.p.fron[i] = sqrt(t(aux) %*% sigma.mat %*% aux)
}

points(sig.p.fron, mu.p.fron, pch = 16, cex = 0.9, col = 'lightgreen')

# (ix Opcional) Acrescentar GMVP ao grafico
# Coloquei após a bala de Markowitz por motivos visuais
sig.gmvp = gmvp.zivot$sd
mu.gmvp = gmvp.zivot$er

points(sig.gmvp, mu.gmvp, pch = 16, col = 'red', cex = 1.2)
text(sig.gmvp, mu.gmvp, 'GMVP', pos = 2, cex = 0.75, col = 'red')


# (xi) Calculando portfolio tangente pelos multiplicadores de Lagrange
one.vec = rep(1, 18)
rf <- 0
tangente.formula = sigma.inv %*% (mu.vec - rf * one.vec) / as.numeric(t(one.vec) %*% sigma.inv %*% (mu.vec - rf * one.vec))

# (xii) Calculando portfolio tangente pelo portfolio.r
rf <- 0
tangente.zivot = tangency.portfolio(er = mu.vec, cov.mat = sigma.mat, risk.free = rf)

# (xiii Opcional) Colocar portfolio tangente no grafico

# (xiv Opcional) Coloca reta de investimentos eficientes

# (xv e xvi) Determinando GMVP onde vendas a descoberto não são permitidas
n <- length(mu.vec) # número de ativos
Dmat <- 2 * sigma.mat # Dmat na função solve.QP() é 2 vezes a matriz de covariância
dvec <- rep(0, n) # dvec é um vetor de zeros
Amat <- cbind(1, diag(n)) # matriz de restrições
bvec <- c(1, rep(0, n)) # vetor de restrições
sol <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1) # resolvendo o problema de otimização quadrática
(weights <- sol$solution) # pesos do portfólio GMVP

# (xix) Determinando o portfólio tangente sem vendas a descoberto
rf <- 0 # retorno do ativo livre de risco
(tg.s.short <- tangency.portfolio(er = mu.vec, cov.mat = sigma.mat, risk.free = rf, shorts = FALSE))
