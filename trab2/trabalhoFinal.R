# Integrantes do Grupo
# Christopher Henry Pinheiro Gonzaga DRE: 122047223
# Victor Wohlers Cardoso DRE: 119157174


# Carregar bibliotecas
library(readxl)
library(quadprog)
source("portfolio.r")
getwd()
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# (i) Retirando informacoes da planilha
dados <- read_excel("dados.xlsx", sheet = "close.df")
dados_codigos <- dados[sapply(dados, is.character)]
dados_numericos <- dados[sapply(dados, is.numeric)]
close.df <- data.frame(dados_numericos)

# (ii, iii, iv, v) Calculando retorno cc, esperanca e volatilidade media e covariancias
ret.df = apply(log(close.df), 2, diff)
mu.vec = apply(ret.df, 2, mean)
sig.vec = apply(ret.df, 2, sd)
sigma.mat = cov(ret.df)

# (vi Opcional) Criando grafico risco x retorno
plot(sig.p, mu.p, pch = 16, col = 'black',
     xlim = c(0, max(sig.p)),
     xlab = expression(sigma[p]), ylab = expression(mu[p]))

# (vii) Determinando GMVP pela formula do slide
sigma.inv = solve(sigma.mat)
gmvp.formula = rowSums(sigma.inv)/sum(sigma.inv)

# (viii) Calculando GMVP pelo portfolio.r
gmvp.zivot = globalMin.portfolio(er = mu.vec ,cov.mat = sigma.mat, shorts = TRUE)

# (ix Opcional) Acrescentar GMVP ao grafico

# (x Opcional) Calcular parametros e acrescentar no grafico

# (xi) Calculando portfolio tangente pelos multiplicadores de Lagrange
one.vec = rep(1, 18)
rf <- 0
tangente.formula = sigma.inv %*% (mu.vec - rf * one.vec) / as.numeric(t(one.vec) %*% sigma.inv %*% (mu.vec - rf * one.vec))

# (xii) Calculando portfolio tangente pelo portfolio.r
rf <- 0
tangente.zivot = tangency.portfolio(er = mu.vec, cov.mat = sigma.mat, risk.free = rf)

# (xiii Opcional) Colocar portfolio tangente no grafico

# (xiv Opcional) Coloca reta de investimentos eficientes

# (xv) 
