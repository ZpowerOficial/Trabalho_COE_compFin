# Carregar bibliotecas
library(readxl)
library(mvtnorm)

# Carregar os dados de preços de fechamento dos ativos
dados <- read_excel("Documents/Comp e Finanças/close.xlsx")

# Selecionar apenas as colunas numéricas (preços de fechamento dos ativos)
dados_numericos <- dados[sapply(dados, is.numeric)]

# Calculando os retornos dos ativos
#retornos <- dados_numericos[-1, ] / head(dados_numericos, -1) - 1
retornos = apply(X = log(dados_numericos), MARGIN = 2, FUN = diff)

# Calculando a matriz de covariância dos retornos
cov_matriz <- cov(retornos)

# Calculando as médias dos retornos
mu.vec <- colMeans(retornos, na.rm = TRUE)
# mu.vec <- apply(retornos, 2, mean)

# Inicializando as contas
conta <- rep(0, 7)

# Inicializar a matriz para armazenar os resultados
num_simulacoes <- 1e5
n_ativos <- ncol(dados_numericos)
n_retornos <- nrow(dados_numericos)
PF <- matrix(0, nrow = n_retornos, ncol = n_ativos)

# Laço de Monte Carlo
P0 <- dados_numericos[nrow(dados_numericos), ]

for (i in 1:num_simulacoes) {
  ret.sim <- rmvnorm(753, mean = mu.vec, sigma = cov_matriz)
  PF <- data.frame(matrix(0, nrow = nrow(ret.sim), ncol = n_ativos))
  for(j in 1:ncol(ret.sim)) {
    PF[1,j] = P0[j] * (1 + ret.sim[1,j])
    if(length(PF) > 0 && nrow(PF) > 1){
      PF[,j] = PF[1,j] * exp(cumsum(ret.sim[,j]))
    }
  }
  
  if (all(PF[126,] > P0)) {
    conta[1] <- conta[1] + 1
  } else if (all(PF[250,] > P0)) {
    conta[2] <- conta[2] + 1
  } else if (all(PF[378,] > P0)) {
    conta[3] <- conta[3] + 1
  } else if (all(PF[502,] > P0)) {
    conta[4] <- conta[4] + 1
  } else if (all(PF[629,] > P0)) {
    conta[5] <- conta[5] + 1
  } else if (all(PF[753,] > P0)) {
    conta[6] <- conta[6] + 1
  } else {
    conta[7] <- conta[7] + 1
  }
}

# Imprimindo as probabilidades para cada conta
(probabilidades <- conta / num_simulacoes)
