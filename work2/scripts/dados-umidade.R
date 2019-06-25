rm(list = ls())

# Especificações ----------------------------------------------------------
setwd("C:/Users/User/Dropbox/5° Série/Controle Estatístico de Qualidade/trabalho2/dados")
source("C:/Users/User/Dropbox/functionsR/getSeason.R")
# setwd("/home/andrefbm/Dropbox/5° Série/Controle Estatístico de Qualidade/trabalho2")
# source("/home/andrefbm/Dropbox/functionsR/getSeason.R")

# Formatação banco --------------------------------------------------------
dados            <- read.table(file = "umidade.txt", sep = ";", header = TRUE, skip = 16, na.strings = "")
names(dados)     <- tolower(names(dados))
dados$data       <- as.Date(dados$data, format = '%d/%m/%Y')
dados00          <- subset(dados, dados$hora == 0, select = c("estacao", "data", "umidade.relativa.media", "velocidade.do.vento.media"))
dados12          <- subset(dados, dados$hora == 1200,  select = c("estacao", "data", "precipitacao"))
todos            <- merge(dados00, dados12, by = c("estacao", "data"))
todos$estacaoano <- getSeason(todos$data)
tail(todos)

# Exportação --------------------------------------------------------------
write.table(x = todos, file = "umidade83767.csv", sep = ";", row.names = F, append = F)

