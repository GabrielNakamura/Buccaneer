

library(ggplot2)
library(dplyr)
#library(buccaneer)
library(tibble)
library(purrr)
library(tidyr)


load(here::here("data",
                "tabela_final_valores_completos_todos_clados_Canidae_e_Felidae.RData"))


# Aplicar o filtro a cada data frame na lista
# para fazer o betwwen apenas entre Canidae mesmo
listas_tabelas_finais_between <-
  lapply(listas_tabelas_finais, function(df) {
  df %>% filter(clade %in% c("Borophaginae", "Hesperocyoninae", "Caninae"))
})


grupo_focal <- "Caninae"

# Função que faz TUDO para cada tabela:
processar_tabela <- function(tabela, grupo_focal) {

  # Calcular a distância Euclidiana baseada em mass e rbl
  atributos <- tabela[, c("mass", "rbl")]
  dist_eucl <- dist(atributos)
  matriz_dist <- as.matrix(dist_eucl)

  # Adicionar a coluna group (focal_group x other)
  tabela <- tabela %>%
    mutate(group = ifelse(clade == grupo_focal, "focal_group", "other"))

  # Retorna uma lista com duas coisas: a tabela com group + a matriz de distância
  return(list(
    tabela_com_group = as_tibble(tabela),
    matriz_distancia = matriz_dist
  ))
}

# Agora aplicar em todas as 50 tabelas:
listas_processadas <- lapply(listas_tabelas_finais_between, processar_tabela, grupo_focal = grupo_focal)



##################################################################################################################
# aqui vou pegar todas as especies de caninae q coexistem na fatia de tempo de 10 myr
dados_caninae <- subset(listas_processadas[[1]][["tabela_com_group"]], clade == "Caninae")


# Defina sua fatia de tempo
inicio_fatia <- 10.0
fim_fatia <- 10.1

# Filtra as espécies que estavam vivas em algum ponto dentro dessa janela
coexistentes <- subset(dados_caninae, TS >= fim_fatia & TE <= inicio_fatia)

# Visualiza as espécies que coexistem
coexistentes$spp



##################################################################################################################
# agora vou calcular a distancia média entre essas 5 sps usando o vizinho mais proximo (1 esecie)


# Criar o data frame
dados <- data.frame(
  spp = c("Leptocyon_matthewi", "Leptocyon_vafer", "Metalopex_macconnelli",
          "Metalopex_merriami", "Vulpes_stenognathus"),
  massa = c(0.738, 0.643, 0.658, 0.852, 0.857),
  rbl = c(0.675, 0.6885297, 0.6330419, 0.6780723, 0.6666442)
)


library(ggplot2)
library(ggrepel)

ggplot(dados, aes(x = massa, y = rbl)) +
  geom_point(size = 3, color = "blue") +
  geom_text_repel(aes(label = spp), size = 3) +
  labs(x = "Massa corporal", y = "RBL") +
  theme_minimal()



# Calcular matriz de distâncias euclidianas
matriz_dist <- as.matrix(dist(dados[, c("massa", "rbl")]))

# Substituir a diagonal por NA (não queremos comparar uma espécie com ela mesma)
diag(matriz_dist) <- NA

# Para cada espécie, encontrar a menor distância (vizinho mais próximo), ignorando os NAs
dist_nn <- apply(matriz_dist, 1, min, na.rm = TRUE)

# Visualizar os valores
dist_nn

mnnd <- mean(dist_nn)
mnnd


#######
# agora vamos conferir se bate!
#######

tabela=listas_processadas[[1]][["tabela_com_group"]]


atributos <- tabela[, c("mass", "rbl")]
dist_eucl <- dist(atributos)
matriz_dist <- as.matrix(dist_eucl)


# Rodar o regional_mnd_between
regional_mnd_between_mnd_2 <- clade_regional_distance(
  df.TS.TE = tabela,
  time.slice = 0.1,
  dist.trait = matriz_dist,
  nearest.taxon = 1,
  round.digits = 1,
  species = "spp",
  TS = "TS",
  TE = "TE",
  group = "group",
  group.focal.compare = c("focal_group", "other"),
  type.comparison = "between"
)

# Rodar o regional_mnd_within
regional_mnd_within_mnnd_2 <- clade_regional_distance(
  df.TS.TE = tabela,
  time.slice = 0.1,
  dist.trait = matriz_dist,
  nearest.taxon = 1,
  round.digits = 1,
  species = "spp",
  TS = "TS",
  TE = "TE",
  group = "group",
  group.focal.compare = c("focal_group", "other"),
  type.comparison = "within"
)






##################################################################################################################
# agora vou calcular a distancia média entre essas 5 sps usando os
# dois vizinhos mais proximos (2 esecie)


# Calcular matriz de distâncias
matriz_dist <- as.matrix(dist(dados[, c("massa", "rbl")]))
diag(matriz_dist) <- NA  # ignorar comparação da espécie com ela mesma

# Para cada espécie, pegar as duas menores distâncias e tirar a média
media_dois_vizinhos <- apply(matriz_dist, 1, function(x) {
  sort(x, na.last = NA)[1:2] |> mean()
})

# Visualizar as médias individuais
media_dois_vizinhos

mnnd_2 <- mean(media_dois_vizinhos)
mnnd_2

#######
# agora vamos conferir se bate!
#######

tabela=listas_processadas[[1]][["tabela_com_group"]]


atributos <- tabela[, c("mass", "rbl")]
dist_eucl <- dist(atributos)
matriz_dist <- as.matrix(dist_eucl)


# Rodar o regional_mnd_between
regional_mnd_between_mnd_2 <- clade_regional_distance(
  df.TS.TE = tabela,
  time.slice = 0.1,
  dist.trait = matriz_dist,
  nearest.taxon = 2,
  round.digits = 1,
  species = "spp",
  TS = "TS",
  TE = "TE",
  group = "group",
  group.focal.compare = c("focal_group", "other"),
  type.comparison = "between"
)

# Rodar o regional_mnd_within
regional_mnd_within_mnnd_2 <- clade_regional_distance(
  df.TS.TE = tabela,
  time.slice = 0.1,
  dist.trait = matriz_dist,
  nearest.taxon = 2,
  round.digits = 1,
  species = "spp",
  TS = "TS",
  TE = "TE",
  group = "group",
  group.focal.compare = c("focal_group", "other"),
  type.comparison = "within"
)




##################################################################################################################
# agora vou calcular a distancia média entre essas 5 sps usando os
# tres vizinhos mais proximos (3 esecie)


# Calcular matriz de distâncias
matriz_dist <- as.matrix(dist(dados[, c("massa", "rbl")]))
diag(matriz_dist) <- NA  # ignorar comparação da espécie com ela mesma

# Para cada espécie, pegar as duas menores distâncias e tirar a média
media_dois_vizinhos <- apply(matriz_dist, 1, function(x) {
  sort(x, na.last = NA)[1:3] |> mean()
})

# Visualizar as médias individuais
media_dois_vizinhos

mnnd_3 <- mean(media_dois_vizinhos)
mnnd_3

#######
# agora vamos conferir se bate!
#######

tabela=listas_processadas[[1]][["tabela_com_group"]]


atributos <- tabela[, c("mass", "rbl")]
dist_eucl <- dist(atributos)
matriz_dist <- as.matrix(dist_eucl)


# Rodar o regional_mnd_between
regional_mnd_between_mnd_2 <- clade_regional_distance(
  df.TS.TE = tabela,
  time.slice = 0.1,
  dist.trait = matriz_dist,
  nearest.taxon = 3,
  round.digits = 1,
  species = "spp",
  TS = "TS",
  TE = "TE",
  group = "group",
  group.focal.compare = c("focal_group", "other"),
  type.comparison = "between"
)

# Rodar o regional_mnd_within
regional_mnd_within_mnnd_2 <- clade_regional_distance(
  df.TS.TE = tabela,
  time.slice = 0.1,
  dist.trait = matriz_dist,
  nearest.taxon = 3,
  round.digits = 1,
  species = "spp",
  TS = "TS",
  TE = "TE",
  group = "group",
  group.focal.compare = c("focal_group", "other"),
  type.comparison = "within"
)




##################################################################################################################
# agora vou calcular a distancia média entre essas 5 sps usando todos
# os vizinhos mais proximos (as 5 especies nesse caso aqui)


# Calcular matriz de distâncias
matriz_dist <- as.matrix(dist(dados[, c("massa", "rbl")]))
diag(matriz_dist) <- NA  # ignorar comparação da espécie com ela mesma

# Para cada espécie, pegar as duas menores distâncias e tirar a média
media_dois_vizinhos <- apply(matriz_dist, 1, function(x) {
  sort(x, na.last = NA)[1:4] |> mean()
})

# Visualizar as médias individuais
media_dois_vizinhos

mnnd_all <- mean(media_dois_vizinhos)
mnnd_all

#######
# agora vamos conferir se bate!
#######

tabela=listas_processadas[[1]][["tabela_com_group"]]


atributos <- tabela[, c("mass", "rbl")]
dist_eucl <- dist(atributos)
matriz_dist <- as.matrix(dist_eucl)


# Rodar o regional_mnd_between
regional_mnd_between_mnd_2 <- clade_regional_distance(
  df.TS.TE = tabela,
  time.slice = 0.1,
  dist.trait = matriz_dist,
  nearest.taxon = "all",
  round.digits = 1,
  species = "spp",
  TS = "TS",
  TE = "TE",
  group = "group",
  group.focal.compare = c("focal_group", "other"),
  type.comparison = "between"
)

# Rodar o regional_mnd_within
regional_mnd_within_mnnd_2 <- clade_regional_distance(
  df.TS.TE = tabela,
  time.slice = 0.1,
  dist.trait = matriz_dist,
  nearest.taxon = "all",
  round.digits = 1,
  species = "spp",
  TS = "TS",
  TE = "TE",
  group = "group",
  group.focal.compare = c("focal_group", "other"),
  type.comparison = "within"
)

