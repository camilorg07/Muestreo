library(sqldf)
library(readr)
library(dplyr)
library(ggplot2)
library(corrplot)

setwd("D:/User/Lenovo/Documents/ProyectoMuestreo")

#1. Tratamiento de datos -------------------------------------------------------


brasil <- read.csv("microdados_ed_basica_2023.csv", sep = ";")


brasil=sqldf("select *
        from brasil
        where QT_MAT_BAS is not NULL")

brasil_filtrado <- sqldf("
  SELECT 
    CO_UF,       -- CC3digo del estado
    NO_UF,       -- Nombre del estado
    CO_MUNICIPIO, -- CC3digo del municipio
    NO_MUNICIPIO, -- Nombre del municipio
    CO_ENTIDADE, -- CC3digo de la escuela
    NO_ENTIDADE, -- Nombre de la escuela
    QT_MAT_BAS,  -- NC:mero de matrC-culas en educaciC3n bC!sica
    TP_DEPENDENCIA, -- Tipo de administraciC3n (PC:blica, Privada, Federal, etc.)
    TP_LOCALIZACAO, -- UbicaciC3n (Urbana/Rural)
    TP_LOCALIZACAO_DIFERENCIADA, -- UbicaciC3n diferenciada (ej. en comunidades indC-genas)
    IN_INTERNET_ALUNOS, -- Internet disponible para los alumnos (SC-/No)
    IN_BIBLIOTECA, -- Presencia de biblioteca (SC-/No)
    QT_SALAS_UTILIZADAS, -- NC:mero de salas de aula en uso
    QT_DOC_BAS,  -- NC:mero de docentes en educaciC3n bC!sica
    QT_DESKTOP_ALUNO,
    QT_TUR_BAS
    
  FROM brasil
")


brasil_filtrado<-brasil_filtrado%>%filter(!is.na(QT_DOC_BAS))

brasil_filtrado<-brasil_filtrado%>%filter(QT_DESKTOP_ALUNO!=88888)

rm(brasil)


#2. Muestreo en 2 etapas -------------------------------------------------------

# ------------------------------------------
# ParC!metros Poblacionales
# ------------------------------------------


# NC:mero total de UFs (conglomerados primarios)
N_I <- brasil_filtrado %>% 
  distinct(CO_UF) %>% 
  nrow()

# Total de municipios (subconglomerados) en la poblaciC3n
M_total <- brasil_filtrado %>% 
  distinct(CO_ENTIDADE) %>% 
  nrow()

# TamaC1o promedio de municipios por UF
d <- brasil_filtrado %>% 
  group_by(CO_UF) %>% 
  summarise(n = n()) %>% 
  summarise(mean(n)) %>% 
  pull()

# ------------------------------------------
# CC!lculo de Varianzas y Rho
# ------------------------------------------

# Varianza intra-UFs (S_w)
S_w <- brasil_filtrado %>%
  group_by(CO_UF) %>%
  summarise(s = var(QT_MAT_BAS), n = n()) %>%
  summarise(S_w = sum((n - 1) * s) / (nrow(brasil_filtrado) - N_I)) %>% 
  pull()

# Varianza total (S_U)
S_U <- var(brasil_filtrado$QT_MAT_BAS)

# Coeficiente de correlaciC3n intraclase
rho <- 1 - S_w / S_U

S_w / S_U


# ------------------------------------------
# TamaC1o de Muestra y DiseC1o
# ------------------------------------------

# ParC!metros de diseC1o
alpha <- 0.05
Z <- qnorm(1 - alpha/2)
e <- 3  # PrecisiC3n deseada

# TamaC1o de muestra MAS
n_mas <- (Z^2 * S_U) / e^2

# Efecto de diseC1o
deff <- 1 + (d - 1) * rho
n <- ceiling(deff * n_mas)

# Conglomerados (UFs) a seleccionar
m <- sqrt(27 * (1 - rho) / rho) %>% ceiling()

# Subconglomerados (municipios) por UF
k <- ceiling(n / m)

# ------------------------------------------
# Etapa 1: Seleccionar UFs (MAS)
# ------------------------------------------

conglomerados_seleccionados <- brasil_filtrado %>%
  distinct(CO_UF) %>%
  slice_sample(n = m) %>%
  pull(CO_UF)

# ------------------------------------------
# Etapa 2: Seleccionar Municipios por UF (MAS)
# ------------------------------------------

muestra_final <- brasil_filtrado %>%
  filter(CO_UF %in% conglomerados_seleccionados) %>%
  group_by(CO_UF) %>%
  group_modify(~ {
    municipios_disponibles <- unique(.x$CO_ENTIDADE)
    k_actual <- min(k, length(municipios_disponibles))
    sample_mun <- sample(municipios_disponibles, k_actual)
    .x %>% filter(CO_ENTIDADE %in% sample_mun)
  }) %>%
  ungroup()

# ------------------------------------------
# CC!lculo de Pesos
# ------------------------------------------

# Total de municipios por UF en la poblaciC3n
M_i_pob <- brasil_filtrado %>%
  group_by(CO_UF) %>%
  summarise(M_i = n_distinct(CO_ENTIDADE))

# Unir datos poblacionales a la muestra
muestra_final <- muestra_final %>%
  left_join(M_i_pob, by = "CO_UF") %>%
  group_by(CO_UF) %>%
  mutate(
    m_i = n_distinct(CO_ENTIDADE),  # Municipios seleccionados
    peso = (N_I / m) * (M_i / m_i)    # Peso final
  )


#3. Funci??n para hacer inferencia ----------------------------------------------

Inferencia <- function(variable, muestra_final, brasil_filtrado, alpha, tipo, cate=1){
  
  temp <- muestra_final
  N <- nrow(brasil_filtrado %>% distinct(CO_UF))  
  n <- length(unique(temp$CO_UF))        
  M <- nrow(brasil_filtrado)                     
  M_bar <- M / N 
  
  if(tipo=="p"){
    
    y_i_data <- temp %>%
      group_by(CO_UF) %>%
      summarise(
        y_i = mean(!!sym(variable) == cate),  # ProporciC3n en el cluster
        M_i = first(M_i),                     # Total de municipios en la UF (poblaciC3n)
        m_i = n()                             # Municipios muestreados en la UF
      )
    
    sum_Mi_yi <- sum(y_i_data$M_i * y_i_data$y_i)
    media_estimada <- (N / M) * (sum_Mi_yi / n)
    
    h <- media_estimada
    
    
    s_b2 <- sum((y_i_data$M_i * y_i_data$y_i - M_bar * h)^2) / (n - 1)
    
    s_i2 <- temp %>%
      group_by(CO_UF) %>%
      summarise(s_i2 = var(!!sym(variable)==cate)) %>%
      pull(s_i2)
    
    term1 <- ((N - n)/N) * (s_b2 / (n * M_bar^2))
    term2 <- (1/(n * N * M_bar^2)) * sum(y_i_data$M_i^2 * ((y_i_data$M_i - y_i_data$m_i)/y_i_data$M_i) * (s_i2 / y_i_data$m_i))
    var_h <- term1 + term2
    
    estimacion<-media_estimada
    total <- mean(brasil_filtrado[[variable]]==cate)
    
  }

  else{
  total_estimado <- sum(temp[[variable]] * temp$peso)
  media_estimada <- total_estimado / nrow(brasil_filtrado)
  h <- media_estimada
  
  y_i_data <- temp %>%
    group_by(CO_UF) %>%
    summarise(
      y_i = mean(!!sym(variable)),
      M_i = first(M_i),      # Total de municipios en la UF (de la poblaciC3n)
      m_i = n()              # Municipios muestreados en la UF
    )
  
  s_b2 <- sum((y_i_data$M_i * y_i_data$y_i - M_bar * h)^2) / (n - 1)
  
  s_i2 <- temp %>%
    group_by(CO_UF) %>%
    summarise(s_i2 = var(!!sym(variable))) %>%
    pull(s_i2)
  
  if (tipo =='m'){
    term1 <- ((N - n)/N) * (s_b2 / (n * M_bar^2))
    term2 <- (1/(n * N * M_bar^2)) * sum(y_i_data$M_i^2 * ((y_i_data$M_i - y_i_data$m_i)/y_i_data$M_i) * (s_i2 / y_i_data$m_i))
    var_h <- term1 + term2
    estimacion<-media_estimada
    total <- mean(brasil_filtrado[[variable]])
  }
  
  if (tipo=="t"){
    
    term1 <- ((N - n)/N)*(N^2/n) * (s_b2 / (n * M_bar^2))
    term2 <- (N/n) * sum(y_i_data$M_i^2 * ((y_i_data$M_i - y_i_data$m_i)/y_i_data$M_i) * (s_i2 / y_i_data$m_i))
    var_h <- term1 + term2
    estimacion <- total_estimado
    total <- sum(brasil_filtrado[[variable]])
    
  }
  }
  CV <- (sqrt(var_h) / estimacion) * 100
  
  IC_inf <- estimacion - qnorm(1-alpha/2) * sqrt(var_h)
  IC_sup <- estimacion + qnorm(1-alpha/2) * sqrt(var_h)
  
  return(list(ValorReal = total,Estimacion=estimacion,
              SE=sqrt(var_h),CV=CV,LI=IC_inf, LS=IC_sup ))
}



Inferencia("IN_BIBLIOTECA", muestra_final = muestra_final,
           brasil_filtrado = brasil_filtrado,
           alpha = alpha, tipo = "p", cate = 0)


#4. Gr??ficos de la presentacion -----------------------------------------------

ggplot(data=brasil_filtrado, aes(x=QT_MAT_BAS))+
  geom_histogram(binwidth = 10, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "NC:mero de docentes",
       x = "Docentes",
       y = "Frecuencia") +
  scale_x_continuous(limits = c(0, 1500), expand = c(0, 0),breaks = seq(0, 1500, by = 250)) +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))


ggplot(data=brasil_filtrado, aes(x=QT_DOC_BAS))+
  geom_histogram(binwidth = 10, fill = "#CD8162", color = "black", alpha = 0.7) +
  labs(title = "NC:mero de docentes",
       x = "Docentes",
       y = "Frecuencia") +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0),breaks = seq(0, 100, by = 20)) +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))


ggplot(data=brasil_filtrado, aes(x=QT_DESKTOP_ALUNO))+
  geom_histogram(binwidth = 10, fill = "#CDB5CD", color = "black", alpha = 0.7) +
  labs(title = "NC:mero de computadores",
       x = "Computadores",
       y = "Frecuencia") +
  scale_x_continuous(limits = c(0, 70), expand = c(0, 0),breaks = seq(0, 70, by = 10)) +
  theme_minimal()+
  scale_y_continuous(limits = c(0,30000))+
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggplot(data=brasil_filtrado, aes(x=QT_TUR_BAS))+
  geom_histogram(binwidth = 10, fill = "#FFFACD", color = "black", alpha = 0.7) +
  labs(title = "NC:mero de clases",
       x = "Clases",
       y = "Frecuencia") +
  scale_x_continuous(limits = c(0, 80), expand = c(0, 0),breaks = seq(0, 80, by = 10)) +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))


brasil_filtrado %>%
  count(TP_LOCALIZACAO) %>%  # Contar ocurrencias de cada categorC-a
  mutate(perc = n / sum(n) * 100) %>%  # Calcular porcentaje
  ggplot(aes(y = as.factor(TP_LOCALIZACAO), x = perc, fill = as.factor(TP_LOCALIZACAO))) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(perc, 1), "%")),  size = 4, vjust=.5, hjust=-.1) +
  labs(title = "LocalizaciC3n", 
       y = "", 
       x = "Porcentaje") +
  scale_fill_brewer(name = "CategorC-as", 
                    labels = c("Urbana", "Rural"),
                    palette = "Set2") +
  theme_minimal() +
  theme(legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))


brasil_filtrado %>%
  count(TP_DEPENDENCIA) %>%  # Contar ocurrencias de cada categorC-a
  mutate(perc = n / sum(n) * 100) %>%  # Calcular porcentaje
  ggplot(aes(y = as.factor(TP_DEPENDENCIA), x = perc, fill = as.factor(TP_DEPENDENCIA))) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(perc, 1), "%")),  size = 4, vjust=.5, hjust=-.1) +
  labs(title = "Tipo de Dependencia", 
       y = "", 
       x = "Porcentaje") +
  scale_fill_brewer(name = "CategorC-as", 
                    labels = c("Federal", "Estadual", "Municipal", "Privada"),
                    palette = "Accent") +
  theme_minimal() +
  theme(legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

brasil_filtrado %>%
  count(IN_INTERNET_ALUNOS) %>%  # Contar ocurrencias de cada categorC-a
  mutate(perc = n / sum(n) * 100) %>%  # Calcular porcentaje
  ggplot(aes(y = as.factor(IN_INTERNET_ALUNOS), x = perc, fill = as.factor(IN_INTERNET_ALUNOS))) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(perc, 1), "%")),  size = 4, vjust=.5, hjust=-.1) +
  labs(title = "B?Acceso a internet?", 
       y = "", 
       x = "Porcentaje") +
  scale_fill_brewer(name = "CategorC-as", 
                    labels = c("SC-", "SC-"),
                    palette = "Dark2") +
  theme_minimal() +
  theme(legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

brasil_filtrado %>%
  count(IN_BIBLIOTECA) %>%  # Contar ocurrencias de cada categorC-a
  mutate(perc = n / sum(n) * 100) %>%  # Calcular porcentaje
  ggplot(aes(y = as.factor(IN_BIBLIOTECA), x = perc, fill = as.factor(IN_BIBLIOTECA))) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(perc, 1), "%")),  size = 4, vjust=.5, hjust=-.1) +
  labs(title = "B?Biblioteca en la escuela?", 
       y = "", 
       x = "Porcentaje") +
  scale_fill_brewer(name = "CategorC-as", 
                    labels = c("No", "SC-"),
                    palette = "Accent") +
  theme_minimal() +
  theme(legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggplot(brasil_filtrado,aes(x=QT_DOC_BAS, y=QT_MAT_BAS, fill=as.factor(TP_LOCALIZACAO)))+
  geom_point(shape = 21, color = "black", alpha = 0.6, size = 1) +
  scale_x_continuous(limits=c(0,300))+
  scale_y_continuous(limits=c(0,10000))+
  scale_fill_brewer(labels=c("Urbana", "Rural"), palette = "Dark2",
                    name="LocalizaciC3n")+
  labs(title = "Docentes vs Matriculas",
       x = "NC:mero de docentes",
       y= "NC:mero de matrC-culas")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))+
  annotate("text", x = 250, y = 9000, 
           label = paste0("r = ", round(cor(brasil_filtrado$QT_DOC_BAS, brasil_filtrado$QT_MAT_BAS, use = "complete.obs"), 2)),
           size = 5, color = "black", fontface = "bold")

corrplot(corr = cor(brasil_filtrado[,c(7,13:16)]), 
         is.corr = FALSE,
         addgrid.col = NA, 
         method = "color",
         cl.lim = c(0, 1),
         tl.pos = "l", 
         tl.cex = 0.8,
         tl.col = 'black')
title(main = 'Matriz de correlaciones', line = 1.8, cex.main = 1.6)





ggplot(brasil_filtrado, aes(x=as.factor(CO_UF), y=QT_MAT_BAS, fill=as.factor(CO_UF)))+
  geom_violin(alpha=.5)+
  geom_boxplot(outlier.color="red", outlier.shape=16, outlier.size=1)+
  scale_fill_viridis_d()+
  scale_y_continuous(limits = c(0,5000))+
  theme_minimal()+
  labs(title = "DistribuciC3n por estado",
       x = "Estado",
       y = "NC:mero de matrC-culas")+
  theme(legend.position = "none",  # Ocultar leyenda
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),  # Rotar etiquetas del eje X
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))


brasil_filtrado%>%group_by(CO_UF)%>%summarise(s=var(QT_MAT_BAS))%>%
  ggplot(aes(y=s, x=CO_UF))+
  geom_line()



brasil_filtrado %>%
  count(CO_UF) %>%  # Contar ocurrencias de cada categorC-a
  mutate(perc = n / sum(n) * 100) %>%  # Calcular porcentaje
  ggplot(aes(y = as.factor(CO_UF), x = perc, fill = as.factor(CO_UF))) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(perc, 1), "%")),  size = 4, vjust=.5, hjust=-.1) +
  labs(title = "ProporciC3n de observaciones por UFs", 
       y = "", 
       x = "Porcentaje") +
  scale_fill_viridis_d(name = "CategorC-as")+
  theme_minimal() +
  theme(legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
