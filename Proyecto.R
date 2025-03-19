##########################################
#####Paquetes y librerias ################

install.packages("survey")
install.packages("rnaturalearth")  
install.packages("sf")             
install.packages("ggplot2")        
install.packages("dplyr")         
install.packages("RColorBrewer")  

library(rnaturalearth)
library(sf)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(readr)
library(survey)


getwd()
setwd("C:/Users/User/Desktop/GENERAL/UNI/MUESTREO ESTADISTICO/PROYECTO")


#-----------------------------------------
#-------arreglo base de datos-------------
#-----------------------------------------



rm(list=ls())

Censocanada <- read_csv("C:/Users/User/Desktop/GENERAL/UNI/MUESTREO ESTADISTICO/PROYECTO/censo_canada_2021.csv")


Provincias_tamaño_edad <- Censocanada %>%
  group_by(PR,HHSIZE,AGEGRP) %>% 
  count() %>%
  as.data.frame()





#########################################
####### datos geoespaciales##############

Provincias <- Censocanada %>%
  group_by(PR) %>% 
  count() %>%
  as.data.frame()

canada_shapefile <- ne_states(country = "Canada", returnclass = "sf")

canada_shapefile <- canada_shapefile %>%
  mutate(name = ifelse(name %in% c("Yukon", "Northwest Territories", "Nunavut"), "Northern Canada", name))


province_codes <- c("10" = "Newfoundland and Labrador", "11" = "Prince Edward Island", 
                    "12" = "Nova Scotia", "13" = "New Brunswick", "24" = "Quebec", 
                    "35" = "Ontario", "46" = "Manitoba", "47" = "Saskatchewan", 
                    "48" = "Alberta", "59" = "British Columbia", "70" = "Northern Canada")


Provincias$region <- province_codes[as.character(Provincias$PR)]


population_data <- Provincias %>%
  rename(population = n)


canada_shapefile <- canada_shapefile %>%
  left_join(population_data, by = c("name" = "region"))

canada_shapefile <- canada_shapefile %>%
  mutate(population = ifelse(name == "Quebec", 224250, population))

##################################################
################### mapa de calor ################


ggplot(data = canada_shapefile) +
  geom_sf(aes(fill = population)) + 
  geom_sf_text(data = canada_shapefile, 
               aes(label = ifelse(population > 30000, name, "")),  
               color = "black", fontface = "bold", size = 2, 
               fun.geometry = function(x) sf::st_centroid(x)) +  
  scale_fill_gradientn(colours = brewer.pal(9, "YlGnBu"), labels = scales::comma) +  
  theme_minimal() + 
  labs(fill = "Población", title = "Mapa de Calor de la Población en Canadá")  


brewer.pal




# ------------------------------
# 1. Tamaño de Muestra y Cálculo de Varianzas
# ------------------------------

# Valor Z para un 95% de confianza
Z <- qnorm(0.975)

# Precisión deseada en la estimación del total
d <- 10

# Calcular varianza entre conglomerados (S1^2) por provincia:
var_conglomerados <- Provincias_tamaño_edad %>%
  group_by(PR, HHSIZE) %>%
  summarise(Total_HH = sum(n), .groups = "drop") %>%   # Suma de n por combinación (PR, HHSIZE)
  group_by(PR) %>%
  summarise(S1_sq = var(Total_HH, na.rm = TRUE), .groups = "drop")

# Calcular varianza entre subconglomerados (S2^2) dentro de cada conglomerado:
var_subconglomerados <- Provincias_tamaño_edad %>%
  group_by(PR, HHSIZE, AGEGRP) %>%
  summarise(Total_AGE = sum(n), .groups = "drop") %>%
  group_by(PR, HHSIZE) %>%
  summarise(S2_sq = var(Total_AGE, na.rm = TRUE), .groups = "drop") %>%
  group_by(PR) %>%
  summarise(S2_sq = mean(S2_sq, na.rm = TRUE), .groups = "drop")

# Calcular varianza total de los subconglomerados (S_T^2) por provincia:
var_totales_subcong <- Provincias_tamaño_edad %>%
  group_by(PR, AGEGRP) %>%
  summarise(Total_AGE = sum(n), .groups = "drop") %>%
  group_by(PR) %>%
  summarise(S_T_sq = var(Total_AGE, na.rm = TRUE), .groups = "drop")

# Unir todas las varianzas en un solo data frame
varianzas <- var_conglomerados %>%
  left_join(var_subconglomerados, by = "PR") %>%
  left_join(var_totales_subcong, by = "PR") %>%
  mutate(
    rho = S1_sq / (S1_sq + S2_sq),   # Correlación intra-conglomerado
    rho_s = S2_sq / (S2_sq + S_T_sq)  # Correlación intra-subconglomerado
  )

# Tamaño de muestra bajo muestreo aleatorio simple (SRS)
n_SRS <- (Z^2 * varianzas$S1_sq) / d^2

# Efecto de diseño (DEFF)
DEFF <- 1 + (varianzas$S1_sq / varianzas$S2_sq) * (varianzas$rho / (1 - varianzas$rho))

# Para la primera etapa se fija el total de conglomerados disponibles en cada provincia
total_conglomerados <- 8

# Número efectivo de conglomerados a seleccionar
n_eff <- ceiling(n_SRS / DEFF)

# Se selecciona el mínimo entre n_eff y total_conglomerados
n_h <- pmin(n_eff, total_conglomerados)

# Para la segunda etapa se fija el total de subconglomerados disponibles (por ejemplo, 22)
total_subconglomerados <- 22
# Número óptimo de subconglomerados por conglomerado
b_opt <- sqrt(varianzas$S1_sq / varianzas$S2_sq)

# Crear data frame de resultados
resultado <- data.frame(
  Provincia = Provincias$PR,
  "Tamaño SRS" = ceiling(n_SRS),
  "Tamaño Ajustado (n_eff)" = ceiling(n_h),
  "Conglomerados totales a Muestrear" = ceiling(n_h),
  "Subconglomerados por Conglomerado (b*)" = floor(round(b_opt, 2))
)

# ------------------------------
# 2. Preparar Data Frame con Tamaños de Muestra por Provincia
# ------------------------------
resultado_mod <- resultado %>%
  rename(
    PR = Provincia,
    n_h = `Conglomerados.totales.a.Muestrear`,  # Se asume que al usar read.csv o similar, los puntos se reemplazan por _
    b  = `Subconglomerados.por.Conglomerado..b..`
  ) %>%
  select(PR, n_h, b)

# ------------------------------
# 3. Primera Etapa: Selección de Conglomerados (UPM: HHSIZE) por Provincia
# ------------------------------

# 3.1 Extraer la lista única de conglomerados (combinación de PR y HHSIZE) de la base original
clusters <- Provincias_tamaño_edad %>%
  distinct(PR, HHSIZE)

# 3.2 Unir el número de conglomerados a seleccionar (n_h) por provincia desde resultado_mod
clusters <- clusters %>%
  left_join(resultado_mod %>% select(PR, n_h), by = "PR")

# 3.3 Realizar el muestreo aleatorio simple (MAS) por provincia:
clusters_sample <- clusters %>%
  group_by(PR) %>%
  group_modify(~ {
    sample_size <- unique(.x$n_h)  # Número de conglomerados a seleccionar en cada provincia
    slice_sample(.x, n = min(sample_size, nrow(.x)))
  }) %>%
  ungroup()

print("Muestra de Conglomerados (UPM) seleccionados:")
print(clusters_sample)

# ------------------------------
# 4. Segunda Etapa: Selección de Subconglomerados (USM: AGEGRP) dentro de cada Conglomerado Seleccionado
# ------------------------------

# 4.1 Filtrar la base original para quedarse solo con registros pertenecientes a los conglomerados seleccionados
subclusters <- Provincias_tamaño_edad %>%
  semi_join(clusters_sample, by = c("PR", "HHSIZE")) %>%
  distinct(PR, HHSIZE, AGEGRP)

# 4.2 Unir el número de subconglomerados a seleccionar (b) por provincia desde resultado_mod
# Evitamos duplicados (si no existe ya 'b', se agrega directamente)
subclusters <- subclusters %>%
  left_join(resultado_mod %>% select(PR, b), by = "PR")

# 4.3 Realizar el MAS dentro de cada conglomerado (agrupado por PR y HHSIZE)
subclusters_sample <- subclusters %>%
  group_by(PR, HHSIZE) %>%
  group_modify(~ {
    b_value <- unique(.x$b)  # Número de subconglomerados a seleccionar en cada conglomerado
    slice_sample(.x, n = min(b_value, nrow(.x)))
  }) %>%
  ungroup()

print("Muestra de Subconglomerados (USM) seleccionados:")
print(subclusters_sample)

# ------------------------------
# 5. Cálculo de Probabilidades de Inclusión y Pesos
# ------------------------------

# 5.1 Calcular N_h: total de conglomerados (HHSIZE) por provincia en la base original
N_h_df <- Provincias_tamaño_edad %>%
  distinct(PR, HHSIZE) %>%
  group_by(PR) %>%
  summarise(N_h = n(), .groups = "drop")

# 5.2 Calcular B_HHSIZE: total de subconglomerados (AGEGRP) en cada conglomerado (HHSIZE) en la base original
B_HHSIZE_df <- Provincias_tamaño_edad %>%
  distinct(PR, HHSIZE, AGEGRP) %>%
  group_by(PR, HHSIZE) %>%
  summarise(B_HHSIZE = n(), .groups = "drop")

# 5.3 Unir la información de N_h y B_HHSIZE a la muestra final de subconglomerados
subclusters_sample <- subclusters_sample %>%
  left_join(N_h_df, by = "PR") %>%                        # Añade N_h por provincia
  left_join(B_HHSIZE_df, by = c("PR", "HHSIZE")) %>%      # Añade B_HHSIZE por conglomerado
  left_join(resultado_mod, by = "PR")                    # Une n_h y b; este join agrega columnas n_h y b
# Si aparecen columnas duplicadas para b, por ejemplo b.x y b.y, las resolvemos:
subclusters_sample <- subclusters_sample %>%
  mutate(b = coalesce(b.x, b.y)) %>%  # Se conserva el primer valor no-NA
  select(-b.x, -b.y)

# 5.4 Calcular las probabilidades de inclusión:
#    - pi_HHSIZE: probabilidad de inclusión de un conglomerado en la primera etapa
#    - pi_AGEGRP_given_HHSIZE: probabilidad de inclusión de un subconglomerado dado el conglomerado
#    - pi_AGEGRP: probabilidad total de inclusión del subconglomerado
subclusters_sample <- subclusters_sample %>%
  mutate(
    pi_HHSIZE = n_h / N_h,
    pi_AGEGRP_given_HHSIZE = b / B_HHSIZE,
    pi_AGEGRP = pi_HHSIZE * pi_AGEGRP_given_HHSIZE,
    w_ij = 1 / pi_AGEGRP  # Peso de Horvitz-Thompson
  )

# Incorporar la variable n (conteo) a partir de la base original, si aún no está
subclusters_sample <- subclusters_sample %>%
  left_join(
    Provincias_tamaño_edad %>% select(PR, HHSIZE, AGEGRP, n),
    by = c("PR", "HHSIZE", "AGEGRP")
  )

# ------------------------------
# 6. Estimación de Objetivos
# ------------------------------

# OBJETIVO PRINCIPAL: Estimar la media del número de personas por hogar en cada provincia.
# Se estima como el ratio:
#   (Total de personas estimado) / (Total de hogares estimado)
# Donde:
#   - total_persons = sum(n * w_ij)
#   - total_households = sum((n / HHSIZE) * w_ij)
mean_size_by_prov <- subclusters_sample %>%
  group_by(PR) %>%
  summarise(
    total_persons = sum(n * w_ij),
    total_households = sum((n / HHSIZE) * w_ij),
    mean_household_size = total_persons / total_households
  )

print("Media del número de personas por hogar en cada provincia:")
print(mean_size_by_prov)

# OBJETIVO ESPECÍFICO 1: Estimar el total de personas que residen en los hogares de cada provincia.
total_persons_by_prov <- subclusters_sample %>%
  group_by(PR) %>%
  summarise(
    total_persons = sum(n * w_ij)
  )

print("Total de personas estimado por provincia:")
print(total_persons_by_prov)

# OBJETIVO ESPECÍFICO 2: Estimar la proporción de personas mayores en cada hogar o provincia.
# En este caso, AGEGRP es un código de grupo de edad. Asumiremos que los códigos >= 17 corresponden a mayores de 65 años.
subclusters_sample <- subclusters_sample %>%
  mutate(
    is_elder = if_else(AGEGRP >= 17, 1, 0)  # Considera mayor de 65 si el código AGEGRP es 17 o más
  )

prop_elder_by_prov <- subclusters_sample %>%
  group_by(PR) %>%
  summarise(
    total_elder = sum(n * w_ij * is_elder),
    total_persons = sum(n * w_ij),
    prop_elder = total_elder / total_persons
  )

print("Proporción de personas mayores estimada por provincia:")
print(prop_elder_by_prov)

# ------------------------------
# OBJETIVO ESPECÍFICO 3:
# Analizar la relación entre la edad de los habitantes y el tamaño del hogar,
# estimando el total (expansionado) de cada conglomerado y comparándolo con el valor nominal de HHSIZE.
# ------------------------------

# Data frame que asocia cada AGEGRP con una edad representativa (mid_age)
age_mapping <- data.frame(
  AGEGRP = c(1:21, 88),
  mid_age = c(
    2,    # 1: 0-4
    5.5,  # 2: 5-6
    8,    # 3: 7-9
    12,   # 4: 10-14
    17,   # 5: 15-19
    22,   # 6: 20-24
    27,   # 7: 25-29
    32,   # 8: 30-34
    37,   # 9: 35-39
    42,   # 10: 40-44
    47,   # 11: 45-49
    52,   # 12: 50-54
    57,   # 13: 55-59
    62,   # 14: 60-64
    67,   # 15: 65-69
    72,   # 16: 70-74
    77,   # 17: 75-79
    82,   # 18: 80-84
    87,   # 19: 85-89
    92,   # 20: 90-94
    90,   # 21: 85+ (aproximado, se solapa con 19 y 20)
    NA    # 88: Not available
  )
)


# 1. Unir la información de edad representativa (mid_age) a la muestra
subclusters_sample <- subclusters_sample %>%
  filter(AGEGRP != 88) %>% 
  left_join(age_mapping, by = "AGEGRP")

# 2. Para cada conglomerado (HHSIZE) en cada provincia (PR), estimar:
#    a) Total de personas (usando n * w_ij).
#    b) Edad promedio ponderada: sum(mid_age * n * w_ij) / (total estimado de personas).
# Se asume que HHSIZE (valor nominal) es el tamaño declarado del hogar.
cluster_summary <- subclusters_sample %>%
  group_by(PR, HHSIZE) %>%
  summarise(
    # Estimación total de personas en el conglomerado (expansión Horvitz-Thompson)
    estimated_total = sum(n * w_ij),
    # Estimación de la suma ponderada de edad
    weighted_age = sum(mid_age * n * w_ij),
    # Edad promedio en el conglomerado
    average_age = weighted_age / estimated_total,
    .groups = "drop"
  ) %>%
  mutate(
    # Diferencia absoluta entre el total estimado y el tamaño nominal declarado (HHSIZE)
    difference = estimated_total - HHSIZE,
    # Diferencia relativa (en porcentaje)
    relative_diff = difference / HHSIZE
  )

# Imprimir la tabla resumen por conglomerado
print("Resumen por conglomerado (por PR y HHSIZE):")
print(cluster_summary)

# 3. Graficar la relación entre el tamaño declarado del hogar (HHSIZE) y el total estimado:
ggplot(cluster_summary, aes(x = HHSIZE, y = estimated_total)) +
  geom_point(color = "darkgreen") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(title = "Comparación: HHSIZE vs. Total Estimado de Personas",
       x = "HHSIZE (Tamaño Nominal del Hogar)",
       y = "Total Estimado de Personas (Expansión HT)") +
  theme_minimal()

# 4. Graficar la relación entre el tamaño del hogar y la edad promedio estimada:
ggplot(cluster_summary, aes(x = HHSIZE, y = average_age)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "orange") +
  labs(title = "Relación entre Tamaño del Hogar y Edad Promedio Estimada",
       x = "HHSIZE (Tamaño Nominal del Hogar)",
       y = "Edad Promedio Estimada") +
  theme_minimal()
