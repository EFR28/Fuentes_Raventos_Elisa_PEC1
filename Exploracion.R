## PEC 1

# Datos: human_cachexia

# Cargamos dataset
human_cachexia <- read.csv("C:/Users/efuen/OneDrive/Escritorio/UOC/Semestre 3/Datos ómicos/human_cachexia.csv",
                           row.names=1, stringsAsFactors=TRUE)

# Guardamos los datos en formato binario y en formato texto para tenerlos en el respositorio
save(human_cachexia, file = "human_cachexia.Rda").
write.table(human_cachexia, file = "human_cachexia.txt", sep = ",", row.names = TRUE)

# Cargamos la librería necesaria de Bioconductor
library(SummarizedExperiment)

# Observamos la estructura de nuestros datos
str(human_cachexia)

# Trasponemos la matriz para tener los datos dispuestos de tal manera que cada columna sea una muestra
datos_met <- t(as.matrix(human_cachexia)[,2:64])
variables <- rownames(datos_met)
datos_met <- apply(datos_met, 2, as.numeric)
rownames(datos_met) <- variables

samples <- as.matrix(t(as.matrix(human_cachexia))[1,])

colnames(samples)[1] <- "Group" 

features <- as.matrix(rownames(datos_met))
colnames(features) <- "features"

se <- SummarizedExperiment(
  assays = list(counts = datos_met),
  rowData = features,
  colData = samples
)

se

# Normalización del conjunto de datos
# Cargamos los paquetes
library(POMA)

senorm <- se %>% 
  PomaNorm(method = "log_pareto")

senorm

# Representamos diagramas de cajas 
PomaBoxplots(senorm)


# Distribución de los sujetos
table(colData(senorm)$Group)

# Calcular estadísticas descriptivas para la matriz de conteos
counts <- t(assay(senorm))

# Cálculo de estadísticas por columnas (por variable)

apply(counts, 2, function(x) {
  c(mean = mean(x), 
    sd = sd(x), 
    median = median(x), 
    min = min(x), 
    max = max(x))
})

# Calcular la matriz de correlación
cor_matrix <- cor(counts, method = "pearson")

# Visualizar la matriz de correlación usando pheatmap
library(pheatmap)
pheatmap(cor_matrix,
         color = colorRampPalette(c("white", "skyblue", "navy"))(50),
         main = "Matriz de Correlación de Variables")


# Crear la anotación de los grupos desde rowData
annotation <- data.frame(Group = colData(senorm)$Group)
rownames(annotation) <- rownames(counts)

# Crear el heatmap con anotación
pheatmap(counts,
         annotation_row = annotation,   # Anotación de grupo
         scale = "column",              
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",      
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Niveles de Metabolitos por sujeto")


# Realizar el PCA
pca <- prcomp((counts), center = TRUE, scale. = FALSE)

# Calcular la varianza explicada
var_exp <- pca$sdev^2 / sum(pca$sdev^2) * 100
var_exp

# Graficar la varianza explicada por cada componente principal
barplot(var_exp[1:10], main = "Varianza Explicada por Componentes Principales", 
        xlab = "Componente Principal", ylab = "Porcentaje de Varianza Explicada")

loadings <- pca$rotation
(loadings[, 1:2])

# Gráfico por grupos
library(ggplot2)

pca_data <- as.data.frame(pca$x)
pca_data$Group <- annotation$Group

ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  labs(title = "PCA de Muestras", x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(legend.position = "bottom")