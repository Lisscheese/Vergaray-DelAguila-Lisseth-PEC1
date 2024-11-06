## ----setup, include=FALSE----------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, encoding = "UTF-8",error = FALSE)

# Librerías usadas en la realización de este informe
library(tinytex)
library(knitr)
library(rmdformats)
library(SummarizedExperiment)
library(readxl)
library(ggplot2)
library(reshape2)
library(corrplot)
library(pheatmap)


## ----echo=FALSE, message=FALSE, warning=FALSE, paged.print=FALSE-------------------------------------------------------------------------------------

url_dataset <- "https://raw.githubusercontent.com/nutrimetabolomics/metaboData/main/Datasets/2023-CIMCBTutorial/GastricCancer_NMR.xlsx"

download.file(url_dataset, destfile = "GastricCancer_NMR.xlsx", mode = "wb")

data <- read_excel("GastricCancer_NMR.xlsx")

# Visualizamos las primeras filas
head(data)


## ----message=FALSE, warning=FALSE, include=FALSE, paged.print=FALSE----------------------------------------------------------------------------------
# Instalamos SummarizedExperiment
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("SummarizedExperiment")


## ----echo=FALSE, message=FALSE, warning=FALSE, paged.print=FALSE-------------------------------------------------------------------------------------

# Leeremos cada sheet del data set por separado
data_sheet <- read_excel("GastricCancer_NMR.xlsx", sheet = "Data")
peak_sheet <- read_excel("GastricCancer_NMR.xlsx", sheet = "Peak")

# Construimos la matriz que contendra los valores de los metabolitos
metabolite_data <- as.matrix(data_sheet[, grep("^M", names(data_sheet))]) 

# Contruimos los row names y col names de ambas sheets
rownames(metabolite_data) <- data_sheet$SampleID  
colnames(metabolite_data) <- peak_sheet$Name      

# Creamos la matrix traspuesta
metabolite_data <- t(metabolite_data)

# Construimos los datos de las columnas
row_data <- DataFrame(
  Label = peak_sheet$Label,
  Perc_missing = peak_sheet$Perc_missing,
  QC_RSD = peak_sheet$QC_RSD,
  row.names = peak_sheet$Name
)
row_data

# Construimos los datos de las filas
col_data <- DataFrame(
  SampleType = data_sheet$SampleType,
  Class = data_sheet$Class,
  row.names = data_sheet$SampleID
)
col_data



## ----echo=FALSE, message=FALSE, warning=FALSE, paged.print=FALSE-------------------------------------------------------------------------------------

# Construimos el objeto
rse <- SummarizedExperiment(
  assays = list(counts = metabolite_data),
  rowData = row_data,
  colData = col_data
)


## ----------------------------------------------------------------------------------------------------------------------------------------------------
# Mostramos el obejto y los metadatos del mismo.
rse
dim(rse)
head(assay(rse),1)



## ----echo=FALSE, message=FALSE, warning=FALSE, paged.print=FALSE-------------------------------------------------------------------------------------

url_description <- "https://raw.githubusercontent.com/nutrimetabolomics/metaboData/main/Datasets/2023-CIMCBTutorial/description.md"
download.file(url_description, destfile = "description.md", mode = "wb")

# Leer el contenido (opcional)
contenido <- readLines("description.md")


dataset_description <- list(
  content = contenido,
  title = "1H-NMR urinary metabolomic profiling for diagnosis of gastric cancer",
  description = "Data from a gastric cancer study in article Chan et al. (2016), in the British Journal of Cancer",
  project_DOI = "ID del proyecto PR000699.",
  pub_year = "2016"
)
# Agregar esta información a los metadatos del objeto SummarizedExperiment
metadata(rse) <- dataset_description

# Verificar los metadatos
metadata(rse)


## ----echo=FALSE, message=FALSE, warning=FALSE, paged.print=FALSE-------------------------------------------------------------------------------------
# Filtramos los metabolitos con Perc_missing <= 20
filtered_indices <- which(rowData(rse)$Perc_missing <= 20)

# Creamps un nuevo objeto SummarizedExperiment con el filtrado aplicado
rse_filtered <- rse[filtered_indices, ]

# Verificamos la estructura del nuevo objeto
rse_filtered




## ----echo=FALSE, message=FALSE, warning=FALSE, paged.print=FALSE-------------------------------------------------------------------------------------

# Aplicamos filtro para QC_RSD
final_filtered_indices <- which(rowData(rse_filtered)$QC_RSD <= 20)

# Creamos un nuevo objeto SummarizedExperiment con ambos filtros aplicados
rse_final <- rse_filtered[final_filtered_indices, ]

# Verificamos la estructura del nuevo objeto filtrado
rse_final



## ----echo=FALSE, message=FALSE, warning=FALSE, paged.print=FALSE-------------------------------------------------------------------------------------
# Guardamos la matriz de ultimo objetivo SummarizedExperiment 
data_matrix <- assay(rse_final, "counts")

# Aplicamos transformación logarítmica
metabolite_data_log <- log2(data_matrix + 1)

# Guardamos la matriz normalizada en el objeto SummarizedExperiment
assay(rse_final, "normalized_counts") <- metabolite_data_log

rse_final



## ----------------------------------------------------------------------------------------------------------------------------------------------------
# Guardar el objeto SummarizedExperiment en un archivo .Rda
save(rse_final, file = "datos_metabolomica.Rda")


## ----echo=FALSE, message=FALSE, warning=FALSE, paged.print=FALSE-------------------------------------------------------------------------------------
par(mfrow = c(1, 2))

# Histograma de los datos originales (counts)
hist(assay(rse_final, "counts"), 
     main = "Distribución Metabolitos", 
     xlab = "Concentración (counts)", 
     ylab = "Frecuencia", 
     col = "purple", 
     breaks = 20)


# Histograma de los datos transformados logarítmicamente (log2)
hist(assay(rse_final, "normalized_counts"), 
     main = "Distribución Metabolitos (Normalized)", 
     xlab = "Concentración (log2)", 
     ylab = "Frecuencia", 
     col = "pink", 
     breaks = 20)

# Restablecer la configuración de gráficos a una única gráfica
par(mfrow = c(1, 1))


## ----echo=FALSE, message=FALSE, warning=FALSE, paged.print=FALSE-------------------------------------------------------------------------------------
sum(is.na(assay(rse_final, "normalized_counts")))

# Eliminamos filas con NAs
data_matrix_clean <- assay(rse_final, "normalized_counts")[complete.cases(assay(rse_final, "normalized_counts")), ]


# Realizamos PCA 
pca <- prcomp(t(data_matrix_clean), scale. = TRUE)


pca_data <- as.data.frame(pca$x)
pca_data$Condition <- colData(rse_final)$Class

# Graficamos PCA
ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 3) +
  labs(title = "Análisis de Componentes Principales (PCA)",
       x = "PC1", y = "PC2") +
  theme_classic()



## ----echo=FALSE, message=FALSE, warning=FALSE, paged.print=FALSE-------------------------------------------------------------------------------------
# Graficamos el heatmap  
pheatmap(assay(rse_final, "normalized_counts"),
         scale = "row", 
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         main = "Heatmap de Metabolitos Normalizados",
         color = colorRampPalette(c("beige", "orange", "red"))(50))



