# Análisis de Datos de Metabolómica - PEC1

Este repositorio contiene los archivos y el código necesarios para el análisis de un dataset de metabolómica, como parte del PEC1. El objetivo del proyecto es la exploración de datos metabolómicos, utilizando el entorno **R** y la creación de un contenedor **SummarizedExperiment** para almacenar los datos y metadatos de manera estructurada.

## Estructura del Proyecto

- **Informe**: El documento PDF con el informe final (`informe_PEC1.pdf`), que describe cada paso del análisis, incluyendo la selección del dataset, creación del contenedor `SummarizedExperiment`, exploración de datos, y las conclusiones.
- **Datos**:
  - `GastricCancer_NMR.xlsx`: Archivo con los datos originales descargados.
  - `description.md`: Archivo Markdown que con detalles del dataset.
- **Código R**:
  - `Vergaray-DelAguila-Lisseth-PEC1.R`: Script en R que contiene el código utilizado para cargar, explorar, y analizar el dataset.
- **Contenedor SummarizedExperiment**:
  - `summarized_experiment.Rda`: Objeto `SummarizedExperiment` en formato `.Rda`, que contiene los datos y metadatos en una estructura integrada.

## Descripción del Análisis

El análisis sigue los siguientes pasos:
1. **Selección del Dataset**: Se elige un dataset de metabolómica del repositorio [MetaboData](https://github.com/nutrimetabolomics/metaboData) o [MetabolomicsWorkbench](https://www.metabolomicsworkbench.org/).
2. **Creación del Contenedor SummarizedExperiment**: Los datos se estructuran en un contenedor `SummarizedExperiment`, una clase especializada en R para el manejo de datos ómicos y sus metadatos.
3. **Exploración de los Datos**: Se realizan estadísticas descriptivas y visualizaciones para comprender mejor las características del dataset.
4. **Documentación**: Se registra cada paso en el informe, proporcionando transparencia y reproducibilidad.

## Requisitos

- **R** y los paquetes **SummarizedExperiment** y **BiocManager**.
- **RStudio**.

## Cómo Ejecutar el Código

1. **Instalar los Paquetes Requeridos**:
   ```r
   if (!requireNamespace("BiocManager", quietly = TRUE)) {
     install.packages("BiocManager")
   }
   BiocManager::install("SummarizedExperiment")
