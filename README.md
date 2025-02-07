# RNAseq-Analisis
#### Escrito por Sonia García Llorens
#### Actualizado el 7 Febrero 2025
### Contenido
1. [Preparación de los datos](#1-Preparación-de-los-datos)
     - 1.1 Descarga de los datos de RNA-seq del repositorio SRA con SRAtools  

2. [Procesamiento de los datos RNA-seq](#2-procesamiento-de-los-datos-rna-seq)  
     - 2.1 Control de calidad, recorte de adaptadores y extremos de mala calidad\
      - 2.2 Alineamiento contra genoma de referencia  
      - 2.2.1 Preparación del genoma de referencia  
            Búsqueda:  UCSC Genome Browser o herramienta HISAT2  
            Si el genoma no está indexado, indexación con hisat2-build  
      - 2.2.2 Alineamiento de las lecturas contra el genoma de referencia con HISAT2  
              HISAT2 usa menos recursos computacionalmente que STAR, pero STAR genera resultados más precisos  
              Elementos que mapean 1 vez  
      - 2.2.3 Modificación y conversión de archivos SAM con SAMtools  
  - 2.3 Identificación y recuento de features o características  
      - 2.3.1 Preparación del archivo de anotaciones gff  
      - 2.3.2 Recuento de características con htseq-count  
      - 2.3.3 Obtención de la matriz de recuentos  

3. [Analisis estadístico de los datos de RNAseq y Genes Diferencialmente Expresados](#3-Analisis-estadístico-de-los-datos-de-RNAseq-y-Genes-Diferencialmente-Expresados) 
  - 3.1 Instalación de edgeR  
              BiocManager::install("edgeR")  
      - 3.1.1 Importación de la matriz de recuentos y metadatos 
              seqdata <- read.csv(file, sep=",", header=T)
      - 3.1.2 Conversión de la matriz de recuentos al objeto DGEList
              y <- DGEList(seqdata)
              Objeto DGElist: 2 apartados: $counts, $sample  +añadimos $genes
      - 3.1.3 Eliminación de genes con recuentos bajos
              keep <- filterbyExpr(y)
              y <- y[keep, keep.lib.sizes=F]
      - 3.1.4 Normalización de librerias y recuentos
              y <- calcNormFactors(y)
      - 3.1.5 Estimacion de la variabilidad biologica entre muestras y réplicas
              plotMDS(y) / PCA
      - 3.1.6 Estudio de la dispersión de los genes
            Creación de la matriz de diseño 
              design <- model.matrix()
            Dispersión
              y <- estimateDisp(y, design, robust=T)
      - 3.1.7 Ajuste de la variabilidad de cada gen según la dispersión
              glmQLfit()
      - 3.1.8 Prueba de significancia o Test de expresión diferencial
              glmQLFTest()
      - 3.1.9 Correción FP
              topTags()
      - 3.1.10 Filtrado de genes segun FDR<0.05 o <0.01 y logFC >= 2

4. [Anotación de los genes](4-Anotación-de-los-genes)
  - 4.1 GeneOntology
      Bioconductor as the R package GO.db.
  - 4.2 KEGG
      Bioconductor , R package KEGG.db
  - 4.3 Reactome (open source and fully open acces)
  - 4.4 MSigDB

## 1 Preparación de los datos

```console 
fastqc read1.fq
```

## 2 Procesamiento de los datos RNA-seq  

## 3 Analisis estadístico de los datos de RNAseq y Genes Diferencialmente Expresados

## 4 Anotación de los genes
