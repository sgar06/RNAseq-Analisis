# RNAseq-Analisis
#### Escrito por Sonia García Llorens
#### Actualizado el 7 Febrero 2025
### Contenido
1. [Preparación de los datos](#1-Preparación-de-los-datos)
     * 1.1 Descarga de los datos de RNA-seq del repositorio SRA con SRAtools  

2. [Procesamiento de los datos RNA-seq](#2-procesamiento-de-los-datos-rna-seq)
     * 2.1 Control de calidad, recorte de adaptadores y extremos de mala calidad
     * 2.2 Alineamiento contra genoma de referencia  
       - 2.2.1 Preparación del genoma de referencia  
       - 2.2.2 Alineamiento de las lecturas contra el genoma de referencia con HISAT2              
       - 2.2.3 Modificación y conversión de archivos SAM con SAMtools
     * 2.3 Identificación y recuento de features o características  
       - 2.3.1 Preparación del archivo de anotaciones gff  
       - 2.3.2 Recuento de características con htseq-count  
       - 2.3.3 Obtención de la matriz de recuentos  

3. [Analisis estadístico de los datos de RNAseq y Genes Diferencialmente Expresados](#3-Analisis-estadístico-de-los-datos-de-RNAseq-y-Genes-Diferencialmente-Expresados) 
     * 3.1 Instalación de edgeR  
       - 3.1.1 Importación de la matriz de recuentos y metadatos 
       - 3.1.2 Conversión de la matriz de recuentos al objeto DGEList
       - 3.1.3 Eliminación de genes con recuentos bajos
       - 3.1.4 Normalización de librerias y recuentos
       - 3.1.5 Estimacion de la variabilidad biologica entre muestras y réplicas
       - 3.1.6 Estudio de la dispersión de los genes
       - 3.1.7 Ajuste de la variabilidad de cada gen según la dispersión
       - 3.1.8 Prueba de significancia o Test de expresión diferencial
       - 3.1.9 Correción FP
       - 3.1.10 Filtrado de genes segun FDR y logFC

4. [Anotación de los genes](4-Anotación-de-los-genes)
     * 4.1 GeneOntology
     * 4.2 KEGG
     * 4.3 Reactome (open source and fully open acces)
     * 4.4 MSigDB

## Estructura de archivos
```console
RNAseq_analysis
|-- Code
|-- Data
|   |-- 1_Raw
|   |-- 2_Processed
|   |   |-- 1_Quality_Control
|   |   |-- 2_Trimming
|   |   `-- 3_Alignment
|   |-- 3_Annotation
|   |-- Reference_genome
|   `-- Supplementary
`-- Results
```

## Instalación de las herramientas a través de conda
Instalación de miniconda
Por defecto el programa se instala en el directorio *home*
```console
# Downloading miniconda
$ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
# Installing miniconda
$ bash miniconda.sh -b -u -p $HOME/miniconda
# Ejecutar conda por defecto en la terminal
$HOME/miniconda/bin/conda init bash
# Updating conda
$ conda update -q conda
```

Establecimiento de los canales de instalacion
```console
conda config --add channels default
conda config --add channels bioconda
conda config --add channels conda-forge
# Mostrar la preferencia de repositorios
conda config --show-sources
```

Creación de un nuevo environment y activación
```console
# Creación
conda create -n genomic_analysis
# Activación
conda activate genomic_analysis
```

Herramientas a usar 
| Programa | Versión | Canal | Comando de instalación | Utilidad | 
|---------|---------|----------|----------|----------|
sratools | 3.1.0 | bioconda | conda install -c bioconda sra-tools |
cutadapt | 3.5 |  | conda install cutadapt=3.5 |
trim-galore | 0.6.10 | | conda install -c bioconda trim-galore=0.6.10 |
hisat2 | 2.2.1 | bioconda | 
samtools | 1.21  | bioconda | conda install -c bioconda samtools |
htseq | 0.13.5  | bioconda | conda install -c bioconda htseq |
fastqc | 0.11.9 | bioconda | conda install -c bioconda fastqc |
multiqc | 1.19 | conda install -c bioconda multiqc |
IGV |   |   |



## 1 Preparación de los datos
### 1.1 Descarga de los datos de RNA-seq del repositorio SRA con SRAtools  

Datos [GSE261866](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261866)  
Bioproject PRJNA1089226  
Secuenciador: Illumina Novaseq 6000  
Lecturas: paired-end  
Información de hebra específica de ARN  (paired-end strand specific RNA)
Longitud: 101  
Alineamiento de secuencias a genoma hg38 con el alineador STAR. 
[Stranded or non-stranded reads](https://eclipsebio.com/eblogs/stranded-libraries/)  
![image](https://github.com/user-attachments/assets/fc97efa9-a336-4203-b60d-3a4602b8c204)  
![image](https://github.com/user-attachments/assets/b77955c8-6c20-4d15-a905-90c5987efe23)



```console 
fastq-dump --gzip --readids --split-3 SRR
```

--gzip: Compress output using gzip.  
--readids or -I: Append read ID after spot ID as ‘accession.spot.readid’. With this flag, one sequence gets appended the ID .1 and the other .2. Without this option, pair-ended reads will have identical IDs.  

--split-3 separates the reads into left and right ends. If there is a left end without a matching right end, or a right end without a matching left end, they will be put in a single file.

## 2 Procesamiento de los datos RNA-seq  
### 2.1 Control de calidad, recorte de adaptadores y extremos de mala calidad

trim-galore               0.6.10 (Recorte Phred Score <20, deteccion de adaptadore y filtrado de lect <20pb

### 2.2 Alineamiento contra genoma de referencia  
![image](https://github.com/user-attachments/assets/6f42b6c5-30d6-41b0-b99a-8e57c317e667)
**2.2.1 Preparación del genoma de referencia**  

Búsqueda:  UCSC Genome Browser o herramienta HISAT2  
Si el genoma no está indexado, indexación con hisat2-build 
El genoma de referencia se puede buscar en la base de Ensembl o en [HISAT2](http://daehwankimlab.github.io/hisat2/)  
Descargar el genoma de referencia hg38_genome.tar.gz (o GRCh38) y descomprimir  
```console
# Dentro de la carpeta Reference_genome
tar -xvf hg38_genome.tar.gz
```

**2.2.2 Alineamiento de las lecturas contra el genoma de referencia con HISAT2**  
HISAT2 usa menos recursos computacionalmente que STAR, pero STAR genera resultados más precisos  
Elementos que mapean 1 vez  
```console
hisat2 -k1 -U ../02.Trimming/SRR1552444_trimmed.fq.gz -x ../../Reference_genome/mm10/genome -S SRR1552444_hisat2.sam
```
-x : prefijo del índice del genoma de referencia [genome]
-U : lista de lecturas para ser alineadas [trimmed]
-S : archivo de salida en formato SAM
-k : define el número máximo de alineamientos por lectura.


**2.2.3 Modificación y conversión de archivos SAM con SAMtools**

Creación del archivo BAM, ordenación e indexado
```console
# Conversión al archivo BAM
samtools view -Sbh SRR1552444_hisat2.sam > SRR1552444_hisat2.bam
# Ordenación del archivo BAM por coordenadas genómicas
samtools sort SRR1552444_hisat2.bam -o SRR1552444_hisat2.sorted.bam
# Indexación del archivo BAM
samtools index SRR1552444_hisat2.sorted.bam
```
Se genera un archivo SRR1552444_hisat2.sorted.bam.bai cuyo alineamiento podemos visualizar en programas de visualización como IGV.
En este caso, podemos elegir el genoma de referencia y
cargar nuestro archivo sorted.bam previamente indexado para observar las lecturas y su alineamiento
sobre el genoma de referencia.

### 2.3 Identificación y recuento de features o características  

**2.3.1 Preparación del archivo de anotaciones GTF**  

El archivo de anatociones empleados es GRCh38.p14 descargado del [GENCODE](https://www.gencodegenes.org/human/)

**2.3.2 Recuento de características con htseq-count**  
```console
htseq-count -t exon -i gene_id --stranded=yes?? -f bam -r pos 2_Processed/3_Alignment/SRRxx_hisat2.sorted.bam 3_Annotation/gencode.xxx,gtf > ../Results/SRRxxx_counts.tsv
```
Con la opción -t exon indicamos que cuente las lecturas alineadas específicamente contra los exones y, con la opción -i gene_id, indicamos que agrupe las diferentes lecturas atendiendo al identificador del gen al que pertenecen. Con la opción --stranded=no y -s no, determinamos que las lecturas no provienen de un experimento de hebra específica y, por tanto, la lectura puede mapear en ambas hebras del genoma.  
Finalmente, las opciones –f y –r , las usamos para indicar el formato del archivo a emplear (bam) y cómo están ordenados los alineamientos, en este caso por posición o coordenadas genómicas (pos).  
Una vez hemos definido todas las opciones, indicamos la ruta de los archivos BAM y el archivo de anotaciones (GTF), así como la ruta de salida y el nombre del nuevo archivo que contendrá la información.

**2.3.3 Obtención de la matriz de recuentos** 

## 3 Analisis estadístico de los datos de RNAseq y Genes Diferencialmente Expresados
### 3.1 Instalación de edgeR  

```console
BiocManager::install("edgeR")
```
**3.1.1 Importación de la matriz de recuentos y metadatos**
```console
seqdata <- read.csv(file, sep=",", header=T)
```
**3.1.2 Conversión de la matriz de recuentos al objeto DGEList**
```console
y <- DGEList(seqdata)
```
Objeto DGElist: 2 apartados: $counts, $sample  +añadimos $genes
**3.1.3 Eliminación de genes con recuentos bajos**
```console
keep <- filterbyExpr(y)
y <- y[keep, keep.lib.sizes=F]
```
**3.1.4 Normalización de librerias y recuentos**
```console
y <- calcNormFactors(y)
```
**3.1.5 Estimacion de la variabilidad biologica entre muestras y réplicas**
```console
plotMDS(y) / PCA
```
**3.1.6 Estudio de la dispersión de los genes**
Creación de la matriz de diseño 
```console
design <- model.matrix()
```
Dispersión
```console
y <- estimateDisp(y, design, robust=T)
```
**3.1.7 Ajuste de la variabilidad de cada gen según la dispersión**
```console
glmQLfit()
```
**3.1.8 Prueba de significancia o Test de expresión diferencial**
```console
glmQLFTest()
```
**3.1.9 Correción FP**
```console
topTags()
```
**3.1.10 Filtrado de genes segun FDR y logFD**
FDR<0.05 o <0.01 y logFC >= 2

## 4 Anotación de los genes
### 4.1 GeneOntology
Bioconductor as the R package GO.db.
### 4.2 KEGG
Bioconductor , R package KEGG.db
### 4.3 Reactome (open source and fully open acces)
### 4.4 MSigDB
Panther? Biocarta? EnrichR, GSEA, DAVID, GOstats
## Interactome analysis: STRING-DB
