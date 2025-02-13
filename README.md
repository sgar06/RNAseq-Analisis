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
| Programa | Versión | Comando de instalación | Utilidad | 
|---------|---------|----------|----------|
sratools | 3.2.0 | conda install -c bioconda sra-tools |
cutadapt | 5.0 | conda install cutadapt |
trim-galore | 0.6.10 | conda install -c bioconda trim-galore |
hisat2 | 2.2.1 | conda install -c bioconda hisat2 |
samtools | 1.21  | conda install -c bioconda samtools |
htseq | 2.0.5 | conda install -c bioconda htseq |
fastqc | 0.12.1 | conda install -c bioconda fastqc |
multiqc | 1.27 | conda install -c bioconda multiqc |
RseQC | 5.0.4 | conda install bioconda::rseqc |
bedops | | conda install -c bioconda bedops |


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

![image](https://github.com/user-attachments/assets/6a04a794-3a5c-46ec-83aa-4a3a5b83413b)
![image](https://github.com/user-attachments/assets/43ba8b94-eae8-454c-810e-326a7d8d7da3)

## Cómo saber la hebra de procedencia de las lecturas
* A subset of 200 000 reads is first made from the input FASTQ files
* Next the reads are aligned against the selected reference genome using hisat2. The alignment is then compared to reference annotation to infer the strandedness of reads.
* For strand specific experiments there are two scenarios:
  * Reads in file 1 are always on the same strand as the gene (sense)
  * Reads in file 2 are always on the same strand as the gene
* Preparation of annotation file in bed format
  Descarga del archivo GTF
![image](https://github.com/user-attachments/assets/d24311bb-f95f-4b54-9a5c-fa0d52745cde)
![image](https://github.com/user-attachments/assets/98d1f57c-9493-4829-904c-8026c8ed7bb7)

  Tras la descarga del archivo de anotaciones GTF, descomprimimos el archivo con gunzip. Es imporante que pinchemos sobre incluir transcrit_id porque si no el programa convert2bed da error!
* ```console
  gunzip Homo_sapiens.gtf.gz
  # empleo del script convert2bed del paquete bedops
  convert2bed --input=gtf < Homo_sapiens.gtf > Homo_sapiens.bed
  ```
* RSeQC script: infer_experiment.py
* ```console
  infer_experiment.py -r Homo_sapiens.bed  -i sample.bam | tee infer_strand.txt
  ```
  * Options:
  * -i : input alignment file SAM or BAM format
  * -r : reference gene model in bed  format
  
Resultados para la muestra 65: 
This is PairEnd Data  (tipo de librería antisentido: Reverse, reverse stranded)  
Fraction of reads failed to determine: 0.1925  
Fraction of reads explained by "1++,1--,2+-,2-+": 0.0158  
Fraction of reads explained by "1+-,1-+,2++,2--": 0.7917  

Resultados para la muestra 65:  
$ infer_experiment.py -r ~/Descargas/Homo_sapiens.bed -i sample66_alignment.bam   
Reading reference gene model /home/sgarciallorens/Descargas/Homo_sapiens.bed ... Done  
Loading SAM/BAM file ...  Total 200000 usable reads were sampled  

This is PairEnd Data  
Fraction of reads failed to determine: 0.2869  
Fraction of reads explained by "1++,1--,2+-,2-+": 0.0054  
Fraction of reads explained by "1+-,1-+,2++,2--": 0.7077  


**Descarga de datos crudos**
```console 
xargs -n1 fastq-dump --gzip --split-3 < SRR_Acc_List.txt
```
xargs -n1 
--gzip: Compress output using gzip.  
--split-3 separates the reads into left and right ends. If there is a left end without a matching right end, or a right end without a matching left end, they will be put in a single file.

## 2 Procesamiento de los datos RNA-seq  
### 2.1 Control de calidad, recorte de adaptadores y extremos de mala calidad

trim-galore               0.6.10 (Recorte Phred Score <20, deteccion de adaptadore y filtrado de lect <20pb
```console
trim_galore --paired SAMPLE_R1.fastq.gz SAMPLE_R2.fastq.gz
```

### 2.2 Alineamiento contra genoma de referencia  
![image](https://github.com/user-attachments/assets/6f42b6c5-30d6-41b0-b99a-8e57c317e667)  

**2.2.1 Preparación del genoma de referencia**  

Búsqueda:  UCSC Genome Browser o herramienta HISAT2  
Si el genoma no está indexado, indexación con hisat2-build 
El genoma de referencia se puede buscar en la base de Ensembl o en [HISAT2](http://daehwankimlab.github.io/hisat2/)  
Descargar el genoma de referencia hg38_genome.tar.gz (o GRCh38) y descomprimir  
```console
# Dentro de la carpeta Reference_genome
tar -xvf grch38_genome.tar.gz
```
Se nos va a generar un carpeta /grch38/ con el genoma de referencia. Va a tener diferentes archivos genome.1 , genome.2 ... genome.8 y también el ejecutable. De esta forma ya lo tenemos indexado. 

**2.2.2 Alineamiento de las lecturas contra el genoma de referencia con HISAT2**  
HISAT2 usa menos recursos computacionalmente que STAR, pero STAR genera resultados más precisos  
Elementos que mapean 1 vez  
```console
#Paired-end reads
hisat2 -k1 -x (/ruta-genoma-ref/grch38/genome) -1 sample_R1.fg.gz -2 sample_R2.fg.gz |\
samtools view -Sbh > sample_alignment.bam |\
tee alignment.txt
```
-x : prefijo del índice del genoma de referencia [genome]
-1 y -2: lecturas a alinear
-k : define el número máximo de alineamientos por lectura


**2.2.3 Modificación y conversión de archivos SAM con SAMtools**

Creación del archivo BAM, ordenación e indexado
```console
# Conversión al archivo BAM
samtools view -Sbh SRR1552444_hisat2.sam > SRR1552444_hisat2.bam
# Ordenación del archivo BAM por coordenadas genómicas
samtools sort SRR1552444_hisat2.bam -o SRR1552444_hisat2.sorted.bam
# Indexación del archivo BAM
samtools index SRR1552444_hisat2.sorted.bam
# Opción 2:
samtools sort --write-index 
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
cargamos todas las librerías necesarias para el análisis y establecemos el directorio de trabajo.  

**3.1.1 Importación de la matriz de recuentos y metadatos**
```console
seqdata <- read.csv(file, sep=",", header=T)
```

With the example data set, we see 61852 rows and 25 columns, meaning 61852 annotated genes and 25 samples. That's a large number of genes, but are they all actually informative? We can firstly check their average expression levels.

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
