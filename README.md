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

## Decarga del genoma de referencia y el archivo de anotaciones de la especie *Homo sapiens*
Primero descarga del genoma re ferencia 
So we need to firstly retrieve the reference genome. Specifically for the example data set, we need a human reference genome.  
Búsqueda:  UCSC Genome Browser o herramienta HISAT2  
Si el genoma no está indexado, indexación con hisat2-build 
El genoma de referencia se puede buscar en la base de Ensembl o en [HISAT2](http://daehwankimlab.github.io/hisat2/)  
In the download page, data are grouped by species. At the Index section, you can see the links of different genome data such as human genome. At the human section, you can see the links of different human genome data, which are further grouped by different human reference genome versions. Here we want to the newest human reference genome (GRCh38/hg38). We need the [GRCh38 genome](https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz) 
Descargar el genoma de referencia hg38_genome.tar.gz (o GRCh38) y descomprimir  
```console
cd Reference_genome
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
# Dentro de la carpeta Reference_genome
tar -xvf grch38_genome.tar.gz
```
Se nos va a generar un carpeta /grch38/ con el genoma de referencia. Va a tener diferentes archivos genome.1 , genome.2 ... genome.8 y también el ejecutable. De esta forma ya lo tenemos indexado.   

Segundo decarga del archivo de anotaciones de referencia
![image](https://github.com/user-attachments/assets/6a04a794-3a5c-46ec-83aa-4a3a5b83413b)
![image](https://github.com/user-attachments/assets/43ba8b94-eae8-454c-810e-326a7d8d7da3)

  Descarga del archivo GTF
![image](https://github.com/user-attachments/assets/d24311bb-f95f-4b54-9a5c-fa0d52745cde)
![image](https://github.com/user-attachments/assets/98d1f57c-9493-4829-904c-8026c8ed7bb7)

  Tras la descarga del archivo de anotaciones GTF, descomprimimos el archivo con gunzip. Es imporante que pinchemos sobre incluir transcrit_id porque si no el programa convert2bed da error!
* ```console
  gunzip Homo_sapiens.gtf.gz
  ```


## 1 Preparación de los datos
### 1.1 Descarga de los datos de RNA-seq del repositorio SRA con SRAtools  

Datos [GSE261866](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261866)  
Bioproject PRJNA1089226  
Secuenciador: Illumina Novaseq 6000  
Lecturas: paired-end  
Información de hebra específica de ARN  (paired-end strand specific RNA)  
Longitud: 101  
Alineamiento de secuencias a genoma hg38 con el alineador STAR.   

**Descarga de datos crudos**
```console 
xargs -n1 fastq-dump --gzip --split-3 < SRR_Acc_List.txt
```
`xargs -n1`  
`--gzip`: Compress output using gzip.  
`--split-3` separates the reads into left and right ends. If there is a left end without a matching right end, or a right end without a matching left end, they will be put in a single file.

## 2 Procesamiento de los datos RNA-seq  

### 2.1.1 Estimation of the strandness
Para determinar si las lecturas RNA-seq son de hebra específica, primero se realizó un subsampling de lecturas a partir de una de las muestras. Para ello, se empleó la herramienta `seqkit` y se seleccionaron de forma aleatoria las lecturas tanto en el archivo `R1.fastq.gz`como `R2.fastq-gz`. Sin embargo, es importante que en el caso del subsampling aleatorio se seleccionen ambos extremos de cada par de lecturas, tanto  Forward como Reverse en el orden correcto.

Para ello se empleó el siguiente comando:
```console
seqkit sample -p 0.1 -s 100 {sample}_1.fastq.gz -o subsampled_{sample}_1.fastq.gz
seqkit sample -p 0.1 -s 100 {sample}_2.fastq.gz -o subsampled_{sample}_2.fastq.gz
```
> NOTA  
> `-p` se emplea para seleccionar la proporción de lecturas a seleccionar. En nuestro caso el 10% de lecturas totales.  
> `-s` se emplea para determinar el random seed.
> Ambos parametros tanto `-p`como `-s` tienen que ser los mismos en ambos archivos R1 y R2.

Posteriormente, las lecturas seleccionadas se alinean contra el genoma de referencia usando HISAT2. Seguidamente el alineamiento se compara contra el archivo de anotación de referencia para la especie Homo sapiens mediante la herramienta infer_experiment.py del paquete RseQC para determinar el tipo de librería empleada. 
En el caso de los experimento de hebra específicos, se pueden dar 2 escenarios:
* Lecturas forward o R1 situadas en la misma hebra del gene
* Lecturas reverse o R2 situadas en la misma hebra del gen




infer_experiment.py samples a few hundred thousand reads from your bam/sam and tells you the portion that would be explained by each type of strandedness, e.g  
[Stranded or non-stranded reads](https://eclipsebio.com/eblogs/stranded-libraries/)  
![image](https://github.com/user-attachments/assets/fc97efa9-a336-4203-b60d-3a4602b8c204)  
![image](https://github.com/user-attachments/assets/b77955c8-6c20-4d15-a905-90c5987efe23)  


| Tipo de Librería | Infer experiment | HISAT2 | htseq-count | 
|---------|---------|----------|----------|
Paired-End (PE) - SF | 1++,1–,2+-,2-+ | Second Strand F/FR | yes |
PE-SR | 1+-,1-+,2++,2– | First Strand R/RF | reverse |
Single-End (SE) - SF |	+,– | Second Strand F/FR | yes
SE - SR |	+-,-+ |	First Strand R/RF |	reverse
PE, SE - U |	undecided |	default	| no

## Cómo saber la hebra de procedencia de las lecturas
* A subset of reads is first made from the input FASTQ files
* Primero se hizo una subselección de lecturas
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

`--rna-strandedness` option in HISAT2  sets how reads are expected to align against genes. With this option being used, every read alignment will have an XS attribute tag: '+' means a read belongs to a transcript on '+' strand of genome. '-' means a read belongs to a transcript on '-' strand of genome.  
Most stranded protocols in use these days follow the dUTP-method, where read #2 in a pair has the same orientation as the transcript from which it arose (2++ or 2--). So either `R` or `RF` would typically be appropriate  
Use 'RF' if the first read in the pair corresponds to a transcript on the reverse strand, and the second read corresponds to the forward strand. When you use the `--rna-strandness` option with either 'FR' or 'RF' for paired-end reads, HISAT2 will assign an XS attribute tag to each read alignment, indicating whether the read belongs to a transcript on the '+' (plus) or '-' (minus) strand of the genome.  

### 2.1 Control de calidad, recorte de adaptadores y extremos de mala calidad

FastQC is a tool providing a simple way to do some quality control checks on the sequencing data. It checks different aspect of data quality and provides a graphical report so that one can intuitively get the idea about the data quality. Outputs an html report and a .zip file with the raw quality data  
```console
cd 1_Raw
mkdir initial_qc
fastqc -o initial_qc  *.fastq.gz
```
MultiQC Aggregates FastQC results of multiple analyses into a single report.  
```console
multiqc ??
```
Artefact removal.  Adapter trimming and quality-based trimming
trim-galore               0.6.10 (Recorte Phred Score <20, deteccion de adaptadore y filtrado de lect <20pb
```console
trim_galore --paired SAMPLE_R1.fastq.gz SAMPLE_R2.fastq.gz -o /2_Processed/2_Trimming/
```
![image](https://github.com/user-attachments/assets/96ca0af9-6caf-4244-872a-5405249788ce)  
Para ver el número de lecturas después 
```console
zcat Data/2_Processed/2_Trimming/SRR155244_trimmed.fq.gz | grep -c "@SRR"
```


### 2.2 Alineamiento contra genoma de referencia  
Once the quality of the data is confirmed, we need to convert those millions of reads per sample into the gene- or transcript-level quantification. This would need the assignment of reads to genes or transcripts.  
![image](https://github.com/user-attachments/assets/6f42b6c5-30d6-41b0-b99a-8e57c317e667)  

**2.2.1 Preparación del genoma de referencia**  
So we need to firstly retrieve the reference genome. Specifically for the example data set, we need a human reference genome.  
Búsqueda:  UCSC Genome Browser o herramienta HISAT2  
Si el genoma no está indexado, indexación con hisat2-build 
El genoma de referencia se puede buscar en la base de Ensembl o en [HISAT2](http://daehwankimlab.github.io/hisat2/)  
In the download page, data are grouped by species. At the Index section, you can see the links of different genome data such as human genome. At the human section, you can see the links of different human genome data, which are further grouped by different human reference genome versions. Here we want to the newest human reference genome (GRCh38/hg38). We need the [GRCh38 genome](https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz) 
Descargar el genoma de referencia hg38_genome.tar.gz (o GRCh38) y descomprimir  
```console
cd Reference_genome
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
# Dentro de la carpeta Reference_genome
tar -xvf grch38_genome.tar.gz
```
Se nos va a generar un carpeta /grch38/ con el genoma de referencia. Va a tener diferentes archivos genome.1 , genome.2 ... genome.8 y también el ejecutable. De esta forma ya lo tenemos indexado.  

**2.2.2 Alineamiento de las lecturas contra el genoma de referencia con HISAT2**  
HISAT2 usa menos recursos computacionalmente que STAR, pero STAR genera resultados más precisos  
Once the genome indexing is done, you are ready to map the reads to the reference genome with HISAT2.
Elementos que mapean 1 vez  
```console
#Paired-end reads
hisat2 -k1 --summary-file {sample}.summary.txt --rna-strandedness {STRING} -x (/ruta-genoma-ref/grch38/genome) -1 {sample_R1.fg.gz} -2 {sample_R2.fg.gz} |\
samtools view -Sbh > sample_alignment.bam 
```
> NOTA
> `-x` : prefijo del índice del genoma de referencia [genome]
> `--summary-file`
> `--rna-strandedness` en nuestro caso RF
> `-1` y `-2`: lecturas a alinear
> `-k`: define el número máximo de alineamientos por lectura


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
htseq-count -t exon -i gene_id --stranded=reverse -f bam -r pos 2_Processed/3_Alignment/SRRxx_hisat2.sorted.bam 3_Annotation/gencode.xxx,gtf > ../Results/SRRxxx_counts.tsv
```
Con la opción -t exon indicamos que cuente las lecturas alineadas específicamente contra los exones y, con la opción -i gene_id, indicamos que agrupe las diferentes lecturas atendiendo al identificador del gen al que pertenecen. Con la opción --stranded=no y -s no, determinamos que las lecturas no provienen de un experimento de hebra específica y, por tanto, la lectura puede mapear en ambas hebras del genoma.  
Finalmente, las opciones `–f` y `–r` , las usamos para indicar el formato del archivo a emplear (bam) y cómo están ordenados los alineamientos, en este caso por posición o coordenadas genómicas (pos).  
Una vez hemos definido todas las opciones, indicamos la ruta de los archivos BAM y el archivo de anotaciones (GTF), así como la ruta de salida y el nombre del nuevo archivo que contendrá la información.

**2.3.3 Obtención de la matriz de recuentos** 

## 3 Analisis estadístico de los datos de RNAseq y Genes Diferencialmente Expresados
### 3.1 Instalación de edgeR  

edgeR is implemented as R packages in Bioconductor. It expects the raw count matrix without normalization.  
```R
BiocManager::install("edgeR")
```
cargamos todas las librerías necesarias para el análisis y establecemos el directorio de trabajo.  
```R
# Libreria para la anotacion del gneoma humano
#no se si esta bien
BiocManager::install("Homo.sapiens")
library(Homo.sapiens)

#Innstalacion de  edgeR
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("edgeR")
library(edgeR)

#Librerias complementarias
```  

**3.1.1 Importación de la matriz de recuentos y metadatos**
```R
seqdata <- read.csv(file, sep=",", header=T)
```

With the example data set, we see 61852 rows and 25 columns, meaning 61852 annotated genes and 25 samples. That's a large number of genes, but are they all actually informative? We can firstly check their average expression levels.

Cambiamos los nombres de las muestras, para que sean los mismos que los especificados en el archivo de metadatos.

**3.1.2 Conversión de la matriz de recuentos al objeto DGEList**
```R
y <- DGEList(seqdata)
```
los nombres de las filas representan los identificadores de los genes en formato **(??)**, correspondiente a la base de  datos **(??)**
Para cada identificador del genoma Homo sapiens,  vamos a buscar los identificadores para estos mismos genes en formato Symbol. 

Objeto DGElist: 2 apartados: $counts, $sample  +añadimos $genes  
modificación de la columna sample$group para especificar grupo ctrl o enfermos (conversion primero a una variable categórica group <- as.factor(group)

**3.1.3 Eliminación de genes con recuentos bajos**
```R
keep <- filterbyExpr(y)
y <- y[keep, keep.lib.sizes=F]
```
Se eliminan los genes sin expresar o con expresión muy baja (0-10)  
**3.1.4 Normalización de librerias y recuentos**
Los recuentos filtrados y obtenidos previamente para cada gen se tienen que normalizar para corregir las diferencias debido a las profundidades de secuenciación irregular en cada muestra
```R
y <- calcNormFactors(y)
```
**3.1.5 Estimacion de la variabilidad biologica entre muestras y réplicas**
Another way to check similarities between samples is to use dimension reduction techniques.  
```R
plotMDS(y) / PCA
```
Para estimar la variabilidad biológica podemos usar un gráfico de escala multidimensional (MDS), con la función plotMDS del paquete limma. Este gráfico nos permite ver las relaciones entre muestras, de forma que las muestras con perfiles de expresión de genes similares estarán más cerca en el gráfico.

PCA: prcomp function considers rows as samples and columns as features.  
**3.1.6 Estudio de la dispersión de los genes**  

Creación de la matriz de diseño 
Para estimar la sobredispersión de los genes, vamos a emplear la función estimateDisp() del paquete edgeR. Esta función necesita que se le proporcione una matriz que contenga el diseño experimental que especifique cómo se asocian o agrupan las muestras.
Por tanto, primero creamos la matriz con la función model.matrix(~0 + group). Vamos a construir una matriz a partir de la variable categórica “group” de nuestro Environment y vamos a establecer las diferentes relaciones entre los grupos con el símbolo (~). Si no queremos que se emplee ningúngrupo de ref escribimos 0 **¿QUEREMOS QUE SE USE UN GRUPO COMO REFERENCIA?**  
```R
design <- model.matrix(~)
```
Dispersión de los genes
```R
y <- estimateDisp(y, design, robust=T)
```
robust=T protege la estimación contra los outliers  
![image](https://github.com/user-attachments/assets/9174f674-7070-4bed-a096-783874e07bb6)  

El resultado de aplicar la función estimateDisp() genera una serie de resultados estadísticos dentro del objeto “y”. La dispersión de los genes se va a evaluar desde 3 puntos de vista diferentes: “common dispersion”, “trended dispersion” y “tagwise dispersion”.
podemos evaluar cómo se ajustan los datos y decidir el tipo de dispersión más apropiada  

**3.1.7 Ajuste de la variabilidad de cada gen según la dispersión**  

Una vez hemos computado los valores de dispersión, edgeR usa estos valores para ajustar el nivel de variabilidad de cada gen. Por defecto, la dispersión empleada para datos de RNA-seq es “trended dispersion”.   Para modelar los datos de conteo de RNAseq y capturar la sobredispersión, usaremos una aproximación Binomial Negativa. Para ello, emplearemos una función específica de edgeR, glmQLfit. Esta función permite ajustar nuestros datos empleando un modelo lineal generalizado y usando un enfoque de QL (quasi-likelihood). De esta manera, se tiene en cuenta la sobredispersión de los genes con una mayor precisión.  
```R
fit <- glmQLfit(y, design, robust =T)
```
El nuevo objeto creado `fit` será un objeto de tipo DGEGLM. Si observamos el objeto “fit” con la función View, vemos que tendrá diferentes parámetros computados para cada uno de nuestros genes, tales como coeficientes, valores ajustados...
![image](https://github.com/user-attachments/assets/90c4e4ff-a3ac-4710-bc71-1f74d4642f6f)  

**3.1.8 Prueba de significancia o Test de expresión diferencial**
Control vs enfermedad  
Una vez computada la dispersión y ajustados los datos, vamos a llevar a cabo una prueba de significancia o test de expresión diferencial, de forma que ahora sí que queremos saber lo genes diferencialmente expresados entre grupos comparados o contraste.  
Primero, con la función makeContrast del paquete limma, vamos a definir el tipo de comparación que vamos a hacer entre los grupos experimentales. Además, la función makeContrast, requiere de la matriz de diseño experimental para saber qué muestras se asocian con estos grupos. La función makeContrast genera una matriz numérica que representa los grupos indicados a contrastar. El valor 1 y –1 corresponderá a los grupos a comparar.  
```R
# Grupos contraste a comparar
CvsL <- makeContrast(control-lupus, leves = design)
```
Test de expresión diferencial con la función glmQLFTest() del paquete de edgeR y guardamos los resultados en la variable. Se genera un objeto de tipo DGELRT. Le tenemos que indicar el objeti `fit` con el modelo ajustado y la variable con los grupos de contraste
```R
res_CvsL <- glmQLFTest(fit, contrast = CvsL)
```
El objeto res_CvsL, contendrá parámetros estadísticos comunes al objeto “fit” con el modelo ajustado. Sin embargo, aparecen unos subapartados nuevos que contendrán los resultados de la comparación. En concreto, los resultados obtenidos tras la comparación se guardarán en el subapartado “table” .
![image](https://github.com/user-attachments/assets/6f75917a-48cc-4b58-a31a-480b85843e59)  
Dentro de esta tabla tendremos diferentes estadísticos como logFC, logCPM, F y Pvalue obtenidos a partir de la comparación. El valor de Pvalue es un parámetro estadístico que se usa para estimar la probabilidad de que un resultado sea obtenido al azar, de forma que cuando la probabilidad de obtener ese mismo resultado al azar es menor al 5% (p < 0.05), decimos que se trata de un resultado estadísticamente significativo. Sin embargo, al comparar un número de genes tan alto, aumenta la probabilidad de tener valores de p significativos (p< 0.05), por puro azar, aunque no sean realmente significativos.
Of course, we shouldn't directly take those estimated p-values. In statistics, there is the multiple testing problem. 
![image](https://github.com/user-attachments/assets/d035687a-c58e-4ee1-b852-1fd980d37239)  


**3.1.9 Correción FP**
multiple testing correction techniques to make statistical tests more stringent in order to counteract the problem of multiple testing. There are quite some different techniques.
 There are other alternative approaches, for instance, the Benjamini–Hochberg (BH) correction, or FDR (False Discovery Rate) correction, which estimates an FDR for each test based on the assumption of 0-1 uniform distribution of the p-values when the null hypothesis holds. In practice, it counts the number of test with p-values no larger than the observed p-value of a test (k), and then it estimates the expected number of tests with no-larger p-value as p×n. The FDR is thus $ \frac{k}{p \times n}$. In R, both methods are implemented, together with others, as the function p.adjust

hay que corregir los valores de P. Para ello usamos la función topTags() del paquete edgeR. Esta función aplicará el método de Benjamini-Hochberg para corregir los valores de P, reduciendo así los falsos positivos.  
```R
res_correct_CvsL <- topTags(res_CvsL, n = Inf)
```

Se genera un objeto de tipo TopTags que contiene los resultados de la prueba de significancia para los grupos comparados. El dataset contenido en el elemento “table” lo almacenamos en una nueva variable, denominada de_data.
![image](https://github.com/user-attachments/assets/c5167d7d-1afe-4f9f-99e7-645e81239917)  
  
El dataset contenido en el elemento “table” lo almacenamos en una nueva variable, denominada de_data.
```R
de_CvsL <- res_correct_CvsL$table
```
Esta función topTags aplica la corrección de los P valores y crea una nueva columna con los valores corregidos (FDR).  

**3.1.10 Filtrado de genes segun FDR y logFD**
 FDR<0.05 o <0.01 y logFC >= 2 
A partir del objeto "de_CvsL", llevamos a cabo la selección de todas aquellas filas que cumplan ambos filtros especificados: valores de FDR <= 0.05 y valor absoluto de logFC >= 1.
Es decir, vamos a seleccionar todas las filas que cumplen los criterios especificados y vamos a reportar los resultados filtrados para todas las columnas. En el caso de la escala de logFC, escribimos el valor absoluto de 1 puesto que al ser una escala simétrica me van a interesar valores positivos y negativos. Valores de log de FC por debajo de –1 los anotaremos como genes infraexpresados, mientras que valores de log FC mayores a 1 los anotaremos como genes sobreexpresados.
```R
# numero de genes DE o bien UP o bien Down
isDE_CvsL <- de_CvsL[ de_CvsL$FDR =< 0.05 & abs(de_CvsL$logFC) >= 1,]

#Creacion de una nueva columna con Up y Down
de_CvsL$DEgene <- "Not Sign"
de_CvsL[de_CvsL$FDR =< 0.05 & de_CvsL$logFC >= 1, de_CvsL$DEgene] <- "Up"  #cambio de anotación de genes sobreexpresados
de_CvsL[de_CvsL$FDR =< 0.05 & de_CvsL$logFC <= -1, de_CvsL$DEgene] <- "Down" #cambio de anotacion de genes infraexpresados

# grafico
ggplot(de_CvsL, aes(x=logFC , y = -log10(FDR), color = DEgene) +  geom_point(size =1)
```

## 4 Anotación de los genes
we need to figure out that those different groups of DEGs mean biologically. This is not a simple task, as we need to know the functions of every gene within the gene set, summarize them together, and then compare with genes outside of the gene set and see whether the gene set is significantly enriched in participating or representing certain functions, processes, pathways, components or other biological features.  
### 4.1 Gene Ontology
Bioconductor as the R package GO.db. databases providing functional annotation of genes.  The Gene Ontology project provides an ontology of defined terms representing gene product properties.
GO ontology information is available in Bioconductor as the R package GO.db  
### 4.2 KEGG
Bioconductor , R package KEGG.db
database resource for understanding high-level functions and utilities of the biological system. It curates huge amount of molecular pathways in different species, providing information of how different genes, gene products and small molecules regulate, interact, or have chemical interactions with each other. It is one of the most comprehensive database of biological pathways, chemicals and drugs. The database is also available in Bioconductor, as the R package KEGG.db.
It is based on a subscription model. (only subscribers can download the complete data, others have limited access via web portal or REST API)  
### 4.3 Reactome (open source and fully open acces)
### 4.4 MSigDB
It is a resource of tens of thousands of annotated gene sets, but is only available for human and mouse. One can download the gene sets from the website, or alternatively, to use the R package msigdbr.
### Enrichment Analysis: frequency comparison
DAVID is probably the most commonly used tool for biologist to check functional enrichment given a gene list
### Enrichment analysis: rank distribution comparison
Among them, GSEA (Gene Set Enrichment Analysis) is the most widely used one. Instead of using a gene list of interest, it takes a ranked list of all the genes in the analysis. 

Panther? Biocarta? EnrichR, GSEA, DAVID, GOstats
## Interactome analysis: STRING-DB
