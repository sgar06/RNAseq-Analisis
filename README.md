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
|   |-- 2_Infer_strandedness
|   |-- 3_Processed
|   |-- 4_Alignment
|   |   `-- Reference_genome
|   |-- 5_Annotation
|   `-- Supplementary
`-- Results
```

## Instalación de las herramientas a través de conda
Instalación de miniconda
Por defecto el programa se instala en el directorio *home*
```console
# Descarga de miniconda desde la terminal
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
# Instalación de miniconda
bash miniconda.sh -b -u -p $HOME/miniconda
# Ejecución de conda por defecto en la terminal
HOME/miniconda/bin/conda init bash
# Actualización de conda
conda update -q conda
```

Establecimiento de los canales de instalación
```console
conda config --add channels default
conda config --add channels bioconda
conda config --add channels conda-forge
# Mostrar la preferencia de repositorios
conda config --show-sources
```

Creación de un nuevo entorno y activación
```console
# Creación de un nuevo entorno denominada genomic_analysis
conda create -n genomic_analysis
# Activación del entorno
conda activate genomic_analysis
```

Herramientas a instalar dentro del entorno:  
| Programa | Versión | Comando de instalación | Utilidad | 
|---------|---------|----------|----------|
sratools | 3.2.0 | conda install -c bioconda sra-tools | Descarga de lecturas desde el repositorio NCBI |
seqkit | 2.9.0 | conda install bioconda::seqkit | Subselección de lecturas en una muestra |
bedops | 2.4.41 | conda install -c bioconda bedops | Conversión de archivos |
RseQC | 5.0.4 | conda install bioconda::rseqc | Análsis de calidad |
fastqc | 0.12.1 | conda install -c bioconda fastqc | Análisis de calidad |
multiqc | 1.27 | conda install -c bioconda multiqc | Análisis de calidad |
cutadapt | 5.0 | conda install cutadapt | Recorte de adaptadores |
trim-galore | 0.6.10 | conda install -c bioconda trim-galore | Procesado de lecturas |
hisat2 | 2.2.1 | conda install -c bioconda hisat2 | Mapeador de lecturas |
samtools | 1.21  | conda install -c bioconda samtools | Manejo de archivos BAM |
htseq | 2.0.5 | conda install -c bioconda htseq | Anotación de características genómicas |



## Decarga del genoma de referencia y el archivo de anotaciones de la especie *Homo sapiens*
**Descarga del genoma de referencia humano (GRcH38) indexado desde la página web HISAT2**   
_HISAT2_ es un programa de alineamiento rápido y eficiente capaz de alinear lecturas obtenidas tras la secuenciación contra diferentes genomas.  Desde su [repositorio online](http://daehwankimlab.github.io/hisat2/), se pueden descargar directamente diversos genomas indexados. En el caso de que el genoma de interés no esté indexado, la herramienta _HISAT2_ permite realizar una indexación manual con la función `hisat2-build`, aunque es un proceso lento y costoso.  

Es por ello, que descargamos directamente el genoma humano indexado desde su página web.
Para ello, primero, es necesario acceder a la sección de descargas y, posteriormente, en la sección _Index_ encontramos diferentes _links_ según el genoma de interés, siendo en nuestro caso el genoma perteneciente a la especie _H.sapiens_.  

![image](https://github.com/user-attachments/assets/ddfb3af8-1e12-4a78-be2b-753cf76d5133)  
  
Dentro de la sección _H.sapiens_, encontramos diferentes genomas de referencia según la versión, y para cada uno de ellos, se disponen de distintos _links_ de descarga según el genoma de interés. En este caso, se selecciona el genoma de referencia humano más reciente GRCh38.   
![image](https://github.com/user-attachments/assets/2a07a735-3247-4f0a-92fd-9e0da1b8994e)

  
![image](https://github.com/user-attachments/assets/d88c15e4-f20d-466d-a3de-b822aa82f16f)

Una vez conocido el _link_ de nuestro genoma de referencia, lo descargamos directamente a través de la terminal. Este nuevo archivo se almacena en una nueva carpeta.  
```console
mkdir ~/RNAseq_analysis/Data/4_Alignment/Reference_genome && cd $_
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
# Descomprisión del archivo descargado
tar -xvf grch38_genome.tar.gz
```
> NOTA  
> Con el comando `mkdir` creamos un nuevo directorio para el genoma de referencia dentro de nuestra carpeta 4_Alignment. El comando `cd $_` permite movernos a esta nueva carpeta dónde descargamos el archivo de interés con el comando `wget`.  
> Tras la descompresión del archivo, se genera el directorio `/grch38/` con el genoma de referencia y los archivos necesarios para la indexación. Esta indexación permitirá al alineador _HISAT2_, llevar a cabo un alineamiento de las lecturas más eficiente.    

**Descarga del archivo de anotaciones de referencia GRCh38 desde el repositorio _ENSEMBLE_**  

Para descargar el archivo de anotaciones de referencia de la especie _H.sapiens_, se emplea el repositorio [ENSEMBL](https://www.ensembl.org/Homo_sapiens/Tools/FileChameleon). Por otro lado, es importante que el archivo de anotaciones tenga ciertas características específicas necesarias para la correcta ejecución de programas posteriores. Es por ello, que se emplea la herramienta _File Chamaleon_ con el fin de formatear el archivo de anotaciones.  

Una vez seleccionado el genoma de interés (GRCh38.p14), se descarga el archivo de anotaciones en formato GTF y se marca la casilla de _transcript_id_ para incluir este campo en el archivo descargado.  
  
![image](https://github.com/user-attachments/assets/5141e4d9-41f4-4362-a1b2-ec78c2e01690)  
  
A continuación descargamos y descomprimimos el archivo de anotaciones en formato GTF en nuestra computadora, y ya podemos almacenarlo en el directorio de interés.  
```console
gunzip Homo_sapiens.gtf.gz
mv Homo_sapiens.gtf ~/RNAseq_analysis/Data/5_Annotation/
```
 
Una vez realizados todos los pasos anteriores, la estructura de directorios dentro de nuestro proyecto _RNAseq_analysis_ resulta de la siguiente manera:  
```console
RNAseq_analysis
|-- Code
|-- Data
|   |-- 1_Raw
|   |-- 2_Infer_strandedness
|   |-- 3_Processed
|   |-- 4_Alignment
|   |   |-- Reference_genome
|   |   |   |-- grch38
|   |   |   |   |-- genome.1.ht2
|   |   |   |   |-- genome.2.ht2
|   |   |   |   |-- genome.3.ht2
|   |   |   |   |-- genome.4.ht2
|   |   |   |   |-- genome.5.ht2
|   |   |   |   |-- genome.6.ht2
|   |   |   |   |-- genome.7.ht2
|   |   |   |   |-- genome.8.ht2
|   |   |   |   `-- make_grch38.sh
|   |   |   `-- grch38_genome.tar.gz
|   |-- 5_Annotation
|   |   `-- Homo_sapiens.gtf
|   `-- Supplementary
`-- Results

```

## 1 Preparación de los datos
### 1.1 Descarga de las lecturas crudas a partir del repositorio _Sequence Read Archive_(SRA) del NCBI  

Para llevar a cabo el análisis, se emplean lecturas crudas depositadas en el repositorio público _Gene Expression Omnibus _(GEO) bajo el índice [GSE261866](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261866). En esta página encontramos distinta información relativa al estudio, tales como el número de muestras, el organismo, el tipo de experimento o la plataforma de secuenciación. Además, en el apartado de _Supplementary file_ se encuentra un enlace hacia el repositorio _SRA Run selector_.   
  
![image](https://github.com/user-attachments/assets/adcc0a0a-690a-487e-a592-f34eb0c3fdfc)  

A partir del repositorio _SRA Run selector_ se pueden descargar tanto las lecturas crudas como los metadatos de las muestras. Además, existe una herramienta implementada por NCBI denominada _SRA-toolkit_ que permite descargar las lecturas directamente en la computadora.  

En el presente tutorial, se emplean las muestras de un estudio de RNAseq de neutrófilos aislados de pacientes con lupus eritematoso sistémico y controles sanos. Para llevar a cabo el flujo de trabajo, solamente se emplean 6 muestras de cada grupo seleccionadas aleatoriamente.  

Secuenciador: Illumina Novaseq 6000   
Librería: lecturas pareadas con información específica de hebra  
Longitud: 101 pb  

| Nombre de la muestra | Identificador de Lecturas | Estado | Lecturas (M) | Bases | Tamaño |
|---------|---------|----------|----------|----------|----------|
GSM8153253 | SRR28380565 | Sano | 28.5 | 5.7G | 1.8 Gb |
GSM8153252 | SRR28380566 | Sano | 50.8 | 10.3G | 3.1Gb |
GSM8153250 | SRR28380568 | Sano | 55.1 | 11.1 G | 3.5 Gb |
GSM8153248 | SRR28380570 | Sano | 52.7 | 10.6G	| 3.3Gb |
GSM8153246 | SRR28380572 | Sano | 33.4 | 6.8G | 2.1Gb |
GSM8153245 | SRR28380573 | Sano | 33.2 | 6.7G |	2.1Gb |
GSM8153238 | SRR28380580 | LES | 31.1M	| 6.3G | 2.0Gb |
GSM8153236 | SRR28380582 | LES | 47.8M  | 9.6G	| 2.9Gb
GSM8153234 | SRR28380584 | LES | 39.8M	| 8.0G	| 2.6Gb
GSM8153232 | SRR28380586 | LES | 39.5M	| 8.0G	| 2.5Gb
GSM8153230 | SRR28380588 | LES | 29.6M	| 6.0G	| 1.9Gb
GSM8153229 | SRR28380589 | LES | 31.5M	| 6.4G	| 2.0Gb
  
Para descargar la lista con los identificadores de las lecturas en formato de texto plano, primero se seleccionan las muestras de interés, se filtran y se descarga la lista con los identificadores o _Accesion list_.   
  
![image](https://github.com/user-attachments/assets/afbcce1d-228e-4994-ba4e-32d585c77471)  

  
**Descarga de lecturas crudas con _SRA-toolkit_**
Una vez descargado el archivo de texto plano _SRR_Acc_List.txt_ con los identificadores, se gurada en la ruta `~/RNAseq_analysis/Data/Supplementary` y se emplea el comando ``fastq-dump` de la herramienta _SRA-toolkit_ desde la terminal para obtener las lecturas crudas en nuestra computadora.  
Las lecturas crudas se almacenan en el directorio `1_Raw` de nustra estructura de archivos.  
```console
cd /home/sgarciallorens/RNAseq_analysis/Data/1_Raw
```  
Una vez situados en la carpeta de interés, empleamos el siguiente comando para descargar las lecturas:  
```console 
xargs -n1 fastq-dump --gzip --split-3 < ../Supplementary/SRR_Acc_List.txt
```  
> NOTA  
> * `xargs -n1` ejecuta el comando `fastq-dump` en las lecturas cuyos nombres figuran en el archivo _SRR_Acc_List.txt_.  Para llegar al archivo, indicamos la ruta relativa desde el directorio dónde nos encontramos.  
> * `fastq-dump` permite la descarga de lecturas desde el repositorio SRA en formato FASTQ y de forma comprimida con la opcción `--gzip`. Al tratarse de lecturas pareadas, la opción `--split-3` separa las lecturas _forward_ y _reverse_ en diferentes archivos, y en el caso de que alguna lectura no esté pareada, estas se anotan en un tercer archivo que contiene las lecturas desparejadas.  

El formato FASTQ contiene las siguientes características:  

## 2 Procesamiento de los datos RNA-seq  

### 2.1.1 Determinación de la direccionalidad de las lecturas
Para determinar si las lecturas RNA-seq contienen información de hebra específica, se emplea inicialmente una sola muestra y se realiza una selección aleatoria de entorno al 10% de lecturas totales incluyendo ambos extremos gracias a la herramienta _Seqkit_. Los reslutados los guardamos en el directorio `2_Infer_strandedness`.  

En este paso solamente nos interesa determinar el tipo de librería empleada en el estudio y la especificidad de hebra, por lo que se realiza un *subsampling* de lecturas a partir de una de las muestras con el fin de aumentar la eficiencia en el paso posterior. Con la herramienta `seqkit` se selccionan de forma aleatoria las lecturas tanto en el archivo `{sample}_1.fastq.gz` como `{sample}_2.fastq.gz`. Es importante que en el caso del *subsampling* aleatorio se seleccionen ambos extremos de cada par de lecturas, tanto  *Forward* como *Reverse* en el orden correcto por lo que se emplean los mismos parámetros.  

```console
cd ~/RNAseq_analysis/Data/2_Infer_strandedness
seqkit sample -p 0.1 -s 100 {sample}_1.fastq.gz -o subsampled_{sample}_1.fastq.gz
seqkit sample -p 0.1 -s 100 {sample}_2.fastq.gz -o subsampled_{sample}_2.fastq.gz
```
> NOTA  
> `-p` permite determinar la proporción de lecturas a seleccionar. En nuestro caso el 10% de lecturas totales.   
> `-s` permite determinar el random _seed_. En este caso, se especifica el valor de 100, permitiendo que las lecturas aleatorias sean las mismas en ambos archivos.  
> Ambos parametros tanto `-p`como `-s` tienen que ser los mismos en ambos archivos R1 y R2 para seleccionar las mismas lecturas aleatorias.  

Para comprobar el número de pares de lecturas seleccionadas en la submuestra, se puede especificar el siguiente comando:   
```console
zcat subsampled_{sample}_1.fastq.gz | grep -c @SRR
zcat subsampled_{sample}_2.fastq.gz | grep -c @SRR
```
En nuestro caso, como la submuestra se obtiene a partir de la muestra SRR28380565, cuyo contenido de lecturas iniciales era de 28.5 millones, al llevar a cabo la selección del 10% de lecturas, obtenemos archivos con 2.8 millones de lecturas.  
Una vez preparada la submuestra, se lleva a cabo un alineamiento rápido contra el genoma de referencia usando el mapeador _HISAT2_ y los resultados del alineamiento se comparan respecto al archivo de anotación de referencia mediante la herramienta `infer_experiment.py` del paquete _RseQC_ con el fin de determinar la direccionalidad de la librería.  

Como se ha comentado previamente, para realizar el mapeo se emplea la herramienta _HISAT2_ y el genoma de referencia humano indexado descargado previamente (_grch38_genome.tar.gz_). El resultado del alineamiento en formato `.bam` lo guardamos en la carpeta `2_Infer_strandedness`.  
```console
hisat2 -k1 -x ../4_Alignment/Reference_genome/grch38/genome \
-1 subsampled_{sample}_1.fastq.gz -2 subsampled_{sample}_2.fastq.gz |\
samtools view -Sbh > subsampled_alignment.bam
```  
> NOTA   
> `-k` permite determinar el número de veces que se permite que una misma lectura pueda alinear en varias ubicaciones dentro del genoma de referencia. En nuestro caso especificamos el valor 1, para impedir los alineamientos múltiples.  
> `-x` sirve para indicar la ruta de directorios hasta llegar al genoma de referencia indexado.  
> `-1` y `-2` permite indicar las lecturas filtradas _forward_ y _reverse_ de la submuestra.  
> La salida del comando `hisat2` genera un archivo SAM con los resultados del alineamiento. Sin embargo, como se trata de un archivo muy pesado, la salida del comando de `hisat2` se concatena directamente a la herramienta `samtools` a través de una tubería `|` con el fin evitar la creación de archivos intermedios. De esta forma el comando `samtools view -Sbh > subsampled_alignment.bam`  permite crear directamente el archivo BAM en código binario, ocupando menos espacio. La opción `-Sbh` permite indicar el  formato de tipo SAM que se usa como entrada (`S`), el archivo de salida en formato BAM (`b`) y el mantenimiento del encabezado (`h`).  

Para visualizar el resultado del alineamiento en formato binario (BAM), se puede emplear la herramienta `samtools`.  
```console
samtools view subsampled_alignment.bam | head
```
Los resultados de las primeras 10 líneas del alineamiento contenidos en el archivo BAM se muestran en la salida por pantalla de la terminal:  
```console
SRR28380565.3	83	1	184595	1	101M	=	14009	-170687	 CTCTCAACCACTTGAG[...]	 FFFFF:FFFFFFFF[...]	AS:i:0	ZS:i:0	XN:i:0	XM:i:0[...]
SRR28380565.3	163	1	14009	1	101M	=	184595	170687	 CACAGCCTTGCCTGGA[...]	 FFFFFFFFFFFFFF[...]	AS:i:0	ZS:i:0	XN:i:0	XM:i:0[...]
SRR28380565.11	99	12	14583	1	101M	=	14583	101	GCGCAGGCTGGGTGGAG[...]   FFFFFFFFFFFFFF[...]	AS:i:0	ZS:i:0	XN:i:0	XM:i:0[...]
SRR28380565.11	147	12	14583	1	101M	=	14583	-101	GCGCAGGCTGGGTGGAG[...] 	 FFFFFFFFFFFFFF[...]	AS:i:0	ZS:i:0	XN:i:0	XM:i:0[...]
SRR28380565.43	97	16	14232	0	3S98M	1	14678	0	GTTAGCCTTCCGCTCCC[...]   FFFFFFFFFFFFFF[...]	AS:i:-18  ZS:i:-18  XN:i:0  XM:i:3[...]
SRR28380565.43	145	1	14678	60	2S99M	16	14232	0	GAAAGGTGTCATGGAGCC[...]  FF:FFFFFFF:FFF[...]	AS:i:-2	XN:i:0	XM:i:0	XO:i:0[...]
SRR28380565.44	99	1	14678	60	4S97M	=	14678	103	GGGAAAGGTGTCATGGAG[...]  FFFFFFFFFFFFFF[...]	AS:i:-4	XN:i:0	XM:i:0	XO:i:0[...]
SRR28380565.44	147	1	14678	60	2S99M	=	14678	-103	GAAAGGTGTCATGGAGCC[...]  F:,FFFFF,F::FF[...]	AS:i:-5	XN:i:0	XM:i:1	XO:i:0[...]
SRR28380565.47	99	1	14630	60	101M	=	14693	164	TGGCTGTGTCCATGTCAG[...]  FFFFFFFFFFFFFF[...]	AS:i:-5	ZS:i:-5	XN:i:0	XM:i:1[...]
SRR28380565.47	147	1	14693	60	101M	=	14630	-164	CCCCTACGATTCCCAGTC[...]  FFFFFFFFF:FFFFF[...]	AS:i:0	ZS:i:-7	XN:i:0	XM:i:0[...]
```
Como podemos ver, al tratarse de una submuestra en la cual las lecturas han sido seleccionadas aleatoriamente, observamos que el número especificado después del ID de las lecturas (SRR28380565) no está ordenado, sino que representa pares de lecturas aleatorias (3, 11, 43, 44, 47...).  

| tabla con los campos |
| ---|
|...| 
|...|

Para conocer el significado de cada FLAG, podemos usar el siguiente [enlace](https://broadinstitute.github.io/picard/explain-flags.html)  

Una vez realizado el alineamiento, se emplea el archivo de anotaciones y la herramienta _RseQC_ para derminar la direccionalidad de las lecturas. El archivo de anotaciones a emplear lo descargamos previamente en formato GTF y activamos la opción para incluir los identificadores de los transcritos (_transcript_id_). Sin embargo, para emplear a herramienta _RseQC_, se necesita el archivo de anotaciones en formato BED, por lo que primero convertimos el archivo de anotaciones al nuevo formato con el paquete de herramientas _Bedops_.  

```console
# Cambio del directorio de trabajo al directorio con el archivo de anotaciones
cd ~/RNAseq_analysis/Data/5_Annotation
# Empleo del script convert2bed del paquete bedops
convert2bed --input=gtf < Homo_sapiens.gtf > Homo_sapiens.bed
```
> NOTA  
> Es importante que el archivo de anotaciones en formato GTF contenga los idenficadores de los transcritos ya que si no el programa de conversión no funciona. Este campo se puede adicionar manualmente, o bien, se puede descargar directamente el archivo de anotaciones como describimos previamente.  

Una vez el archivo de anotaciones presenta el formato correcto, se emplea el resultado del alineamiento de la submuestra (`subsampled_alignment.bam`) y la herramienta _infer_experiment.py_ de _RseQC_.  En este caso, la herramienta infer_experiment.py selecciona unos pocos miles de lecturas del archivo de alineamiento BAM/SAM y representa la proporción de la dirección de las hebras en el experimento.   
```console
cd ~/RNAseq_analysis/Data/2_Infer_strandedness
infer_experiment.py -r ../5_Annotation/Homo_sapiens.bed  -i subsampled_alignment.bam | tee infer_strand.txt
```
> NOTA  
> `-r` : archivo de anotaciones de referencia en formato BED.  
> `-i` : resultado del archivo de alineammiento en formato SAM o BAM.   
> La salida del comando se concatena con el comando `tee` que permite mostrar los resultados en pantalla y a la vez generar un nuevo archivo `infer_strand.txt`con la información mostrada.  
  
Como resultado se pueden dar diferentes escenarios según la forma de preparación de la librería: lecturas pareadas o de extremo único, lecturas con información de hebra específica o no, y a su vez, dentro de los experimentos de hebra específica, se pueden diferenciar:  
* Librerias sentido: Lecturas *forward* o R1 situadas en la misma dirección que el transcrito.  
* Librerias antisentido: Lecturas *reverse* o R2 situadas en la misma dirección que el transcrito.  

![image](https://github.com/user-attachments/assets/6a04a794-3a5c-46ec-83aa-4a3a5b83413b)
![image](https://github.com/user-attachments/assets/43ba8b94-eae8-454c-810e-326a7d8d7da3)


infer_experiment.py samples a few hundred thousand reads from your bam/sam and tells you the portion that would be explained by each type of strandedness, e.g  
[Stranded or non-stranded reads](https://eclipsebio.com/eblogs/stranded-libraries/)  
![image](https://github.com/user-attachments/assets/fc97efa9-a336-4203-b60d-3a4602b8c204)  

Según los resultados obtenidos, se deben especificar diferentes paramétros en las herramientas posteriores.   

| Tipo de Librería | Infer experiment | HISAT2 | htseq-count | 
|---------|---------|----------|----------|
Paired-End (PE) - SF (sense Forward) | 1++,1–,2+-,2-+ | Second Strand F/FR | yes |
PE-SR (Sense-Reverse) | 1+-,1-+,2++,2– | First Strand R/RF | reverse |
Single-End (SE) - SF |	+,– | Second Strand F/FR | yes
SE - SR |	+-,-+ |	First Strand R/RF |	reverse
PE, SE - U |	undecided |	default	| no  




Resultados de la ejecución de la herramienta `infer_experiment.py` para la submuestra SRR28380565:  
```console
This is PairEnd Data
Fraction of reads failed to determine: 0.3406
Fraction of reads explained by "1++,1--,2+-,2-+": 0.0185
Fraction of reads explained by "1+-,1-+,2++,2--": 0.6409
```

Como resultado se obtuvo un experimento de tipo **LIBRERIA ANTISENTIDO**   

### 2.1 Control de calidad, recorte de adaptadores y extremos de mala calidad  
FastQC es una herramienta que permite hacer un control de calidad de las lecturas generando un reporte gráfico con diferentes aspectos del alineamiento. Genera una salida html y un archivo .zip con los datos de calidad.  
Para lanzar la herramienta, primero nos situamos en el directorio con las lecturas crudas y creamos una nueva carpeta para almacenar los resultados de FastQC.  
```console
cd ~/RNAseq_analysis/Data/1_Raw
mkdir qc_raw_reads
fastqc -o qc_raw_reads *.fastq.gz
```
> NOTA  
> La opción `-o` permite indicar el directorio de salido dónde almacenar los reportes generados.
  
Al tratarse de muchas muestras, la herramienta MultiQC permite aunar los reportes generados por FastQC y comparar los resultados para todas las muestras en un mismo informe html.  
```console
cd ~/RNAseq_analysis/Data/1_Raw/qc_raw_reads
multiqc .
```
> NOTA  
> Con el (`.`) indicamos el  directorio actual dónde MultiQC busca los archivos resultantes de FastQC o de otros programas para aunarlos y generar un informe `multiqc_report.html`. Además, genera una carpeta `multiqc_data` con archivos con información complementaria del análisis.  

**Recorte de adaptadores y de bases anotadas de mal calidad**
Seguidamente se lleva a cabo el recorte de adaptadores y filtrado por calidad, con el fin de eliminar artefactos técnicos. Para ello, empleamos la herramienta TrimGalore y procesamos las lecturas crudas almacenadas en la ruta (`~/RNAseq_analysis/Data/1_Raw`). Los resultaos los almacenamos en el direcctorio (`3_Processed`).    
  
```console
cd ~/RNAseq_analysis/Data/1_Raw
trim_galore --paired {sample}_1.fastq.gz {sample}_2.fastq.gz -o ../3_Processed/
```
> NOTA  
> Con la opción `--paired` indicamos que se trata de lecturas pareadas y seguidamente indicamos los archivos con las lecturas _forward_ y _reverse_.  

El programa TrimGalore recorta los extremos de las lecturas en función a su calidad (Phred Q<20), detecta adaptadores y filtra las lecturas con longitud menor a 20pb. Al tratarse de lecturas pareadas, el programa procesa primero cada archivo con las lecturas _forward_ y _reverse_ de forma independiente generando archivos intermedios `*_trimmed.fq.gz` y, una vez completado el recorte, se genera una etapa de validación de los archivos intermedios que serán eliminados, generando los archivos finales validados `*_val_1.fq.gz` y `*_val_2.fq.gz`.  Además, para cada archivo FASTQ, se genera un reporte `trimming_report.txt`con los resultados del recorte y filtrado.  

En cuanto al procesamiento, TrimGalore primero recorta las bases con mala calidad del extremo 3' eliminado las partes de las lecturas con poca calidad (Phred-score < 20). Seguidamente, la herramienta Cutadapt contenida en TrimGalore, encuentra y elimina las secuencias de adaptadores de los extremos 3' de  las lecturas y, finalmente, aquellas lecturas que tras el recorte de bases y adaptadores contengan una longitud inferior a 20pb por defecto son eliminadas.  

Se puede encontrar más información del funcionamiento de TrimGalore en la página de sus creadores, pinchando [aquí](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md).  

Para procesar todos los archivos `.fastq.gz` con las lecturas crudas simultaneamente, se puede emplear el script siguiente, denominado `trimgalore.sh` para ejecutar el comando anterior sobre varios archivos.  
```console
#! /usr/bin/bash
# Bash version 4.2.46(2)
# Date 2025
# Project Name: RNAseq_analysis

# Bash script for trimgalore; quality and adapter trimming for all .fastq.gz files
# Storage of results in the path ~/RNAseq_analysis/Data/3_Processed/

cd ~/RNAseq_analysis/Data/1_Raw

SAMPLES="SRR28380566 SRR28380565 SRR28380570 SRR28380572 SRR28380573 SRR28380568 SRR28380580 SRR28380582 SRR28380584 SRR28380586 SRR28380588 SRR28380589"

for SAMPLE in $SAMPLES; do
	trim_galore --paired ${SAMPLE}_1.fastq.gz ${SAMPLE}_2.fastq.gz -o ../3_Processed/
done
```
El script `trimgalore.sh` se almacena en la ruta  ~/RNAseq_analysis/Code/ y para ejecutarlo se usa el comando `. trimgalore.sh`.  
  
Para observar los nuevos resultados, se emplea FastQC y MultiQC nuevamente sobre las lecturas filtradas, de forma que podemos comparar los resultados antes y después del empleo de la herramienta Trimgalore.  
``` console
cd ~/RNAseq_analysis/Data/3_Processed
mkdir qc_processed_reads
fastqc -o qc_processed_reads *.fq.gz
cd ./qc_processed_reads
mutliqc .
```
Tras el filtrado nos quedan los siguientes resultados:   
![image](https://github.com/user-attachments/assets/ca03a63f-90e0-43ce-84ce-af3c6ac2e049)
  
  
### 2.2 Alineamiento de las lecturas procesadas contra el genoma de referencia con HISAT2
Once the quality of the data is confirmed, we need to convert those millions of reads per sample into the gene- or transcript-level quantification. This would need the assignment of reads to genes or transcripts.  
![image](https://github.com/user-attachments/assets/6f42b6c5-30d6-41b0-b99a-8e57c317e667)  

**2.2.1 Preparación del genoma de referencia**  
Como explicamos previamente, el genoma de referencia indexado se descarga desde el repositorio de _HISAT2_ y, una vez tenemos el genoma en nuestra computadora, se lleva a cabo el alineamiento.  

**2.2.2 Alineamiento de las lecturas**  
Hoy en día, existen diferentes alineadores o mapeadores que se pueden emplear. Entre los más conocidos se encuentran HISAT2 y STAR. Algunas de las diferencias entre ellos son la cantidad  de recursos necesarios, siendo mennos en el caso de HISAT2; o la precisión de los resultados,  siendo más precisos en el caso de STAR. 
En el presente trabajo, se emplea el alineador HISAT2 para mapear las lecturas procesadas al genoma de referencia indexado.  

Para llevar a cabo el alineamiento de las lecturas pareadas con HISAT2, empleamos el siguiente comando:  
```console
#Paired-end reads
cd ~/RNAseq_analysis/Data/4_Alignment
hisat2 -k1 --summary-file {sample}.summary.txt  --rna-strandness {STRING} 
-x (./Reference_genome/grch38/genome) \
-1 {SAMPLE_1_val_1.fq.gz} -2 {SAMPLE_2_val_2.fq.gz} |\
samtools view -Sbh > {SAMPLE}.bam 
````
> NOTA  
> * `-k`: permite establecer el número máximo de alineamientos permitidos por lectura. El valor 1 impide la existencia de alineamientos múltiples.  
> * `--summary-file`: crea un nuevo archivo con los resultados del alineamiento.  
> * `-x` : ruta relativa y nombre principal de los archivos indexados para el genoma de referencia. El nombre principal incluye el nombre de los archivos indexados sin incluir el final (.1.ht2, .2.ht2, etc).  
> `-1` y `-2`: lecturas a alinear  
> ` --rna-strandness` {STRING} option in HISAT2  sets how reads are expected to align against genes. With this option being used, every read alignment will have an XS attribute tag: '+' means a read belongs to a transcript on '+' strand of genome. '-' means a read belongs to a transcript on '-' strand of genome.  
Most stranded protocols in use these days follow the dUTP-method, where read #2 in a pair has the same orientation as the transcript from which it arose (2++ or 2--). So either `R` or `RF` would typically be appropriate   
Use 'RF' if the first read in the pair corresponds to a transcript on the reverse strand, and the second read corresponds to the forward strand. When you use the `--rna-strandness` option with either 'FR' or 'RF' for paired-end reads, HISAT2 will assign an XS attribute tag to each read alignment, indicating whether the read belongs to a transcript on the '+' (plus) or '-' (minus) strand of the genome.
> `samtools`  creamos un BAM directamente. A BAM file is the binary version of a SAM file
> con la opción -S indicamos el archivo de entrada en formato SAM; con la opción -b, indicamos el formato de salida de tipo BAM; y con la opción -h, indicamos que queremos mantener el encabezado

Here is a bash script for the above HISAT2 command called hisat2.sh that will run all the .fastq.gz files for you simultaneously.

```console
#!/usr/bin/bash

#bash script for hisat2; align all .fastq.gz files to indexed reference genome to generate .bam files

SAMPLES="SRR28380566 SRR28380565 SRR28380570 SRR28380572 SRR28380573 SRR28380568 SRR28380580 SRR28380582 SRR28380584 SRR28380586 SRR28380588 SRR28380589"

cd  ~/RNAseq_analysis/Data/4_Alignment/

for SAMPLE in $SAMPLES; do
	hisat2 -k1 --summary-file ${SAMPLE}_align.summary.txt  --rna-strandness RF \
	-x ./Reference_genome/grch38/genome \
	-1 ../3_Processed/${SAMPLE}_1_val_1.fq.gz -2 ../3_Processed/${SAMPLE}_2_val_2.fq.gz |\
 	samtools view -Sbh > ${SAMPLE}.bam
done
```
Los resultados se guardarán en la carpeta `~/RNAseq_analysis/Data/4_Alignment/`

Posteriormente, We will use the samtools command with the options: ‘sort’ to sort the alignments by the leftmost coordinates
```console
samtools sort {SAMPLE}.bam -o {SAMPLE}.sorted.bam | samtools index
```

Here is a bash script for the above sort command called bam.sh
```console
#!/usr/bin/bash
#bash script for samtools; sort and index the .bam files to obtain .bam.bai files

SAMPLES="SRR28380566 SRR28380565 SRR28380570 SRR28380572 SRR28380573 SRR28380568 SRR28380580 SRR28380582 SRR28380584 SRR28380586 SRR28380588 SRR28380589"
cd  ~/RNAseq_analysis/Data/4_Alignment/

for SAMPLE in $SAMPLES; do
	samtools sort ./${SAMPLE}.bam -o ./${SAMPLE}.sorted.bam | samtools index
done
```

paquete **RseQC**, tenemos distintas funciones para analizar el alineamiento.
Vamos a usar el comando `bam_stat.py`, para obtener un resumen de las estadísticas del mapeo del archivo BAM. En este caso, se determina una calidad de mapeo para cada lectura y se calcula la probabilidad de que esa lectura esté mal posicionada en función de un umbral mínimo. El comando a emplear es el siguiente:
```console
bam_stat.py -i ~/RNAseq_analysis/Data/4_Alignment/{SAMPLE}.sorted.bam
```


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
cd ~RNAseq_analysis/Data/

htseq-count -t exon -i gene_id --stranded=reverse -f bam -r pos 4_Alignment/{SAMPLE}.sorted.bam \
5_Annotation/Homo_sapiens.gtf > ../Results/{SAMPLE}_counts.tsv
```
> NOTA Con la opción -t exon indicamos que cuente las lecturas alineadas específicamente contra los exones y, con la opción -i gene_id, indicamos que agrupe las diferentes > lecturas atendiendo al identificador del gen al que pertenecen. Con la opción --stranded=no y -s no, determinamos que las lecturas no provienen de un experimento de hebra > específica y, por tanto, la lectura puede mapear en ambas hebras del genoma.  
> Finalmente, las opciones `–f` y `–r` , las usamos para indicar el formato del archivo a emplear (bam) y cómo están ordenados los alineamientos, en este caso por posición > o coordenadas genómicas (pos).  
> Una vez hemos definido todas las opciones, indicamos la ruta de los archivos BAM y el archivo de anotaciones (GTF), así como la ruta de salida y el nombre del nuevo > archivo que contendrá la información.

```console
SAMPLES="SRR28380566 SRR28380565 SRR28380570 SRR28380572 SRR28380573 SRR28380568 SRR28380580 SRR28380582 SRR28380584 SRR28380586 SRR28380588 SRR28380589"
cd  ~RNAseq_analysis/Data/
for SAMPLE in $SAMPLES; do
	htseq-count -t exon -i gene_id --stranded=reverse -f bam -r pos \
	4_Alignment/${SAMPLE}.sorted.bam \
	5_Annotation/Homo_sapines.gtf > ../Results/${SAMPLE}_counts.tsv
	
done

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
