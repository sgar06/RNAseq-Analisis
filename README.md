# RNAseq-Analisis
#### Escrito por Sonia García Llorens
#### Actualizado el 7 Febrero 2025
### Contenido
1. [Preparación](#1-Preparación)
     * 1.1 Linux y Bash
     * 1.2 Conda e instalación de herramientas para el análisis RNA-seq
     * 1.3 Estructura de directorios
     * 1.4 Obtención del genoma de referencia y el archivo de anotaciones de la especie _Homo sapiens_
     * 1.5 Descarga de los datos de RNA-seq del repositorio público SRA

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

## 1. Preparación 
### 1.1 Linux y Bash  
  
Con el fin de realizar el presente tutorial se emplea un escritorio virutal basado en la nube, denominado Amazon WrokSpaces proporcionado por Amazon Web Services, con el sistema operativo Linux. Para llevar a cabo la ejecucion de los comandos siguientes, se emplea el shell de linea de comandos Bash (Bourne Again Shell).   
  
### 1.2 Conda e instalación de herramientas para el análisis RNA-seq  
  
[Conda](https://docs.conda.io/en/latest/) es un gestor de paquetes y entornos de código abierto. Además, cuenta con grandes repositorios de _software_ o canales, que permiten la instalación y actualización de paquetes o herramientas, así como sus dependencias. Para emplear Conda, descargaremos la distribución de [Miniconda](https://docs.anaconda.com/miniconda/), un instalador gratuito de conda más ligero que incluye solamente conda, Python, y un pequeño número de paquetes o _softwares_.  
  
Para la instalación empleamos el siguiente bloque de comandos.  
```console
# Creación de la carpeta miniconda3 para almacenar el contenido del software
mkdir ~/miniconda3
# Descarga de miniconda a través de la terminal
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
# Instalación de miniconda
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
# Ejecución de conda por defecto en la terminal
~/miniconda3/bin/conda init bash
```
> NOTA  
> El comando `conda init` permite inicializar conda por defecto en la terminal de _bash_. Seguidamente, para que el cambio tenga efecto, se cierra y abre una nueva terminal de _bash_, y nos encontraremos en el entorno 'base' de conda.
  
Una vez tenemos conda instalado por defecto en la terminal, podemos instalar nuevas herramientas o paquetes con el comando `conda install`, de forma que conda verifica si el programa está disponible dentro sus repositorios o canales, y si está  disponible, descargará el _software_ y las dependencias necesarias para su buen funcionamiento. 
Para la correcta instalación de herramientas, debemos especificar a conda los canales o repositorios donde buscar, asi cómo el orden de preferencia. De esta forma, nos aseguramos de la instalación adecuada de los programas y de sus dependendencias con las versiones necesarias. 

Para configurar los canales o repositorios de preferencia y el orden jerárquico de búsqueda, empleamos los siguientes comandos, de forma que el último repositorio añadido tiene mayor prioridad que el anterior.  El orden de preferencia para la búsqueda de paquetes será, por tanto, _conda-forge_, seguido de _bioconda_ y finalmente _default_.  
  
```console
conda config --add channels default
conda config --add channels bioconda
conda config --add channels conda-forge
```
  
Una vez establecidos los parámetros de conda, creamos un entorno específico dónde trabajaremos e instalaremos los paquetes o herramientas para el análisis RNA-seq.  
```console
# Creación de un nuevo entorno denominado genomic_analysis
conda create -n genomic_analysis
# Activación del entorno
conda activate genomic_analysis
```

Posteriormente, instalamos la lista de herramientas mostrada a continuación, necesarias para el análisis transcriptómico dentro de nuestro entorno:  

| Programa | Versión | Comando de instalación | Función | 
|---------|---------|----------|----------|
SRA-Toolkit | 3.2.0 | conda install -c bioconda sra-tools | Descarga de lecturas desde el repositorio SRA. |
SeqKit | 2.9.0 | conda install bioconda::seqkit | Subselección de lecturas en una muestra. |
Bedops | 2.4.41 | conda install -c bioconda bedops | Conversión de archivos. |
RseQC | 5.0.4 | conda install bioconda::rseqc | Evualuación de datos de RNAseq. |
FastQC | 0.12.1 | conda install -c bioconda fastqc | Análisis de calidad de archivos FASTQ. |
MultiQC | 1.27 | conda install -c bioconda multiqc | Creación de reportes. Capacidad para aglomerar los resultados de otras herramientas en un único informe interactivo. |
Trim-Galore | 0.6.10 | conda install -c bioconda trim-galore | Recorte de adaptadores con la herramienta Cutadapt (v5.0) y de bases de mala calidad en archivos FASTQ. |
HISAT2 | 2.2.1 | conda install -c bioconda hisat2 | Mapeador de lecturas. |
Samtools | 1.21  | conda install -c bioconda samtools | Manejo de archivos SAM/BAM. |
HTSeq | 2.0.5 | conda install -c bioconda htseq | Anotación de características genómicas. |

  
### 1.3 Estructura de directorios  
  
Para llevar a cabo el análisis transcriptómico a partir de las lecturas crudas hasta la obtención de genes diferencialmente expresados, empleamos la siguiente estructura de directorios. Siguiendo esta misma estructura, se podrán ejecutar todos los comandos especificados más adelante en el manual.  
  
Dentro del directorio _home_ de nuestro ordenador creamos la siguiente estructura de archivos:   
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
> NOTA  
> La carpeta **Code**, contendrá los _scripts_ necesarios para la ejecución de herramientas.  
> La carpeta **Data**, contendrá las lecturas crudas y procesadas, los archivos necesarios para la deducción de la direccionalidad de la librería y los resultados tras el alineamiento, entre otros.  
> La carpeta **Results**, se empleará para almacenar los recuentos de características genómicas para cada muestra y la matriz de recuentos final.  
  
### 1.4 Obtención del genoma de referencia y el archivo de anotaciones de la especie *Homo sapiens*  
  
* **Descarga del genoma de referencia GRCh38 desde la página web de _HISAT2_**  
  
_HISAT2_ es un programa de alineamiento rápido y eficiente capaz de alinear lecturas resultantes de la secuenciación contra diferentes genomas.  Desde su [repositorio online](http://daehwankimlab.github.io/hisat2/), se cuenta con diversos genomas indexados disponibles para su descarga. La indexación del genoma permite al alineador llevar a cabo su función de manera mucho más dinámica. En el caso de que el genoma de interés no esté indexado, la herramienta _HISAT2_ permite realizar una indexación manual con la función `hisat2-build`, aunque es un proceso lento y costoso. Es por ello, que descargaremos directamente el genoma de referencia humano (GRcH38) indexado desde su página web.  
  
Para descargar el genoma, primero, es necesario acceder a la sección de descargas y, posteriormente, en la sección _Index_ encontramos diferentes _links_ según el genoma de interés, siendo en nuestro caso el genoma perteneciente a la especie _H.sapiens_.   

![image](https://github.com/user-attachments/assets/ddfb3af8-1e12-4a78-be2b-753cf76d5133)  
  
Dentro de la sección _H.sapiens_, encontramos diferentes genomas de referencia según la versión, y para cada uno de ellos, se disponen de distintos _links_ de descarga según el genoma de interés. En este caso, se selecciona el genoma de referencia humano más reciente GRCh38.   
  
![image](https://github.com/user-attachments/assets/2a07a735-3247-4f0a-92fd-9e0da1b8994e)
  
Una vez conocido el _link_ de nuestro genoma de referencia, lo descargamos directamente a través de la terminal y lo almacenamos en una nueva carpeta dentro de la ruta `~RNAseq_analysis/Data/4_Alignment`.   
```console
# Creación de la carpeta Reference_genome y cambio de directorio a la nueva carpeta
mkdir ~/RNAseq_analysis/Data/4_Alignment/Reference_genome && cd $_
# Descarga del genoma de referencia de interés
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
# Descompresión del archivo descargado
tar -xvf grch38_genome.tar.gz
```
> NOTA   
> Con el comando `mkdir` creamos un nuevo directorio para el genoma de referencia dentro de nuestra carpeta 4_Alignment. El comando `cd $_` permite movernos a esta nueva carpeta dónde descargamos el archivo de interés con el comando `wget`.  
> Tras la descompresión del archivo, se genera el directorio `/grch38/` con el genoma de referencia, el script `make_grch38.sh` y los archivos necesarios para la indexación identificados por la palabra `genome` seguidos de la terminación `.1.ht2, .2.ht2`, etc.  
  
* **Descarga del archivo de anotaciones de referencia GRCh38 desde el repositorio _ENSEMBL_**   
  
Para descargar el archivo de anotaciones de referencia de la especie _H.sapiens_, se emplea el repositorio [ENSEMBL](https://www.ensembl.org/Homo_sapiens/Tools/FileChameleon). Además, es importante que el archivo de anotaciones tenga ciertas características específicas necesarias para la correcta ejecución de programas posteriores. Es por ello, que se emplea la herramienta _File Chamaleon_ con el fin de formatear el archivo de anotaciones.  
  
Una vez seleccionado el genoma de interés (GRCh38.p14), se descarga el archivo de anotaciones en formato GTF y se marca la casilla de _transcript_id_ para incluir este campo en el archivo descargado.  
  
![image](https://github.com/user-attachments/assets/5141e4d9-41f4-4362-a1b2-ec78c2e01690)  
  
A continuación descargamos y descomprimimos el archivo de anotaciones en formato GTF en nuestra computadora, y lo almacenamos en la ruta `~/RNAseq_analysis/Data/5_Annotation/`.  
```console
# Descompresión del archivo descargado
gunzip Homo_sapiens.gtf.gz
# Almacenaje del archivo de anotaciones en la ruta de interés
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

  
### 1.5  Descarga de los datos de RNA-seq del repositorio público SRA  
  
Para llevar a cabo el análisis, se emplean lecturas crudas depositadas en el repositorio público _Gene Expression Omnibus_ (GEO) bajo el índice [GSE261866](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261866). En esta página encontramos distinta información relativa al estudio, tales como el número de muestras, el organismo, el tipo de experimento o la plataforma de secuenciación. Además, en el apartado de _Supplementary file_ se encuentra un enlace hacia el repositorio _SRA Run selector_.   
  
![image](https://github.com/user-attachments/assets/d6a8b328-8015-441d-9369-3c9876b0c39f)
  
A partir del repositorio _SRA Run selector_ se pueden descargar tanto las lecturas crudas como los metadatos de las muestras. Además, existe una herramienta implementada por NCBI denominada _SRA-toolkit_ que permite descargar las lecturas directamente en la computadora.  
  
En el presente tutorial, se emplean las muestras de un estudio de RNAseq de neutrófilos aislados de pacientes con lupus eritematoso sistémico y controles sanos. Sin embargo, en lugar de emplear las 27 muestras totales de la investigación, solamente se emplean 6 muestras seleccionadas aleatoriamente para cada grupo.     

**Diseño del experimento**  
Secuenciador: Illumina Novaseq 6000   
Librería: lecturas pareadas con información específica de hebra   
Longitud de las lecturas: 101 pb   
  
**Muestras seleccionadas para el análisis de genes diferencialmente expresados**    
| Nombre de la muestra | ID de las Lecturas | Estado | Número de Lecturas por Millón | Número de Bases | Tamaño |
|---------|---------|----------|----------|----------|----------|
GSM8153253 | SRR28380565 | Sano | 28.5 | 5.7G | 1.8 Gb |
GSM8153252 | SRR28380566 | Sano | 50.8 | 10.3G | 3.1Gb |
GSM8153250 | SRR28380568 | Sano | 55.1 | 11.1 G | 3.5 Gb |
GSM8153248 | SRR28380570 | Sano | 52.7 | 10.6G	| 3.3Gb |
GSM8153246 | SRR28380572 | Sano | 33.4 | 6.8G | 2.1Gb |
GSM8153245 | SRR28380573 | Sano | 33.2 | 6.7G |	2.1Gb |
GSM8153238 | SRR28380580 | LES | 31.1M	| 6.3G | 2.0Gb |
GSM8153236 | SRR28380582 | LES | 47.8M  | 9.6G	| 2.9Gb |
GSM8153234 | SRR28380584 | LES | 39.8M	| 8.0G	| 2.6Gb |
GSM8153232 | SRR28380586 | LES | 39.5M	| 8.0G	| 2.5Gb |
GSM8153230 | SRR28380588 | LES | 29.6M	| 6.0G	| 1.9Gb |
GSM8153229 | SRR28380589 | LES | 31.5M	| 6.4G	| 2.0Gb |
  
Para descargar la lista con los identificadores de las lecturas en formato de texto plano, primero se seleccionan las muestras de interés, se filtran y se descarga la lista con los identificadores o _Accesion list_.   
  
![image](https://github.com/user-attachments/assets/afbcce1d-228e-4994-ba4e-32d585c77471)  
  
Una vez descargado el archivo de texto plano _SRR_Acc_List.txt_ con los identificadores, se guarda en la ruta `~/RNAseq_analysis/Data/Supplementary` y se emplea el comando `fastq-dump` de la herramienta _SRA-toolkit_ desde la terminal para obtener las lecturas crudas en nuestra computadora, que se almacenarán en la carpeta `1_Raw`.  
  
```console
# Cambio de directorio
cd ~/RNAseq_analysis/Data/1_Raw
# Descarga de las lecturas crudas a partir del archivo SRR_Acc_List.txt
xargs -n1 fastq-dump --gzip --split-3 < ../Supplementary/SRR_Acc_List.txt
```  
> NOTA  
> * `xargs -n1` ejecuta el comando `fastq-dump` en las lecturas cuyos nombres figuran en el archivo _SRR_Acc_List.txt_.  Para llegar al archivo, indicamos la ruta relativa desde el directorio dónde nos encontramos.  
> * `fastq-dump` permite la descarga de lecturas desde el repositorio _SRA_ en formato FASTQ y de forma comprimida con la opcción `--gzip`. Además, al tratarse de lecturas pareadas, la opción `--split-3` separa las lecturas _forward_ y _reverse_ en diferentes archivos, y en el caso de que alguna lectura no esté pareada, estas se anotan en un tercer archivo que contiene las lecturas desparejadas.  
   
El **formato FASTQ** contiene la secuencia de las lecturas, así como la puntuación de calidad asociada para cada una de las bases anotadas. Estas puntuaciones de calidad se conocen como "_phred quality score"_ y surgieron en los años 90 para la secuenciación de Sanger, aunque posteriormente se expandieron para su uso en las tecnologías de Secuenciación de Nueva Generación o NGS.    
  
Por ejemplo, si nos fijamos en uno de los archivos **FASTQ** correspondiente a las lecturas obtenidas con el identificador SRR28380565 (muestra GSM8153253), y analizamos la estructura de sus primeras 4 líneas, observamos lo siguiente:  
```console
# Comando empleado para visualizar las 4 primeras líneas
zcat SRR28380565_1.fastq.gz | head -n 4
# Resultado
@SRR28380565.1 A00177:787:HG2KYDMXY:1:1488:15980:33927 length=101
CCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACACTAACCCTAACCCTAACCCTAACCCTAACCCTAACACAAAACAAAAGCAAAACACATA
+SRR28380565.1 A00177:787:HG2KYDMXY:1:1488:15980:33927 length=101
FF:FFFF::FF,FFFFFF,F:FFF,F,FFFF:FFFFFFFFFF,F,FF,FF,FFF,F,FFF:FFFFF,F,FF:FF:FFF::,::,,FFFF,F,:F:,::,,F
```
  
Los resultados de la secuenciación para cada lectura vienen dados por 4 líneas: la primera, representa el símbolo @ seguido del identificador de la lectura; la segunda, contiene la secuencia de la lectura; la tercera, comienza por el símbolo '+' seguido del identificador nuevamente; y la cuarta, contiene la calidad asociada a cada base anotada para dicha lectura.  
  
## 2 Procesamiento de los datos RNA-seq  
  
### 2.1 Determinación de la direccionalidad de las lecturas  
  
Para determinar si las lecturas RNA-seq contienen información de hebra específica, se emplea inicialmente una sola muestra y se realiza una selección aleatoria de entorno al 10% de lecturas totales incluyendo ambos extremos gracias a la herramienta **_Seqkit_**, y los resultados se almancenan en la ruta `~/RNAseq_analysis/Data/2_Infer_strandedness`.  
  
En este paso solamente nos interesa determinar el tipo de librería empleada en el estudio y la especificidad de hebra, por lo que se realiza un selección de lecturas a partir de una de las muestras con el fin de aumentar la eficiencia y rapidez en los pasos posteriores. Con la herramienta `seqkit` se seleccionan de forma aleatoria las lecturas tanto _forward_ como _reverse_ contenidas en los archivos `{sample}_1.fastq.gz` y `{sample}_2.fastq.gz`, respectivamente. Además, es importante que a la hora de realizar el filtrado aleatorio, se incluyan ambos extremos de cada par de lecturas, tanto  *forward* como *reverse* en el orden correcto por lo que se deben emplear los mismos parámetros en ambos procesos.   

Para la realizar la submuestra, se emplean las lecturas correspondientes al identificador SRR28380565.  
```console
# Cambio de directorio
cd ~/RNAseq_analysis/Data/2_Infer_strandedness
# Creación de submuestras
seqkit sample -p 0.1 -s 100 SRR28380565_1.fastq.gz -o subsampled_SRR28380565_1.fastq.gz
seqkit sample -p 0.1 -s 100 SRR28380565_2.fastq.gz -o subsampled_SRR28380565_2.fastq.gz
```
> NOTA  
> * `-p` permite determinar la proporción de lecturas a seleccionar. En nuestro caso, se selecciona el 10% de lecturas totales.   
> * `-s` permite determinar una semilla aleatoria específica para generar resultados idénticos. En este caso, se especifica el valor de 100, permitiendo que las lecturas aleatorias sean las mismas en ambos archivos.
  
Posteriomente, para comprobar el número de pares de lecturas seleccionadas en cada submuestra, se puede especificar el siguiente comando:  
```console
# Comando para determinar el número de lecturas totales en la submuestra 1 con las lecturas forward
zcat subsampled_SRR28380565_1.fastq.gz | grep -c @SRR
# Resultado
2846460
# Comando para determinar el número de lecturas totales en la submuestra 2 con las lecturas reverse
zcat subsampled_SRR28380565_2.fastq.gz | grep -c @SRR
# Resultado
2846460
```
  
En nuestro caso, como la submuestra se obtiene a partir del identificador SRR28380565, cuyo contenido de lecturas iniciales era de 28.5 millones, al llevar a cabo la selección del 10% de lecturas, obtenemos archivos con 2.8 millones de lecturas.  
  
Una vez preparada la submuestra, se lleva a cabo un alineamiento rápido contra el genoma de referencia usando el mapeador _HISAT2_ y los resultados del alineamiento se comparan respecto al archivo de anotación de referencia mediante la herramienta `infer_experiment.py` del paquete _RseQC_ con el fin de determinar la direccionalidad de la librería.  
  
Para realizar el mapeo se emplea la herramienta **_HISAT2_**, las lecturas de la submuestra, y el genoma de referencia humano indexado descargado previamente (_grch38_genome.tar.gz_), y como resultado, se obtiene un archivo con la información del alineamiento de las lecturas en formato BAM.  
```console
# Alineamiento de las lecturas filtradas con HISAT2
hisat2 -k1 -x ../4_Alignment/Reference_genome/grch38/genome \
-1 subsampled_SRR28380565_1.fastq.gz -2 subsampled_SRR28380565_2.fastq.gz |\
samtools view -Sbh > subsampled_alignment.bam
```  
> NOTA   
> * `-k` permite determinar el número de ubicaciones diferentes en las que una misma lectura puede alinearse dentro del genoma de referencia. En nuestro caso especificamos el valor 1, para impedir los alineamientos múltiples.  
> * `-x` indica la ruta de directorios hasta llegar al genoma de referencia indexado, así como la palabra común de todos los archivos de indexación sin la extensión final (.1.ht2, .2.ht2, etc).  
> * `-1` y `-2` permite indicar las lecturas filtradas _forward_ y _reverse_ de la submuestra.  
> * Como resultado del alineamiento, el comando `hisat2` genera un archivo SAM. Sin embargo, como se trata de un archivo muy pesado, la salida del comando `hisat2` se concatena directamente a la herramienta `samtools` a través de una tubería `|` con el fin evitar la creación de archivos intermedios. De esta forma el comando `samtools view -Sbh > subsampled_alignment.bam`  permite crear directamente el archivo BAM en código binario, ocupando menos espacio. La opción `-Sbh` permite indicar el  formato de tipo SAM que se usa como entrada (`S`), el archivo de salida en formato BAM (`b`) y el mantenimiento de la cabecera (`h`).  

Para visualizar el resultado del alineamiento en formato binario (BAM), se puede emplear la herramienta `samtools`.  
```console
samtools view subsampled_alignment.bam | head
```
Como resultado, se muestran las primeras 10 líneas del alineamiento contenidos en el archivo **BAM** en la salida por pantalla de la terminal:  
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
  
Al tratarse de una submuestra, en la cual las lecturas han sido seleccionadas aleatoriamente, observamos que el número especificado después del ID de las lecturas (SRR28380565) no está ordenado, sino que representa pares de lecturas aleatorias (3, 11, 43, 44, 47...).  
  
El **formato BAM** incluye una parte de encabezado opcional, y otra parte con los resultados del alineamiento. La parte del encabezado se diferencia porque incluye el símbolo de '@' al comienzo de cada línea, e incluye informaciones como el nombre de la muestra o el método de alineamiento, entre otras. En este caso, al usar el comando `samtools view` se muestra directamente la sección con el alineamiento, sin mostrar el encabezado.  
  
A su vez, la sección con el alineamiento, se divide en 11 campos obligatorios y, en general, cada línea representa el alineamiento para un segmento dentro del genoma de referencia. Sin embargo, algunas lecturas pueden ocupar varias líneas del alineamiento, cuando se permiten los alineamientos múltiples contra el genoma de referencia.  

| Campo | Descripción |
| ---| ---|
|1.QNAME | Representa un identificador único para cada una de las lecturas. |
|2.FLAG |  Indica un valor de puntuación en base a las características del alineamiento de cada lectura.  |
|3.RNAME | Nombre de la secuencia de referencia dónde mapea la lectura. En nuestro caso, representa el número de cromosoma. Si no hay información del mapeo, se representa el símbolo "*". |
|4.POS | Posición de la primera base de la lectura que mapea en el genoma de referencia. |
|5.MAPQ | Calidad del mapeo. Indica cómo de correcto o incorrecto es el alineamiento. |
|6.CIGAR | Secuencia de números y letras que describen el alineamiento, inserciones y deleciones. |
| 7.RNEXT | Análogo al campo 3, pero con la información relativa al otro par de la lectura en el caso de lecturas pareadas. |
| 8.PNEXT | Análogo al campo 4, pero con la información relativa al otro par de la lectura en el caso de lecturas pareadas. |
| 9.TLEN | Longitud  de la secuencia de referencia dónde mapean las lecturas. |
| 10.SEQ | Secuencia de la lectura actual, corresponde a la misma secuencia especificada en el archivo FASTQ. | 
| 11.QUAL | Calidad de la lectura actual, corresponde a la calidad asociada a cada una de las bases en formato _Phred score_ y también esta especificada en el archivo FASTQ. |
| 12.TAGS | Etiquetas específicas con información adicional del alineamiento o de la lectura. |
   
Además, existen diferentes recursos online que permiten conocer el significado de los distintos identificadores. Por ejemplo, para conocer el significado de los valores del campo FLAG, se puede usar el siguiente [enlace](https://broadinstitute.github.io/picard/explain-flags.html).   
  
Una vez realizado el alineamiento, se emplea el archivo de anotaciones y la herramienta **_RseQC_** para determinar la direccionalidad de las lecturas. El archivo de anotaciones a emplear lo descargamos previamente en formato GTF y activamos la opción para incluir los identificadores de los transcritos (_transcript_id_). Sin embargo, para emplear la herramienta _RseQC_, se necesita el archivo de anotaciones en formato BED, por lo que primero convertimos el archivo de anotaciones al nuevo formato con el paquete de herramientas _Bedops_.  

```console
# Cambio del directorio de trabajo al directorio con el archivo de anotaciones
cd ~/RNAseq_analysis/Data/5_Annotation
# Empleo de la herramienta convert2bed del paquete Bedops
convert2bed --input=gtf < Homo_sapiens.gtf > Homo_sapiens.bed
```
> NOTA   
> Para llevar a cabo la conversión, es importante que el archivo de anotaciones, en formato GTF, contenga los identificadores de los transcritos, ya que si no el programa de conversión no funciona. Este campo se puede adicionar manualmente, o bien, se puede descargar directamente el archivo de anotaciones como describimos previamente.    
  
El **formato BED** representa regiones específicas del genoma y puede tener diferentes elementos, aunque en su forma más básica incluye solamente 3 columnas para indicar el cromosoma, la posición inicial y la posición final de una región concreta.  Sin embargo, pueden contar con hasta 9 campos adicionales y opcionales, tales como los identificadores, los nombres de los elementos, las puntuaciones, la hebra del genoma y el número de exones, entre otros.  
  
```console
head Homo_sapiens.bed
1	11120	11211	ENSG00000290825	.	+	havana_tagene	exon	.	transcript_version "1"; transcript_biotype "lncRNA"; exon_version "1"; gene_id "ENSG00000290825"; transcript_name "DDX11L16-260"; gene_source "havana"; exon_id "ENSE00004248723"; gene_version "2"; gene_name "DDX11L16"; transcript_id "ENST00000832824"; gene_biotype "lncRNA"; exon_number "1"; transcript_source "havana_tagene"
1	11120	14413	ENSG00000290825	.	+	havana_tagene	transcript	.	gene_biotype "lncRNA"; gene_source "havana"; transcript_source "havana_tagene"; transcript_version "1"; transcript_biotype "lncRNA"; gene_id "ENSG00000290825"; gene_version "2"; transcript_id "ENST00000832824"; transcript_name "DDX11L16-260"; gene_name "DDX11L16"
1	11120	24894	ENSG00000290825	.	+	havana	gene	.	gene_name "DDX11L16"; transcript_id "ENSG00000290825"; gene_source "havana"; gene_biotype "lncRNA"; gene_version "2"; gene_id "ENSG00000290825"
1	11124	11211	ENSG00000290825	.	+	havana_tagene	exon	.	gene_biotype "lncRNA"; transcript_source "havana_tagene"; exon_number "1"; gene_version "2"; exon_id "ENSE00004248721"; gene_name "DDX11L16"; transcript_id "ENST00000832825"; gene_source "havana"; exon_version "1"; transcript_biotype "lncRNA"; gene_id "ENSG00000290825"; transcript_version "1"; transcript_name "DDX11L16-261"
1	11124	14405	ENSG00000290825	.	+	havana_tagene	transcript	.	transcript_biotype "lncRNA"; gene_version "2"; gene_id "ENSG00000290825"; transcript_version "1"; gene_name "DDX11L16"; transcript_name "DDX11L16-261"; transcript_id "ENST00000832825"; gene_source "havana"; gene_biotype "lncRNA"; transcript_source "havana_tagene"
1	11409	11671	ENSG00000290825	.	+	havana_tagene	exon	.	exon_version "1"; transcript_biotype "lncRNA"; gene_id "ENSG00000290825"; transcript_version "1"; transcript_name "DDX11L16-262"; gene_source "havana"; exon_id "ENSE00004248726"; gene_version "2"; transcript_id "ENST00000832826"; gene_name "DDX11L16"; gene_biotype "lncRNA"; transcript_source "havana_tagene"; exon_number "1"
1	11409	14413	ENSG00000290825	.	+	havana_tagene	transcript	.	gene_source "havana"; gene_biotype "lncRNA"; transcript_source "havana_tagene"; transcript_biotype "lncRNA"; gene_id "ENSG00000290825"; gene_version "2"; transcript_version "1"; gene_name "DDX11L16"; transcript_name "DDX11L16-262"; transcript_id "ENST00000832826"
1	11410	11671	ENSG00000290825	.	+	havana_tagene	exon	.	transcript_name "DDX11L16-263"; transcript_version "1"; transcript_biotype "lncRNA"; exon_version "1"; gene_id "ENSG00000290825"; gene_source "havana"; gene_name "DDX11L16"; transcript_id "ENST00000832827"; gene_version "2"; exon_id "ENSE00004248710"; exon_number "1"; transcript_source "havana_tagene"; gene_biotype "lncRNA"
1	11410	14413	ENSG00000290825	.	+	havana_tagene	transcript	.	transcript_source "havana_tagene"; gene_source "havana"; gene_biotype "lncRNA"; gene_name "DDX11L16"; transcript_name "DDX11L16-263"; transcript_id "ENST00000832827"; gene_version "2"; transcript_biotype "lncRNA"; gene_id "ENSG00000290825"; transcript_version "1"
1	11425	11671	ENSG00000290825	.	+	havana_tagene	exon	.	transcript_id "ENST00000832828"; gene_name "DDX11L16"; exon_id "ENSE00004248702"; gene_version "2"; exon_number "1"; transcript_source "havana_tagene"; gene_biotype "lncRNA"; transcript_name "DDX11L16-264"; transcript_version "1"; transcript_biotype "lncRNA"; exon_version "1"; gene_id "ENSG00000290825"; tag "gencode_basic"; gene_source "havana"
```
   
Una vez el archivo de anotaciones presenta el formato correcto, se emplea el resultado del alineamiento de la submuestra (`subsampled_alignment.bam`) y la herramienta `infer_experiment.py` de **_RseQC_**.  En este caso, la herramienta `infer_experiment.py` selecciona unos pocos miles de lecturas del archivo de alineamiento BAM/SAM y representa la proporción de la dirección de las hebras en el experimento.   
  
```console
# Cambio de directorio
cd ~/RNAseq_analysis/Data/2_Infer_strandedness
# Deducción de la direccionalidad de las hebras
infer_experiment.py -r ../5_Annotation/Homo_sapiens.bed  -i subsampled_alignment.bam | tee infer_strand.txt
```
> NOTA  
> * `-r` : indica el archivo de anotaciones de referencia en formato BED.  
> * `-i` : indica el archivo resultante tras el alineamiento en formato SAM o BAM.    
> La salida del comando se concatena con el comando `tee` que permite mostrar los resultados en la salida por pantalla, y a su vez generar un nuevo archivo `infer_strand.txt` con la información visualizada.    
   
Como resultado se pueden dar diferentes escenarios según la forma de preparación de la librería: lecturas pareadas o de extremo único, lecturas con información de hebra específica o no, y a su vez, dentro de los experimentos de hebra específica, se pueden diferenciar:  
* Librerias sentido: Lecturas *forward* o R1 situadas en la misma dirección que el transcrito.  
* Librerias antisentido: Lecturas *reverse* o R2 situadas en la misma dirección que el transcrito.

Se puede encontrar más información en la página [Eclipsebio](https://eclipsebio.com/eblogs/stranded-libraries/).  
  
Según los resultados obtenidos, existen diferentes parámetros que se deben especificar en los análisis posteriores.  
  
| Tipo de Librería | Dirección de las lecturas | Paramétros para la herramienta HISAT2 | Parámentros para la herramienta HTseq | 
|---------|---------|----------|----------|
Lecturas pareadas y con Sentido Forward | 1++,1–-,2+-,2-+ | -rna-strandness FR | --stranded=yes |
Lecturas pareadas y con Sentido Reverse o Antisentido | 1+-,1-+,2++,2-- | --rna-strandness RF |  --stranded=reverse |
Lecturas únicas y con sentido Forward |	+,– | -rna-strandness F | --stranded=yes |
Lecturas únicas y con Sentido Reverse o Antisentido |	+-,-+ |	--rna-strandness R | --stranded=reverse |
Lecturas pareas o únicas sin información de heabra | Sin determinar | parámetros por defecto | --stranded=no  |
  
Tras la ejecución de la herramienta `infer_experiment.py`, se obtienen los siguientes resultados para la submuestra SRR28380565, indicando que la librería del experimento se trata de una librería con **lecturas pareadas antisentido** (1+-,1-+,2++,2--).  
```console
This is PairEnd Data
Fraction of reads failed to determine: 0.3406
Fraction of reads explained by "1++,1--,2+-,2-+": 0.0185
Fraction of reads explained by "1+-,1-+,2++,2--": 0.6409
```
Determinar la direccionalidad de las lecturas, es importante  especialmente en genomas complejos con genes solapantes en direcciones opuestas, ya que de esta forma es posible conocer la hebra de ADN con la información codificante para cada transcrito.  
  
### 2.2 Control de calidad, recorte de adaptadores y extremos de baja calidad  
#### 2.2.1 Control de calidad 
Para comprobar la calidad de la secuenciación, se empea la herramienta **FastQC** que permite hacer un control de calidad de las lecturas generando un reporte gráfico con diferentes aspectos del alineamiento, además de generar una salida html y un archivo .zip con los datos de calidad.  
  
Para lanzar la herramienta, primero nos situamos en el directorio con las lecturas crudas y creamos una nueva carpeta para almacenar los resultados generados por FastQC.   
```console
# Cambio de directorio
cd ~/RNAseq_analysis/Data/1_Raw
# Creación de una nueva carpeta
mkdir qc_raw_reads
# Análisis de calidad de las lecturas crudas con FastQC
fastqc -o qc_raw_reads *.fastq.gz
```
> NOTA  
> * La opción `-o` de fastqc permite indicar el directorio de salida dónde almacenar los reportes generados.  
  
Además, al tratarse de varias muestras, la herramienta **MultiQC** permite aunar los reportes generados por FastQC y comparar los resultados para todas las muestras en un mismo informe html.  
```console
# Cambio de directorio
cd ./qc_raw_reads
# Generación de un único reporte con MultiQC
multiqc .
```
> NOTA  
> * Con el símbolo '.' indicamos el directorio actual, donde MultiQC busca los ficheros resultantes de FastQC u otros programas para unirlos y generar tanto un informe `multiqc_report.html`, como una nueva carpeta `multiqc_data` que contiene información complementaria del análisis.  
  
#### 2.2.2 Recorte de adaptadores y de bases anotadas de mala calidad  
A continuación, se lleva a cabo el recorte de adaptadores y filtrado por calidad, con el fin de eliminar artefactos técnicos. Para ello, empleamos la herramienta **TrimGalore** para procesar las lecturas crudas, y las lecturas procesadas las almacenamos en el directorio `~/RNAseq_analysis/Data/3_Processed`.  
  
```console
# Cambio de directorio a la ruta con las lecturas crudas 
cd ~/RNAseq_analysis/Data/1_Raw
# Procesamiento de lecturas con la herramienta TrimGalore
trim_galore --paired {sample}_1.fastq.gz {sample}_2.fastq.gz -o ../3_Processed/
```
> NOTA  
> * La opción `--paired` de TrimGalore permite indicar que se trata de lecturas pareadas, y seguidamente, se indican los archivos con las lecturas _forward_ y _reverse_.  
  
El programa TrimGalore recorta los extremos de las lecturas en función a su calidad (Phred Q<20), detecta adaptadores y filtra las lecturas con longitud menor a 20pb. Al tratarse de lecturas pareadas, el programa procesa primero cada archivo con las lecturas _forward_ y _reverse_ de forma independiente generando archivos intermedios `*_trimmed.fq.gz`, y una vez completado el recorte, se genera una etapa de validación de los archivos intermedios, que serán eliminados, generando los archivos finales validados `*_val_1.fq.gz` y `*_val_2.fq.gz`.  Además, para cada archivo FASTQ, se genera un reporte `trimming_report.txt` con los resultados del recorte y filtrado.  
  
En cuanto al procesamiento, TrimGalore primero recorta las bases con mala calidad del extremo 3' eliminado las partes de las lecturas con baja calidad (Phred-score < 20). Seguidamente, la herramienta Cutadapt, empaquetada dentro del programa TrimGalore, encuentra y elimina las secuencias de adaptadores en los extremos 3' de  las lecturas y, finalmente, aquellas lecturas que tras el recorte de bases y adaptadores contengan una longitud inferior a 20pb por defecto, son eliminadas.  
  
Se puede encontrar más información del funcionamiento de TrimGalore en la página de sus creadores, pinchando [aquí](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md).  
  
En el caso de querer optimizar el proceso o contar con varias muestras, se pueden procesar todos los archivos FASTQ con las lecturas crudas empleando el siguiente script, denominado `trimgalore.sh`, desde el directorio  ~/RNAseq_analysis/Code.  
```console
#! /usr/bin/bash
# Versión de Bash: 4.2.46(2)
# Fecha: 2025
# Nombre del proyecto: RNAseq_analysis

# Script de bash para la herramienta TrimGalore; recorte de bases de mala calidad y adaptadores de las lecturas crudas
# Almacenamiento de los resultados en el directorio ~/RNAseq_analysis/Data/3_Processed/

cd ~/RNAseq_analysis/Data/1_Raw

SAMPLES="SRR28380565 SRR28380566 SRR28380568 SRR28380570 SRR28380572 SRR28380573
	SRR28380580 SRR28380582	SRR28380584 SRR28380586	SRR28380588 SRR28380589"

for SAMPLE in $SAMPLES; do
	trim_galore --paired ${SAMPLE}_1.fastq.gz ${SAMPLE}_2.fastq.gz -o ../3_Processed/
done
```
> NOTA  
> * Dentro del script, se crea una nueva variable denominada `SAMPLES` con la lista de los identificadores de todas las lecturas a procesar, y con el bucle `for` se recorre cada elemento de la lista y se ejecuta la herramienta Trimgalore para cada uno de ellos.    
  
Todos los script `.sh` generados se almacenan en la ruta  `~/RNAseq_analysis/Code/`, y para ejecutarlo se usa el comando `. {nombre_script}.sh` desde el directorio `Code`.  
    
Para analizar los resultados tras el procesamiento de las lecturas, se realiza de nuevo un análisis de calidad con FastQC y MultiQC, de forma que comparamos los resultados antes y después del empleo de la herramienta TrimGalore.  
``` console
# Cambio de directorio
cd ~/RNAseq_analysis/Data/3_Processed
# Creación de una nueva carpeta
mkdir qc_processed_reads
# Análisis de calidad de las lecturas procesadas con FastQC
fastqc -o qc_processed_reads *.fq.gz
cd ./qc_processed_reads
# Generación de un único reporte con MultiQC
mutliqc .
```
Tras el procesamiento y filtrado de las lecturas se obtienen los siguientes resultados:   
![image](https://github.com/user-attachments/assets/307556a3-8a03-4141-a91c-c5f10e02e243)  

Además, los informes generados por la herramienta MultiQC, están disponibles para su visualización en la carpeta **MultiQC** del presente manual.  
      
### 2.3 Alineamiento de las lecturas procesadas contra el genoma de referencia con HISAT2   
Una vez que se tienen las lecturas procesadas, se alinean contra el genoma de referencia con el fin de asignar la posición y coordenadas cromosómicas para cada una de ellas.  
  
#### 2.3.1 Preparación del genoma de referencia  
Como explicamos previamente, el genoma de referencia indexado se descarga desde el repositorio de _HISAT2_ y, una vez tenemos el genoma en nuestra computadora, se lleva a cabo el alineamiento.  
  
#### 2.3.2 Alineamiento de las lecturas
Hoy en día, existen diferentes alineadores o mapeadores que se pueden emplear. Entre los más conocidos se encuentran TopHat2, HISAT2 o STAR, y algunas de las diferencias entre ellos son la cantidad de recursos necesarios o la precisión de los resultados. TopHat2 es el alineador menos eficiente de los tres, ya que puede tomar varios días en procesar los experimentos de ARNseq, mientras que HIAST2 y STAR son alineadores mucho más eficaces. Sin embargo, HISAT2 es un alineador menos preciso en comparación con STAR, aunque es adecuado cuando hay menos recursos disponibles ya que necesita menos menoria RAM disponible.  
  
En el presente manual, se emplea el alineador HISAT2, y para llevar a cabo el alineamiento de las lecturas procesadas, se utiliza el siguiente comando:  
```console
# Cambio de directorio
cd ~/RNAseq_analysis/Data/4_Alignment
# Alineamiento de las lecturas con HISAT2
hisat2 -k1 --summary-file {sample}.summary.txt  --rna-strandness {STRING} 
-x ./Reference_genome/grch38/genome \
-1 ../3_Processed/{sample}_1_val_1.fq.gz -2 ../3_Processed/{sample}_2_val_2.fq.gz |\
samtools view -Sbh > {sample}.bam 
````
> NOTA  
> * `-k`: permite establecer el número máximo de alineamientos permitidos por lectura. El valor 1 impide la existencia de alineamientos múltiples.  
> * `--summary-file`: crea un nuevo archivo de texto plano con los resultados del alineamiento.  
> * ` --rna-strandness {STRING}`: permite indicar cómo se espera que alineen las lecturas en función de la direccionalidad de la librería. En nuestro caso, al tratarse de una librería pareada antisentido, se usa la opción `RF` para indicar que la lectura 2 o _reverse_ tiene la misma orientación que el transcrito a partir del cual se genera (1+-,1-+,2++,2--).  
> * `-x` : ruta relativa y nombre común de los archivos indexados para el genoma de referencia. El nombre común incluye el nombre de los archivos indexados sin incluir la extensión final (.1.ht2, .2.ht2, .3.ht2 etc).  
> * `-1` y `-2`: lecturas _forward_ y _reverse_, respectivamente.  
> *`| samtools view -Sbh > {sample}.bam`: permite crear el archivo BAM con los resultados del alineamiento de forma directa, sin necesidad de crear el archivo SAM intermedio.  
  
Para mapear todos los archivos FASTQ con las lecturas procesadas, se puede emplear el script siguiente, denominado `hisat2_alignment.sh` para ejecutar el comando anterior sobre varios archivos.  
  
```console
#! /usr/bin/bash
# Versión de Bash: 4.2.46(2)
# Fecha: 2025
# Nombre del proyecto: RNAseq_analysis

# Script de bash para el programa hisat2; alineamiento de las lecturas procesadas al genoma de referencia indexado
# Creación de archivos BAM como resultado

SAMPLES="SRR28380565 SRR28380566 SRR28380568 SRR28380570 SRR28380572 SRR28380573
	SRR28380580 SRR28380582	SRR28380584 SRR28380586	SRR28380588 SRR28380589"

cd  ~/RNAseq_analysis/Data/4_Alignment/

for SAMPLE in $SAMPLES; do
	hisat2 -k1 --summary-file ${SAMPLE}_alignment.summary.txt  --rna-strandness RF \
	-x ./Reference_genome/grch38/genome \
	-1 ../3_Processed/${SAMPLE}_1_val_1.fq.gz -2 ../3_Processed/${SAMPLE}_2_val_2.fq.gz |\
 	samtools view -Sbh > ${SAMPLE}.bam
done
```  
Los archivos resultantes con el alineamiento en formato BAM se almacenan en la carpeta `~/RNAseq_analysis/Data/4_Alignment/`, junto con los archivos de texto plano `${SAMPLE}_alignment.summary.txt`, que contienen la información resumen derivada del mapeo. Para ilustrar esta información, podemos observar alguno de los archivos de texto plano generados, tal como `SRR28380565_alignment.summary.txt`, de forma que se muestran los siguientes resultados:
  
```console
28427141 reads; of these:      				   # Lecturas totales (x2= 56854282)
  28427141 (100.00%) were paired; of these:
    1123549 (3.95%) aligned concordantly 0 times           # Lecturas no conformes (x2 = 2247098)
    27303592 (96.05%) aligned concordantly exactly 1 time  # Lecuturas con mapeo conforme y con una ubicación (x2 = 54607184)
    0 (0.00%) aligned concordantly >1 times		   # Lecturas con mapeo conforme y con más de una ubicación
    ----
    1123549 pairs aligned concordantly 0 times; of these:
      88122 (7.84%) aligned discordantly 1 time    	   # Lecturas mapeadas pero discordantes (x2 = 176244)
    ----
    1035427 pairs aligned 0 times concordantly or discordantly; of these:
      2070854 mates make up the pairs; of these:   
        1094218 (52.84%) aligned 0 times	   	   # Lecturas sin mapear
        976636 (47.16%) aligned exactly 1 time     	   # Lecturas con solo un extremo alineado 
        0 (0.00%) aligned >1 times
98.08% overall alignment rate
```
Como se puede observar, se genera un resumen del alineamiento con las lecturas mapeadas y no mapeadas. Los comentarios anexados a la derecha de los parámetros describen qué representan cada uno. En el caso de las lecturas mapeadas de manera discordante o no conforme (7.84%), hacen referencia a pares de lecturas con mapeos únicos pero que no cumplen las expectitivas para un alineamiento de lecturas pareadas, como puede ser la orientación esperada de las lecturas o el rango de distancia esperado entre ambas lecturas.  
El total de lecturas alineadas se calcula a partir de la suma de las lecturas con mapeo conforme y una ubicación (54607184), lecturas con mapeo discordantes (176244) y lecturas con un solo extremo alineado (976636), haciendo un total de 55.760.064 lecturas alineadas y 98.08% de tasa total de alineamiento.  
   
Una vez realizado el alineamiento y obtenidos los archivos BAM, se emplea la herramienta **Samtools** para ordenar las lecturas según sus coordenadas genómicas.  
  
```console
# Ordenación de los archivos BAM con la herramienta samtools
samtools sort {sample}.bam -o {sample}.sorted.bam 
```
  
Para ordenar todos los archivos resultantes del alineamiento en formato BAM, se puede emplear el siguiente script `bam_order.sh` desde el directorio  ~/RNAseq_analysis/Code.   
```console
#! /usr/bin/bash
# Versión de Bash: 4.2.46(2)
# Fecha: 2025
# Nombre del proyecto: RNAseq_analysis

# Script de bash para el programa samtools; ordenación por posición genómica de los archivos .bam

SAMPLES="SRR28380565 SRR28380566 SRR28380568 SRR28380570 SRR28380572 SRR28380573
	SRR28380580 SRR28380582	SRR28380584 SRR28380586	SRR28380588 SRR28380589"

cd  ~/RNAseq_analysis/Data/4_Alignment/

for SAMPLE in $SAMPLES; do
	samtools sort ./${SAMPLE}.bam -o ./${SAMPLE}.sorted.bam
done
```
   
Finalmente, se lleva a cabo la indexación de los archivos BAM ordenados y almacenados en el directorio `~/RNAseq_analysis/Data/4_Alignment/`.   
```console
# Indexación de los archivos BAM ordenados
find . -name "*sorted.bam" | xargs -n1 samtools index
```
> NOTA  
> * El comando `find` busca y lista en el directorio actual ('.') todos los archivos en cuyo nombre se contiene la terminación `*sorted.bam`, y el resultado se concatena al comando `xargs -n1 samtools index`, que para cada archivo identificado, ejecuta el comando de indexación.  
> * Como resultado, el comando `samtools index` genera un archivo `sample.sorted.bam.bai`, dónde la terminación BAI hace referencia al archivo índice para cada BAM.  
  
La indexación de los archivos BAM, genera un archivo índice complementario con un tamaño mucho menor, que actúa como una tabla de contenidos y permite identificar dónde se localizan las lecturas. Esto es importante en la ejecución de algunos programas, ya que permite localizar partes específicas del archivo BAM de manera más rápida.  
  
#### 2.3.3 Calidad del alineamiento   
Una vez alineadas las lecturas procesadas contra el genoma de referencia, existen herramientas para comprobar la calidad de dicho alineamiento. Para ello, vamos a emplear el programa _RseQC_, dentro del cual existen diferentes módulos y scripts creados para analizar el alineamiento.  
En concreto, vamos a emplear el script `bam_stat.py`, para obtener un resumen de las estadísticas del mapeo del archivo BAM. Primero, se determina la calidad del mapeo para cada lectura y seguidamente se calcula la probabilidad de que esa lectura esté mal posicionada en función de un umbral mínimo.   
  
Para emplear la herramienta usamos el siguiente comando:  
```console
# Creación de un nuevo directorio dentro de la carpeta 4_Alignment
mkdir ./statistics
# Calidad del alineamiento con RseQC
bam_stat.py -i ./{sample}.sorted.bam > statistics/{sample}.bamstats.txt
```
  
Para comprobar la calidad del alineamiento para cada uno de los archivos BAM resultantes, se puede emplear el siguiente script `bam_stats.sh` desde el directorio  ~/RNAseq_analysis/Code.  
```console
#! /usr/bin/bash
# Versión de Bash: 4.2.46(2)
# Fecha: 2025
# Nombre del proyecto: RNAseq_analysis

# Script de bash para el programa RseQC; calidad del alineamiento con la herramienta bam_stats.py

SAMPLES="SRR28380565 SRR28380566 SRR28380568 SRR28380570 SRR28380572 SRR28380573
	SRR28380580 SRR28380582	SRR28380584 SRR28380586	SRR28380588 SRR28380589"

cd ~/RNAseq_analysis/Data/4_Alignment/

for SAMPLE in $SAMPLES; do
	bam_stat.py -i ${SAMPLE}.sorted.bam > ./statistics/${SAMPLE}.bamstats.txt
done
```
Como resultado se obtiene un archivo de texto plano con diferentes estadísticas del alineamiento. Por ejemplo, en el caso del archivo SRR28380565 se obtienen los siguientes resultados:   
```console
#==================================================
#All numbers are READ count
#==================================================

Total records:                          56854282  # Lecturas totales

QC failed:                              0
Optical/PCR duplicate:                  0
Non primary hits                        0
Unmapped reads:                         1094218   # Lecturas no mapeadas
mapq < mapq_cut (non-unique):           1663991   # Lecturas con mapeo de baja calidad

mapq >= mapq_cut (unique):              54096073  # Lecturas con mapeo de alta calidad
Read-1:                                 27138442  # Lecturas foward
Read-2:                                 26957631  # Lecturas reverse
Reads map to '+':                       27053055 
Reads map to '-':                       27043018  
Non-splice reads:                       40899546  # Lecturas mapeadas en un único exón
Splice reads:                           13196527  # Lecturas mapeadas en diferentes exones
Reads mapped in proper pairs:           53037578
Proper-paired reads map to different chrom:0
```

Dentro del resumen anterior, podemos observar diferentes estadísticas del mapeo. Además, cabe destacar que el genoma de referencia humano contiene intrones y exones dentro de sus genes, y que los transcritos de ARNm se pueden obtener a partir del empalme de varios exones, es por ello que es importante que el alineador empleado sea capaz de tener en cuenta la presencia de intrones para alinear las lecturas correctamente.   
  
Finalmente, los resultados obtenidos tras el alineamiento de todas las muestras son los siguientes:   
  
![image](https://github.com/user-attachments/assets/0b49ef3a-0d83-4844-8ac0-034608e37af9)
  
  
### 2.4 Identificación y recuento de features o características  

#### 2.4.1 Preparación del archivo de anotaciones GTF

Como explicamos al principio del manual, el archivo de anotaciones de la especie humana (GRCh38.p14) en formato GTF se descargó desde el repositorio [ENSEMBL](https://www.ensembl.org/Homo_sapiens/Tools/FileChameleon), y su estructura es la siguiente:   
   
```console
#!genome-build GRCh38.p14
#!genome-version GRCh38
#!genome-date 2013-12
#!genome-build-accession GCA_000001405.29
#!genebuild-last-updated 2024-07
1	ensembl_havana	gene		3069168	3438621	.	+	.	gene_id "ENSG00000142611"; gene_version "17"; gene_biotype "protein_coding"; gene_source "ensembl_havana"; transcript_id "ENSG00000142611"; gene_name "PRDM16"
1	havana		transcript	3069168	3434342	.	+	.	tag "gencode_primary"; transcript_support_level "5"; gene_source "ensembl_havana"; transcript_name "PRDM16-206"; transcript_biotype "protein_coding"; gene_id "ENSG00000142611"; transcript_version "5"; transcript_source "havana"; gene_biotype "protein_coding"; transcript_id "ENST00000511072"; gene_name "PRDM16"; gene_version "17"
1	havana		exon		3069168	3069296	.	+	.	gene_biotype "protein_coding"; transcript_source "havana"; exon_number "1"; gene_version "17"; exon_id "ENSE00002048533"; gene_name "PRDM16"; transcript_id "ENST00000511072"; gene_source "ensembl_havana"; transcript_support_level "5"; tag "gencode_primary"; transcript_biotype "protein_coding"; exon_version "1"; gene_id "ENSG00000142611"; transcript_version "5"; transcript_name "PRDM16-206"
1	havana		CDS		3069260	3069296	.	+	.	protein_version "1"; protein_id "ENSP00000426975"; transcript_name "PRDM16-206"; transcript_biotype "protein_coding"; gene_id "ENSG00000142611"; transcript_version "5"; tag "gencode_primary"; gene_source "ensembl_havana"; transcript_support_level "5"; gene_name "PRDM16"; transcript_id "ENST00000511072"; gene_version "17"; transcript_source "havana"; exon_number "1"; gene_biotype "protein_coding"
1	havana		start_codon	3069260	3069262	.	+	.	gene_version "17"; transcript_id "ENST00000511072"; gene_name "PRDM16"; gene_biotype "protein_coding"; transcript_source "havana"; exon_number "1"; gene_id "ENSG00000142611"; transcript_biotype "protein_coding"; transcript_version "5"; transcript_name "PRDM16-206"; transcript_support_level "5"; gene_source "ensembl_havana"; tag "gencode_primary"
```
  
El **formato GTF** representa las características y anotaciones específicas para un genoma de referencia.  Cada línea representa una región particular del genoma, y el número de columnas puede variar. Se pueden tener hasta 10 columnas, las cuales representan los siguientes indicadores:   
```
[cromosoma] [fuente] [característica] [posición.inicio] [posición.final] [puntuación] [hebra] [marco] [atributos][comentarios]
```
   
| Columna | Significado |
|-------|-------|  
| 1. Cromosoma | Indica el cromosoma donde se localiza la anotación. |
| 2. Fuente | Fuente de la anotación. Hace referencia al programa de predicción empleado o a la base de datos de anotaciones. |
| 3. Característica | Indica el tipo de característica: gen, transcrito, exón, CDS, etc. |
| 4. Posición de inicio | Posición de inicio de la característica en el genoma de referencia. |
| 5. Posición final |  Posición final de la característica en el genoma de referencia. |
| 6. Puntuación | Puntuación asociada a la característica. Si no hay puntuación, se presenta solamente un '.' |
| 7. Hebra | Permite indicar la hebra de procedencia de la característica.  Puede ser '+' , '-' o '.' (si la hebra es desconocida o no aplicable). |
| 8. Marco | Puede presentarse como '0', '1' o '2', e indica la fase o el marco donde empieza un codón. |
| 9. Atributos | Incluyen identificadores y otra información suplementaria. |
| 10. Comentarios | Es una columna con comentarios adicionales. |

  
#### 2.4.2 Recuento de características con htseq-count

Una vez visualizada la estructura del archivo de anotaciones (_Homo_sapiens.gtf_), este se emplea para anotar las características de las lecturas según la posición del genoma dónde alineen. Para ello, se debe comprobar que la forma de nombrar los cromosomas, tanto en los archivos BAM como en el archivo de anotaciones GTF, es la misma, ya que si no esto puede conducir a errores.  
  
```console
# Cambio de directorio
cd ~RNAseq_analysis/Data/
# Anotación de las características de las lecturas mapeadas con HTseq
htseq-count -t exon -i gene_id --stranded=reverse -f bam -r pos \
4_Alignment/{sample}.sorted.bam 5_Annotation/Homo_sapiens.gtf > ../Results/{saple}_counts.tsv
```
> NOTA  
> * Con la opción `-t exon` indicamos que anote aquellas lecturas alineadas específicamente contra la característica exón y, con la opción `-i gene_id`, indicamos que agrupe las diferentes lecturas atendiendo al identificador del gen al que pertenecen.  
> * Con la opción `--stranded=reverse`, determinamos que las lecturas provienen de una librería de hebra específica antisentido y, por tanto, las lecturas _forward_ mapean en la dirección contraria a la característica, mientras que las lecturas _reverse_ mapean en la misma dirección que la hebra con la caracterísitca (1+-,1-+,2++,2--).   
> * Finalmente, las opciones `–f` y `–r` , las usamos para indicar el formato del archivo de entrada (`bam`) y el orden de los alineamientos, en este caso por posición o coordenadas genómicas (`pos`).  
> * Una vez especificadas todas las opciones, indicamos la rutas relativas para llegar a los archivos BAM y al archivo de anotaciones (GTF), así como la ruta de salida y el nombre del nuevo archivo que contendrá la información.  

Para anotar las características de cada uno de los archivos BAM, se puede emplear el siguiente script `htseqcount_features.sh` desde el directorio  ~/RNAseq_analysis/Code.   
```console
#! /usr/bin/bash
# Versión de Bash: 4.2.46(2)
# Fecha: 2025
# Nombre del proyecto: RNAseq_analysis

# Script de bash para el programa HTseq; anotación de las características genómicas para las lecturas alineadas

SAMPLES="SRR28380565 SRR28380566 SRR28380568 SRR28380570 SRR28380572 SRR28380573
	SRR28380580 SRR28380582	SRR28380584 SRR28380586	SRR28380588 SRR28380589"

cd  ~RNAseq_analysis/Data/

for SAMPLE in $SAMPLES; do
	htseq-count -t exon -i gene_id --stranded=reverse -f bam -r pos \
	4_Alignment/${SAMPLE}.sorted.bam \
	5_Annotation/Homo_sapines.gtf > ../Results/${SAMPLE}_counts.tsv
	
done
```

En el caso de lecturas superpuestas en varias características genómicas, existen distintos modos que se pueden especificar en el comando `htseq-count`, según como queramos que resuelva esta información. Por defecto, se emplea el modo `union`, aunque existen otros modos diferentes según se muestra en la imágen siguiente. Para cambiar el modo se emplea la opción `-m <mode>` y se pueden especificar cualquiera de los 3 modos: `union`, `intersection-strict` o `intersection-nonempty`, aunque se recomienda trabajar con el modo por defecto (`union`).  

![image](https://github.com/user-attachments/assets/b75eb615-535e-47ae-9fe1-c2cabb82a8e8)

Finalmente, podemos observar los resultados obtenidos tras la anotación para cada una de las muestras. En el caso del archivo `SRR28380565_counts.tsv`, se obtiene la siguiente  estructura:    
```console
ENSG00000000003	0
ENSG00000000005	0
ENSG00000000419	1207
ENSG00000000457	702
[...]
ENSG00000310556	0
ENSG00000310557	0
__no_feature	3102207
__ambiguous	1097442
__too_low_aQual	1595129
__not_aligned	165990
__alignment_not_unique	0
```
> NOTA  
> El archivo resultante muestra una tabla con los conteos para cada característica, seguido por 5 grupos que contienen las lecturas que no han sido anotadas a ninguna característica.   
> `__no_feature`, hace referencia a los pares de lecturas que no se han asignado a ninguna característica.   
> `__ambiguous`, lecturas que pueden ser asignadas a más de una característica y que por tanto en el modo `union` no son asignadas.  
> `__too_low_aQual, lecturas con un mapeo de baja calidad y que no se tienen en cuenta en el análisis.   
> `__not_aligned`, lecturas no alineadas.   
> `__alignment_not_unique`, lecturas con mapeos múltiples. En este caso, el valor es 0 puesto que restringimos el alineamiento a posiciones únicas con HISAT2.

Las lecturas anotadas como ambiguas, son aquellas que mapean en regiones del genoma con genes solapantes, ya sea en la misma hebra o en hebras opuestas. Sin embargo, en nuestro caso, al tratarse de una librería con información específica de hebra, somos capacez de distinguir los genes solapantes de hebras contrarias, y por tanto, las lecturas anotadas como ambiguas pertenecen a genes solapamentes dentro de la misma hebra que tienen un sitio de comienzo y final de transcripción diferente.   
  
#### 2.4.3 Obtención de la matriz de recuentos  
  
Finalmente, juntamos todas las anotaciones de cada muestra, en un archivo común para obtener la matriz de recuentos final.  
```console
#!/usr/bin/bash

SAMPLES="SRR28380566 SRR28380565 SRR28380570 SRR28380572 SRR28380573"
cd ~/RNAseq_analysis/Results


for SAMPLE in $SAMPLES; do
	sed -i "1s/^/feature\t${SAMPLE}\n/" "${SAMPLE}_counts.tsv" 
done

```



## 3 Analisis estadístico de los datos de RNAseq y Genes Diferencialmente Expresados
### 3.1 Instalación y carga en memoria de las librerías empleadas. 
  
En primer lugar, se cargan todas las librerías necesarias para el análisis y se establece el directorio de trabajo.   
```R

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("edgeR") # Análisis de genes diferencialmente expresados
BiocManager::install("biomaRt") # Conversión de identificadores de los genes
BiocManager::install("limma") #Análisis de datos, modelos lineales y expresión diferencial
BiocManager::install("clusterProfiler") #Análisis de enriquecimiento génico
BiocManager::install("org.Hs.eg.db") # Base de datos para la anotación humana

# Librerías complementarias
## Librerías del paquete Tidyverse: ggplot2, dplyr, tidyr, ggrepel, tibble

# Creación de una variable con todas las librerías necesarias
list_packages <- c("edgeR",  "ggplot2", "dplyr", "tidyr", "ggrepel", "tibble",
                   "biomaRt", "limma", "clusterProfiler", "org.Hs.eg.db", "enrichplot")

# Comprobación de las librerías instaladas e instalación de las librerías faltantes
new_packages <- list_packages[!(list_packages %in% installed.packages())]
if(length(new_packages)>0) 
  install.packages(new_packages)

# Carga en memoria de todas las librerías especificadas
invisible(lapply(list_packages, FUN=library, character.only =TRUE))

# Establecimiento del directorio de trabajo
setwd("C:/Users/sonni/OneDrive/Escritorio/2025/RNAseq_analysis/")
```
  
**3.2 Importación de la matriz de recuentos y metadatos**  
  
Para realizar el análisis, primero se importan los datos de interés y se ajustan los nombres de las muestras.   
```R
metadatos <- read.csv(file ="Results/Matadata.csv") #carga de los metadatos

seqdata = read.delim(file = "Results/matriz_conteos.tsv", sep = " ") #importación de la matriz de conteos
colnames(seqdata) <- c("feature", metadatos$Sample.Name) #Cambio de los nombres de las columnas
```
Si observan las dimensiones de la matriz de conteos, se tienen 78937 filas y 13 columnas.  
```R
> dim(seqdata)
[1] 78937    13
> head(seqdata)
          feature GSM8153253 GSM8153252 GSM8153250 GSM8153248 GSM8153246 GSM8153245 GSM8153238 GSM8153236 GSM8153234 GSM8153232 GSM8153230 GSM8153229
1 ENSG00000000003          0          0          5          0          0          1          0          2          2          1          0          1
2 ENSG00000000005          0          0          0          0          0          0          0          0          0          0          0          0
3 ENSG00000000419       1207       1885       1501       1535        832       1002       1410       1370       1122       1234        954        488
4 ENSG00000000457        702        726        823       1175        620        690        656        729       1018       1345        495        773
5 ENSG00000000460         70         64        103        128         91         80         55         70        121        138         55         46
6 ENSG00000000938      20342      45609      81968      44188      25032      28710      17354      36541      19381      18132      18277      22895
> tail(seqdata)
                     feature GSM8153253 GSM8153252 GSM8153250 GSM8153248 GSM8153246 GSM8153245 GSM8153238 GSM8153236 GSM8153234 GSM8153232 GSM8153230 GSM8153229
78932        ENSG00000310557          0          0          0          0          0          0          0          0          0          0          0          0
78933           __no_feature    3102207    4248945    5368135    5778559    6082434    5070486    4308555    7094467    5194665    4659772    3686613    3506258
78934            __ambiguous    1097442    2447154    2567359    2496826    1193122    1173681    1100280    2363146    1611862    1701720    1244623    1296521
78935        __too_low_aQual    1595129    3221199    4157612    3383771    1924811    1925120    1683095    3203467    2386385    2391636    1854569    1919572
78936          __not_aligned     165990     239233     258154     221530     178293     122302     163687     245773     189997     251794     176381     173047
78937 __alignment_not_unique          0          0          0          0          0          0          0          0          0          0          0          0
```
  
**3.3 Representación de las características genómicas anotadas**  
  
A partir de la matriz de conteos final, se calcula el porcentaje total de lecturas asignadas, no_feature, ambiguous, too_low_aQual y not_aligned; y el resultado se muestra en un gráfico de barras.  
```R
annot <- seqdata  
annot[1:78932,1] <- "gene"  #conversión de los identificadores de los genes a una categoría común

annot <- annot %>% group_by(feature) %>% summarise(across(where(is.numeric), sum)) #recuento de lecturas para cada categoria
```
Tras la ejecución se observa el número de recuentos para cada categoría:  
```R
> head(annot)
# A tibble: 6 × 13
  feature                GSM8153253 GSM8153252 GSM8153250 GSM8153248 GSM8153246 GSM8153245 GSM8153238 GSM8153236 GSM8153234 GSM8153232 GSM8153230 GSM8153229
  <chr>                       <int>      <int>      <int>      <int>      <int>      <int>      <int>      <int>      <int>      <int>      <int>      <int>
1 __alignment_not_unique          0          0          0          0          0          0          0          0          0          0          0          0
2 __ambiguous               1097442    2447154    2567359    2496826    1193122    1173681    1100280    2363146    1611862    1701720    1244623    1296521
3 __no_feature              3102207    4248945    5368135    5778559    6082434    5070486    4308555    7094467    5194665    4659772    3686613    3506258
4 __not_aligned              165990     239233     258154     221530     178293     122302     163687     245773     189997     251794     176381     173047
5 __too_low_aQual           1595129    3221199    4157612    3383771    1924811    1925120    1683095    3203467    2386385    2391636    1854569    1919572
6 gene                     22466408   40491330   42397848   40694703   24028088   24842924   23833401   34717738   30330906   30500265   22604489   24532419
```
Posteriormente, se eliminan aquellas categorías para las cuales no se tiene información y se cambia la estructura de la tabla a un formato largo con una función del paquete tidyr.
```R
annot <- annot[-1,] #eliminación de las categorías no informativas
annot_long <- pivot_longer(annot, 
                           cols = 2:ncol(annot), 
                           names_to = "IDLecturas", 
                           values_to = "Contaje")
```
Como resultado se obtiene la siguiente estructura:   
```R
> head(annot_long)
# A tibble: 6 × 3
  feature     IDLecturas Contaje
  <chr>       <chr>        <int>
1 __ambiguous GSM8153253 1097442
2 __ambiguous GSM8153252 2447154
3 __ambiguous GSM8153250 2567359
4 __ambiguous GSM8153248 2496826
5 __ambiguous GSM8153246 1193122
6 __ambiguous GSM8153245 1173681
```
A continuación, se calcula el total de lecturas anotadas para cada muestra:  
```R
total <- annot_long %>% group_by(IDLecturas) %>% summarise(across(where(is.numeric), sum))
```
```R
> head(total, n=12)
# A tibble: 12 × 2
   IDLecturas  Contaje
   <chr>         <int>
 1 GSM8153229 31427817
 2 GSM8153230 29566675
 3 GSM8153232 39505187
 4 GSM8153234 39713815
 5 GSM8153236 47624591
 6 GSM8153238 31089018
 7 GSM8153245 33134513
 8 GSM8153246 33406748
 9 GSM8153248 52575389
10 GSM8153250 54749108
11 GSM8153252 50647861
12 GSM8153253 28427176
```
Finalmente, se lleva a cabo la creación de una nueva columna con el porcentaje de lecturas correspondiente a cada categoría en cada una de las muestras:    
```R
annot_long <- annot_long[order(annot_long$IDLecturas),] #ordenación de las muestras

annot_long$Total <- rep(total$Contaje, each = 5)  #creación de una nueva columna con el total de lecturas por muestra
annot_long <- annot_long %>% mutate(Porcentaje = Contaje/Total*100) #creación de una columna con el porcentaje para cada cetegoria
```
Tras la ejecución de los comandos, se obtiene el siguiente resultado:  
```R
> head(annot_long)
# A tibble: 6 × 5
  feature         IDLecturas  Contaje    Total Porcentaje
  <chr>           <chr>         <int>    <int>      <dbl>
1 __ambiguous     GSM8153229  1296521 31427817      4.13 
2 __no_feature    GSM8153229  3506258 31427817     11.2  
3 __not_aligned   GSM8153229   173047 31427817      0.551
4 __too_low_aQual GSM8153229  1919572 31427817      6.11 
5 gene            GSM8153229 24532419 31427817     78.1  
6 __ambiguous     GSM8153230  1244623 29566675      4.21 
```
Una vez se tienen los porcentajes de los recuentos para cada categoría, se crea un gráfico de frecuencias con la librería ggplot2.  
```R
# Gráfico con la anotación de lecturas
ggplot(annot_long, aes(x=IDLecturas, y=Porcentaje, fill=feature)) + 
  geom_bar(stat="identity") + #selección de la geometría
  labs(fill="", x="") + #cambio de nombre del eje X y leyenda
  scale_fill_manual(values=c("lightblue","blue3","lightgreen","green4","lightsalmon"),
                    labels=c('Lecturas ambiguas', 'Lecturas sin anotar', 'Lecuras sin alinear', 'Lecturas con mapeo de baja calidad', 'Lecturas anotadas como genes')) + #personalización de los colores y etiquetas
  guides(x=guide_axis(angle=90 )) + #cambio de orientación de los ejes 
  scale_y_continuous(expand = expansion(mult = 0)) + #ajuste del grafico
  theme_classic(base_size=20)  #selección del tema
```
![image](https://github.com/user-attachments/assets/2442616f-f5d7-423c-b0dd-028887e17054)
  
**3.4 Análisis de los genes diferencialmente expresados con el paquete edgeR**. 
  
Para el análisis de Genes Diferencialmente Expresados se utiliza el paquete edgeR de Bioconductor y las funciones recogidas en la tabla siguiente:    

| Etapa | Función de edgeR |
|-------|-----------------|
| Importacion de la matriz de conteos | DGEList |
| Filtración de las lecturas con recuentos bajos o nulos | filterbyExpr | 
| Normalización de las librerías | calcNormFactors |
| Estimacion de la dispersion y visualización | estimateDisp / plotBCV |
| juste de los datos y test de expresión diferencial | enfoque _quasi-likelihood_ (QL): glmQLFit / glmQLFTest |
| Resultados | decideTest / topTags |  
  
**Conversión de la matriz de recuentos al objeto DGEList**  
  
A partir de la matriz con los recuentos crudos, se genera un objeto DGEList. Este objeto está formado por 2 componentes esenciales: la matriz de conteos y la información de las muestras.  
  
```R
seqdata <- seqdata[1:78932,] # Selección de las lecturas anotadas como genes
seqdata <- seqdata %>% tibble::column_to_rownames(var = "feature") #conversión de la columna feature a los nombres de las filas

y <- DGEList(seqdata) #creación del objeto DGEList
```
  
Una vez creado el objeto DGEList, se especifica el grupo experimental de cada muestra: sano o enfermo.   
```R
disease <- rep(c("Sano","LES"),each=6)  #creación de una variable con las condiciones experimentales
disease <- factor(disease, levels =c('Sano','LES')) #conversión a una variable categórica

y$samples$group <- disease #adición de los grupos experimentales al objeto DGEList
```
  
Una vez modificado el objeto DGEList, se verifican que los cambios se hayan efectuado.  
```R
> head(y)
An object of class "DGEList"
$counts
                GSM8153253 GSM8153252 GSM8153250 GSM8153248 GSM8153246 GSM8153245 GSM8153238 GSM8153236 GSM8153234 GSM8153232 GSM8153230 GSM8153229
ENSG00000000003          0          0          5          0          0          1          0          2          2          1          0    	  1
ENSG00000000005          0          0          0          0          0          0          0          0          0          0          0   	  0
ENSG00000000419       1207       1885       1501       1535        832       1002       1410       1370       1122       1234        954	488
ENSG00000000457        702        726        823       1175        620        690        656        729       1018       1345        495	773
ENSG00000000460         70         64        103        128         91         80         55         70        121        138         55	 46
ENSG00000000938      20342      45609      81968      44188      25032      28710      17354      36541      19381      18132      18277      22895

$samples
           group lib.size norm.factors
GSM8153253  Sano 22466408            1
GSM8153252  Sano 40491330            1
GSM8153250  Sano 42397848            1
GSM8153248  Sano 40694703            1
GSM8153246  Sano 24028088            1
7 more rows ...
```
Además, se puede comprobar el número de genes totales incluidos en la matriz de conteos.  
```R
> nrow(y$counts)
[1] 78932
```
  
**3.1.3 Eliminación de genes con recuentos bajos**  
  
Sin embargo, pese a tener tantos genes no todos se expresan por lo que se filtran los genes con recuentos nulos o inferiores a 10. 
```R
keep.genes <- filterbyExpr(y)  #filtrado de genes con expresión baja o nula
y <- y[keep.genes, keep.lib.sizes=F] #seleccion de genes y cálculo de nuevo del tamaño de librería
```
Para comprobar el número de genes mantenidos tras el filtrado, se emplea el siguiente comando:  
```R
> nrow(y$counts)
[1] 14617
> nrow(y$counts) / nrow(seqdata)*100 #Porcentaje de genes mantenidos respecto al número de genes iniciales
[1] 18.51847
```
De todos los genes iniciales, solo el 18.5% se expresa en las muestras estudiadas.  
  
**3.1.4 Cálculo de los factores de normalización**  
  
Los recuentos filtrados obtenidos previamente para cada gen se normalizan mediante el cálculo de factores de normalización para corregir las diferencias de composición en las librerías. Para ello, se usa la función calcNormFactors() del paquete edgeR, que permite la normalización mediante el método de la Media Recortada de los valores M (TMM, del inglés Trimmed Mean of M-values).     
```R
y <- calcNormFactors(y)
```
Los factores de normalización calculados se añaden al objeto ‘DGEList’  dentro del apartado con la información para las muestras.  
```R
> y$samples
 group lib.size norm.factors
GSM8153253  Sano 22443109    0.9462243
GSM8153252  Sano 40411595    1.0450250
GSM8153250  Sano 42316106    1.0706549
GSM8153248  Sano 40636040    1.0167785
GSM8153246  Sano 23994900    0.9328915
GSM8153245  Sano 24810656    0.9682306
GSM8153238   LES 23787855    0.9871341
GSM8153236   LES 34659877    1.0637101
GSM8153234   LES 30280095    0.9769371
GSM8153232   LES 30454394    0.9879841
GSM8153230   LES 22574501    0.9503878
GSM8153229   LES 24480946    1.0677673
```
El producto entre el tamaño de librería real y el factor de normalización computado para cada muestra, resulta en el tamaño de librería efectivo, reduciendo así las diferencias de composición.  
  
**3.1.5 Estimacion de la variabilidad biologica**   

Para estimar la variabilidad y similitud entre las réplicas biológicas se emplea un gráfico de escala multidimensional MDS con la función plotMDS del paquete _limma_.     
```R
colors <- c("steelblue2","indianred2") 

limma::plotMDS(y,
               col = colors[y$samples$group],
               pch = 16,  #forma
               cex = 1.4) #tamaño

legend("bottomright", as.character(unique(y$samples$group)), 
       pch = 16,
       col = colors,
       ncol = 2 , cex = 0.9)
```
Este gráfico permite hacerse una idea de las relaciones entre las muestras, de forma que las muestras con perfiles de expresión de genes similares estarán más cerca en el gráfico.  
  
**3.1.6 Dispersión de los genes**  

Para estimar la sobredispersión de los genes, primero se crea una matriz con el diseño experimental que especifique cómo se asocian o agrupan las muestras.    
```R
> mdesign <- model.matrix(~0 + disease)           # Creación de la matriz con el diseño experimental
> colnames(mdesign) <- gsub("disease","", colnames(mdesign)) # Modificación de la cabecera 
> rownames(mdesign) <- colnames(seqdata)          # Modificación de los nombres de las filas
> mdesign
           Sano LES
GSM8153253    1   0
GSM8153252    1   0
GSM8153250    1   0
GSM8153248    1   0
GSM8153246    1   0
GSM8153245    1   0
GSM8153238    0   1
GSM8153236    0   1
GSM8153234    0   1
GSM8153232    0   1
GSM8153230    0   1
GSM8153229    0   1
attr(,"assign")
[1] 1 1
attr(,"contrasts")
attr(,"contrasts")$disease
[1] "contr.treatment"
```

Seguidamente se estima la dispersión de los genes con la función estimateDisp, y se especifica tanto la matriz de diseño como el parámetro robust=T para proteger la estimación contra los outliers.  
```R
y <- estimateDisp(y, mdesign, robust = TRUE) #Dispersión por defecto con reducción bayesiana
```

El resultado de aplicar la función anterior genera una serie de resultados estadísticos dentro del objeto ‘DGEList’. La dispersión de los genes se va a evaluar desde 3 puntos de vista diferentes: “common dispersion”, “trended dispersion” y “tagwise dispersion”. 
 
```R
> y
An object of class "DGEList"
$counts
                GSM8153253 GSM8153252 GSM8153250 GSM8153248 GSM8153246 GSM8153245 GSM8153238 GSM8153236 GSM8153234 GSM8153232 GSM8153230 GSM8153229
ENSG00000000419       1207       1885       1501       1535        832       1002       1410       1370       1122       1234        954        488
ENSG00000000457        702        726        823       1175        620        690        656        729       1018       1345        495        773
14615 more rows ...
$samples
           group lib.size norm.factors
GSM8153253  Sano 22443109    0.9462243
GSM8153252  Sano 40411595    1.0450250
10 more rows ...
$design
           Sano LES
GSM8153253    1   0
GSM8153252    1   0
10 more rows ...

$common.dispersion
[1] 0.1727263
$trended.dispersion
[1] 0.1261033 0.1350513 0.2065515 0.1105549 0.1858998
14612 more elements ...
$tagwise.dispersion
[1] 0.09745349 0.09531198 0.15104251 0.06725113 0.15657217
14612 more elements ...
$AveLogCPM
[1] 5.353583 4.796460 1.566854 9.959342 2.391398
14612 more elements ...
$trend.method
[1] "locfit"
$prior.df
[1] 5.286141 5.286141 5.286141 5.286141 5.286141
14612 more elements ...
$prior.n
[1] 0.5286141 0.5286141 0.5286141 0.5286141 0.5286141
14612 more elements ...
$span
[1] 0.2938649
```
En la estimación anterior, la dispersión individual de los genes o _tagwise_, se corrige mediante un modelo de reducción bayesiana para aproximar las dispersiones individuales a la línea de tendencia, con el fin de modelar los genes de manera más confiable y precisa.    
  
Para observar las diferencias entre los valores de dispersión individuales corregidos y sin corregir, se puede emplear el código siguiente:
```R
y.priordf0 <- estimateDisp(y, mdesign, prior.df = 0, robust = TRUE) #Dispersión sin reducción bayesiana

#Representación de los resutlados 
par(mfrow = c(1,2))
plotBCV(y.priordf0)
title("Dispersión tagwise sin ajuste Bayesiano")
plotBCV(y)
title("Dispersión tagwise con ajuste Bayesiano")
```
  
**3.1.7 Ajuste de los datos de conteo con el enfoque _Quasi-likelihood_**  

Una vez computado los valores de dispersión, se ajustan los datos mediante el enfoque de cuasidispersión o quasi-likelihood (QL), que permite obtener resultados más fiables.    
```R
fit <- glmQLFit(y, mdesign, robust = TRUE) # Ajuste del modelo con la función glmQLFit
```
```R
> class(fit)  #creación de un objeto DGEGLM
[1] "DGEGLM"
attr(,"package")
[1] "edgeR"
```
El nuevo objeto creado `fit` es un objeto de tipo DGEGLM y cuenta con diferentes parámetros computados para cada uno de los genes, tales como coeficientes, valores ajustados...  
  
**3.1.8 Prueba de significancia o Test de expresión diferencial**  
  
Una vez computada la dispersión y ajustados los datos, se lleva a cabo la prueba de significancia o test de expresión diferencial.  
Primero, con la función makeContrast del paquete limma, se define el tipo de comparación a realizar entre los grupos experimentales.   
```R
mcontrast.LvsS <- makeContrasts(lupusVSsano = LES-Sano, levels = mdesign) # Creación de la matriz con los grupos contraste a comparar
```
La función makeContrast genera una matriz numérica que representa los grupos indicados a contrastar. El valor 1 y –1 corresponde a los grupos a comparar.  
```R
> mcontrast.LvsS
      Contrasts
Levels lupusVSsano
  Sano          -1
  LES            1
```

A continuación, se realiza el test de expresión diferencial con la función glmQLFTest() del paquete de edgeR y se genera un objeto de tipo DGELRT.
```R
test.LvsS <- glmQLFTest(fit, contrast = mcontrast.LvsS)  #test de expresión diferencial
```
```R
> class(test.LvsS)  #creación de un objeto DGELRT
[1] "DGELRT"
attr(,"package")
[1] "edgeR"
```
El objeto DGELRT contiene algunos parámetros estadísticos comunes al objeto DGEGLM con el modelo ajustado. Sin embargo, aparecen unos subapartados nuevos con los resultados de la comparación. Los resultados obtenidos tras la comparación se guardan en el subapartado “table” , que contiene diferentes estadísticos como logFC, logCPM, F y P-valores.  
  
```R
> head(test.LvsS$table)
                       logFC   logCPM            F      PValue
ENSG00000000419 -0.066630968 5.353583 7.994665e-02 0.781124674
ENSG00000000457  0.229701188 4.796460 9.960150e-01 0.333659564
ENSG00000000460 -0.007569388 1.566854 6.056389e-04 0.980681644
ENSG00000000938 -0.616529724 9.959342 1.068466e+01 0.005016953
ENSG00000001036  0.297974957 2.391398 9.517483e-01 0.344313726
ENSG00000001084  0.430007747 1.975376 4.473338e+00 0.051073759
```

El valor de P es un parámetro estadístico que se usa para estimar la probabilidad de que un resultado sea obtenido al azar, de forma que cuando la probabilidad de obtener ese mismo resultado al azar es menor al 5% (p < 0.05), decimos que se trata de un resultado estadísticamente significativo. Sin embargo, al comparar un número de genes tan alto, aumenta la probabilidad de tener valores de p significativos (p< 0.05), simplemente por azar, aunque no lo sean realmente. Por lo que los valores de P se deben corregir.  
  
**3.1.9 Corrección de P-valores y filtrado de genes diferencialmente expresados**  
  
Para resolver el problema de testeado múltiple, los valores de P calculados se ajustan mediante la técnica de corrección _Benjamini–Hochberg_, que calcula los valores de FDR (_False Discovery Rate_) con el objetivo de reducir la tasa de falsos positivos. La función decideTest, permite llevar a cabo esta corrección e identificar los genes diferencialmente expresados en función al umbral de FDR y logFC establecido.  

```R
dTest_LvsS <- decideTests(test.LvsS, adjust.method = "BH", 
                             p.value = 0.01,   # Umbral de p-valor corregido 
			     lfc = 2)          # Umbral de log2-Fold-Change
```
La función `summary` permite observar los resultados obtenidos:  
```R
> summary(dTest_LvsS)
       -1*Sano 1*LES
Down               1
NotSig         14558
Up                58
```
  
Seguidamente, los genes se ordenan en función al valor de FDR con la función topTags(), que crea un objeto específico de edgeR con el mismo nombre ‘topTags’.  

```R 
toptag.LvsS <- topTags(test.LvsS, n = Inf)  #Creación de un objeto de clase topTags
```
```R
> names(toptag.LvsS)
[1] "table"         "adjust.method" "comparison"    "test"
```
El objeto TopTags contiene diferentes elementos, sin embargo, solamente nos interesa el elemento "table" con los resultados para los genes, por lo que se modifica la variable `toptag.LvsS`.  
```R
toptag.LvsS <- toptag.LvsS$table       
```
La variable `toptag.LvsS` contiene un data frame con todos los genes y los resultados para la comparación entre los grupos experimentales.  
```R
> head(toptag.LvsS)
                    logFC   logCPM        F       PValue         FDR
ENSG00000115155 10.983012 5.905536 62.94245 1.088771e-07 0.001591456
ENSG00000299483  4.081678 2.248540 75.21752 2.668500e-07 0.001922562
ENSG00000137965  4.880786 8.563740 69.78641 5.928582e-07 0.001922562
ENSG00000305126  5.045554 1.783654 62.22386 7.416759e-07 0.001922562
ENSG00000304354  2.430758 3.421995 61.26940 9.308613e-07 0.001922562
ENSG00000126709  3.069027 7.076709 60.70911 1.028534e-06 0.001922562
```

Con el fin de conocer los genes diferencialmente expresados, se lleva cabo la selección de todas aquellas filas que cumplan los filtros especificados: valores de FDR menor a 0.01 y valor absoluto de logFC mayor o igual al valor absoluto de 2. 

```R
# Anotación de los genes sobreexpresados e infraexpresados en una nueva columna
toptag.LvsS <- toptag.LvsS %>%  mutate(DE = case_when(logFC >= 2 & FDR < 0.01 ~ "Up",
                                                      logFC <= -2 & FDR < 0.01 ~ "Down",
                                                      .default="NotSign"))
```
En el caso del logFC, es una escala simétrica por lo que nos interesan cambios positivos y negativos. Valores de log de FC  mayores a 2 los anotamos como genes sobreexpresados, mientras que valores de log FC por debajo de –2 los anotamos como genes infraexpresados.  
  
```R
> nrow(toptag.LvsS[toptag.LvsS$DE == "Up",]) # Número de genes anotados como UP
[1] 58
> nrow(toptag.LvsS[toptag.LvsS$DE == "Down",]) # Número de genes anotados como DOWN
[1] 1
```

## 4 Conversión de los identificadores de los genes diferencialmente expresados

Para la anotación de los genes, se emplea el paquete biomaRt, que emplea información de bases de datos en línea y convierte los identificadores (IDs) de los genes entre distintos sistemas (Ensembl, Entrez, Symbol, etc.). Además, también permite obtener las descripciones de los genes.

Para ello, primero se crea una tabla que contenga únicamente los genes diferencialmente expresados.  
```R
de_table <- toptag.LvsS[toptag.LvsS$DE %in% c("Up","Down"),] #selección de genes sobre e infraexpresados
de_table <- de_table %>% arrange(row.names(de_table))  #ordenación de los identificadores Ensembl ID de los genes
```
```R
> dim(de_table)
[1] 59  6
```
  
A continuación, se lleva a cabo la conversión de los identificadores en formato _Ensembl ID_ a otros sistemas de anotación.  
```R
# Conexión con la base de datos Ensembl
ensembl_113 <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",
			  version = 113) #misma versión que el archivo de anotaciones descargado desde la base de datos Ensembl, empleado para la anotación de lecturas

#Anotación de genes con diferentes sistemas
geneID <- getBM(attributes = c('ensembl_gene_id', 'entrezgene_id', 'external_gene_name', 'description'),   
                  filters = 'ensembl_gene_id', 
                  values = row.names(de_table), 
                  mart = ensembl_113)
```
Como resultado, se obtiene un objeto con las anotaciones correspondientes encontradas para cada gen.  
```R
> head(geneID, n = 10)
   ensembl_gene_id entrezgene_id external_gene_name                                                                                    description
1  ENSG00000059378         64761             PARP12               poly(ADP-ribose) polymerase family member 12 [Source:HGNC Symbol;Acc:HGNC:21919]
2  ENSG00000078081         27074              LAMP3                    lysosomal associated membrane protein 3 [Source:HGNC Symbol;Acc:HGNC:14582]
3  ENSG00000088827          6614            SIGLEC1                       sialic acid binding Ig like lectin 1 [Source:HGNC Symbol;Acc:HGNC:11127]
4  ENSG00000089127          4938               OAS1                           2'-5'-oligoadenylate synthetase 1 [Source:HGNC Symbol;Acc:HGNC:8086]
5  ENSG00000106785          9830             TRIM14                             tripartite motif containing 14 [Source:HGNC Symbol;Acc:HGNC:16283]
6  ENSG00000108771         79132              DHX58                                       DExH-box helicase 58 [Source:HGNC Symbol;Acc:HGNC:29517]
7  ENSG00000111331          4940               OAS3                           2'-5'-oligoadenylate synthetase 3 [Source:HGNC Symbol;Acc:HGNC:8088]
8  ENSG00000111335          4939               OAS2                           2'-5'-oligoadenylate synthetase 2 [Source:HGNC Symbol;Acc:HGNC:8087]
9  ENSG00000115155          9381               OTOF                                                   otoferlin [Source:HGNC Symbol;Acc:HGNC:8515]
10 ENSG00000119917          3437              IFIT3 interferon induced protein with tetratricopeptide repeats 3 [Source:HGNC Symbol;Acc:HGNC:5411]
```
Sin embargo, existen algunos identificadores que el programa no es capaz de encontrar, de forma que se lleva a cabo una búsqueda manual y se completan los identificadores faltantes.   
```R
geneID[geneID$ensembl_gene_id == "ENSG00000225964",2] <- "104326052"
geneID[geneID$ensembl_gene_id == "ENSG00000228318",2] <- "130890644"
geneID[geneID$ensembl_gene_id == "ENSG00000233975",2] <- "111216282"
geneID[geneID$ensembl_gene_id == "ENSG00000279296",2] <- "109245082"
geneID[geneID$ensembl_gene_id == "ENSG00000289234",2] <- "107985115"
geneID[geneID$ensembl_gene_id == "ENSG00000275676",3] <- "lnc-FOXD4L5-97"
```
Finalmente, se combinan los identificadores encontrados con los resultados estadísticos obtenidos para los genes diferencialmente expresados.  
```R
de_results <- cbind(geneID, de_table$logFC, de_table$PValue, de_table$FDR, de_table$DE)
colnames(de_results) <- c("Ensemble ID","Entrez ID", "Nombre del Gen", "Descripción", "logFC", "P-valor", "FDR", "DE" )  
```
```R
> head(de_results)
      Ensemble ID Entrez ID Nombre del Gen                                                                      Descripción    logFC      P-valor         FDR DE
1 ENSG00000059378     64761         PARP12 poly(ADP-ribose) polymerase family member 12 [Source:HGNC Symbol;Acc:HGNC:21919] 2.181316 1.301791e-06 0.001922562 Up
2 ENSG00000078081     27074          LAMP3      lysosomal associated membrane protein 3 [Source:HGNC Symbol;Acc:HGNC:14582] 3.104834 3.494549e-05 0.007663624 Up
3 ENSG00000088827      6614        SIGLEC1         sialic acid binding Ig like lectin 1 [Source:HGNC Symbol;Acc:HGNC:11127] 6.430987 3.054285e-06 0.002125928 Up
4 ENSG00000089127      4938           OAS1             2'-5'-oligoadenylate synthetase 1 [Source:HGNC Symbol;Acc:HGNC:8086] 5.504131 8.135505e-06 0.003129386 Up
5 ENSG00000106785      9830         TRIM14               tripartite motif containing 14 [Source:HGNC Symbol;Acc:HGNC:16283] 2.219036 1.994167e-05 0.005397914 Up
6 ENSG00000108771     79132          DHX58                         DExH-box helicase 58 [Source:HGNC Symbol;Acc:HGNC:29517] 2.742959 1.270625e-06 0.001922562 Up
```

## 5 Visualización de los genes diferencialmente expresados en un gráfico de Volcán   

Una vez convertidos los identificadores _Ensembl_ de los genes diferencialmente expresados a los nombres comunes, se lleva a cabo un gráfico de tipo Volcano Plot para mostrar los resultados, y se etiquetan los 5 genes con mayor cambio de expresión resultantes tras la comparación.  Para ello, se emplea el objeto `topTags.LvsS` con todos los genes testeados en el análisis diferencial.  
```R
> head(toptag.LvsS)
                    logFC   logCPM        F       PValue         FDR DE
ENSG00000115155 10.983012 5.905536 62.94245 1.088771e-07 0.001591456 Up
ENSG00000299483  4.081678 2.248540 75.21752 2.668500e-07 0.001922562 Up
ENSG00000137965  4.880786 8.563740 69.78641 5.928582e-07 0.001922562 Up
ENSG00000305126  5.045554 1.783654 62.22386 7.416759e-07 0.001922562 Up
ENSG00000304354  2.430758 3.421995 61.26940 9.308613e-07 0.001922562 Up
ENSG00000126709  3.069027 7.076709 60.70911 1.028534e-06 0.001922562 Up
> dim(toptag.LvsS)
[1] 14617     6
```
Primero, se crea una nueva columna y se anotan solamente los nombres comunes para 5 genes sobreexpresados con mayor valor de logFC. 
```R
toptag.LvsS$Annot <- NA  #creación de una nueva columna
toptag.LvsS <- toptag.LvsS %>% arrange(desc(abs(logFC)))  #rrden descendiente del valor absoluto del logFC

ensemblID <- head(toptag.LvsS,5) # selección de las primeras 5 filas

top5 <- geneID[(geneID$ensembl_gene_id %in% rownames(ensemblID)),c(-2,-4)] #búsqueda de los nombres comunes para los genes en el objeto geneID
```
```R
> head(top5)
   ensembl_gene_id external_gene_name
3  ENSG00000088827            SIGLEC1
9  ENSG00000115155               OTOF
22 ENSG00000137959             IFI44L
33 ENSG00000160932               LY6E
41 ENSG00000184979              USP18
```
Una vez identificados los nombres de los genes se agregan al _data frame_ `toptag.LvsS`  
```R
toptag.LvsS <- toptag.LvsS %>% rownames_to_column(var = "ensemblID")
toptag.LvsS[toptag.LvsS$ensemblID == "ENSG00000088827", "Annot"] <- "SIGLEC1"
toptag.LvsS[toptag.LvsS$ensemblID == "ENSG00000115155", "Annot"] <- "OTOF"
toptag.LvsS[toptag.LvsS$ensemblID == "ENSG00000137959", "Annot"] <- "IFI44L"
toptag.LvsS[toptag.LvsS$ensemblID == "ENSG00000160932", "Annot"] <- "LY6E"
toptag.LvsS[toptag.LvsS$ensemblID== "ENSG00000184979", "Annot"] <- "USP18"
```
```R
> head(toptag.LvsS)
        ensemblID     logFC   logCPM        F       PValue         FDR DE   Annot
1 ENSG00000115155 10.983012 5.905536 62.94245 1.088771e-07 0.001591456 Up    OTOF
2 ENSG00000184979  6.506025 4.393928 51.75408 4.219425e-06 0.002422494 Up   USP18
3 ENSG00000088827  6.430987 4.775736 54.65455 3.054285e-06 0.002125928 Up SIGLEC1
4 ENSG00000160932  6.195741 9.043884 48.97921 5.487287e-06 0.002506490 Up    LY6E
5 ENSG00000137959  5.853624 8.842057 44.92089 8.900121e-06 0.003215684 Up  IFI44L
6 ENSG00000111335  5.755085 9.088278 55.40027 2.713972e-06 0.001983506 Up    <NA>
```
Finalmente, se representa los resultados en el gráfico Volcano Plot.  
```R
ggplot(toptag.LvsS, aes(x=logFC, y=-log10(FDR), color=DE, label=Annot)) + 
  geom_point(size = 2) + 
  scale_color_manual(values = c("blue","grey","indianred2"))+
  geom_vline(xintercept=c(2,-2), linetype = 3) +
  geom_label_repel(color="black", size = 5, box.padding = 0.3, nudge_y = 0.6) +
  labs(color="") +
  theme_classic(base_size = 20) 
```
  
### 4.1 Análisis de enriquecimiento funcional

A pesar de determinar los genes diferencialmente expresados, en general, es difícil interpretar su significado biológico. Es por ello, que existen herramientas informáticas que permiten establecer las rutas o funciones metabólicas más representadas en los genes diferencialmente expresados. En este caso solamente se analizarán los genes sobreexpresados.  

```R
# Filtrado de genes sobreexpresados
res_up <- de_results[de_results$DE == "Up",]

#Análisis de enriquecimiento funcional de genes en términos de Gene Ontology
# Funciones moleculares sobrerrepresentadas
ego_MF <- enrichGO(gene = res_up$`Nombre del Gen`,  # Lista de genes
                   OrgDb = org.Hs.eg.db,            # Base de datos de anotación
                   keyType = 'SYMBOL',              # Tipo de identificador
                   ont="MF",                        # Ontología: BP (Biological Process)
                   pAdjustMethod = "BH",            # Método de correción de p-valores
                   pvalueCutoff = 0.01,             # Umbral de significancia
                   qvalueCutoff = 0.05)             # Umbral de FDR

# Procesos biológicos sobrerrepresentados
ego_BP <- enrichGO(gene = res_up$`Nombre del Gen`, 
                   OrgDb = org.Hs.eg.db,
                   keyType = 'SYMBOL',
                   ont="BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05)

#Simplificación de los procesos biológicos y eliminación de términos redundantes
ego_BP2 <- simplify(ego_BP, cutoff=0.7, by="p.adjust", select_fun=min)
```
Finalmente, se lleva a cabo la visualización con el paquete 'enrichplot'.    
```R
# Gráfico de barras para la visualización de Funciones Moleculares sobrerrepresentadas
barplot(ego_MF, font.size = 20)  

# Gráfico con las relaciones entre las Funciones Moleculares sobrerrepresentadas
cnetplot(ego_MF, color.params = list(category = "blue", gene = "grey"),
         cex.params = list(category_node = 3, gene_node = 1, category_label = 1.8, gene_label = 1.5))

# Gráfico con las relaciones jerárquicas entre los Procesos Biológicos sobrerrepresentados
treeBP <- pairwise_termsim(ego_BP2)   #Cálculo de la similitud entre los términos
treeplot(treeBP, label_format=10, fontsize=5) #Representación gráfica
```
