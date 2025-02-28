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
> La carpeta **Data**, contendrá las lecturas crudas y procesadas, losarchivos necesarios para la deducción de la direccionalidad de la librería y los resultados tras el alineamiento, entre otros.  
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
# Descomprisión del archivo descargado
tar -xvf grch38_genome.tar.gz
```
> NOTA   
> Con el comando `mkdir` creamos un nuevo directorio para el genoma de referencia dentro de nuestra carpeta 4_Alignment. El comando `cd $_` permite movernos a esta nueva carpeta dónde descargamos el archivo de interés con el comando `wget`.  
> Tras la descompresión del archivo, se genera el directorio `/grch38/` con el genoma de referencia, el script `make_grch38.sh` y los archivos necesarios para la indexación identificados por la palabra `genome` seguidos de la terminación `.1.ht2, .2.ht2`, etc.  
  
* **Descarga del archivo de anotaciones de referencia GRCh38 desde el repositorio _ENSEMBLE_**   
  
Para descargar el archivo de anotaciones de referencia de la especie _H.sapiens_, se emplea el repositorio [ENSEMBL](https://www.ensembl.org/Homo_sapiens/Tools/FileChameleon). Además, es importante que el archivo de anotaciones tenga ciertas características específicas necesarias para la correcta ejecución de programas posteriores. Es por ello, que se emplea la herramienta _File Chamaleon_ con el fin de formatear el archivo de anotaciones.  
  
Una vez seleccionado el genoma de interés (GRCh38.p14), se descarga el archivo de anotaciones en formato GTF y se marca la casilla de _transcript_id_ para incluir este campo en el archivo descargado.  
  
![image](https://github.com/user-attachments/assets/5141e4d9-41f4-4362-a1b2-ec78c2e01690)  
  
A continuación descargamos y descomprimimos el archivo de anotaciones en formato GTF en nuestra computadora, y lo almacenamos en la ruta `~/RNAseq_analysis/Data/5_Annotation/`.  
```console
# Descomprisión del archivo descargado
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
  
![image](https://github.com/user-attachments/assets/a7046b80-04b2-4e10-aca4-9e97654c851d)  
  
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
Hoy en día, existen diferentes alineadores o mapeadores que se pueden emplear. Entre los más conocidos se encuentran HISAT2 y STAR. Algunas de las diferencias entre ellos son la cantidad de recursos necesarios o la precisión de los resultados. HISAT2 es un alineador menos preciso en comparación con STAR, aunque es adecuado cuando hay menos recursos disponibles.  
  
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
28427141 reads; of these:      # Lecturas totales (x2= 56854282)
  28427141 (100.00%) were paired; of these:
    1123549 (3.95%) aligned concordantly 0 times           # Lecturas no conformes (x2 = 2247098)
    27303592 (96.05%) aligned concordantly exactly 1 time  # Lecuturas con mapeo conforme y con una ubicación (x2 = 54607184)
    0 (0.00%) aligned concordantly >1 times		   # Lecturas con mapeo conforme y con más de una ubicación
    ----
    1123549 pairs aligned concordantly 0 times; of these:
      88122 (7.84%) aligned discordantly 1 time     # Lecturas mapeadas pero discordantes (x2 = 176244)
    ----
    1035427 pairs aligned 0 times concordantly or discordantly; of these:
      2070854 mates make up the pairs; of these:   
        1094218 (52.84%) aligned 0 times	    # Lecturas sin mapear
        976636 (47.16%) aligned exactly 1 time      # Lecturas con solo un extremo alineado 
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
  
Para ordenar todos los archivos resultantes del alineamiento formato BAM, se puede emplear el siguiente script `bam_order.sh` desde el directorio  ~/RNAseq_analysis/Code.   
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
  
La indexación de los archivos BAM, genera un archivo índice complementario con un tamaño mucho menor, que actúa como una tabla de contenidos y permite identificar dónde se localizan las lecturas. Esto es importante en la ejecución de algunas programas, ya que permite localizar partes específicas del archivo BAM manera más rápida.  

#### 2.3.3 Calidad del alineamiento   
Una vez alineadas las lecturas procesadas contra el genoma de referencia, existen herramientas para comprobar la calidad de dicho alineamiento. Para ello, vamos a emplear el programa _RseQC_, dentro del cual existen diferentes funciones y scripts creados para analizar el alineamiento.  
En concreto, vamos a emplear el script `bam_stat.py`, para obtener un resumen de las estadísticas del mapeo del archivo BAM. Primero, se determina una calidad de mapeo para cada lectura y seguidamente se calcula la probabilidad de que esa lectura esté mal posicionada en función de un umbral mínimo.   

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

# Script de bash para el programa RseQC; estadísticas de la calidad del alineamiento con la herramienta bam_stats.py

SAMPLES="SRR28380565 SRR28380566 SRR28380568 SRR28380570 SRR28380572 SRR28380573
	SRR28380580 SRR28380582	SRR28380584 SRR28380586	SRR28380588 SRR28380589"

cd ~/RNAseq_analysis/Data/4_Alignment/

for SAMPLE in $SAMPLES; do
	bam_stat.py -i ${SAMPLE}.sorted.bam > ./statistics/${SAMPLE}.bamstats.txt
done
```
Como resultado se obtiene un archivo de texto plano con diferentes estadísticos del alineamiento. Por ejemplo, en el caso del archivo SRR28380565 se obtienen los siguientes resultados:   
```console
#==================================================
#All numbers are READ count
#==================================================

Total records:                          56854282

QC failed:                              0
Optical/PCR duplicate:                  0
Non primary hits                        0
Unmapped reads:                         1094218
mapq < mapq_cut (non-unique):           1663991

mapq >= mapq_cut (unique):              54096073
Read-1:                                 27138442
Read-2:                                 26957631
Reads map to '+':                       27053055
Reads map to '-':                       27043018
Non-splice reads:                       40899546
Splice reads:                           13196527
Reads mapped in proper pairs:           53037578
Proper-paired reads map to different chrom:0
```

Finalmente, los resultados obtenidos tras el alineamiento de todas las muestras son los siguientes:  
  
![image](https://github.com/user-attachments/assets/0b49ef3a-0d83-4844-8ac0-034608e37af9)
  
  
### 2.4 Identificación y recuento de features o características  

#### 2.4.1 Preparación del archivo de anotaciones GTF

Como explicamos previamente, el genoma de anotaciones de la especie humana (GRCh38.p14) en formato GTF se descargó desde el repositorio [ENSEMBL](https://www.ensembl.org/Homo_sapiens/Tools/FileChameleon).  

La estructura del archivo es la siguiente:   
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
  
El formato GTF representa las características y anotaciones específicas para un genoma de referencia.  Cada línea representa una región particular del genoma, y el número de columnas puede variar. Se pueden tener hasta 10 columnas, las cuales representan los siguientes indicadores:  
```
[cromosoma] [fuente] [característica] [posición.inicio] [posición.final] [puntuación] [hebra] [marco] [atributos] [comentarios]
```
   
| Columna | Significado |
| -------| ------------|
1. Cromosoma | Indica el cromosoma donde se localiza la anotación. |
2. Fuente | Fuente de la anotación. Hace referencia al programa de predicción empleado o a la base de datos de anotaciones. |
3. Característica | Indica el tipo de característica: gen, transcrito, exón, CDS, etc. |
4. Posición de inicio | Posición de inicio de la característica en el genoma de referencia. |
5. Posición final |  Posición final de la característica en el genoma de referencia. |
6. Puntuación | Puntuación asociada a la característica. Si no hay puntuación, se presenta solamente un '.' |
7. Hebra | Permite indicar la hebra de procedencia de la característica.  Puede ser '+' , '-' o '.' (si  la hebra es desconocida o no aplicable |
8. Marco | Puede presentarse como '0', '1' o '2' |
9. Atributos | Incluyen identificadores y otras información suplementaria |
10. Comentarios | Es una columna con información complementaria |
  

#### 2.4.2 Recuento de características con htseq-count

Una vez visualizada la estructura del archivo de anotaciones (_Homo_sapiens.gtf_), este se emplea para anotar las características de las lecturas según la posición del genoma dónde alineen. Para ello, se debe comprobar que la forma de nombrar los cromosomas en ambos archivos es la misma, ya que si no esto puede conducir a errores.

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
```

Estructura ejemplo SRR28380565:
```console
ENSG00000000003	0
ENSG00000000005	0
ENSG00000000419	1207
ENSG00000000457	702
.
.
.
ENSG00000310556	0
ENSG00000310557	0
__no_feature	3102207
__ambiguous	1097442
__too_low_aQual	1595129
__not_aligned	165990
__alignment_not_unique	0
```


#### 2.4.3 Obtención de la matriz de recuentos

Finalmente 
```console
#!/usr/bin/bash

SAMPLES="SRR28380566 SRR28380565 SRR28380570 SRR28380572 SRR28380573"
cd ~/RNAseq_analysis/Results


for SAMPLE in $SAMPLES; do
	sed -i "1s/^/feature\t${SAMPLE}\n/" "${SAMPLE}_counts.tsv" 
done

```



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
