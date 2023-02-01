

<img src="https://user-images.githubusercontent.com/52743495/173834089-526c540a-df4b-452f-964e-26104bf6f261.png" width="350" />

**JRC_seeker** is a Snakemake pipeline for the genome-wide discovery of jointly regulated CpGs (JRCs) in pooled whole genome bisulfite sequencing (WGBS) data. 

#### Tested on:
```
- Operating system: Ubuntu 18.04.6 LTS (Bionic Beaver)
- R version: 3.6.1 (2019-07-05) -- "Action of the Toes" and 4.2.2 Patched (2022-11-10 r83330) -- "Innocent and Trusting"
- Python version: Python 3.9.12
- Organism: Homo sapiens (3.1 GB reference genome). For a 98 GB BAM file, it takes approximately 1 day (using 20 cores). 

Do not attempt to run without at least 50 GB of RAM.

Test data available at Zenodo ($URL$)
```

## Dependencies

#### [conda](https://www.anaconda.com/products/individual)

Download and install conda by:
```
wget https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh
bash Anaconda3-2021.11-Linux-x86_64.sh
```

#### [mamba](https://github.com/mamba-org/mamba)

```
conda install -n base -c conda-forge mamba
```
Activate the Conda base environment (which now includes Mamba).
```
conda activate base
```
Installation sometimes fails; nevertheless, it is possible to install all dependencies without Mamba (more instructions below).

#### [snakemake](https://snakemake.readthedocs.io/) (at least v4.3.1)

Create a new conda environment called ```jrc_seeker``` with python 3.9 in it.
```
mamba create -c conda-forge -c bioconda -n jrc_seeker snakemake python=3.9
```

**Option:**
If Mamba failed to install, simply create a conda environment as follows and install snakemake with the pip command:
```
conda create -n jrc_seeker python=3.9
conda activate jrc_seeker
pip install snakemake
```

Activate the ```jrc_seeker``` conda environment.
```
conda activate jrc_seeker
```
Check whether Snakemake is succesfully installed by running the following command:
```
snakemake --help
```

#### [biopython-1.79](https://biopython.org/docs/1.79/api/Bio.html)

```
conda install -c conda-forge biopython=1.79
```
or
```
pip install "biopython==1.79"
```

#### [tabix](https://github.com/samtools/htslib)

```
conda install -c bioconda tabix
```
or
```
pip install tabix
```

#### [pandas](https://pandas.pydata.org/)

```
conda install -c anaconda pandas
```
or
```
pip install pandas
```

#### [bgzip](https://github.com/xbrianh/bgzip.git)

```
conda install -c bioconda bgzip
```
or
```
pip install bgzip
```

#### [ChromHMM](http://compbio.mit.edu/ChromHMM/)

Quick instructions on downloading ChromHMM:

1. Install Java 1.5 or later if not already installed.
2. Download and unzip the ChromHMM.zip file using the following code:

```
wget http://compbio.mit.edu/ChromHMM/ChromHMM.zip
unzip ChromHMM.zip
```

#### [bedtools](https://bedtools.readthedocs.io/en/latest/)

```
conda install -c bioconda bedtools
```
or
```
pip install bedtools
```

#### [samtools](http://www.htslib.org/doc/samtools.html)

Recommended version: 1.14

```
conda install -c bioconda samtools
```
or
```
pip install "samtools==1.14"
```

#### [samblaster](https://github.com/GregoryFaust/samblaster)

```
conda install -c bioconda samblaster
```
or
```
pip install samblaster
```

#### [BISCUIT](https://huishenlab.github.io/biscuit/)

```
conda install -c bioconda biscuit
```

#### [R](https://cran.r-project.org/)

To install R, click on the hyperlink above and follow instructions. To download dependency R-packages, run the following:

Open R by running:
```bash
R
```
And run the following R commands:

```r
# To install
install.packages('data.table')
install.packages('parallel')
install.packages('MASS')
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.16")
BiocManager::install('GenomicRanges')

# To verify that installation was succesful
library(data.table)
library(parallel)
library(MASS)
library(GenomicRanges)
```

#### JRC_Seeker

Clone the repository with the following command:

```
mkdir jrc_seeker
cd jrc_seeker
git clone https://github.com/b-kolar/jrc_seeker.git
```

# Test run

Download test data from [here](http://compbio.mit.edu/ChromHMM/) and uncompress. This data is not included in the Github repository since its size exceed the 100 MB limit. This includes:

* ```sample_data.bam``` - Alignment of WGBS reads from a pooled whole blood experiment. Only reads at the beginning of chr12 have been included (for the sake of timely debugging). This bam file has been sorted and indexed (.csi file).
* ```reference_genome``` directory - includes .fa sequence for chromosome 12, already indexed by BISCUIT.
* ```chromosomes.txt``` - specifies on which chromosomes to run JRC_seeker (in this case, only chr12).
* ```test_config.json``` - contains the parameters and directory paths for the JRC_seeker Snakemake pipeline.
* ```expected_output.zip``` - contains the expected output of JRC_seeker on this data for debugging purposes.

To run the example, edit ```test_config.json``` with absolute paths. Specifically, adjust:

* ```output_folder``` - where you want the results to be stored. This directory must exist already before running JRC_seeker.
* ```path_to_jrc_seeker``` - path to git cloned jrc_seeker directory.
* ```chromhmm``` - path to the ChromHMM.jar script (as part of the ChromHMM software).
* ```path_to_config_file``` - path to the test_config.json file from ```sample_data```.
* ```path_to_reference_genome``` - path to the chr12.fa file from ```sample_data```.
* ```path_to_bam``` - path to the sample_data.bam from ```sample_data```.
* ```chromosomes_file``` - path to the chromosomes.txt from ```sample_data```.
* ```path_to_chrom_length_file``` - path to chromosome lengths for your organism and assembly of choice. See ```/ChromHMM/CHROMSIZES/``` for examples (as part of the ChromHMM software).
* ```path_to_mappability_file``` - BED files for low mappability exclusion. We provide examples for humans (hg19, hg38) at ```JRC_seeker/assets/mappability_files/```.
* ```path_to_blacklist``` - - BED files for regions with anomalous coverage. We provide examples for humans (hg19, hg38) at ```JRC_seeker/assets/blacklist_regions/```.
* ```PARAMETERS``` - Make sure to adapt ```binokulars_cores``` to whatever number of threads are available in your machine.


To run the pipeline on the example, activate the jrc_seeker conda environment and run Snakemake as follows:

```
conda activate jrc_seeker
snakemake -s <path_to_JRC_seeker>/Snakefile --cores 1 --configfile <path_to_sample_data>/test_config.json
```

This pipeline should deliver the same results as stored in ```expected_output.zip``` from the ```sample_data``` and in the absence of errors.


## Pipeline overview Overview

The pipeline pre-selects for intermediately-methylated regions (IMRs) and runs a read-level randomization test (i.e. Binokulars) to identify jointly-regulated CpGs (JRCs). The pipeline runs the following steps:

1. **Methylation calling** - It uses BISCUIT (biscuit pileup, vcf2bed, mergecg) to call CpG methylation (both strands are combined).
2. **Binarization of methylation data** - Custom python scripts are used to create a methylated/unmethylated binned binary signals.
3. **Genome segmentation** - A hidden Markov model (HMM) representation of the data is obtained with 4 states (no data, unmethylated, methylated, intermediately methylated) via ChromHMM.
4. **Genome segmentation polishing** - A custom R-script is employed to remove small intermediately-methylated regions (IMRs) and regions collocating with low-mappability or blacklisted regions.
5. **Pile-up parsing** - The whole methylation data is processed to count number of methylated/unmethylated Cs per read. This is done with BISCUIT and bash commands.
5. **JRC statistical test (Binokulars)** - The permutation test implemented in R gives rise to permutation p-values. When $\text{p-value} < 1/N_\text{iter}$, a parametric approximation is employed.

For more information on Binokulars (Binomial likelihood function-based bootstrap hypothesis test for co-methylation within reads), visit [here](https://github.com/BenjaminPlanterose/Binokulars).


## How to run on your own data?

Please make sure that you have succesfully installed all dependencies and run the example on the ```test_data```. To run JRC_seeker on your own data, follow the next steps:

### 1. Prepare input files

#### BAM 

A sorted, indexed BAM with duplicate marked reads is required to run this analysis. We used BISCUIT to align our data and produced an index (.csi) for the BAM file using:

```samtools sort --write-index -o my_output.bam -O BAM -```

#### Reference Genome

A reference genome is required to run JRC_seeker and JRC_seeker has asset files supporting the use of the hg19 and hg38 reference genomes.

The reference genome you use must be indexed using the BISCUIT ```index``` command. To download and index a reference genome, you can run the following command (for hg19):

```
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/latest/hg19.fa.gz
gunzip hg19.fa.gz
biscuit index hg19.fa
```

Be aware that indexing a reference genome can take time, but only needs to be done once for each reference.

#### Chromosomes file

The Snakemake pipeline requires a text file that contains a list of chromosomes found in the BAM file. To generate it, use the following Linux shell commands:

```
samtools idxstats your_file.bam | grep 'chr' |  cut -f 1 | tee temp.txt
awk 'NR==FNR{A[$1];next}$1 in A' /ChromHMM/CHROMSIZES/hg19.txt temp.txt > chromosomes.txt
rm temp.txt
```

##### Mappability Files

Regions with low mappability can result in false positive JRC regions, which is why regions with low mappability are removed during BinPolish. In the ```/assets/mappability_files``` directory, two files are included (for hg19 and hg38) that contain a list of regions of the genome that are uniquely mappable by at least one k-mer (in our case, k=100).

These files are the Bismap individual k-mer files for Human hg19 and Human hg38 genomes (k100 Single-read), which were downloaded from: https://bismap.hoffmanlab.org/

##### Blacklist Regions

Regions that overlap with regions have been labelled as "blacklist" regions are also removed during BinPolish. For hg19, a blacklist file was created by combining regions from the following two blacklist region datasets:

1) DAC blacklisted regions ([DBR - ENCFF001TDO.bed](https://www.encodeproject.org/annotations/ENCSR636HFF/))

2) Duke Excluded Regions ([DER - ENCFF001THR.bed](https://www.encodeproject.org/annotations/ENCSR797MUY/))

For hg38, the hg19 file was lifted over using [UCSC LiftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver).

Further details on the LiftOver settings and the blacklist regions files are listed in the ```/assets/blacklist_regions/README.md```.

### Edit the configuration file

Snakemake uses a configuration file to locate external files and parameter values. Make a copy of the sample configuration file ```/sample_data/test_config.json``` and edit the file and directory paths to point towards the respective input files and directories on your machine (use absolute paths). Furthermore, be sure to edit the parameters in the configuration file. 

The configuration file is organized into three sections, for readability:
_DIRECTORIES_: Directories to add your own path to.
_FILES_: Files to add your own path to.
_PARAMETERS_: Parameters to change. The listed default values for parameters that are listed at the end of this tutorial are the suggested values per parameter. 

Here are a couple of important reminders:
- If running JRC_seeker on your entire BAM file, be sure to change the region value to none:
```"region" : "none"```
- Ensure you change the number of binokulars cores to an amount your machine can handle, such as 4:
```"binokulars_cores" : 4```

A detailed list of the parameters in the configuration file are found at the end of this tutorial.

### Run JRC_seeker

Go into the JRC_seeker directory by:
```
cd <Path to JRC_seeker>
```

Make sure to activate the conda environment by: 
```
conda activate jrc_seeker
```

To make sure there are no problems with your configuration file, do a dry run of the snakemake pipeline:

```
snakemake -n --configfile [path to config file]
```

To run the JRC_Seeker pipeline, simply run the command below, where the number of cores and the path to your edited configuration file is specified.
```
snakemake --cores [amount of cores] --configfile [path to config file]
```

*Note:* It is recommended to run with 1 core and specify in the ```config.json``` file the number of cores to run the binokulars step with. If the number of snakemake cores is above 1, the number of cores the final binokulars step uses will be the product of the snakemake cores and the value of the ```binokulars_cores``` field in the ```config.json``` file.

## A note on time 

The final step of the snakemake pipeline is the actual binokulars run. From a 98 GB BAM file of pooled blood samples of 12 individuals and a 3.1 GB reference genome (entire hg19 reference genome), ~400,000 regions are found to be intermediately methylated regions after BinPolish. It takes approximately 2 seconds for binokulars to process each region. As a result, for 4 cores we expect a binokulars run to take ~55 hours. Using 20 cores, this takes ~11 hours. 

The binokulars step is the bottleneck of this pipeline, but it is also the most parallelizable. The remainder of the pipeline takes under 10 hours for the above mentioned dataset. We recommend using ~20 cores to run the entire pipeline within a day.

## Overview of sample data

In the ```/sample_data``` folder you can find some sample files to test if JRC Seeker is running properly on your machine. Below you can find an explanation of where these files come from and how we created them, just in case you were interested:

**sample_data.bam - Sample BAM file**

This BAM file is of data from chromosome 20 of the pooled blood data from old and young men from the following study: https://www.ebi.ac.uk/ena/browser/view/PRJEB28044?show=reads.

```
samtools view -b pooled_blood_hg38.bam chr20 > chr_20.bam
```

To reduce the file size, The first million aligned reads were selected using the following command:

```
samtools view -h chr_20.bam \
    | head -n 1000000 \
    | samtools view -b -o sample_data.bam
```

**sample_data.bam.csi - BAM index**

The BAM file was indexed using the following command:

```
samtools index -c sample_data.bam
```

**chromosomes.txt - chromosomes list**

A list of chromosmes, as outlined earlier. As this comprised of only data from chromosome 20, the chromosomes.txt file is simply one line:

```
chr20
```

**test_config.json - sample config file**

The sample config file contains the settings for this run. As we are only looking at chromosome 20, the region is set to ```"chr20"```.

**Reference genome**

In the ```/reference_genome``` directory, you can find ```chr20.fa```, a FASTA file of the reference genome of hg38 for chromosomome 20, which was downloaded [from this Genome Reference Consortium source](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/) using the following command: 

```
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz
gunzip chr20.fa.gz
```

The FASTA file was indexed by BISCUIT using the following command, which generated the additional files:

```
biscuit index chr20.fa
```

## Config file parameters

DIRECTORIES:
```
output_folder : folder all files will be outputted in (absolute path)

path_to_scripts : path to JRC_seeker scripts folder

temp_folder : some computations might be too large for default temp folders. Put the default temp folder path here to add a new path (absolute path)

path_to_jrc_seeker : absolute path to JRC_seeker (e.g. /home/opt/jrc_seeker). This is the folder where the Snakefile is located.
```

FILES:
```
path_to_config_file : path to Snakemake config.json file (absolute path)

path_to_reference_genome : path to FASTA reference genome (absolute path)

path_to_bam : path to BAM file (absolute path)

path_to_chrom_length_file : path to chromosome length text file used by ChromHMM. Found in the CHROMSIZES folder of ChromHMM. Use the file corresponding the the genome version your are using (e.g. hg19.txt) (absolute path)

chromhmm : path to ChromHMM jar file (absolute path)

path_to_mappability_file : path to mappability file. For hg19 or hg38, use the ones in the assets/mappability_files folder. Otherwise, add the absolute path to the version for your genome.

chromosomes_file : path to chromosomes.txt file outlined in the "Input Files" step of this tutorial (absolute path)

path_to_blacklist : path to blacklist regions file. For hg19 or hg38, use the ones in the assets/blacklist_regions folder. Otherwise, add the absolute path to the version for your genome.
```

PARAMETERS:
```
sample_name : name of your sample (don't use spaces or the following characters: "/" "," "." "\")

bin_size : size of bins for binzarization and ChromHMM segmentation. Measured in base pairs. Recommended not to change. (default: 200)

binokulars_output_directory : name of directory for binokulars output (don't use spaces or the following characters: "/" "," "." "\")

region : if a specific region of the BAM file is to be investigated OR if the BAM file contains one chromosome, specify this here (e.g. "chr1" or "chr20:0-1000"). Only one region is permitted and if changing this setting, ensure that the chromosomes.txt file is manually created (see instructions above). In most circumstances, the entire BAM file is to be processed and this should thus be set to "none". (default: "none")

lower_im_methylation_bound : intermediately methylated methylation value for lower boundary for binarization. Recommended not to change (default: 0.2)

upper_im_methylation_bound : intermediately methylated methylation value for upper boundary for binarization. Recommended not to change (default: 0.8)

data_assignment_setting : if (un)methylated counts are on the bin boundary, they are added to the earlier or "left" bin. Recommended not to change (default: "left", "right" is other option)

k_threshold : threshold for number of reads needed for bin to not be set to a no-data state (see "Binarize methylation data" section). (default: 3)

n_states : number of states for ChromHMM to predict. Binarize subroutine depends on 4 states. Recommended not to change. (default: "4")

chromhmm_it : number of ChromHMM iterations. Recommend 500 to ensure convergence, but ChromHMM default for maxiterations in the LearnModel command is 200. (default: "500")

map_threshold : threshold overlap coverage (as a percent) of intermediately methylated regions overlapping with low mappability regions for them to be discarded. Recommend not to change (default: 0.95)

segment_min_sz : BinPolish discards regions equal or smaller to this threshold size. (default: 200)

permutation_iterations : number of permutations that binokulars runs per region. (default: 1000)

seed_binokulars : binokulars seed value. (default: 4)

binokulars_cores : number of cores that the binokulars subroutine uses. (default: 4)

flank_length : number of base pairs that binokulars flanks regions by. Recommend not to change (default: 500)
```

## Sources

ChromHMM: automating chromatin-state discovery and characterization, Nature methods, 2012, Jason Ernst & Manolis Kellis

Faust, G.G. and Hall, I.M., “SAMBLASTER: fast duplicate marking and structural variant read extraction,” Bioinformatics Sept. 2014; 30(17): 2503-2505.

Amemiya, H.M., Kundaje, A. & Boyle, A.P. The ENCODE Blacklist: Identification of Problematic Regions of the Genome. Sci Rep 9, 9354 (2019). https://doi.org/10.1038/s41598-019-45839-z

Karimzadeh M, Ernst C, Kundaje A, Hoffman MM. 2018. Umap and Bismap: quantifying genome and methylome mappability. doi: https://doi.org/10.1093/nar/gky677 Nucleic Acids Research, Volume 46, Issue 20, 16 November 2018, Page e120.
