<img src="https://user-images.githubusercontent.com/52743495/173834089-526c540a-df4b-452f-964e-26104bf6f261.png" width="350" />

## A Snakemake pipeline for the genome-wide discovery of jointly regulated CpGs (JRCs) in pooled whole genome bisulfite sequencing (WGBS) data. 
#### Benjamin Planterose Jiménez, Brontë Kolar, Manfred Kayser, Athina Vidaki
Department of Genetic Identification, Erasmus MC University Medical Center Rotterdam, Rotterdam, the Netherlands.

This repository together with [Binokulars](https://github.com/BenjaminPlanterose/Binokulars), [JRC_sorter](https://github.com/BenjaminPlanterose/JJRC_sorter) 
and [JRC_downstream_analysis ](https://github.com/BenjaminPlanterose/JRC_downstream_analysis) are part of the "Simultaneous mapping of epigenetic inter-haplotype, inter-cell and inter-individual 
variation via the discovery of jointly regulated CpGs in pooled sequencing data" publication (under review).


#### Tested on:
```
- Operating system: Ubuntu 18.04.6 LTS (Bionic Beaver)
- R version: 3.6.1 (2019-07-05) -- "Action of the Toes" and 4.2.2 Patched (2022-11-10 r83330) -- "Innocent and Trusting"
- Python version: Python 3.9.12
- Organism: Homo sapiens (3.1 GB reference genome). For a 98 GB BAM file, it takes approximately 1 day (using 20 cores). 

Do not attempt to run without at least 50 GB of RAM.
```

**IMPORTANT**: Test data available at [Zenodo](https://zenodo.org/record/7625657#.Y-UOORzMJu0).


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

The test data consists of pooled whole blood WGBS data from old and young men from the study of [Laurentino et al](https://onlinelibrary.wiley.com/doi/full/10.1111/acel.13242). Download test data from [Zenodo](https://zenodo.org/record/7625657#.Y-UOORzMJu0) and uncompress. This data is not included in the Github repository since its size exceeds the 100 MB limit. This includes:

* ```sample_data.bam``` - Alignment of WGBS reads from a pooled whole blood experiment. Only reads at the beginning of chr12 have been included (for the sake of timely debugging). This bam file has been sorted and indexed (.csi file).
* ```reference_genome``` directory - includes .fa sequence for chromosome 12, already indexed by BISCUIT.
* ```chromosomes.txt``` - specifies on which chromosomes to run JRC_seeker (in this case, only chr12).
* ```test_config.json``` - contains the parameters and directory paths for the JRC_seeker Snakemake pipeline.
* ```expected_output.zip``` - contains the expected output of JRC_seeker on this data for debugging purposes.

To run the example, edit ```test_config.json``` with absolute paths. Specifically, adjust:

* DIRECTORIES
	* ```output_folder``` - where you want the results to be stored. This directory must exist already before running JRC_seeker.
	* ```path_to_jrc_seeker``` - path to git cloned jrc_seeker directory.
* PROGRAMMES
	* ```chromhmm``` - path to the ChromHMM.jar (as part of the ChromHMM software).
* FILES
	* ```path_to_config_file``` - path to the test_config.json file from ```sample_data```.
	* ```path_to_reference_genome``` - path to chr12.fa from ```sample_data```.
	* ```path_to_bam``` - path to sample_data.bam from ```sample_data```.
	* ```chromosomes_file``` - path to chromosomes.txt from ```sample_data```.
	* ```path_to_chrom_length_file``` - path to chromosome lengths for your organism and assembly of choice. See ```/ChromHMM/CHROMSIZES/``` for examples (as part of the ChromHMM software).
	* ```path_to_mappability_file``` - BED files for low mappability exclusion. We provide examples for humans (hg19, hg38) at ```JRC_seeker/assets/mappability_files/```.
	* ```path_to_blacklist``` - BED files for regions with anomalous coverage. We provide examples for humans (hg19, hg38) at ```JRC_seeker/assets/blacklist_regions/```.
* PARAMETERS - Make sure to adapt ```binokulars_cores``` to whatever number of threads are available in your machine. All other parameters must remain as default for this test.


To run the pipeline on the example, activate the jrc_seeker conda environment and run Snakemake as follows:

```
conda activate jrc_seeker
snakemake -s <path_to_JRC_seeker>/Snakefile --cores 1 --configfile <path_to_sample_data>/test_config.json
```

This pipeline should deliver the same results as stored in ```expected_output.zip``` from the ```sample_data``` and in the absence of errors.


## Pipeline Overview

The pipeline pre-selects for intermediately-methylated regions (IMRs) and runs a read-level randomization test (i.e. Binokulars) to identify jointly-regulated CpGs (JRCs). The pipeline runs the following steps:

1. **Methylation calling** - It uses BISCUIT (biscuit pileup, vcf2bed, mergecg) to call CpG methylation (both strands are combined).
2. **Binarization of methylation data** - Custom python scripts are used to create a methylated/unmethylated binned binary signals.
3. **Genome segmentation** - A hidden Markov model (HMM) representation of the data is obtained with 4 states (no data, unmethylated, methylated, intermediately methylated) via ChromHMM.
4. **Genome segmentation polishing** - A custom R-script is employed to remove small intermediately-methylated regions (IMRs) and regions collocating with low-mappability or blacklisted regions.
5. **Pile-up parsing** - The whole methylation data is processed to count number of methylated/unmethylated Cs per read. This is done with BISCUIT and bash commands.
5. **JRC statistical test (Binokulars)** - Our custom permutation test, implemented in R, computes permutation p-values. When $\text{p-value} < 1/N_\text{iter}$, a parametric approximation is employed.

For more information on Binokulars (Binomial likelihood function-based bootstrap hypothesis test for co-methylation within reads), visit [here](https://github.com/BenjaminPlanterose/Binokulars).


## How to run on your own data?

Please make sure that you have succesfully installed all dependencies and ran the example on the ```test_data```.

### 1. Prepare input files

#### Reference Genome & BAM files

Starting from fastq files, it will be necessary to perform sequence alignment. We used the BISCUIT aligner as described [here](https://huishenlab.github.io/biscuit/docs/alignment). 
Briefly, we make an index of the reference genome and align as:
```
biscuit index my_reference.fa
biscuit align /path/to/my_reference.fa read1.fq.gz read2.fq.gz | \
    samblaster | samtools sort --write-index -o my_output.bam -O BAM -
```

**IMPORTANT**: The output bam must be sorted and indexed. The command above does this already.


#### Chromosomes file

JRC_seeker requires a text file containing chromosomes found in the BAM file. To generate it, you may use the following Linux shell commands:

```
samtools idxstats your_file.bam | grep 'chr' |  cut -f 1 | tee temp.txt
awk 'NR==FNR{A[$1];next}$1 in A' /ChromHMM/CHROMSIZES/hg19.txt temp.txt > chromosomes.txt # here we filter out chromosomes not included in the chromosome size file from ChromHMM.
rm temp.txt
```

#### Mappability Files and Blacklisted Regions

In JRC_seeker, we remove potentially artifactual regions to avoid false positive JRCs.

A list of regions that are uniquely mappable by at least one 100-mer can be found under ```/assets/mappability_files``` as part of the JRC_seeker software (only for *Homo sapiens*). 
These were downloaded from [Hoffman's lab website](https://bismap.hoffmanlab.org/).

A list of "blacklisted" regions with anomalous coverage can be found under ```/assets/blacklist_regions```. For hg19, the blacklist file was created by combining regions from the following repositories:

* DAC blacklisted regions ([DBR - ENCFF001TDO.bed](https://www.encodeproject.org/annotations/ENCSR636HFF/))

* Duke Excluded Regions ([DER - ENCFF001THR.bed](https://www.encodeproject.org/annotations/ENCSR797MUY/))

The hg38 version is simply the liftover of the hg19 file, obtained with [UCSC LiftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver).


### 2. Edit the configuration file

As for the example, edit the configuration file```test_config.json``` with absolute paths. Specifically, adjust:

* DIRECTORIES
	* ```output_folder``` - where you want the results to be stored. This directory must exist already before running JRC_seeker.
	* ```path_to_jrc_seeker``` - path to git cloned jrc_seeker directory.
* PROGRAMMES
	* ```chromhmm``` - path to the ChromHMM.jar (as part of the ChromHMM software).
* FILES
	* ```path_to_config_file``` - path to the test_config.json file.
	* ```path_to_reference_genome``` - path to the reference genome (.fa).
	* ```path_to_bam``` - path to the alignment file (.bam).
	* ```chromosomes_file``` - path to the chromosomes.txt file.
	* ```path_to_chrom_length_file``` - path to chromosome lengths. See ```/ChromHMM/CHROMSIZES/``` for examples (as part of the ChromHMM software).
	* ```path_to_mappability_file``` - BED files for low mappability exclusion. See ```JRC_seeker/assets/mappability_files/```.
	* ```path_to_blacklist``` - BED files for regions with anomalous coverage. See ```JRC_seeker/assets/blacklist_regions/```.
* PARAMETERS
	* ```sample_name```: name of your sample (don't use spaces or the following characters: "/" "," "." "\").
	* Binarization and binning.
		* ```bin_size```: size of bins for binarization and ChromHMM segmentation. Measured in base pairs. (default: 200).
		* ```lower_im_methylation_bound```: intermediately methylated methylation value for lower boundary for binarization. Recommended not to change (default: 0.2).
		* ```upper_im_methylation_bound```: intermediately methylated methylation value for upper boundary for binarization. Recommended not to change (default: 0.8).
		* ```data_assignment_setting``` : if (un)methylated counts are on the bin boundary, they are added to the earlier or "left" bin. Recommended not to change (default: "left", "right" is other option).
		* ```k_threshold```: threshold for number of reads needed for bin to not be set to a no-data state (see "Binarize methylation data" section). (default: 3)

	* ChromHMM
		* ```chromhmm_it```: number of ChromHMM iterations. Recommend 500 to ensure convergence, but ChromHMM default for maxiterations in the LearnModel command is 200. (default: "500").
		* ```n_states```: number of states for ChromHMM to predict. Binarize subroutine depends on 4 states. Recommended not to change. (default: "4").
		* 
	* BinPolish
		* ```binpolishflank```: .
		* ```segment_min_sz``` : BinPolish discards regions equal or smaller to this threshold size. (default: 200)
	* Binokulars
		* ```permutation_iterations```: number of permutations that binokulars runs per region. (default: 1000).
		* ```seed_binokulars```: binokulars seed value. (default: 4).
		* ```binokulars_cores```: number of cores that the binokulars subroutine uses. (default: 4).
		* ```flank_length```: number of base pairs that binokulars flanks regions by. Recommend not to change (default: 500).



### 3. Run JRC_seeker


To make sure there are no problems with your configuration file, do a dry run of the snakemake pipeline:

```
conda activate jrc_seeker
snakemake -n -s <path_to_JRC_seeker>/Snakefile --configfile <path_to_sample_data>/test_config.json
```

To run the JRC_seeker pipeline, simply run the command below, where the number of cores and the path to your edited configuration file is specified.
```
conda activate jrc_seeker
snakemake -s <path_to_JRC_seeker>/Snakefile --cores [number of cores] --configfile <path_to_sample_data>/test_config.json
```

*Note:* We recommended to run with ```--cores 1``` and specify ```binokulars_cores``` in the ```config.json``` file (the most intensive step). If the number of snakemake cores is above 1, the number of cores the final binokulars step uses will be the product of the snakemake cores and the value of the ```binokulars_cores``` field in the ```config.json``` file.



## References


Ernst J., Kellis M. (**2012**). ChromHMM: automating chromatin-state discovery and characterization. *Nat Methods*. 9(3):215-6. [doi: 10.1038/nmeth.1906](https://www.nature.com/articles/nmeth.1906). 

Faust G.G., Hall I.M. (**2014**). SAMBLASTER: fast duplicate marking and structural variant read extraction. *Bioinformatics*. 30(17):2503-5. [doi: 10.1093/bioinformatics/btu314](https://academic.oup.com/bioinformatics/article/30/17/2503/2748175).

Karimzadeh M., Ernst C., Kundaje A., Hoffman M.M. (**2018**). Umap and Bismap: quantifying genome and methylome mappability. *Nucleic Acids Res*. 46(20):e120. [doi: 10.1093/nar/gky677](https://academic.oup.com/nar/article/46/20/e120/5086676).

Laurentino S. *et al* (**2020**). A germ cell-specific ageing pattern in otherwise healthy men. *Aging Cell*. 19(10):e13242. [doi: 10.1111/acel.13242](https://onlinelibrary.wiley.com/doi/full/10.1111/acel.13242).


