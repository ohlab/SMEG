# SMEG
**Strain-level Metagenomic Estimation of Growth rate (SMEG)** measures growth rates of bacterial subspecies or strains from complex metagenomic
samples. SMEG is capable of identifying novel or uncharacterized strains in a given sample prior to growth estimation. 

SMEG pipeline consists of two modules;

1. < build_species > - generates a database for a specie of interest using its member strains
2. < growth_est > - measure strain-specific growth rates in your dataset using either a de novo or reference-based approach

# INSTALLATION

   OPTION 1 - **Singularity**
   
    - # Retrieve the image file
      wget ftp://ftp.jax.org/ohlab/SMEG_installation/smeg.sif
      
    - # Run a quick test
      singularity exec smeg.sif smeg -h
      
   OPTION 2 - **Bioconda**
   
      Please follow the EXACT installation guidelines provided in this section.
       - # Ensure you have gcc compiler >=4.8.5. 
       
          Please set up channels in the following order. NOTE that conda-forge has the highest priority. 
          conda config --add channels defaults
          conda config --add channels bioconda
          conda config --add channels conda-forge

          To avoid compatibility issues with dependencies, we recommend creating a new conda environment for SMEG 
          prior to installation e.g. "conda create --name SMEG" and "source activate SMEG"
          
       - # Install SMEG
          conda install smeg=1.1.5 prokka=1.11 r-base=3.5.1
          
       - # Also, "rename" is not default on some Linux-machines and needs to be installed
          sudo apt-get install -y rename
          
**It is highly recommended you run the example test to ensure proper installation before running SMEG on your dataset. 
Also, make sure you read the "coverage requirements" section before interpreting your results**. 

# USAGE

    Usage:
    smeg build_species <options>   Build species database
    smeg growth_est <options>      Estimate strain-specific growth rate
    smeg -v                        Version
    smeg -h                        Display this help message
    
### build_species module #

The species database is built using strains of a species of interest. Strains are typically downloaded from NCBI Genbank but custom strains can be used. **Downloaded strains MUST contain at least one COMPLETE reference genome**. In the absence of a complete reference genome, the draft genome with the least fragmentation should be reordered (e.g. with Mauve software) using a strain from a closely related species. For species having > 700 strains (e.g. *E. coli*), it is advisable to build the database using only strains with a complete genome. **Strains must have .fna, .fa, or .fasta extensions**. 

For convenience, we provided a script `download_genomes.sh` to retrieve and rename genomes from NCBI Genbank. Simply edit lines 3 and 4 to specify the output directory and species name, respectively, and run the script.

    smeg build_species <options>
        <options>
        -g        Genomes directory
        -o        Output directory
        -l        File listing a subset of genomes for database building
                    [default = use all genomes in 'Genomes directory']
        -p INT    Number of threads [default 4]
        -s FLOAT  SNP assignment threshold (range 0.1 - 1) [default 0.6]
        -t INT    Cluster SNPs threshold for iterative clustering [default 50]
        -i        Ignore iterative clustering
        -a        Activate auto-mode
        -r        Representative genome [default = auto select Rep genome]
        -k        Keep Roary output [default = false]
        -e        Create database ONLY applicable with Reference-based SMEG method
                    [default = generate database suitable for both de novo and ref-based methods]
        -h        Display this message

A core-genome phylogeny is constructed and used to assign strains into clusters. The underlying biological assumption is that strains constituting a cluster have high phylogenetic relatedness and will have similar phenotypic properties like growth rate in a sample. Outlier strains, defined as having pairwise distances 30 times above the median are excluded, as these genomes may have been misclassified taxonomically, or may contain contaminant contigs. 

For each cluster, SMEG identifies cluster-specific unique SNPs, i.e. SNPs that are shared between a given proportion of cluster members, but absent in all strains from other clusters. This proportion is referred to as the SNP assignment threshold **(-s flag)**. For instance, setting the SNP assignment threshold to 0.8 indicates that a unique SNP will only be identified if it is shared between at least 80% of the cluster members and absent in strains from other clusters. If the total number of unique SNPs of a given cluster is below a user-specified threshold **(-t flag)**, SMEG iteratively subclassifies the cluster and infers unique SNPs for each sub cluster. SMEG repeats the sub-classification step for a maximum of 3 iterations or until the threshold is met.  Next,
SMEG retrieves the coordinates of the cluster-specific SNPs in a representative genome, which is randomly selected from the cluster (after favoring for genome completeness) or can be user-defined **(-r flag)**. If the representative genome provided with the -r flag is absent or a draft-genome, SMEG defaults to auto-select. 

**----- "auto" (-a) option -----**

While the optimal value of the SNP assignment threshold will vary depending on the species being analyzed, SMEG has an “auto” option **(-a flag)**, in which different threshold values are tested in parallel and output, giving the user the flexibility to select desired parameters and the associated database. Here, output databases are stored in folders named T.{num} or F.{num} where T and F represent "ignore iterative clustering" and "do not ignore iterative clustering" respectively. {num} is the SNP assignment threshold. The summary statistics is stored in `log.txt`. Specifying the -a flag overrides user-defined options for -s -t -i and -e flags.  Also, strains and their corresponding cluster identity will be located in `clusterOutput.txt`

Below is an example of contents of `log.txt` when the auto option is activated. Ideally, we want a database with a high SNP assignment threshold and sufficient unique SNPs for most clusters. In this example, the preferred database will be that created with a SNP assignment threshold of 0.7 with iterative clustering because SMEG could generate sufficient unique SNPs for all clusters. 


![alt text](https://github.com/ohlab/SMEG/blob/master/smeg_logfile_example.png)
      

***Pre-compiled databases for some bacterial species commonly encountered in human-associated metagenomes can be retrieved from ftp://ftp.jax.org/ohlab/SMEG_database/***

### growth_est module #

    smeg growth_est <options>
        <options>

    MAIN OPTIONS
      -r         Reads directory (single-end reads)
      -x         Sample filename extension (fq, fastq, fastq.gz) [default fastq]
      -o         Output directory
      -s         Species database directory
      -m  INT    SMEG method (0 = de novo-based method, 1 = reference-based method) [default = 0]
      -c  FLOAT  Coverage cutoff (>= 0.5) [default 0.5]
      -u  INT    Minimum number of SNPs to estimate growth rate [default = 100]
      -l         Path to file listing a subset of reads for analysis
                 [default = analyze all samples in Reads directory]

    DE-NOVO BASED APPROACH OPTIONS
      -d  FLOAT  Cluster detection threshold (range 0.1 - 1) [default = 0.2]
      -t  FLOAT  Sample-specific SNP assignment threshold (range 0.1 - 1) [default = 0.6]

    REFERENCE BASED APPROACH OPTIONS
      -g         File listing reference genomes for growth rate estimation
      -a         FIle listing FULL PATH to DESMAN-resolved strains in fasta format (core-genes)
      -n  INT    Max number of mismatch [default = use default bowtie2 threshold]

    OTHER OPTIONS
      -e         merge output tables into a single matrix file and generate heatmap
      -p  INT    Number of threads [default 4]
      -h         Display this message


Growth rate is estimated using either a de novo or reference-based estimation approach. The de novo approach assumes no prior knowledge of strain composition in a sample. Here, we assume that uncharacterized strains in a sample can be assigned to a cluster using information on cluster-specific SNPs. In a given sample, a cluster is deemed present if the proportion of unique SNPs with coverage > 0 exceeds the ‘cluster detection threshold’ **(-d flag)**. In scenarios where a sample does not contain all clusters in
the species database, SMEG further generates a sample-specific SNP profile based solely on the identified clusters in a sample (threshold controlled by **-t flag**). This step increases the number of SNP sites for growth estimation by reducing the number of clusters compared and is especially useful for clusters lacking sufficient unique SNPs in the species database. 

The reference-based approach assumes prior knowledge of strain composition which may have been determined using other tools. Here, SMEG requires the user to provide a file of genome names (if genomes already exist in the species database) **(-g flag)** or file listing **full path** to DESMAN-reconstituted genome sequences **(-a flag)**. 

An example of an input file using the -g flag is as follows:
      
      E_coli_O157_H7.fna
      E_coli_K-12.fna
whereas, a typical file input using the -a flag would be: 
      
      /path/to/folder/haplotype_0.fa
      /path/to/folder/haplotype_1.fa
      /path/to/folder/haplotype_2.fa

**SMEG assumes DESMAN-reconstituted haplotypes are core genes and of the same length and order [this is usually the default output for DESMAN haplotypes anyway]**. Note that incorrect strain identification will impact SMEG’s accuracy using this option, because only user-supplied strains are used to estimate growth rate.

**------ COVERAGE REQUIREMENTS ------**

SMEG can accurately detect clusters at up to 0.5x coverage. However, it requires cluster coverage of 5x and 0.5x for microbes with high- and low- within-species genetic diversity, respectively, to accurately estimate growth rate. We recommend a 5x cutoff without a priori knowledge of the genomic characteristics of the species of interest.

# Output

SMEG outputs four statistics for a given sample; (i) the median SMEG score from 3 imputations, (ii) the coverage of phylogenetic clusters (in the de novo approach) or the user provided strains (in the reference-based approach), (iii) the number of non-zero SNP sites used for the analysis, and (iv) the range of SMEG with different imputations. strains/clusters with SNPs below the minimum cutoff **(-u flag)** will be assigned a SMEG score of 1. Finally, if -e flag is set, all output will be merged into a single matrix file called "merged_table.txt" and a heatmap (.pdf) displaying growth rates (SMEG) across all samples with hierachical clustering is generated.

**Please note that strain coverage is calculated using only available unique SNPs and thus, the output coverage may not be the most accurate**.  

# Example test

**Example 1 - de novo estimation**

This example involves 12 strains of *Faecalibacterium prausnitzii* and the image below shows the resultant phylogenetic tree.

![alt text](https://github.com/ohlab/SMEG/blob/master/smeg_tree.png)

We have provided a mock metagenomic sample containing 2 strains (indicated with red-arrows) which we will exclude from database generation. The excluded strains will act as hypothetical uncharacterized or novel strains. We will thus, create the database using 10 strains. The expected *ori/ter* ratios for CNCM_I_4543.fna and AF10-13.fna are ~1.1 and ~1.8 respectively in the sample. To save runtime, we reoredered majority of the draft genomes.

    wget https://github.com/ohlab/SMEG/archive/1.1.1.tar.gz
    tar xvf 1.1.1.tar.gz
    cd SMEG-1.1.1/test
    
    # First, let's create the species database
    smeg build_species -g . -o test_database -a -p 16
    
    # The 'auto' option is activated and you should have different database directories created using different 
    # parameters. You will note that in `test_database/log.txt`, all parameters resulted in the generation of 
    # sufficient unique SNPs for all clusters. Thus, we will select the database generated with the highest SNP 
    # assignment threshold (e.g. `test_database/T.0.9`). You can also evaluate strains and their corresponding 
    # cluster identity from your selected database e.g. `test_database/T.0.9/clusterOutput.txt`

    # Next, estimate growth rate
    smeg growth_est -o Result_denovo -r . -s test_database/T.0.9 -p 8 -x fastq.gz
**Example 2 - estimation using DESMAN-reconstituted haplotypes**

Using the same excluded strains as example 1 above, we extracted and ordered core genes to represent hypothetical DESMAN-reconstituted haplotypes.
      
    smeg growth_est -m 1 -o Result_DESMAN -r . -a desman/desman_haplotypes.txt -p 8 -x fastq.gz -n 3

### DEPENDENCIES #
    gcc, GNU parallel, Mauve, roary, prokka, bowtie2, samtools, bamtools, bedtools, blastn
    
    R libraries 
    (dplyr, getopt, ggplot2, gsubfn, gplots, ape, dynamicTreeCut, seqinr, data.table)
