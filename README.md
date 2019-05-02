# SMEG
**Strain-level MEtagenomic Growth estimation (SMEG)** measures growth rates of bacterial subspecies or strains from complex metagenomic
samples. SMEG is capable of identifying novel or uncharacterized strains in a given sample prior to growth estimation. 

SMEG pipeline consists of two modules;

1. < build_species > - generates a database for a specie of interest using its member strains
2. < growth_est > - measure strain-specific growth rates in your dataset using either a de novo or reference-based approach

# INSTALLATION
Installation of SMEG is through anaconda/miniconda. Please follow the exact installation guidelines provided in this section.

1.    Ensure you have **gcc compiler >=4.8.5**
              
2.    **Install SMEG**

          Please set up channels in the following order. NOTE that conda-forge has the highest priority.
          
          conda config --add channels defaults
          conda config --add channels bioconda
          conda config --add channels conda-forge

          Install SMEG
          conda install smeg=1.1
          
          Reload .bashrc environment `source ~/.bashrc`

# USAGE

    Usage:
    smeg build_species <options>   Build species database
    smeg growth_est <options>      Estimate strain-specific growth rate
    smeg -v                        Version
    smeg -h                        Display this help message
    
### build_species module #

The species database is built using strains of a species of interest. Strains are typically downloaded from NCBI Genbank but custom strains can be used. **Downloaded strains MUST contain at least one COMPLETE reference genome**. In the absence of a complete reference genome, the draft genome with the least fragmentation should be reordered (e.g. with Mauve software) using a strain from a closely related species. Also, genome names should not exceed 35 characters. For species having > 700 strains (e.g. *E. coli*), it is advisable to build the database using only strains with a complete genome. Strains must have .fna, .fa, or .fasta extensions. 

For convenience, we provided a script `download_genomes.sh` to retrieve and rename genomes from NCBI Genbank. Simply edit lines 2 and 3 to specify the output directory and species name, respectively, and run the script.

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

A core-genome phylogeny is constructed using the provided strains which is used to assign strains into clusters. Outlier strains, defined as having pairwise distances 30 times above the median are excluded, as these genomes may have been misclassified taxonomically, or may contain contaminant contigs. Generally, the underlying biological assumption is that strains constituting a cluster have high phylogenetic relatedness and will have similar phenotypic properties like growth rate in a sample. 

For each cluster, SMEG identifies cluster-specific unique SNPs, i.e. SNPs that are shared between a given proportion of cluster members, but absent in all strains from other clusters. This proportion is referred to as the SNP assignment threshold **(-s flag)**. For instance, setting the SNP assignment threshold to 0.8 indicates that a unique SNP will only be identified if it is shared between at least 80% of the cluster members and absent in strains from other clusters. If the total number of unique SNPs of a given cluster is below a user-specified threshold **(-t flag)**, SMEG iteratively subclassifies the cluster and infers unique SNPs for each sub cluster. SMEG repeats the sub-classification step for a maximum of 3 iterations or until the threshold is met.  Next,
SMEG retrieves the coordinates of the cluster-specific SNPs in a representative genome, which is randomly selected from the cluster (after favoring for genome completeness) or can be user-defined **(-r flag)**. If the value provided by the -r flag is absent from the genomes or a draft-genome, SMEG defaults to auto-select. 

**NOTE:** While the optimal value of the SNP assignment threshold will vary depending on the species being analyzed, SMEG has an “auto” option **(-a flag)**, in which different threshold values are tested in parallel and output, giving the user the flexibility to select desired parameters and the associated database. Here, output databases are stored in folders named either T.<num> or F.<num> where T and F represent "ignore iterative clustering" and "do not ignore iterative clustering" respectively. <num> is the SNP assignment threshold. The summary statistics is stored in "log.txt". Typically, a threshold yielding the highest number of unique SNPs with a high SNP assignment threshold is prefarable. Specifying the -a flag overrides user-defined options from -s -t -i and -e flags.  
      
Pre-compiled species database are available from **ftp://ftp.jax.org/ohlab/SMEG_species_database/**

### growth_est module #

    smeg growth_est <options>
        <options>

    MAIN OPTIONS
      -r         Reads directory (single-end reads)
      -x         Sample filename extension (fq, fastq, fastq.gz) [default fastq]
      -o         Output directory
      -s         Species database directory
      -m  INT    SMEG method (0 = de novo-based method, 1 = reference-based method) [default = 0]
      -c  FLOAT  Coverage cutoff (>= 1) [default 1]
      -u  INT    Minimum number of SNPs to estimate growth rate [default = 50]
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

The reference-based approach assumes prior knowledge of strain composition which may have been determined using other tools. Here, SMEG requires the user to provide a file of genome names (if genomes already exist in the species database) **(-g flag)** or file listing **full path** to DESMAN-reconstituted genome sequences **(-a flag)**. **SMEG assumes DESMAN-reconstituted haplotypes are core genes and of the same length and order**. Note that incorrect strain identification will impact SMEG’s accuracy, because only user-supplied strains are used to estimate growth rate.

# Output

SMEG outputs four statistics for a given sample; (i) the median SMEG score from 3 imputations, (ii) the coverage of phylogenetic clusters (in the de novo approach) or the user provided strains (in the reference-based approach), (iii) the number of non-zero SNP sites used for the analysis, and (iv) the range of SMEG with different imputations. strains/clusters with SNPs below the minimum cutoff **(-u flag)** will be assigned a SMEG score of 1. Finally, if -e flag is set, all output will be merged into a single matrix file called "merged_table.txt" and a heatmap (.pdf) displaying growth rates (SMEG) across all samples with hierachical clustering is generated. 

# Example usage

Assuming our downloaded strains (e.g *E. coli* strains) are present in a directory called `mygenomes`, we can build our species database by running the following command.

`smeg build_species -g mygenomes -p 16 -o E_coli_species_database`

To estimate microbial growth rate from a given dataset of interest, we first build a representative-strains database prior to growth estimation. In this example, we assume all our paired-end reads are present in a directory called `Reads_folder` and growth predictions would be stored in a directory called `Results`.

`smeg build_rep -r Reads_folder -s E_coli_species_database -o E_coli_rep_database -d 4 -p 16` # SNP-method

`smeg growth_est -r Reads_folder -s E_coli_species_database -d E_coli_rep_database -o Results -p 16 -e`

OR

`smeg build_rep -r Reads_folder -s E_coli_species_database -o E_coli_rep_database -d 0 -p 16` # gene-based method

`smeg growth_est -r Reads_folder -s E_coli_species_database -d E_coli_rep_database -o Results -m 1 -p 16 -e`


# DEPENDENCIES
- R 
    - Required R libraries - 
    dplyr,
    getopt,
    ggplot2,
    gsubfn,
    gplots,
    ape,
    dynamicTreeCut,
    WGCNA,
    data.table
    
- gcc 
- OrthoANIu
- usearch 
- parallel 
- Mauve 
- Roary 
- Prokka 
- Bowtie2 
- Pathoscope2 
- samtools >= 1.5 
- bamtools
- featureCounts
