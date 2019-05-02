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
      
      `conda config --add channels defaults`
      
      `conda config --add channels bioconda`
      
      `conda config --add channels conda-forge`
      
      Install SMEG
      
      `conda install smeg=1.1`

      Reload .bashrc environment `source ~/.bashrc`

# USAGE

    Usage:
    smeg build_species <options>   Build species database
    smeg growth_est <options>      Estimate strain-specific growth rate
    smeg -v                        Version
    smeg -h                        Display this help message
    
### build_species module #

The species database is built using strains of a species of interest. Strains are typically downloaded from NCBI Genbank but custom strains can be used. **Downloaded strains MUST contain at least one COMPLETE reference genome**. In the absence of a complete reference genome, the draft genome with the least fragmentation should be reordered (e.g. with Mauve software) using a strain from a closely related species. Also, genome names should not exceed 35 characters. For species having > 700 strains (e.g. *E. coli*), it is advisable to build the database using only strains with a complete genome. 

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

**NOTE:** While the optimal value of the SNP assignment threshold will vary depending on the species being analyzed, SMEG has an “auto” option **(-a flag)**, in which different threshold values are tested in parallel and output, giving the user the flexibility to select desired parameters and the associated database. Here, output databases are stored in folders named either T.<num> or F.<num> where T and F represent "ignore iterative clustering" and "do not ignore iterative clustering" respectively. <num> is the SNP assignment threshold. The summary statistics is stored in "log.txt". Typically, a threshold yielding the highest number of unique SNPs with a high SNP assignment threshold is prefarable.   
      
Pre-compiled species database are available from **ftp://ftp.jax.org/ohlab/SMEG_species_database/**


### build_rep module # 

SMEG first identifies the strains (using pathoscope) present in your dataset (i.e collection of metagenomic reads) to create a dataset-specific representative database. The strain of a cluster with the highest median relative abundance is chosen as the representative strain for that cluster. SMEG assumes strains in a cluster would have similar growth rates and can be represented by an individual member. The -d flag is used to specify the clustering method to select representative strains; deepSplit = 0 would create **fewer representative strains** while deepSplit = 4 generates **many representative strains**.  These representative strains are then used to create the dataset-specific database. **NOTE: delimiter seperating paired reads must be the underscore ( _ ) symbol.  Reads must also have the .fastq (and not .fq) extension.**  

Additionally, if you have *a priori* knowledge of the strains present in your samples, you can simply use the -c flag to specify a file listing your strains which are directly used to create the database.    

### growth_est module #

For the final step, SMEG estimates growth rate of representatives strains using gene-based or SNP-based method. The choice of method (-m flag) depends on the deepSplit method in the previous step (i.e build_rep module). **If deepSplit = 4 option (default) was used to generate representative strains, SNP-based method (default) should be used.** Gene-based method, which is much faster and requires fewer algorithmic steps, is only suitable for deepSplit method = 0. Nonetheless, SNP-based method is applicable to ALL deepSplit options. Strains with coverage below the cutoff (-c flag) are considered as non-replicating and would have a SMEG score of 1.     

*Summary*

*DeepSplit = 0 ---> Gene-based or SNP-method*

*DeepSplit = 4 ---> SNP-method*

*Gene-based method --> Fast method, estimates growth rates for fewer strains*

*SNP-method --> Slower method, estimates growth rates for more strains*

# Output

For every sample, a table of results (.txt) displaying strain-specific growth rate (SMEG) and genome coverage is generated. If -e flag is set, all tables will be merged into a single matrix file called "merged_table.txt" and a heatmap (.pdf) displaying growth rates (SMEG) across all samples with hierachical clustering is generated.

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
