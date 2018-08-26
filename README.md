# SMEG
**Strain-level MEtagenomic Growth estimation (SMEG)** measures growth rates of microbial strains from complex metagenomic dataset. SMEG is able to accurately identify strain(s) associated with different disease conditions even at low coverage (1X). Its high sensitivity to complex admixtures of closely related strains make it applicable to study outbreak epidemiology, infection dynamics, and understanding commensal strain diversity to guide future probiotics development.

SMEG pipeline consists of three modules;

1. < build_species > - generates a database for a specie of interest using its member strains.
2. < build_rep > - generates a dataset-specific representative strains database using strains that are present in a given dataset of interest   
3. < growth_est > - measure strain-specific growth rates in your dataset

# INSTALLATION
Installation of SMEG is through anaconda/miniconda. Please follow the exact installation guidelines provided in this section.

1.    Ensure you have **gcc compiler**

2.    **Install usearch**  - If you do not have usearch already installed, download and install usearch from https://drive5.com/usearch/download.html and follow the installation steps detailed in this link https://www.drive5.com/usearch/manual/install.html

3.    **Install pathoscope**
 
      `wget https://github.com/PathoScope/PathoScope/archive/v2.0.6.tar.gz`
      
      `tar xvf v2.0.6.tar.gz`
      
      `cd PathoScope-2.0.6/`
      
      `python setup.py install`
      
      -- Test the installation --
      
      `pathoscope -h`
      
      **NOTE: do NOT install pathoscope via conda as its samtools dependency conflicts with that of SMEG**
                
4.    **Install SMEG**

      `conda install smeg`

It is highly recommended to run the example test to ensure proper installation before running SMEG on your dataset.

# USAGE

    smeg build_species <options>   Build species database
    smeg build_rep <options>       Build dataset-specific representative strains database
    smeg growth_est <options>      Estimate strain-specific growth rate
    smeg -v                        Version
    smeg -h                        Display this help message


    smeg build_species <options>
    <options>
    -g      Genomes directory
    -o      Output directory
    -l      Path to file listing a subset of genomes for
            database building [default = use all genomes in 'Genomes directory']
    -p INT  Number of threads [default 1]
    -h      Display this message


    smeg build_rep <options>
    <options>
    -r         Reads directory (paired-end reads)
    -o         Output directory
    -s         Species database directory
    -d  INT    deepSplit for grouping strains into clusters (0 or 4) [default = 4]
    -l         Path to file listing a subset of reads to estimate representative
               strains [default = use all samples in Reads directory]
    -c         Path to file with customized list of strains for representative database (ONLY use this option IF
               you are certain your chosen strains belong to different clusters and are present in your dataset)
    -p  INT    Number of threads [default 1]
    -h         Display this message


    smeg growth_est <options>
    <options>
    -r         Reads directory (paired-end reads)
    -o         Output directory
    -s         Species database directory
    -d         Representative strains database directory
    -c  FLOAT  Coverage cutoff (>= 1) [default 1]
    -m  INT    SMEG method (0 = SNP-method, 1 = gene-dosage method) [default = 0]
    -t  FLOAT  Theta value for bin size (bin size = no of unique SNPs x theta)
               Not compatible with gene-dosage method (i.e. -m 1)  [default 0.06]
    -l         Path to file listing a subset of reads for analysis
               [default = analyze all samples in Reads directory]
    -e         merge output tables into a single matrix file and generate heatmap
    -p  INT    Number of threads [default 1]
    -h         Display this message


### build_species module #

Pre-compiled species database are available from xxxxxx. 

The species database is built using strains of a species of interest. Strains are typically downloaded from NCBI but custom strains can be used. For species having > 700 strains (e.g. *E. coli*), it is advisable to build the database using strains with a complete genome. Also, a subset of strains can be used for database building by specifying a file listing specific strains via the -l flag. SMEG aligns the strains to generate a phylogenetic tree which is then used to group strains into clusters. Two different types of clustering output are generated. In the first output (clusters_deepSplit0), strains are grouped in the absence of cluster splitting sensitivity (deepSplit = 0) which generates **fewer and bigger clusters**. In the second output (clusters_deepSplit4), **many smaller clusters** are created in the presence of cluster splitting sensitivity (deepSplit = 4).  


### build_rep module # 

SMEG first identifies the strains (using pathoscope) present in your dataset to create a dataset-specific representative database. The strain of a cluster with the highest median relative abundance is chosen as the representative strain for that cluster. SMEG assumes strains in a cluster would have similar growth rates and can be represented by an individual member. The -d flag is used to specify the clustering method to select representative strains; deepSplit = 0 would create **fewer representative strains** while deepSplit = 4 generates **many representative strains**.  These representative strains are then used to create the dataset-specific database. **NOTE: delimiter seperating paired reads must be the underscore ( _ ) symbol.  Reads must also have the .fastq (and not .fq) extension.**  

Additionally, if you have *a priori* knowledge of the strains present in your samples, you can use the -c flag to specify a file listing your strains which are directly used to create the database.    

### growth_est module #

For the final step, SMEG estimates growth rate of representatives strains using gene dosage or SNP-based method. The choice of method (-m flag) depends on the deepSplit method in the previous step (i.e build_rep module). **If deepSplit = 4 option was used to generate representative strains, SNP-based method (default) should be used.** Gene dosage method, which is much faster and requires fewer algorithmic steps, is only suitable for deepSplit method = 0. Nonetheless, SNP-based method is applicable to ALL deepSplit options.   

*Summary*

*DeepSplit = 0 ---> Gene dosage or SNP-method*

*DeepSplit = 4 ---> SNP-method*

*Gene dosage method --> Fast method, estimates growth rates for fewer strains*

*SNP-method --> Slower method, estimates growth rates for more strains*

# Output

For every sample, a table of results (.txt) displaying strain-specific growth rate (SMEG) and genome coverage is generated. If -e flag is set, all tables will be merged into a single matrix file called "merged_table.txt" and a heatmap (.pdf) displaying growth rates (SMEG) across all samples with hierachical clustering is generated.


# Example test

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
