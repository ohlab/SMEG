# SMEG
**Strain-level MEtagenomic Growth estimation (SMEG)** measures growth rates of microbial strains from complex metagenomic dataset. SMEG is able to accurately identify strain(s) associated with different disease conditions even at low coverage (1X). Its high sensitivity to complex admixtures of closely related strains make it applicable to study outbreak epidemiology, infection dynamics, and understanding commensal strain diversity to guide future probiotics development.

SMEG pipeline consists of three modules;

1. < build_species > - generates a database for a specie of interest using its member strains.
2. < build_rep > - generates a dataset-specific representative strains database using strains that are present in a given dataset of interest   
3. < growth_est > - measure strain-specific growth rates in your dataset

# INSTALLATION
Installation of SMEG is through anaconda/miniconda. Please follow the exact installation guidelines provided in this section.

1.    Ensure you have **gcc compiler >=4.8.5**

2.    **Download OrthoANIu** from https://www.ezbiocloud.net/tools/orthoaniu and add the folder of the downloaded file (i.e. OAU.jar) to your path (echo 'export PATH=$PATH:/path/to/folder' >> ~/.bashrc)  

3.    **Install usearch**  - If you do not have usearch already installed, 

      download usearch from `https://drive5.com/usearch/download.html` 
      
      Rename the downloaded filename to simply "usearch". For example, to rename usearch11.0.667_i86linux32, simply run the following command  `rename usearch11.0.667_i86linux32 usearch usearch11.0.667_i86linux32`
      
      Add execute permissions `chmod +x usearch`
      
      Add folder to path `echo 'export PATH=$PATH:/path/to/folder' >> ~/.bashrc`
      
      -- Test installation --
      
      `usearch --version`
      
4.    **Install pathoscope** **(requires python)**

      `wget https://github.com/PathoScope/PathoScope/archive/v2.0.6.tar.gz`

      `tar xvf v2.0.6.tar.gz`

      `cd PathoScope-2.0.6/`

      `python2.7 setup.py install`

      -- Test the installation --

      `pathoscope -h`

      **NOTE: do NOT install pathoscope via conda as its samtools dependency conflicts with that of SMEG**

                      
4.    **Install SMEG**

      Please set up channels in the following order. NOTE that conda-forge has the highest priority.
      
      `conda config --add channels defaults`
      
      `conda config --add channels bioconda`
      
      `conda config --add channels conda-forge`
      
      Install SMEG
      
      `conda install smeg=1.0.1=h2d50403_2`

      Reload .bashrc environment `source ~/.bashrc`

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
    -m  INT    SMEG method (0 = SNP-method, 1 = gene-based method) [default = 0]
    -t  FLOAT  Theta value for bin size (bin size = no of unique SNPs x theta)
               Not compatible with gene-based method (i.e. -m 1)  [default 0.06]
    -l         Path to file listing a subset of reads for analysis
               [default = analyze all samples in Reads directory]
    -e         merge output tables into a single matrix file and generate heatmap
    -p  INT    Number of threads [default 1]
    -h         Display this message


### build_species module #

The species database is built using strains of a species of interest. Strains are typically downloaded from NCBI Genbank but custom strains can be used. **Downloaded strains MUST contain at least one COMPLETE reference genome. It is advisable to rename your strains from the conventional names accompanying NCBI genomes. For instance, names such as `GCA_000160335.2_ASM16033v2_genomic.fna` should be renamed**.  For species having > 700 strains (e.g. *E. coli*), it is advisable to build the database using only strains with a complete genome. 

For convenience, we provided a script `download_genomes.sh` to retrieve and rename genomes from NCBI Genbank. Simply edit lines 2 and 3 to specify the output directory and species name, respectively, and run the script.

In addition, a subset of strains can be used for database building by specifying a file listing specific strains via the -l flag. SMEG aligns the strains to generate a phylogenetic tree which is then used to group strains into clusters. Outlier strains, defined as having pairwise distances 30 times above the median, are excluded as these may be misclassified genomes or they may contain contaminant contigs. Two different types of clustering output are generated. In the first output (clusters_deepSplit0), strains are grouped in the absence of cluster splitting sensitivity (deepSplit = 0) which generates **fewer and bigger clusters**. In the second output (clusters_deepSplit4), **many smaller clusters** are created in the presence of cluster splitting sensitivity (deepSplit = 4).   

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

To estimate microial growth rate from a given dataset of interest, we first build a representative-strains database prior to growth estimation. In this example, we assume all our paired-end reads are present in a directory called `Reads_folder` and growth predictions would be stored in a directory called `Results`.

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
