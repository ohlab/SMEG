# SMEG
**Strain-level MEtagenomic Growth estimation (SMEG)** measures growth rates of microbial strains from complex metagenomic dataset. SMEG is able to accurately identify strain(s) associated with different disease conditions even at low coverage (1X). Its high sensitivity to complex admixtures of closely related strains make it applicable to study outbreak epidemiology, infection dynamics, and understanding commensal strain diversity to guide future probiotics development.

SMEG pipeline consists of three modules;

1. < build_species > - generates a database for a specie of interest using its member strains.
2. < build_rep > - generates a dataset-specific representative strains database using strains that are present in a given dataset of interest   
3. < growth_est > - measure strain-specific growth rates in your dataset

# INSTALLATION
Installation of SMEG is through anaconda/miniconda. Please follow the exact installation guidelines provided in this section.

1.    Ensure you have **gcc compiler**

2.    **Install usearch**  - If you do not have usearch already installed, download and install usearch from https://drive5.com/usearch/download.html and follow the installation steps detailed here https://www.drive5.com/usearch/manual/install.html

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


**build_species module** 

Pre-compiled species database are available from xxxxxx. 

The species database is built using strains of a species of interest. Strains are typically downloaded from NCBI but custom strains can be used. For species having > 700 strains (e.g. E. coli), it is advisable to build the database using strains with a complete genome. Also, a subset of strains can be used for database building by specifying a file listing specific strains via the -l flag. SMEG aligns the strains to generate a phylogenetic tree which is then used to group strains into clusters. Two different types of clustering output are generated. In the first output (clusters_deepSplit0), strains are grouped in the absence of cluster splitting sensitivity (deepSplit = 0) which generates **fewer and bigger clusters**. In the second output (clusters_deepSplit4), **many smaller clusters** are created in the presence of cluster splitting sensitivity (deepSplit = 4).  


**build_rep module** 

SMEG first identifies the strains present in your dataset to create a dataset-specific representative database. SMEG uses pathoscope to determine the   




If reads are only available in paired-end format, use either of the mate pairs, or concatenate both pairs into a single fastq file. In addition, reads must have the .fastq extension (and not .fq). Using either modules, the default is to analyze all samples present in the reads directory. However, analysis can be restricted to a subset of samples by using the -l flag and specifying a file that lists the subset of samples.    

For the 'multiplex' module, reads mapping to multiple genomes are reassigned using Pathoscope 2 when the -p flag is set. The degree to which reads are reassigned is set by the -t (theta prior) flag. The theta prior value represents the number of non-unique reads that are not subject to reassignment. Finally, when the coverage cutoff (-c flag) is set below 1, only genomes with fragmentation levels below 90 fragments/Mbp are analyzed (see xxx et al. for more details). **Note that to use the 'multiplex' module, you must have downloaded the GRiD database from ftp://ftp.jax.org/ohlab/Index/. However, you do not need the database to run the example test**.

# OUTPUT
`single module` - two output files are generated
- A plot (.pdf) showing coverage information across the genome 
- A table of results (.txt) displaying growth rate (GRiD), 95% confidence interval, unrefined GRiD value, species heterogeneity, genome coverage, *dnaA/ori* ratio, and *ter/dif* ratio. Species heterogeneity is a metric estimating the degree to which closely related strains/species contribute to variance in growth predictions (range between 0 - 1 where 0 indicate no heterogeneity). In most bacteria genomes, *dnaA* is located in close proximity to the *ori* whereas replication typically terminates at/near *dif* sequence. Thus, the closer *dnaA/ori* and *ter/dif* ratios are to one, the more likely the accuracy of GRiD scores.  

`multiplex module` - two output files are generated per sample
- A heatmap (.pdf), displaying growth rate (GRiD) from genomes above the coverage cutoff with hierachical clustering. 
- A table of results (.txt) displaying growth rate (GRiD) of genomes above the coverage cutoff, unrefined GRiD value, species heterogeneity, and genome coverage. If -m flag is set, all tables will be merged into a single matrix file called "merged_table.txt".


# Example test
The test sample contain reads from *Staphylococcus epidermids*, *Lactobacillus gasseri*, and *Campylobacter upsaliensis*, each with coverage of ~ 0.5. Download the GRiD folder and run the test as shown below. 

**NOTE: Ensure GRiD-1.0.6 was installed via conda. Please see installation instructions**  

`wget https://github.com/ohlab/GRiD/archive/1.0.6.tar.gz`

`tar xvf 1.0.6.tar.gz`

`cd GRiD-1.0.6/test`

`chmod +x ../grid`

`../grid single -r . -g S_epidermidis.LRKNS118.fna`

`../grid multiplex -r . -d . -p -c 0.2`

For each module, output files (a pdf and a text file) are generated in the test folder.


# Updating GRiD database 
The database can be updated with metagenomic bins or newly sequenced bacterial genomes by running the update_database script (requires bowtie2).
 

    update_database <options>
    <options>
    -d      GRiD database directory
    -g      Bacterial genomes directory
    -n      Name for new database
    -l      Path to file listing specific genomes
            for inclusion in database [default = include all genomes in directory]
    -h      Display this message

NOTE: bacteria genomes must be in fasta format and must have either .fasta, .fa or .fna extensions.


# DEPENDENCIES
- R (tested using v 3.2.3) 
    - Required R libraries - 
    dplyr,
    getopt,
    ggplot2,
    gsubfn,
    gplots
    
- bowtie2 (tested using v 2.3.1)
- seqtk (tested using v 1.0-r31)
- samtools (tested using v 1.5)
- bedtools (tested using v 2.26.0)
- bamtools (tested using v 1.0.2)
- blast (tested using v 2.6.0)
- Pathoscope2    
