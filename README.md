# HSP-Protein_functions_rUTI-ECCMID2024

### Protein detection and function assingment

The nextflow script performs the analysis from Raw reads to annotated tables

The Phageome analysis has been performed using the nextflow available at Shell_analysis.nf in this repository. To perform the analysis install in your linux shell the following program:

[Nextflow](https://github.com/nextflow-io/nextflow)

[BBduk](https://github.com/BioInfoTools/BBMap) from BBmap tools.

[Bowtie2](https://github.com/BenLangmead/bowtie2)

[diamond](https://github.com/bbuchfink/diamond)

[famli](https://github.com/DerrickWood/kraken2](https://github.com/FredHutch/FAMLI)

[Upimapi](https://github.com/iquasere/UPIMAPI)

### Database needed

[Human genome (GRch38)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/)

[Uniref90](https://www.uniprot.org/help/uniref)

### To run in conda environtment

    ##### Create the conda environtment
    conda env create Proteome_pipeline
    conda activate Proteome_pipeline
    ###### Install all the utils
    conda install -c bioconda nextflow
    conda install -c bioconda bbmap
    conda install -c bioconda bowtie2
    conda install -c bioconda diamond
    conda install -c bioconda upimapi
    #### Install from github famli https://github.com/FredHutch/FAMLI

### Statistical analysis


All the estatistical analysis it is available at Protein_function.R
