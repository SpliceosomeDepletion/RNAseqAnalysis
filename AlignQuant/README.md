# Scripts to run RNAseq analysis of Spliceosome inhibited vs. WT Cal51 cells

## General Info

The following scripts are written to run on the ETH HPC "Euler", for running them in any other environment they will need to be adapted


## Alignment and Transcript quantification

Genome information was downloaded from Gencode v29 (date: 2019/01/30)
```sh
wget -b 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz'
wget -b 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh38.primary_assembly.genome.fa.gz'
```
