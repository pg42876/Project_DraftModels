# Development of metabolic models

## Installation of tools

### MERLIN
merlin is available for download in a single multi-platform version in https://merlin-sysbio.org/download/.

### CARVEME
conda install -c conda-forge mamba
mamba install -c bioconda carveme

### AUREME
docker pull dyliss/aureme-img:2.4

```

## Run tools


```
### CARVEME
carve genome.faa --output model.xml

### AUREME
docker run -it -v C:\Users\Desktop:/outside dyliss/aureme-img:2.4 bash
docker start 40a155090560au
docker exec -it 40a155090560 bash

aureme --init=Formicicum

cp outside/species.faa /shared/Formicicum/genomic_data/

mkdir /shared/Formicicum/orthology_based_reconstruction/mhun

awk '{print $1}' /shared/Formicicum/genomic_data/prot_sequences.fasta > /shared/Formicicum/genomic_data/Formicicum.faa

cp outside/species.faa /shared/Formicicum/orthology_based_reconstruction/mhun/mhun
.faa

cp outside/species.sbml /shared/Formicicum/orthology_based_reconstruction/mhun/mhun
.sbml

(grep '>' /shared/Formicicum/orthology_based_reconstruction/mhun/mhun.faa | awk '{corda=substr($0, 2, length($0)); print corda"\t"corda}' > /shared/Formicicum/orthology_based_reconstruction/mhun/dict_genes.txt)

cp /shared/Formicicum/orthology_based_reconstruction/mhun/orthology (...).sbml outside/aureme_model.sbml
```

## Analysis

Some explanation on methodology of analysis
