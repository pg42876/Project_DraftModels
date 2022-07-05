# Development of metabolic models

### Installation of tools

To develop the models for Methanobacterium formicicum the tools merlin, CarveMe, KBase and AuReMe were used. 

## MERLIN
merlin is available for download in a single multi-platform version in https://merlin-sysbio.org/download/.

## CARVEME
conda install -c conda-forge mamba python=3.6 -n carveme
conda activate carveme
mamba install -c bioconda diamond
mamba install -c bioconda carveme

## AUREME
docker pull dyliss/aureme-img:2.4

### Run tools

## CARVEME
carve genome.faa --output model.xml

## AUREME
docker run -it -v C:\Users\Desktop:/outside dyliss/aureme-img:2.4 bash
docker start 40a155090560au
docker exec -it 40a155090560 bash

aureme --init=Formicicum

cp outside/species.faa /shared/Formicicum/genomic_data/

mkdir /shared/Formicicum/orthology_based_reconstruction/mhun

awk '{print $1}' /shared/Formicicum/genomic_data/prot_sequences.fasta > /shared/Formicicum/genomic_data/Formicicum.faa

cp outside/species.faa /shared/Formicicum/orthology_based_reconstruction/mhun/mhun.faa

cp outside/species.sbml /shared/Formicicum/orthology_based_reconstruction/mhun/mhun.sbml

(grep '>' /shared/Formicicum/orthology_based_reconstruction/mhun/mhun.faa | awk '{corda=substr($0, 2, length($0)); print corda"\t"corda}' > /shared/Formicicum/orthology_based_reconstruction/mhun/dict_genes.txt)

cp /shared/Formicicum/orthology_based_reconstruction/mhun/orthology (...).sbml outside/aureme_model.sbml

### Analysis

The models will be compared based on the number of reactions, genes and metabolites. Taking these results into account, the tools were evaluated for their efficiency and accuracy.
Regarding the reactions/enzymes chosen, the reactions associated with methanogenesis that were marked in the KEGG Pathways were taken into account and those that seemed most relevant were chosen. 
For the metrics, the precision, recall, F1 score and Jaccard Distance will be calculated. 

### Scripts
The sequences used to generate the templates are in the "Sequences" folder.
The draft models are in the "Models" folder.
The scripts for draft models assessment are stored in the "Scripts" folder.
The results of the assessment are stored in the "Results" folder.
