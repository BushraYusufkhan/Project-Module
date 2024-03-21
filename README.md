# Project-Module
Project Module 2024
# Commands used during the project are listed below.
The project was to use mitohifi to assemble mitogenome of some nematode species.
# Plectus Sambesii
Command to find and download the fasta and genbank file of the closely related species that will be used as reference.
```
 docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master findMitoReference.py --species "Plectus sambesii" --outfolder /home/bkhan/Data/Plectus_sambesii --min_length 12000
```
