# Project-Module
Project Module 2024
# Commands used during the project.
The project was to use mitohifi to assemble mitogenome of some nematode species.
# Plectus Sambesii
Command to find and download the fasta and genbank file of the closely related species that will be used as a reference.
```
 docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master findMitoReference.py --species "Plectus sambesii" --outfolder /home/bkhan/Data/Plectus_sambesii --min_length 12000
```
Running mitohifi on the contigs data to assemble and annotate the mitogenome 
```
 docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master mitohifi.py -c /home/bkhan/Data/Plectus_sambesii/psambesii.hifiasm_v019_l3.bp.p_ctg.fasta -f /home/bkhan/Data/Plectus_sambesii/KX017524.1.fasta -g /home/bkhan/Data/Plectus_sambesii/KX017524.1.gb -t 12 -o 5

```
The default annotater in mitohifi pipeline Mitofinder could not predict nad2 and nad6 which are predicted by mitos.
```
 docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master mitohifi.py -c /home/bkhan/Data/Plectus_sambesii/psambesii.hifiasm_v019_l3.bp.p_ctg.fasta -f /home/bkhan/Data/Plectus_sambesii/KX017524.1.fasta -g /home/bkhan/Data/Plectus_sambesii/KX017524.1.gb -t 12 -o 5 --mitos

```
# Plectus PKLO1
Command to search and download reference mitogenome.
```
  docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master findMitoReference.py --species "Plectus" --outfolder /home/bkhan/Data/Plectus_PKL01 --min_length 12000
```
Running Mitohifi on the raw reads.
```
  docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master mitohifi.py -r /home/bkhan/Data/Plectus_PKL01/m64093_230430_124203.hifi_reads.fastq.gz -f /home/bkhan/Data/Plectus_PKL01/KX017524.1.fasta -g /home/bkhan/Data/Plectus_PKL01/KX017524.1.gb -t 12 -o 5
```



