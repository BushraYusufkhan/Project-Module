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
# Diploscapter oocerea
Command to search and download reference mitogenome
```
 docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master findMitoReference.py --species "Diploscapter" --outfolder /home/bkhan/Data/Diploscapter_oocerea/incomplete --min_length 1
2000
```
Running mitohifi on raw sequence data
```
 docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master mitohifi.py -r /home/bkhan/Data/Diploscapter_oocerea/incomplete/m54274Ue_220811_201236.hifi_reads.fastq.gz -f /home/bkhan/Data/Diploscapter_oocerea/incomplete/NC_035106.1.fasta -g /home/bkhan/Data/Diploscapter_oocerea/incomplete/NC_035106.1.gb -t 12 -o 5
```
There was some key error when I used this command, it had something to do with the reference file, and the pipeline couldn't finish. However, the final annotation was done. 36 genes were predicted.
So, I used a different reference and it worked. But this reference was not that closely related to our species so, this time it predicted 35 genes. So, the previous command is still important to have a look at the final annotation.
```
 docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master mitohifi.py -r /home/bkhan/Data/Diploscapter_oocerea/m54274Ue_220811_201236.hifi_reads.fastq.gz -f /home/bkhan/Data/Diploscapter_ooc
erea/NC_009885.1.fasta -g /home/bkhan/Data/Diploscapter_oocerea/NC_009885.1.gb -t 12 -o 5
```
# Tripylella_sp
Command to search and download the reference mitogenome
```
 docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master findMitoReference.py --species "Tripylella" --outfolder /home/bkhan/Data/Tripylella_sp --min_length 12000
```
Running mitohifi on the contigs data
```
 docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master mitohifi.py -c /home/bkhan/Data/Tripylella_sp/assembly.fasta -f /home/bkhan/Data/Tripylella_sp/NC_056391.1.fasta -g /home/bkhan/Data/Tripylella_sp/NC_056391.1.gb -t 12 -o 5 --mitos
```




