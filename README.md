# Project-Module
Project Module 2024
# Commands used during the project.
The project was to use mitohifi to assemble and annotate mitochondrial genomes of some nematode species.
## Circular mitochondrial genome:
### Plectus Sambesii
Command to find and download the fasta and genbank file of the closely related species that will be used as a reference.
```
 docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master findMitoReference.py --species "Plectus sambesii" --outfolder /home/bkhan/Data/Plectus_sambesii --min_length 12000
```
Running mitohifi on the contigs data
```
 docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master mitohifi.py -c /home/bkhan/Data/Plectus_sambesii/psambesii.hifiasm_v019_l3.bp.p_ctg.fasta -f /home/bkhan/Data/Plectus_sambesii/KX017524.1.fasta -g /home/bkhan/Data/Plectus_sambesii/KX017524.1.gb -t 12 -o 5

```
The default annotater in mitohifi pipeline Mitofinder could not predict nad2 and nad6 which are predicted by mitos.
```
 docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master mitohifi.py -c /home/bkhan/Data/Plectus_sambesii/psambesii.hifiasm_v019_l3.bp.p_ctg.fasta -f /home/bkhan/Data/Plectus_sambesii/KX017524.1.fasta -g /home/bkhan/Data/Plectus_sambesii/KX017524.1.gb -t 12 -o 5 --mitos

```
### Plectus PKLO1
Command to search and download reference mitogenome.
```
  docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master findMitoReference.py --species "Plectus" --outfolder /home/bkhan/Data/Plectus_PKL01 --min_length 12000
```
Running Mitohifi on the raw reads.
```
  docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master mitohifi.py -r /home/bkhan/Data/Plectus_PKL01/m64093_230430_124203.hifi_reads.fastq.gz -f /home/bkhan/Data/Plectus_PKL01/KX017524.1.fasta -g /home/bkhan/Data/Plectus_PKL01/KX017524.1.gb -t 12 -o 5
```
### Diploscapter oocerea
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
### Tripylella_sp
Command to search and download reference mitogenome
```
 docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master findMitoReference.py --species "Tripylella" --outfolder /home/bkhan/Data/Tripylella_sp --min_length 12000
```
Running mitohifi on the contigs data
```
 docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master mitohifi.py -c /home/bkhan/Data/Tripylella_sp/assembly.fasta -f /home/bkhan/Data/Tripylella_sp/NC_056391.1.fasta -g /home/bkhan/Data/Tripylella_sp/NC_056391.1.gb -t 12 -o 5 --mitos
```

## linear mitochondrial genomes
### Acrobeloides nanus
Search and Download reference Mitogenome
```
 docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master findMitoReference.py --species "Acrobeloides nanus" --outfolder /home/bkhan/Data/Acrobeloides_nanus --min_length 12000
```
Running mitohifi on the contigs data
```
docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master mitohifi.py -c /home/bkhan/Data/Acrobeloides_nanus/ananus.flye_v29_default.fasta -f /home/bkhan/Data/Acrobeloides_nanus/MK559448.1.fasta -g /home/bkhan/Data/Acrobeloides_nanus/MK559448.1.gb -t 12 -o 5
```
### Acrobeloides_PAP2217
Search and download the reference mitogenome
```
 docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master findMitoReference.py --species "Acrobeloides" --outfolder /home/bkhan/Data/Acrobeloides_PAP2217 --min_length 12000
```
Running mitohifi on the reads
```
docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master mitohifi.py -r /home/bkhan/Data/Acrobeloides_PAP2217/m54274Ue_230602_080202.hifi_reads.fastq.gz -f /home/bkhan/Data/Acrobeloides_PAP2217/MK559448.1.fasta -g /home/bkhan/Data/Acrobeloides_PAP2217/MK559448.1.gb -t 12 -o 5
```
### Acrobeloides_buchneri
Search and download the reference mitogenome
```
docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master findMitoReference.py --species "Acrobeloides buchneri" --outfolder /home/bkhan/Data/Acrobeloides_buchneri --min_length 12000
```
Running mitohifi on the reads
```
docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master mitohifi.py -r /home/bkhan/Data/Acrobeloides_buchneri/m54274Ue_230501_185711.hifi_reads.fastq.gz -f /home/bkhan/Data/Acrobeloides_buchneri/MK559448.1.fasta -g /home/bkhan/Data/Acrobeloides_buchneri/MK559448.1.gb -t 12 -o 5
```
###  Acrobeloides_ARO2205
Command to search and download the reference mitogenome.
```
 docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master findMitoReference.py --species "Acrobeloides" --outfolder /home/bkhan/Data/Acrobeloides_ARO2205 --min_length 12000
```
Running mitohifi on the reads
```
 docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master mitohifi.py -r /home/bkhan/Data/Acrobeloides_ARO2205/m54274Ue_230531_224848.hifi_reads.fastq.gz -f /home/bkhan/Data/Acrobeloides_ARO2205/MK559448.1.fasta -g  /home/bkhan/Data/Acrobeloides_ARO2205/MK559448.1.gb -t 30 -o 5
```

## Could not find mitochondrial genomes
###  Anaplectus_granulosus
Command to find and download the reference mitogenome
```
docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master findMitoReference.py --species "Anaplectus granulosus" --outfolder /home/bkhan/Data/Anaplectus_granulosus --min_length 12000
```
Running mitohifi on the reads
```
 docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master mitohifi.py -r /home/bkhan/Data/Anaplectus_granulosus/m54274Ue_230530_134758.hifi_reads.fa -f  /home/bkhan/Data/Anaplectus_granulosu
s/KX017524.1.fasta -g /home/bkhan/Data/Anaplectus_granulosus/KX017524.1.gb -t 12 -o 5
```
Running mitohifi on the contigs data
```
docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master mitohifi.py -c /home/bkhan/Data/Anaplectus_granulosus/agranulosus.flye_v29_default.fasta -f  /home/bkhan/Data/Anaplectus_granulosus/
KX017524.1.fasta -g /home/bkhan/Data/Anaplectus_granulosus/KX017524.1.gb -t 12 -o 5
```
### Mononchus_laminatus
Running mitohifi on the reads
```
 docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master mitohifi.py -r /home/bkhan/Data/Mononchus_laminatus/m64093_230811_043738.hifi_reads.fastq.gz -f /home/bkhan/Data/Mononchus_laminatus/NC_056391.1.fasta -g /home/bkhan/Data/Mononchus_laminatus/NC_056391.1.gb -t 12 -o 5
```
Running mitohifi on the contigs data
```
docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master mitohifi.py -c /home/bkhan/Data/Mononchus_laminatus/mlaminatus.flye_v29_default.fasta -f /home/bkhan/Data/Mononchus_laminatus/NC_05639
1.1.fasta -g /home/bkhan/Data/Mononchus_laminatus/NC_056391.1.gb -t 12 -o 5
```










