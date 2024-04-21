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
2000 -n 6
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
 docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master findMitoReference.py --species "Tripylella" --outfolder /home/bkhan/Data/Tripylella_sp --min_length 12000 -n 6
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
Command to find and download the reference mitogenome
```
docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master findMitoReference.py --species "Mononchus" --outfolder /home/bkhan/Data/Mononchus_laminatus --min_length 12000
```
Running mitohifi on the reads
```
 docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master mitohifi.py -r /home/bkhan/Data/Mononchus_laminatus/m64093_230811_043738.hifi_reads.fastq.gz -f /home/bkhan/Data/Mononchus_laminatus/NC_056391.1.fasta -g /home/bkhan/Data/Mononchus_laminatus/NC_056391.1.gb -t 12 -o 5
```
Running mitohifi on the contigs data
```
docker run -u $(id -u) -v $(pwd):$PWD ghcr.io/marcelauliano/mitohifi:master mitohifi.py -c /home/bkhan/Data/Mononchus_laminatus/mlaminatus.flye_v29_default.fasta -f /home/bkhan/Data/Mononchus_laminatus/NC_05639
1.1.fasta -g /home/bkhan/Data/Mononchus_laminatus/NC_056391.1.gb -t 12 -o 5
```
## Commands used for the phylogenetic tree:
### Extracting genes from Genbank files
Below is a Python script that you can use to extract genes from Genbank files, make sure to make changes to the script according to your needs. And Biopython should also be installed in your system.
```
import os
from Bio import SeqIO

def extract_features(genbank_file, feature_type):
    extracted_sequences = []
    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type == feature_type:
                if 'gene' in feature.qualifiers:
                    feature_name = feature.qualifiers['gene'][0]
                elif 'product' in feature.qualifiers:
                    feature_name = feature.qualifiers['product'][0]
                else:
                    feature_name = 'Unknown'
                extracted_sequences.append((feature.location.extract(record).seq, feature_name, feature_type))
    return extracted_sequences

def write_sequences_to_fasta(sequences, output_file):
    with open(output_file, "w") as f:
        for idx, (seq, feature_name, feature_type) in enumerate(sequences):
            f.write(f">{feature_type} - {feature_name}\n{seq}\n")

def process_genbank_file(input_file):
    output_file = os.path.splitext(input_file)[0] + ".genes.fasta"
    # Extract genes
    genes = extract_features(input_file, "gene")
    print(f"Number of genes extracted: {len(genes)}")

    # Extract tRNAs
    tRNAs = extract_features(input_file, "tRNA")
    print(f"Number of tRNAs extracted: {len(tRNAs)}")

    # Extract rRNAs
    rRNAs = extract_features(input_file, "rRNA")
    print(f"Number of rRNAs extracted: {len(rRNAs)}")

    # Write all sequences to the same output file
    write_sequences_to_fasta(genes + tRNAs + rRNAs, output_file)

    print(f"All sequences saved to {output_file}")

# Folder containing GenBank files
folder_path = "/home/bushra/geneextract"

# Iterate through each file in the folder
for file_name in os.listdir(folder_path):
    if file_name.endswith(".gb"):
        genbank_file = os.path.join(folder_path, file_name)
        print(f"Processing file: {genbank_file}")
        process_genbank_file(genbank_file)
```
In some cases you will see that the script has extracted more genes than usual, this is because of the logic of the code, the GenBank files sometimes have slightly different formats, and they have used different identifiers for genes, tRNAs, and rrns. For example, in some cases, the identifier "gene" is used for all of them and sometimes the genes are identified by "gene" and the tRNAs and rRNAs are determined by the identifier "product". So, you can use the script below to delete the duplicates, make sure to specify the input and output folder.
```
import os
from Bio import SeqIO

def keep_unique_sequences(input_folder, output_folder):
    for file_name in os.listdir(input_folder):
        if file_name.endswith(".edited.fasta"):
            input_fasta_file = os.path.join(input_folder, file_name)
            output_fasta_file = os.path.join(output_folder, file_name.replace(".edited.fasta", ".unique_sequences.fasta"
))
            sequences = set()

            with open(output_fasta_file, 'w') as output_fasta:
                for record in SeqIO.parse(input_fasta_file, 'fasta'):
                    sequence = str(record.seq)
                    if sequence not in sequences:
                        sequences.add(sequence)
                        SeqIO.write(record, output_fasta, 'fasta')

# Example usage:
input_folder = "/home/bushra/geneextract/editedfiles"
output_folder = "/home/bushra/geneextract/uniquesequences"
keep_unique_sequences(input_folder, output_folder)
```

### Extracting genes from gff files
if you have gff files, you can use this python script to extract the genes, you also need the fasta files for this, make a folder called fastafiles to save your fasta files, it will read the coordinates from the gff files and extract the sequences from the fasta files.
```
from Bio import SeqIO
import os

def extract_sequences_from_gff_folder(input_folder):
    fasta_folder = os.path.join(input_folder, "fastafiles")
    for file_name in os.listdir(input_folder):
        if file_name.endswith(".gff"):
            input_gff_file = os.path.join(input_folder, file_name)
            fasta_file = os.path.join(fasta_folder, file_name.replace(".gff", ".fasta"))
            print(f"GFF file: {input_gff_file}")
            print(f"FASTA file: {fasta_file}")
            if os.path.exists(fasta_file):
                extract_sequences_from_gff(input_gff_file, fasta_file)
            else:
                print(f"FASTA file {fasta_file} not found.")

def extract_sequences_from_gff(input_gff_file, fasta_file):
    # Dictionary to store sequences
    sequences = {}

    # Construct output file name based on input GFF file name
    output_file = os.path.splitext(input_gff_file)[0] + ".genes.fasta"

    # Open output file in write mode
    with open(output_file, 'w') as f_out:
        # Open GFF file
        with open(input_gff_file, 'r') as gff_file:
            for line in gff_file:
                # Skip comment lines
                if line.startswith('#'):
                    continue
                fields = line.split('\t')
                if len(fields) < 9:
                    continue
                seq_id = fields[0]
                feature_type = fields[2]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                attributes = fields[8]

                # Extracting gene ID from attributes
                gene_id = ""
                for attribute in attributes.split(';'):
                    if attribute.startswith('ID='):
                        gene_id = attribute.split('=')[1]
                        break

                # Reading sequence from FASTA file if not already read
                if seq_id not in sequences:
                    seq_record = SeqIO.read(fasta_file, "fasta")
                    sequences[seq_id] = seq_record.seq

                # Extracting sequence based on feature type and coordinates
                if feature_type in ['gene', 'tRNA', 'rRNA']:
                    feature_seq = sequences[seq_id][start - 1:end]
                    # Reverse complement if feature is on the negative strand
                    if strand == '-':
                        feature_seq = feature_seq.reverse_complement()
                    # Write the extracted sequence to the output file with modified header
                    gene_name = gene_id.split(':')[0]  # Extracting gene name
                    f_out.write(f">{seq_id}_{gene_name}_{feature_type}\n")
                    f_out.write(f"{feature_seq}\n")

# Example usage:
input_folder = "/home/bushra/geneextract/more"
extract_sequences_from_gff_folder(input_folder)
```
I had two types of file formats Genbank and gff, after extracting the genes at this point, I highly recommend making changes to the heading of the fasta files. Add the organism name, keep the gene name, and ensure the abbreviations used for the genes in all the files are the same. Genbank files and gff files have used different abbreviations for some genes, which can cause problems.  
Now make a folder with all of the extracted gene files, I hope you have changed the headings. Now is the time to sort the genes into separate files. 
Before that, I looked for the most common genes in all of the species and provided the list of those genes in the script below. I do face some issues with this script, it is not perfect, for nad4 it was adding the gene sequences of nad4l too, so I had to delete nad4l sequences from the nad4 file. Also for tRNA-alanine, it caused some problems, it was extracting all of the gene sequences into the tRNA-ala file, which is why I had to make a few changes to this script for tRNA-ala. So, after running this script make sure that all of the the genes are sorted correctly into their respective files. Provide your input and output folder.

```
import os
from Bio import SeqIO

# Function to extract gene sequences from fasta files
def extract_gene_sequences(folder_path, gene_aliases):
    gene_sequences_dict = {}

    # Iterate through each file in the folder
    for filename in os.listdir(folder_path):
        if filename.endswith(".fasta") or filename.endswith(".fa"):
            file_path = os.path.join(folder_path, filename)

            # Debugging: Print the current file being checked
            print("Checking file:", file_path)

            # Open fasta file and extract gene sequences
            with open(file_path, "r") as fasta_file:
                for record in SeqIO.parse(fasta_file, "fasta"):
                    # Debugging: Print the current record being checked
                    print("Checking record:", record.description)

                    for gene_name, aliases in gene_aliases.items():
                        # Debugging: Print the current gene name being checked
                        print("Gene Name:", gene_name)

                        for alias_set in aliases:
                            for alias in alias_set:
                                if alias.lower() in record.description.lower():
                                    # Debugging: Print the alias found in the record description
                                    print("Found alias:", alias)

                                    if gene_name in gene_sequences_dict:
                                        # Append the sequence to existing list
                                        gene_sequences_dict[gene_name].append((record.description, str(record.seq)))
                                    else:
                                        # Initialize the list with the sequence
                                        gene_sequences_dict[gene_name] = [(record.description, str(record.seq))]

                                    # Debugging: Print when trnA sequence is found
                                    if gene_name == "trnA":
                                        print("Found trnA sequence:", record.description)

                                    break  # Stop searching once the gene is found
                            else:
                                continue
                            break

    return gene_sequences_dict

# Function to save gene sequences to separate files
def save_gene_sequences(gene_sequences_dict, output_folder):
    # Create the output directory if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    for gene_name, sequences in gene_sequences_dict.items():
        # Construct filename
        output_file = os.path.join(output_folder, f"{gene_name}.fasta")
        with open(output_file, "w") as output:
            for description, sequence in sequences:
                output.write(f">{description}\n{sequence}\n")

# Path to the folder containing fasta files
folder_path = "/home/bushra/geneextract/uniquesequences"

# Gene names and their aliases to extract
gene_aliases = {
    "cox2": [["cox2", "COX2"]],
    "nad4l": [["nad4l", "ND4L"]],
    "nad6": [["nad6", "ND6"]],
    "nad2": [["nad2", "ND2"]],
    "cox3": [["cox3", "COX3"]],
    "nad4": [["nad4", "ND4"]],
    "cob": [["cytb", "CYTB", "cob"]],
    "atp6": [["atp6", "ATP6"]],
    "nad5": [["nad5", "ND5"]],
    "cox1": [["cox1", "COX1"]],
    "nad3": [["nad3", "ND3"]],
    "nad1": [["nad1", "ND1"]],
    "rrnS": [["s-rRNA", "rrnS"]],
    "rrnL": [["l-rRNA", "rrnL"]],
    "trnP": [["tRNA-Pro", "trnP"]],
    "trnY": [["tRNA-Tyr", "trnY"]],
    "trnQ": [["tRNA-Gln", "trnQ"]],
    "trnR": [["tRNA-Arg", "trnR"]],
    "trnI": [["tRNA-Ile", "trnI"]],
    "trnE": [["tRNA-Glu", "trnE"]],
    "trnN": [["tRNA-Asn", "trnN"]],
    "trnT": [["tRNA-Thr", "trnT"]],
    "trnA": [["tRNA-trnA"], ["tRNA-Ala"]],  # Updated alias for trnA
    "trnC": [["tRNA-Cys", "trnC"]],
    "trnV": [["tRNA-Val", "trnV"]],
    "trnM": [["tRNA-Met", "trnM"]],
    "trnD": [["tRNA-Asp", "trnD"]],
    "trnG": [["tRNA-Gly", "trnG"]],
    "trnF": [["tRNA-Phe", "trnF"]],
    "trnK": [["tRNA-Lys", "trnK"]],
    "trnH": [["tRNA-His", "trnH"]],
    "trnW": [["tRNA-Trp", "trnW"]],
    "trnL": [["tRNA-Leu", "trnL"]],
    "trnS": [["tRNA-Ser", "trnS"]]
}

# Output folder to save extracted gene sequences
output_folder = "/home/bushra/geneextract/sorted_genes/firstseparate"

# Extract gene sequences
gene_sequences_dict = extract_gene_sequences(folder_path, gene_aliases)

# Save gene sequences to separate files
save_gene_sequences(gene_sequences_dict, output_folder)

print("Gene sequences extracted and saved to", output_folder)
```











