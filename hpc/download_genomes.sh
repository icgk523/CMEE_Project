#!/bin/bash
module load anaconda3/personal

mapfile -t genome_accessions < ../Data/genome_accessions.txt # Open a txt file with a list of genome accession numbers

for acc in "${genome_accessions[@]}"; do
    acc="${acc//\"/}"
    path=$($HOME/edirect/./efetch -db assembly -id $acc -format docsum | $HOME/edirect/./xtract -pattern DocumentSummary -element FtpPath_GenBank | tr -d '"')
    path="https${path#ftp}"
    last_part=$(basename "$path")
    new_path="${path}/${last_part}"
    filename="${new_path}_genomic.fna.gz"
    wget "$filename" -P ../Data/Hymenoptera_Genomes/
    gzip -d ../Data/Hymenoptera_Genomes/${last_part}_genomic.fna.gz
    sed -E -i 's/(^>.*?)[.]1/\1/' ../Data/Hymenoptera_Genomes/${last_part}_genomic.fna
    faToTwoBit ../Data/Hymenoptera_Genomes/${last_part}_genomic.fna ../Data/Hymenoptera_Genomes/${last_part}_genomic.2bit
    rm ../Data/Hymenoptera_Genomes/${last_part}_genomic.fna
done

