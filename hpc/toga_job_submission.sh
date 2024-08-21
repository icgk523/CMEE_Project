#!/bin/bash
#PBS -l select=1:ncpus=256:mem=4000gb
#PBS -l walltime=72:00:00

module load anaconda3/personal

cp -r $HOME/Code/TOGA/* $TMPDIR

mkdir $TMPDIR/Genomes/
mkdir $TMPDIR/lastz_output/

cp -r $HOME/Data/Lepidoptera_Genomes/Genome_1/* $TMPDIR/Genomes/
cp -r $HOME/Results/lastz_output_files_helicoverpa/* $TMPDIR/lastz_output/

cp $HOME/Data/GCF_030705265.1_ASM3070526v1_genomic_filtered.bed $TMPDIR
cp $HOME/Data/GCF_030705265.1_ASM3070526v1_genomic_isoforms.tsv $TMPDIR
cp $HOME/Data/GCF_030705265.1_ASM3070526v1_genomic.2bit $TMPDIR

for twobitfile in "$TMPDIR/Genomes"/*.2bit; do # Iterate through all 2bit files
    genome_name=$(basename "$twobitfile" .2bit) # Take the genome accession name to find the lastz_output file
    chain_file="$TMPDIR/lastz_output/output_$genome_name/helicoverpa.$genome_name.final.chain.gz"
    output_directory="$TMPDIR/TOGA_output_$genome_name"
    mkdir "$output_directory"

    $TMPDIR/./toga.py "$chain_file" $TMPDIR/GCF_030705265.1_ASM3070526v1_genomic_filtered.bed $TMPDIR/GCF_030705265.1_ASM3070526v1_genomic.2bit $TMPDIR/Genomes/"$genome_name".2bit --kt --pn "$output_directory" -i $TMPDIR/GCF_030705265.1_ASM3070526v1_genomic_isoforms.tsv --cb 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50 --cjn 500 --ms
    mv "$output_directory" $HOME/Results/TOGA_Output_Helicoverpa/
    echo "Files for the genome have been moved to $HOME/TOGA_Output/"$output_directory""
done
