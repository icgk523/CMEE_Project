#!/bin/bash
#PBS -l select=1:ncpus=256:mem=100gb
#PBS -l walltime=72:00:00

module load anaconda3/personal # Load Anaconda

echo "Anaconda has been loaded"

cp -r $HOME/Code/make_lastz_chains/ $TMPDIR # Move make_lastz_chains to TMP
cp -r $HOME/Data/Hymenoptera_Genomes/Genome_1/ $TMPDIR/Genome # Move your genomes to TMP
cp -r $HOME/Data/GCF_029169275.1_AcerK_1.0_genomic.2bit $TMPDIR/ # Move reference genome to TMP

echo "Files have been moved - Pipeline is about to run"

for twobitfile in "$TMPDIR/Genome"/*.2bit; do # Iterate through all 2bit files
    genome_name=$(basename "$twobitfile" .2bit) # Take the genome accession name
    output_directory="$TMPDIR/output_$genome_name" # Create an output directory with genome accession name
    mkdir "$output_directory"
    $TMPDIR/make_lastz_chains/./make_chains.py apis "$genome_name" $TMPDIR/GCF_029169275.1_AcerK_1.0_genomic.2bit $TMPDIR/Genome/"$genome_name".2bit --pd "$output_directory" -f --chaining_memory 32 --cluster_executor local
    echo "make_lastz_chains has successfully run for genome:" "$genome_name"
    mv "$output_directory" $HOME/Results/lastz_output_files_apis/
done

echo "The pipeline has finished"
