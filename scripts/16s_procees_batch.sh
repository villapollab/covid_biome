#!/bin/bash

#SBATCH --time=04-00:00:00
#SBATCH --partition=defq
#SBATCH --mail-user=myemail@email.org
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --ntasks-per-node=35
#SBATCH --nodes=1
#SBATCH --job-name=16s_r2
#SBATCH --comment=16s_r2

module load mamba/23.1.0-1
module load dorado/0.5.0

echo "Complete analysis of ONT Nanopore 16S sequencing data"

echo "Created by Austin Marshall, Dept. of Biology @ Clarkson University"

echo "There is one manual part of this script, it is the naming of the files in the next two lines"

# Write what barcodes you used here in the format "SQK-16S024_barcodeXX.fastq.gz"
original_names=("SQK-16S024_barcode01.fastq.gz" "SQK-16S024_barcode02.fastq.gz" "SQK-16S024_barcode03.fastq.gz" "SQK-16S024_barcode04.fastq.gz" "SQK-16S024_barcode05.fastq.gz" "SQK-16S024_barcode06.fastq.gz" "SQK-16S024_barcode07.fastq.gz" "SQK-16S024_barcode08.fastq.gz" "SQK-16S024_barcode09.fastq.gz" "SQK-16S024_barcode10.fastq.gz" "SQK-16S024_barcode11.fastq.gz" "SQK-16S024_barcode12.fastq.gz" "SQK-16S024_barcode13.fastq.gz" "SQK-16S024_barcode14.fastq.gz" "SQK-16S024_barcode15.fastq.gz" "SQK-16S024_barcode16.fastq.gz" "SQK-16S024_barcode17.fastq.gz" "SQK-16S024_barcode18.fastq.gz" "SQK-16S024_barcode19.fastq.gz" "SQK-16S024_barcode20.fastq.gz" "SQK-16S024_barcode21.fastq.gz" "SQK-16S024_barcode22.fastq.gz" "SQK-16S024_barcode23.fastq.gz" "SQK-16S024_barcode24.fastq.gz")

# Write the new names substituting "SQK-16S024_barcodeXX.fastq.gz" for "samplename.fastq.gz"
new_names=("r2_b1.fastq.gz" "r2_b2.fastq.gz" "r2_b3.fastq.gz" "r2_b4.fastq.gz" "r2_b5.fastq.gz" "r2_b6.fastq.gz" "r2_b7.fastq.gz" "r2_b8.fastq.gz" "r2_b9.fastq.gz" "r2_b10.fastq.gz" "r2_b11.fastq.gz" "r2_b12.fastq.gz" "r2_b13.fastq.gz" "r2_b14.fastq.gz" "r2_b15.fastq.gz" "r2_b16.fastq.gz" "r2_b17.fastq.gz" "r2_b18.fastq.gz" "r2_b19.fastq.gz" "r2_b20.fastq.gz" "r2_b21.fastq.gz" "r2_b22.fastq.gz" "r2_b23.fastq.gz" "r2_b24.fastq.gz")

# Check if the new_names array is empty
if [ ${#new_names[@]} -eq 0 ]; then
    echo "No new names have been provided. Exiting..."
    exit 1
fi

# Make sure to remove this line if reusing script
rm -r /home/tmhagm8/scratch/round2
#######################################################

mkdir /home/tmhagm8/scratch/round2
mkdir /home/tmhagm8/scratch/round2/fastqs
mkdir /home/tmhagm8/scratch/round2/emu_output

# Set paths
DIR="/home/tmhagm8/scratch/round2"
INPUT="/home/tmhagm8/scratch/"
DEMUX_OUTDIR="${DIR}/fastqs"
EMU_ABUND_DIR="${DIR}/emu_output"
EMU_DB_DIR="/home/tmhagm8/scratch/emu/emu_database"

mamba activate dorado

# Demultiplex and convert the output.bam to fastqs
dorado demux --kit-name SQK-16S024 --emit-fastq --output-dir "${DEMUX_OUTDIR}" "${INPUT}"/round2_covid.bam

mamba deactivate

mamba activate 16s

# Once a large fastq file for each sample is made we will select only the reads 1350 >= x =< 1650 bps
for f in "${DEMUX_OUTDIR}"/*.fastq; do filename="${f%*.fastq}";
    echo $f "Start Time: `date`";
    cat $f | chopper --minlength 1350 --maxlength 1650 --threads 70 | pigz > $filename".fastq.gz"
    echo $(date);
done

# Loop through the files and rename the actual samples
for ((i=0; i<${#original_names[@]}; i++)); do
    original_name="${original_names[$i]}"
    new_name="${new_names[$i]}"

    # Check if the original file exists in the source directory
    if [ -e "${DEMUX_OUTDIR}/${original_name}" ]; then
        # Rename and move the file to the target directory
        mv "${DEMUX_OUTDIR}/${original_name}" "${DIR}/${new_name}"
        echo "Renamed: ${original_name} to ${new_name}"
    else
        echo "Original file not found: ${original_name}. Skipping..."
    fi
done

# Gunzip the files to run emu on
for f in "${DIR}"/*.fastq.gz; do filename="${f%*.fastq.gz}";
    echo $f "Start Time: `date`";
    pigz -d $f ;
    echo $(date);
done

#Compute emu relative abundances for all fastqs with the EMU database
for f in "${DIR}"/r*.fastq; do filename="${f%*.fastq}";
    echo $f "Start Time: `date`";
    emu abundance $f --db "${EMU_DB_DIR}" --threads 70 --keep-counts --keep-read-assignments --output-dir "${EMU_ABUND_DIR}" ;
    echo $(date);
done

#Combine all samples into one otu type table for species for the EMU output
for f in *rel-abundance.tsv; do filename="${f%*rel-abundance.tsv}";
    echo $f "Start Time: `date`"; emu combine-outputs --split-tables "${EMU_ABUND_DIR}" species ;
    echo $(date);
done

#Combine all samples into one otu type table for species counts for the EMU output
for f in *rel-abundance.tsv; do filename="${f%*rel-abundance.tsv}";
    echo $f "Start Time: `date`"; emu combine-outputs --counts --split-tables "${EMU_ABUND_DIR}" species ;
    echo $(date);
done

#Combine all samples into one otu type table for genus for the EMU output
for f in *rel-abundance.tsv; do filename="${f%*rel-abundance.tsv}";
    echo $f "Start Time: `date`"; emu combine-outputs --split-tables "${EMU_ABUND_DIR}" genus ;
    echo $(date);
done

#Combine all samples into one otu type table for phylum for the EMU output
for f in *rel-abundance.tsv; do filename="${f%*rel-abundance.tsv}";
    echo $f "Start Time: `date`"; emu combine-outputs --split-tables "${EMU_ABUND_DIR}" phylum ;
    echo $(date);
done

echo "congrats, now move to R for visualization"
