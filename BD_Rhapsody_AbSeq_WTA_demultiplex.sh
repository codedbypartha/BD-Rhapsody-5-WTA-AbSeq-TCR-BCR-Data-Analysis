#!/bin/bash

# Get the directory where the script is located
script_dir=$(dirname "$(realpath "$0")")

# Set the input and output directories to the script's directory
input_dir="$script_dir"
output_dir="$script_dir/output"
mkdir -p "$output_dir"

# Define the sample tags
declare -A sample_tags
sample_tags["1"]="ATTCAAGGGCAGCCGCGTCACGATTGGATACGACTGTTGGACCGG"
sample_tags["2"]="TGGATGGGATAAGTGCGTGATGGACCGAAGGGACCTCGTGGCCGG"
sample_tags["3"]="CGGCTCGTGCTGCGTCGTCTCAAGTCCAGAAACTCCGTGTATCCT"
sample_tags["4"]="ATTGGGAGGCTTTCGTACCGCTGCCGCCACCAGGTGATACCCGCT"
sample_tags["5"]="CTCCCTGGTGTTCAATACCCGATGTGGTGGGCAGAATGTGGCTGG"
sample_tags["6"]="TTACCCGCAGGAAGACGTATACCCCTCGTGCCAGGCGACCAATGC"
sample_tags["7"]="TGTCTACGTCGGACCGCAAGAAGTGAGTCAGAGGCTGCACGCTGT"
sample_tags["8"]="CCCCACCAGGTTGCTTTGTCGGACGAGCCCGCACAGCGCTAGGAT"
sample_tags["9"]="GTGATCCGCGCAGGCACACATACCGACTCAGATGGGTTGTCCAGG"
sample_tags["10"]="GCAGCCGGCGTCGTACGAGGCACAGCGGAGACTAGATGAGGCCCC"
sample_tags["11"]="CGCGTCCAATTTCCGAAGCCCCGCCCTAGGAGTTCCCCTGCGTGC"
sample_tags["12"]="GCCCATTCATTGCACCCGCCAGTGATCGACCCTAGTGGAGCTAAG"

# Get the expected number of cells as the first argument and the sample tag numbers as the remaining arguments
if [ "$#" -lt 2 ]; then
    echo "Please provide the expected number of cells and sample tag numbers as arguments."
    exit 1
fi

# remove the diversity inserts from both the Abseq and WTA R1 fastq files
python process_fastq.py

expected_cells=$1
shift 1
sample_tag_numbers=("$@")

# Selectively unzip Abseq files (both R1 and R2)
for file in "$input_dir"/*Abseq*R[12]*.fastq.gz; do
    echo "Unzipping $file..."
    gunzip -c "$file" > "${file%.gz}"
done

# Extract sample tag reads for the specified sample tags
for file in "$input_dir"/*Abseq*R2*.fastq; do
    project_name=$(basename "$file" | cut -d'_' -f1 | tr '[:upper:]' '[:lower:]')

    for tag_num in "${sample_tag_numbers[@]}"; do
        seq="${sample_tags[$tag_num]}"
        if [ -z "$seq" ]; then
            echo "Sample tag number $tag_num is not defined."
            continue
        fi
        output_file="$output_dir/${project_name}st${tag_num}_smk_R2.fastq"
        echo "Extracting reads for sample tag $tag_num from $file..."
        grep -A 3 --no-group-separator "$seq" "$file" > "$output_file"

        # Pair the reads with fastq-pair
        source activate
		conda activate fastq-pair
        smk_reads=$(grep "@" "$output_file" -c)
        r1_file="${input_dir}/$(basename "$file" | sed 's/R2/R1/')"
        fastq_pair "$output_file" "$r1_file" -t "$smk_reads"
        rm *.single.fq "$output_file".*
        mv "${r1_file}.paired.fq" "$output_dir/${project_name}st${tag_num}_smk_R1.fastq"

        # Activate umi_tools environment
        conda activate umi_tools
        umi_file="$output_dir/${project_name}st${tag_num}_smk_R1.fastq"
        whitelist_file="$output_dir/whitelist_${project_name}st${tag_num}.txt"
        
        # Run umi_tools whitelist command
        umi_tools whitelist --stdin "$umi_file" --extract-method=regex --bc-pattern="(?P<cell_1>.{9,12})(?P<discard_1>GTGA)(?P<cell_2>.{9})(?P<discard_1>GACA)(?P<cell_3>.{9})(?P<umi_1>.{6,9})T{3}.*" --expect-cells="$expected_cells" --knee-method=density --log2stderr > "$whitelist_file"
        
        # Check number of lines in whitelist file
        line_count=$(wc -l < "$whitelist_file")
        if [ "$line_count" -lt 100 ]; then
            echo "Whitelist file has less than 100 lines, running alternative command..."
            umi_tools whitelist --stdin "$umi_file" --extract-method=regex --bc-pattern="(?P<cell_1>.{9,12})(?P<discard_1>GTGA)(?P<cell_2>.{9})(?P<discard_1>GACA)(?P<cell_3>.{9})(?P<umi_1>.{6,9})T{3}.*" --set-cell-number="$expected_cells" --log2stderr > "$whitelist_file"
        fi
    done
done

# Process WTA and Abseq reads using whitelist files

whitelist_files=($(ls "$output_dir"/whitelist_*.txt))

# Process each whitelist file one by one
for whitelist_file in "${whitelist_files[@]}"; do
    project_name=$(basename "$whitelist_file" | cut -d'_' -f2 | sed 's/st[0-9]*//')
    tag_num=$(basename "$whitelist_file" | grep -o 'st[0-9]*' | sed 's/st//')
    
    # Remove the extension before further processing
    base_name=$(basename "$whitelist_file" .txt)
    file_name=$(echo "$base_name" | cut -d'_' -f2 | sed 's/st[0-9]*//' | tr '[:lower:]' '[:upper:]')
	
	input_dir="$script_dir"
	for file in "$input_dir"/*Abseq*R2*.fastq.gz; do
		abseq_r2_file="${input_dir}/$(basename "$file" | sed 's/R2/R1/')"
	done
	for file in "$input_dir"/*Abseq*R1*.fastq.gz; do
		abseq_r1_file="${input_dir}/$(basename "$file" | sed 's/R2/R1/')"
	done
	for file in "$input_dir"/*WTA*R2*.fastq.gz; do
		wta_r2_file="${input_dir}/$(basename "$file" | sed 's/R2/R1/')"
	done
	for file in "$input_dir"/*WTA*R1*.fastq.gz; do
		wta_r1_file="${input_dir}/$(basename "$file" | sed 's/R2/R1/')"
	done
	
	echo "file_name: $file_name"
	echo "wta_r1_file: $wta_r1_file"
	echo "wta_r2_file: $wta_r2_file"
	echo "abseq_r1_file: $abseq_r1_file"
	echo "abseq_r2_file: $abseq_r2_file"
	
    # Process WTA files
    (
        echo "Processing WTA files with $whitelist_file..."
        umi_tools extract --extract-method=regex \
                          --bc-pattern="(?P<cell_1>.{9,12})(?P<discard_1>GTGA)(?P<cell_2>.{9})(?P<discard_1>GACA)(?P<cell_3>.{9})(?P<umi_1>.{6,9})T{3}.*" \
                          --stdin "$wta_r1_file" --stdout "${output_dir}/${project_name}st${tag_num}_wta_R1.fastq" \
                          --read2-in "$wta_r2_file" --read2-out="${output_dir}/${project_name}st${tag_num}_wta_R2.fastq" \
                          --whitelist="$whitelist_file"
        
        cat "${output_dir}/${project_name}st${tag_num}_wta_R2.fastq" | sed -e 's/_.* / /' > "${output_dir}/${project_name}st${tag_num}_wta_cleaned_R2.fastq"
        conda deactivate
    ) &

    # Process Abseq files
    (
        echo "Processing Abseq files with $whitelist_file..."
        conda activate umi_tools
        umi_tools extract --extract-method=regex \
                          --bc-pattern="(?P<cell_1>.{9,12})(?P<discard_1>GTGA)(?P<cell_2>.{9})(?P<discard_1>GACA)(?P<cell_3>.{9})(?P<umi_1>.{6,9})T{3}.*" \
                          --stdin "$abseq_r1_file" --stdout "${output_dir}/${project_name}st${tag_num}_abseq_R1.fastq" \
                          --read2-in "$abseq_r2_file" --read2-out="${output_dir}/${project_name}st${tag_num}_abseq_R2.fastq" \
                          --whitelist="$whitelist_file"
        
        cat "${output_dir}/${project_name}st${tag_num}_abseq_R2.fastq" | sed -e 's/_.* / /' > "${output_dir}/${project_name}st${tag_num}_abseq_cleaned_R2.fastq"
        conda deactivate
    ) &

    # Wait for both processes to finish
    wait

    # Pairing reads with fastq_pair and clean up
    (
        conda activate fastq-pair
        gunzip "$wta_r1_file"
		wta_reads=$(grep "@" "${output_dir}/${project_name}st${tag_num}_wta_cleaned_R2.fastq" -c)
        abseq_reads=$(grep "@" "${output_dir}/${project_name}st${tag_num}_abseq_cleaned_R2.fastq" -c)
        
        fastq_pair "${output_dir}/${project_name}st${tag_num}_wta_cleaned_R2.fastq" "${wta_r1_file%.gz}" -t $wta_reads
        rm *.single.fq "${output_dir}/${project_name}st${tag_num}_wta_cleaned_R2.fastq.*" "${output_dir}/${project_name}st${tag_num}_wta_R1.fastq"
        mv "${wta_r1_file%.gz}.paired.fq" "${output_dir}/${project_name}st${tag_num}_wta_cleaned_R1.fastq"
        
        fastq_pair "${output_dir}/${project_name}st${tag_num}_abseq_cleaned_R2.fastq" "$abseq_r1_file" -t $abseq_reads
        rm *.single.fq "${output_dir}/${project_name}st${tag_num}_abseq_cleaned_R2.fastq.*" "${output_dir}/${project_name}st${tag_num}_abseq_R1.fastq"
        mv "${abseq_r1_file%.gz}.paired.fq" "${output_dir}/${project_name}st${tag_num}_abseq_cleaned_R1.fastq"
        
        conda deactivate
    )
done

echo "All processes completed!"
