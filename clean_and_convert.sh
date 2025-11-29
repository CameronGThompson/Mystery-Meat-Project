#!/bin/bash

#######################################
# SAFETY SETTINGS
#######################################

# -e: exit immediately if any command fails
# -u: treat unset variables as errors
# -o pipefail: catch errors in piped commands
set -euo pipefail

#######################################
# ERROR HANDLER
#######################################

# If any command fails, this trap prints the line number
# and the command that caused the failure.

trap 'echo "ERROR on line $LINENO: Command failed -> $BASH_COMMAND" >&2; exit 1' ERR

#######################################
# FUNCTIONS
#######################################

# Check if a file exists, otherwise stop the script
check_file() {
    if [[ ! -f "$1" ]]; 
        then
            echo "Missing input file: $1" >&2
            exit 1

    fi
}

# Print a "starting progress" message
step() {
    echo -e "\n -> $1"
}

# Print a "finished successfully" message
done_step() {
    echo "$1"
}

# Replace headers in the FASTQ files using a temporary file

fix_header() {
    local search="$1"
    local replace="$2"
    local file="$3"
    # Make a temporary file
    tmp=$(mktemp)
    sed "s/${search}/${replace}/" "$file" > "$tmp"
    # Replace the file
    mv "$tmp" "$file"
}

# Convert FASTQ to FASTA and mask low quality bases using 
# seqtk GitHub package within a function of my own making
# including error checking

convert_fastq() {
    local input="$1"
    local output="$2"

    step "Converting $input to $output"
    seqtk seq -a -q20 -n N "$input" > "$output"
    
    # Ensure the output file is not empty
    if [[ ! -s "$output" ]];
        then
            echo "Conversion failed for $input" >&2
            exit 1
    fi

    done_step "$input converted"
}

# Concatenate converted FASTA parts
concat_files() {
    local output="$1"
    shift

    step "Concatenating to $output"
    cat "$@" > "$output"

    # Ensure concatenated file exists and isn't empty
    if [[ ! -s "$output" ]];
        then
            echo "Concatenation failed for $output" >&2
            exit 1
    fi

    done_step "$output created"
}

# Convert concatenated FASTA files into single-line sequences with new headers
create_seq() {
    local infile="$1"
    local header="$2"
    local outfile="$3"

    if [[ ! -f "$infile" ]]; 
        then
            echo "Error: Input file '$infile' not found"
            return 1
    fi

    grep -v "^>" "$infile" | tr -d '\n' | sed "1i >$header" > "$outfile"
}

######################################
# CHECK DEPENDENCIES
######################################

step "Checking dependencies"

# Ensure seqtk is installed and in PATH
if ! command -v seqtk &>/dev/null;
    then
        echo "seqtk is not installed or not in PATH"
        exit 1
fi

done_step "Dependencies OK"

########################################
# CHECK ALL REQUIRED INPUT FASTQ FILES
########################################

step "Checking input FASTQ files"

# Ensure all FASTQ files exist before doing anything else
for file in sample{A,B,C,D}_part{1,2,3}.FASTQ;
    do
        check_file "$file"
done

done_step "All input files present"

#######################################
# FIX HEADERS
#######################################

step "Normalising FASTQ headers"

# Normalise header formats across all samples
fix_header "@sampleB_part2" "@sampleB_part_2" sampleB_part2.FASTQ
fix_header "@sampleB_part3" "@sampleB_part_3" sampleB_part3.FASTQ
fix_header "@sampleC_part1" "@sampleC_part_1" sampleC_part1.FASTQ
fix_header "@sampleC_part2" "@sampleC_part_2" sampleC_part2.FASTQ
fix_header "@sampleC_part3" "@sampleC_part_3" sampleC_part3.FASTQ
fix_header "@sampleD_part1" "@sampleD_part_1" sampleD_part1.FASTQ
fix_header "@sampleD_part2" "@sampleD_part_2" sampleD_part2.FASTQ
fix_header "@sampleD_part3" "@sampleD_part_3" sampleD_part3.FASTQ

done_step "FASTQ headers normalised"

##########################################
# FIX DUPLICATED HEADER IN SAMPLEB_PART1
##########################################

step "Removing duplicate header from sampleB_part1"

# Use awk to remove duplicated header, move the result to a temporary file,
# and normalise the header format
awk '!(($0 == prev) && /^@/) {print} {prev = $0}' sampleB_part1.FASTQ > tmp && mv tmp sampleB_part1.FASTQ && fix_header "@sampleB_part1" "@sampleB_part_1" sampleB_part1.FASTQ

done_step "Duplicate header removed"

#########################################
# FASTQ to FASTA CONVERSION
#########################################

convert_fastq sampleA_part1.FASTQ sampleA1_conv.fa 
convert_fastq sampleA_part2.FASTQ sampleA2_conv.fa
convert_fastq sampleA_part3.FASTQ sampleA3_conv.fa

convert_fastq sampleB_part1.FASTQ sampleB1_conv.fa
convert_fastq sampleB_part2.FASTQ sampleB2_conv.fa
convert_fastq sampleB_part3.FASTQ sampleB3_conv.fa

convert_fastq sampleC_part1.FASTQ sampleC1_conv.fa
convert_fastq sampleC_part2.FASTQ sampleC2_conv.fa
convert_fastq sampleC_part3.FASTQ sampleC3_conv.fa

convert_fastq sampleD_part1.FASTQ sampleD1_conv.fa
convert_fastq sampleD_part2.FASTQ sampleD2_conv.fa
convert_fastq sampleD_part3.FASTQ sampleD3_conv.fa

#################################################
# CONCATENATE PARTS INTO SINGLE FASTA SEQUENCES
#################################################

concat_files sampleA_cat.fa sampleA1_conv.fa sampleA2_conv.fa sampleA3_conv.fa
concat_files sampleB_cat.fa sampleB1_conv.fa sampleB2_conv.fa sampleB3_conv.fa
concat_files sampleC_cat.fa sampleC1_conv.fa sampleC2_conv.fa sampleC3_conv.fa
concat_files sampleD_cat.fa sampleD1_conv.fa sampleD2_conv.fa sampleD3_conv.fa

##################################################
# COMBINE PARTS INTO SINGLE SEQUENCE PER SAMPLE
##################################################

# Combine sequence parts into one singular sequence per sample by removing all headers, 
#joining the sequences together on a single line, and then adding a new header at the top

step "Joining parts into sample sequences"

create_seq sampleA_cat.fa sample_A_full sampleA_final.fas
create_seq sampleB_cat.fa sample_B_full sampleB_final.fas
create_seq sampleC_cat.fa sample_C_full sampleC_final.fas
create_seq sampleD_cat.fa sample_D_full sampleD_final.fas

done_step "Final sample sequences created"

# Finally, remove all transitional files (e.g. anything ending with _cat.fa or _conv.fa)
#so you are only left with the original files and the fully cleaned files

step "removing transitional files"

rm -R *_cat.fa *_conv.fa

done_step "Transitional files removed"

# Print a final message stating that the pipeline is complete
echo -e "\n Pipeline completed successfully."