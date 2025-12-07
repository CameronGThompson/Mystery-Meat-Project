# Mystery-Meat-Project

## Scenario
A suitcase containing several packages of unidentified frozen meat has been intercepted at Heathrow Airport by UK Border Force officers during routine inspections. The shipment had no import documentation, raising concerns that the meat may originate from protected or endangered species. 
Given the serious implications for public health, food safety, and wildlife conservation, the case has been referred to the Animal and Plant Health Agency (APHA). As a bioinformatician assisting with the inquiry, you have been tasked with determining the species of origin for each meat sample. DNA barcoding has already been performed in the laboratory, and you have now received the resulting FASTA sequence files for analysis. 

## Project Goals
The main goal of this project was to take the results of DNA barcoding and prepare the sequences for use in a phylogenetic analysis which would identify what species each sample belonged to.

## Key Features
This project includes:

* A Bash script used to clean up FASTQ files and convert them to FASTA
* A Python script which takes each fasta file and translates the sequences in each file into amino acid sequences
* A secondary Bash script which combines all of the translated sequences into a single FASTA file.


## Repository Structure & Usage Instructions

* clean_and_convert.sh:
  * This shell script takes the raw FASTQ sequence files, normalises their headers, converts the files from FASTQ to FASTA and masks bases with a quality score below 20,        and combines the part sequences in each sample file into a single sequence per sample.
  * You will need to download the seqtk package from GitHub in order to run this script locally.
  * The only requirement for filenames is that the original FASTQ sequence files follow the naming format: e.g. sampleA_part1
  * You will need to save this script into the same directory that the original FASTQ files are saved.

* DNA_translate.py:
  * This Python script loops over FASTA files in a specified input directory, then loops over the sequences in those files and searches for the longest ORF in each sequence.
  * It then takes those ORFs, translates them using the vertebrate mitochondrial genetic code, and deposits them in a specified output directory.
  * This script requires BioPython, so ensure you have this installed before running.
  * In order to run this script, you will need to store FASTA files you want to translate in a directory called "homologs".
  * You will also have to create a directory called "protein_sequences" in which the translated sequences will be deposited.
  * You should save this file into the same directory which contains the "homologs" and "protein_sequences" directories.

* cat_and_clean.sh:
  * This shell script concatenates the translated FASTA files into a single FASTA file for alignment and then replaces all instances of "X" with "-" to allow the alignment      to run correctly.
  * You will need to save this script in the "protein_sequences" directory as well.

## Dependencies
You will need to install the seqtk GitHub package to run the clean_and_convert.sh script. This package can be downloaded [here](https://github.com/lh3/seqtk). Instructions for installation can also be found at that link. Ensure that the seqtk executable is in your PATH before trying to run the script.
  



