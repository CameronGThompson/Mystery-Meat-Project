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

## Repository Structure




