# Load in necessary packages
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pathlib import Path

# Specify input and output directories
input_dir = Path("homologs/")
output_dir = Path("protein_sequences/")

# Loop over fasta files in the input directory
for fasta_file in input_dir.glob("*.fas"):
    # Specifyu the name of the output file
    output_file = output_dir / f"{fasta_file.stem}_prot.fas"
    # Create an empty string to add amino acids to
    aa_sequences = []
    # Loop through sequences in each fasta file
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Nucleotide sequence from the FASTA record
        dna_seq = record.seq
        # Will store the longest ORF found across all frames
        best_orf = ""
        # Will store which frame (0,1,2) contains the longest ORF
        best_frame = None
        # Translate sequence in all three forward reading frames
        for frame in range(3):
            # Translate starting at this frame; keep stop codons (to_stop=False)
            # and specify genetic code 2 for vertebrate mitochondrial
            protein = dna_seq[frame:].translate(table=2, to_stop=False)
            # Split translation into fragments separated by stop codons ("*")
            # Each fragment is a potential ORF (no internal stops)
            fragments = str(protein).split("*")
            # Find the longest continuous amino-acid stretch (longest ORF) in this frame
            longest_in_frame = max(fragments, key=len)
            # If this ORF is longer than the best one we've seen so far, store it
            if len(longest_in_frame) > len(best_orf):
                best_orf = longest_in_frame
                best_frame = frame
        # Convert the longest ORF string back into a Seq object
        best_protein = Seq(best_orf)
        # Create a new SeqRecord containing the longest ORF for this sequence
        aa_record = SeqRecord(
            best_protein,
            # Keep the same FASTA ID
            id=record.id,
            # Preserve full original header
            description=record.description,
            annotations={
                "type": "longest ORF",
                # Which reading frame the ORF came from
                "frame": best_frame
            }
        )
        # Add to output list
        aa_sequences.append(aa_record)
    # Write all translated longest-ORF sequences to a FASTA file
    SeqIO.write(aa_sequences, output_file, "fasta")
    # Print a message stating that each file has been created
    print(f"Created: {output_file.name}")
