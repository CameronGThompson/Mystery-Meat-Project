from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pathlib import Path

input_dir = Path("homologs/")
output_dir = Path("protein_sequences/")


for fasta_file in input_dir.glob("*.fas"):

    output_file = output_dir / f"{fasta_file.stem}_prot.fas"
    aa_sequences = []

    for record in SeqIO.parse(fasta_file, "fasta"):

        dna_seq = record.seq

        best_orf = ""

        best_frame = None

        for frame in range(3):
            protein = dna_seq[frame:].translate(table=2, to_stop=False)
            fragments = str(protein).split("*")
            longest_in_frame = max(fragments, key=len)
            if len(longest_in_frame) > len(best_orf):
                best_orf = longest_in_frame
                best_frame = frame

        best_protein = Seq(best_orf)

        aa_record = SeqRecord(
            best_protein,
            id=record.id,
            description=record.description,
            annotations={
                "type": "longest ORF",
                "frame": best_frame
            }
        )

        aa_sequences.append(aa_record)

    SeqIO.write(aa_sequences, output_file, "fasta")

    print(f"Created: {output_file.name}")
