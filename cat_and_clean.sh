# Concatenate all of the translated files together to use to make the aliginments
# and save it to a new file that has all of the sequences across all samples
cat *_prot.fas > protein_alignment_seqs.fas

# Replace all instances of 'X' in the sequences with '-' so that the alignment software can 
# read the sequences properly
sed -i '/^>/! s/X/-/g' protein_alignment_seqs.fas