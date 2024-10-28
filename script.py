# Import libraries
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline

# Parse the human sequences
seqs = SeqIO.to_dict(SeqIO.parse(open("human.fa"), "fasta"))

# Setup output file
output = open("blast_results.txt", "w")
output.write("HumanSeqID\tMouseSeqID\tAlignment\tE-value\tBitscore\n")

# Iterate through each sequence
for human_id, human_seq in seqs.items():
    """
    Prints the following to a separate file:
        - Human sequence ID
        - Mouse ID of the most similar homolog
        - Corresponding alignment
        - E-value
        - Bitscore
    """
    with open("temp_human.fa", "w") as temp_file:
        SeqIO.write(human_seq, temp_file, "fasta")

    # Searches the mouse database with BLAST
    blast = NcbiblastpCommandline(query="temp_human.fa",
                                  db="mouse_db",
                                  evalue=0.001,
                                  outfmt=5,
                                  out="temp_blast.xml")
    
    stdout, stderr = blast()

    with open("temp_blast.xml") as result_handle:
        blast_record = NCBIXML.read(result_handle)

        if blast_record.alignments:
            best_alignment = blast_record.alignments[0]
            hsp = best_alignment.hsps[0]

            output.write(f"{human_id}\t{best_alignment.hit_id}\t{hsp.sbjct}\t{hsp.expect}\t{hsp.bits}\n")

output.close()