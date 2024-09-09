from tkinter import Tk, filedialog
from Bio import SeqIO

def get_sequences_from_file():
    root = Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename(title="Select Protein Sequences File", filetypes=[("FASTA Files", "*.fa")])

    sequences = []
    protein_names = []

    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(str(record.seq))
        description_parts = record.description.split(" ")
        protein_name = description_parts[-1][1:-1]
        protein_names.append(protein_name)

    return sequences, protein_names
