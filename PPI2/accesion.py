from Bio import Entrez, SeqIO
import tkinter as tk
from tkinter import Entry, Label, Button, StringVar

def get_sequence_from_accession(accession):
    Entrez.email = "nazia.2109018@bau.edu.bd"  # Replace with your email
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    return record.seq

def fetch_sequence():
    user_input = entry.get()
    
    if user_input.startswith("http"):
        pass
    else:
        
        sequence = get_sequence_from_accession(user_input)
        result_var.set("Sequence obtained from accession:\n" + str(sequence))


window = tk.Tk()
window.title("Sequence Retrieval")

# Create and pack widgets
label = Label(window, text="Enter the accession number or link:")
label.pack(pady=10)

entry = Entry(window)
entry.pack(pady=10)

result_var = StringVar()
result_label = Label(window, textvariable=result_var)
result_label.pack(pady=10)

fetch_button = Button(window, text="Fetch Sequence", command=fetch_sequence)
fetch_button.pack(pady=10)

# Start the GUI event loop
window.mainloop()
