from tkinter import Tk, Label, Button
from collections import Counter

def show_min_degree_protein(distance_matrix, protein_names):
    # Calculate the degree of each protein
    degree_counts = Counter((distance_matrix > 0).sum(axis=1))

    # Find the protein with the minimum degree
    min_degree = min(degree_counts.values())
    min_degree_proteins = [protein for protein, degree in degree_counts.items() if degree == min_degree]

    # Create a new window to display the protein with the minimum degree
    root = Tk()
    root.title("Protein with Minimum Degree")

    label = Label(root, text=f"Protein(s) with Minimum Degree ({min_degree}):")
    label.pack(pady=10)

    for protein in min_degree_proteins:
        label = Label(root, text=protein)
        label.pack()

    root.mainloop()
