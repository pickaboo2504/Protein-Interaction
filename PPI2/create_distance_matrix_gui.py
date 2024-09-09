from tkinter import Tk, Label, Entry, Button, messagebox
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import edit_distance
from target_protein import show_min_degree_protein

class DistanceMatrixGUI:
    def __init__(self, root, sequences):
        self.root = root
        self.root.title("Connected Proteins")

        self.sequences = sequences
        self.connected_proteins = []

        # Label and Entry for protein pairs
        self.label = Label(root, text="Enter connected protein pairs (comma-separated):")
        self.label.pack(pady=10)

        self.entry = Entry(root)
        self.entry.pack(pady=10)

        # Button to calculate distance matrix
        self.calculate_button = Button(root, text="Calculate Distance Matrix", command=self.calculate_matrix)
        self.calculate_button.pack(pady=10)

        # Button to show protein with minimum degree
        self.next_button = Button(root, text="Next", command=self.show_min_degree)
        self.next_button.pack(pady=10)

    def calculate_matrix(self):
        # Get user input for connected protein pairs
        input_text = self.entry.get()
        pairs = [pair.strip() for pair in input_text.split(',')]

        # Validate input
        if not all(',' in pair for pair in pairs):
            messagebox.showerror("Error", "Invalid input. Please use a comma-separated format.")
            return

        # Process connected proteins
        self.connected_proteins = [pair.split(',') for pair in pairs]
        self.root.destroy()  # Close the GUI

    def show_min_degree(self):
        if self.connected_proteins:
            distance_matrix, protein_names = calculate_distance_matrix(self.connected_proteins, self.sequences)
            show_min_degree_protein(distance_matrix, protein_names)
        else:
            messagebox.showerror("Error", "Please calculate the distance matrix first.")

def calculate_distance_matrix(connected_proteins, sequences):
    num_sequences = len(sequences)
    distance_matrix = np.zeros((num_sequences, num_sequences))

    for pair in connected_proteins:
        i, j = pair
        i, j = int(i) - 1, int(j) - 1  # Adjust for 0-based indexing

        if 0 <= i < num_sequences and 0 <= j < num_sequences:
            distance = edit_distance.edit_distance(sequences[i], sequences[j])
            distance_matrix[i][j] = distance
            distance_matrix[j][i] = distance

    return distance_matrix, [seq[0] for seq in sequences]

def create_distance_matrix_gui(sequences):
    root = Tk()
    gui = DistanceMatrixGUI(root, sequences)
    root.mainloop()

# Example usage:
# Uncomment and run the following line to test the GUI
# create_distance_matrix_gui(["sequence1", "sequence2", "sequence3"])
