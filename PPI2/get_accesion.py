from tkinter import Tk, filedialog, Button, Entry, Label
from Bio import SeqIO, Entrez
from networkx import get_node_attributes
import create_graph
import visualize_graph
import create_distance_matrix
import get_nodes
class FileSelectionApp:
    def __init__(self, root):
        self.root = root
        self.root.title("File Selection App")

        self.file_path = None
        self.accession_number = None 

       
        self.accession_entry = Entry(root)
        self.accession_entry.grid(row=0, column=0, padx=10, pady=10)
        accession_label = Label(root, text="Accession Number:")
        accession_label.grid(row=0, column=1, padx=10, pady=10)
        self.select_button = Button(root, text="Select File", command=self.select_file)
        self.select_button.grid(row=1, column=0, columnspan=2, pady=10)
        self.next_button = Button(root, text="Next", command=self.proceed)
        self.next_button.grid(row=2, column=0, columnspan=2, pady=10)

    def select_file(self):
        self.file_path = filedialog.askopenfilename(title="Select a File", filetypes=[("FASTA Files", "*.fa")])
        if self.file_path:
            print(f"File selected: {self.file_path}")
            self.proceed() 
        else:
            print("No file selected.")

    def proceed(self):
        accession_number = self.accession_entry.get()
        if accession_number:
            sequences = [(accession_number, self.download_sequence(accession_number))]
            protein_names = [accession_number]
            self.accession_number = accession_number
        elif self.file_path:
            sequences = []
            protein_names = []

            for record in SeqIO.parse(self.file_path, "fasta"):
                sequence = str(record.seq)
                protein_name = record.description.split(" ")[-1][1:-1]
                sequences.append((protein_name, sequence))
                protein_names.append(protein_name)
        else:
            print("No accession number or file selected.")
            return None

        print("Proceeding with sequences:", sequences)
        print("Protein Names:", protein_names)

        return sequences, protein_names


    def download_sequence(self, accession_number):
        try:
            Entrez.email = "nazia.2109018@bau.edu.bd"  # Provide your email for NCBI Entrez

            handle = Entrez.efetch(db="protein", id=accession_number, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            handle.close()

            return str(record.seq)
        except Exception as e:
            print(f"Error downloading sequence: {e}")
            return None

    
def main():
    root = Tk()
    app = FileSelectionApp(root)
    root.mainloop()
    accession_number = app.accession_number

    if accession_number:
        print("Accession Number:", accession_number)
        sequences = [(f"Accession_{accession_number}", app.download_sequence(accession_number))]
        protein_names = [f"Accession_{accession_number}"]
    elif app.file_path:
        print(f"File selected: {app.file_path}")
        sequences = []
        protein_names = []

        for record in SeqIO.parse(app.file_path, "fasta"):
            sequence = str(record.seq)
            protein_name = record.description.split(" ")[-1][1:-1]
            sequences.append((protein_name, sequence))
            protein_names.append(protein_name)

        print("Protein Names:", protein_names)
        print("Sequences:")
        for protein_name, seq in sequences:
            print(f"{protein_name}: {seq}")
        gui = get_nodes.ConnectedNodesGUI(protein_names)
        user_input = gui.run()
        if user_input:
            print("Protein Names and Connected Nodes:")
            for protein_name, connected_nodes in user_input.items():
                print(f"{protein_name}: {connected_nodes}")
            distance_matrix = create_distance_matrix.create_distance_matrix(sequences, user_input)
            G = create_graph.create_graph(distance_matrix, protein_names, connected_nodes)
            visualize_graph.visualize_graph(G, distance_matrix, protein_names)


    else:
        print("No accession number or file selected.")
        return

if __name__ == "__main__":
    main()

