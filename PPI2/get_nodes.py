from tkinter import Tk, Label, Button, filedialog

class ConnectedNodesGUI:
    def __init__(self, protein_names):
        self.protein_names = protein_names
        self.connected_nodes = {}
        self.current_index = 0

        self.root = Tk()
        self.root.title("Connected Nodes GUI")
        self.file_path = None

        self.label = Label(self.root, text="Select file containing protein names and their connected nodes:")
        self.label.grid(row=0, column=0)

        self.select_button = Button(self.root, text="Select File", command=self.select_file)
        self.select_button.grid(row=1, column=0)

    def select_file(self):
        self.file_path = filedialog.askopenfilename(title="Select a File")
        if self.file_path:
            print(f"File selected: {self.file_path}")
            self.submit()  # Call submit method after selecting the file
   
    def read_connected_nodes_from_file(self, file_path):
        connected_nodes = {}

        with open(file_path, 'r') as file:
            for line in file:
                protein, *nodes = line.strip().split(',')
                connected_nodes[protein] = nodes

        return connected_nodes

    def run(self):
        self.root.mainloop()
        return self.connected_nodes

    def submit(self):
        self.connected_nodes = self.read_connected_nodes_from_file(self.file_path)
        self.root.destroy()  

