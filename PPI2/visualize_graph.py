import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
from tkinter import Tk
from matplotlib.widgets import Button

def visualize_graph(G, distance_matrix, protein_names):
    def find_drug_target(event):
        min_degree_nodes = find_drug_target_nodes()
        label.set_text(f"Potential Drug Target: {', '.join(min_degree_nodes)}")
        plt.draw()
    
    def find_drug_target_nodes():
        degrees = {node: 0 for node in G.nodes()} 
        for i in range(len(protein_names)):
            for j in range(i + 1, len(protein_names)):
                if distance_matrix[i][j] != 0:  
                    degrees[protein_names[i]] += 1
                    degrees[protein_names[j]] += 1
        
        min_degree = min(degrees.values())
        min_degree_nodes = [node for node, degree in degrees.items() if degree == min_degree]
        return min_degree_nodes

    plt.figure(figsize=(12, 6))
    plt.subplot(1, 2, 1)
    pos = nx.spring_layout(G)
    labels = {node: G.nodes[node]['label'] for node in G.nodes}
    nx.draw(G, pos, with_labels=True, labels=labels, font_weight='bold')
    plt.title('Protein Interaction Network')

    plt.subplot(1, 2, 2)
    sns.heatmap(distance_matrix, annot=True, cmap='viridis', fmt='g', xticklabels=protein_names, yticklabels=protein_names)
    plt.title('Distance Matrix')

    ax_button = plt.axes([0.1, 0.05, 0.1, 0.05])
    button = Button(ax_button, 'Find Drug Target')
    button.on_clicked(find_drug_target)

    label = plt.figtext(0.1, 0.02, "")

    plt.show()

