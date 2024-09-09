import networkx as nx

def create_graph(distance_matrix, protein_names, connected_nodes):
    G = nx.Graph()
    print("started")
    print(G.nodes(data=True))

    for protein in protein_names:
        G.add_node(protein, label=protein)

    num_nodes = len(distance_matrix)

    for i in range(num_nodes):
        for j in range(i + 1, num_nodes):
            G.add_edge(protein_names[i], protein_names[j], weight=distance_matrix[i][j])

    return G
