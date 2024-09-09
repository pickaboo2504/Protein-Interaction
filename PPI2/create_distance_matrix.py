import numpy as np
from edit_distance import edit_distance

def create_distance_matrix(sequences, connected_nodes):
    num_sequences = len(sequences)
    distance_matrix = np.zeros((num_sequences, num_sequences))

    for i in range(num_sequences):
        for j in range(i + 1, num_sequences):
            if sequences[i][0] in connected_nodes and sequences[j][0] in connected_nodes[sequences[i][0]]:
                distance = edit_distance(sequences[i][1], sequences[j][1])
                distance_matrix[i][j] = distance
                distance_matrix[j][i] = distance

    return distance_matrix
