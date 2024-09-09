import numpy as np

def edit_distance(s1, s2):
    len_s1 = len(s1) + 1
    len_s2 = len(s2) + 1

    matrix = [[0] * len_s2 for _ in range(len_s1)]

    for i in range(len_s1):
        matrix[i][0] = i

    for j in range(len_s2):
        matrix[0][j] = j

    for i in range(1, len_s1):
        for j in range(1, len_s2):
            cost = 0 if s1[i - 1] == s2[j - 1] else 1
            matrix[i][j] = min(
                matrix[i - 1][j] + 1,
                matrix[i][j - 1] + 1,
                matrix[i - 1][j - 1] + cost
            )

    return matrix[len_s1 - 1][len_s2 - 1]
