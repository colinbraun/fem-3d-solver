# File to test out little pieces of code before using them in main.py
import numpy as np

x1, y1, z1 = 1, 2, 3
x2, y2, z2 = 1, 2, 3
x3, y3, z3 = 1, 2, 3
x4, y4, z4 = 1, 2, 3
init_array = np.array([[1, x1, y1, z1],
                       [1, x2, y2, z2],
                       [1, x3, y3, z3],
                       [1, x4, y4, z4]], dtype=float)

all_cofactors = np.zeros([4, 4])
# Iterate over each row
for row in range(4):
    cofactors = np.zeros([4])
    # Iterate over each column, computing the cofactor determinant of the row + column combination
    for col in range(4):
        # Compute the cofactor (remove the proper row and column and compute the determinant)
        cofactors[col] = np.linalg.det(np.delete(np.delete(init_array, row, axis=0), col, axis=1))
    all_cofactors[row] = cofactors
