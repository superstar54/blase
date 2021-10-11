import numpy as np

a = np.array([[1, 2, 3], [4, 5, 6]])
b = np.ones((21, 1))
# print(b)
a = np.tile(a, (21, 1))
print(a*b)
