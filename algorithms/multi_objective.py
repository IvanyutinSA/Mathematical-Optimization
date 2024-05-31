import numpy as np

import sys 

# if "/home/timofey/Desktop/Maths_optimization/Mathematical-Optimization/algorithms" not in sys.path:
#     sys.path.append("/home/timofey/Desktop/Maths_optimization/Mathematical-Optimization/algorithms")
from multi_dimentional import gauss_seidel

def additive_convolution(functions: list, initial_x):
    dim = len(functions)
    alpha = 1/dim
    F = lambda x: sum(alpha*f(x) for f in functions)
    result = gauss_seidel(F, initial_x)
    x_new = result['solution']

    return x_new
