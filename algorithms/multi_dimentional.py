import numpy as np
import numdifftools as nd
from algorithms.one_dimentional import fibonacci

def euclidean_norm(x):
    return np.sqrt(sum(x**2))

def gauss_seidel(f, x, epsilon=.0001):
    norm = euclidean_norm
    one_opt = fibonacci
    inf = False
    n = len(x)
    y = np.array(x)
    x = np.array([z+5 for z in x])
    e = lambda j: np.eye(1, n, j)[0]

    while not (norm(y - x) < epsilon) and not inf:
        x = y
        j = 0
        while j < n:
            lambda_j = one_opt(
                lambda z: f(y + z*e(j)),
                -1,
                1,
                epsilon,
                solution=True,
            )
            y = y + lambda_j * e(j)
            j += 1
        if abs(f(y)) > 10**15:
            solution = y
            inf = True
            

    solution = y
    if not inf:
        val = f(solution)
    else:
        val = np.sign(f(solution))*float('inf')
    return {
        'value': val,
        'solution': solution,
    }

# Прикол с Коши
def gradient_descent(f, x, epsilon=.0001):
    inf = False
    norm = euclidean_norm
    one_opt = fibonacci
    norm_of_grad = float('inf')
    while not norm_of_grad < epsilon and not inf:
        grad = nd.Gradient(f)(x)
        norm_of_grad = norm(grad)
        s = -grad/norm_of_grad
        l = one_opt(
            lambda y: f(x + y*s),
            -1,
            1,
            solution=True,
        )
        x += l*s
        if abs(f(x)) > 10**15:
            val = np.sign(f(x))*float('inf')
            inf = True

    solution = x
    if not inf:
        val = f(solution)

    return {
        'value': val,
        'solution': solution,
    }


