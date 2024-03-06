import numpy as np
from one_dimentional import fibonacci

def euclidean_norm(x, y):
    return sum((x-y)**2)

def gauss_seidel(f, x, e=.0001):
    y = np.array(x)
    x = [y.copy(), y.copy()+5]
    d = len(y)
    z = np.zeros(d)
    while not euclidean_norm(x[0], x[1]) < e:
        for j in range(len(y)):
            z[j] = fibonacci(
                lambda alpha: f(y+alpha*np.eye(1, len(y), j))
                    )

    return f(*x[-1])


