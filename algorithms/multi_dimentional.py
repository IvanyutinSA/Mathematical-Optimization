# import sys 
# if "/home/timofey/Desktop/Maths_optimization/Mathematical-Optimization" not in sys.path:
#     sys.path.append("/home/timofey/Desktop/Maths_optimization/Mathematical-Optimization")

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


import numpy as np
from sympy import Symbol, simplify
import sympy
import numdifftools as nd 

from algorithms import one_dimentional

def onedim_gs_min(f):
    values = dict()
    for l in {2, 5, 10, 15, 100, 1_000, 10_000}: 
        pair = one_dimentional.golden_section(f, -l, l)
        values[pair[0]] = pair[1]
    return min(values, key=values.get)

def simplif_minimiz(currX, S, symb_f, f):
    X2 = currX + Symbol('l')*S
    simplfd = simplify(symb_f(X2))
    final_f = lambda x: simplfd.subs({Symbol('l'): x[0]}) if not str(x)[-1].isdigit() else simplfd.subs({Symbol('l'): x})
    minX = onedim_gs_min(final_f)
    xNew = [float(round((X2[i].subs({Symbol('l'): minX})), 3)) for i in range(len(X2))]
    return xNew, f(xNew) 


def Hook_Jeeves(f, symb_f, deltaX, currX: np.array):
    grad = [1]
    eps = 1e-02
    while not all(abs(component) < eps for component in grad):
        initialX = [comp for comp in currX ]
        for i in range(len(currX)):
            L = [x for x in currX]
            R = [x for x in currX]
            L[i] -= deltaX[i]
            R[i] +=  deltaX[i]
            d = {point:f(point) for point in {tuple(currX), tuple(R), tuple(L)}} 
            minX = list(min(d, key=d.get))
            currX = minX
        if initialX != currX:
            S = np.array(currX) - np.array(initialX)
            point_value = simplif_minimiz(currX, S, symb_f, f)
            xNew = point_value[0]
            grad = nd.Gradient(f)(xNew)
            currX = xNew
            initialX = [comp for comp in currX]
        else: 
            deltaX /= 2

    return xNew, f(xNew) 
    
symb_f = lambda x: 10 - (x[0] - 3)*sympy.exp(-(x[0] - 3)) - (x[1] - 2)*sympy.exp(-(x[1] - 2))

f = lambda x: 10 - (x[0] - 3)*np.exp(-(x[0] - 3)) - (x[1] - 2)*np.exp(-(x[1] - 2)) 




def GS_orthog(old: np.array):
    old = np.array(old).astype(float)
    new = np.array([comp for comp in old])
    for i in range(len(new)):
        new[i] = [float(x) for x in old[i]]
        for j in range(i):
            value = np.dot(old[i], new[j])/np.dot(new[j], new[j]) * new[j]
            new[i] = new[i] - value
    return new


def Rosenbrock(x1, f, symb_f, eps=1e-5):
    dim = len(x1)
    s = np.zeros(shape=(dim, dim))
    for i in range(dim): s[i] = np.eye(1, dim, i)

    k = 0
    y = x1
    xOld = np.array([comp for comp in x1])
    xNew = []
    vNew = float('inf')
    while (abs(vNew - f(xOld)) >= eps):
        k += 1
        if k > 1:         
            left_side = xNew - xOld
            alpha = np.linalg.solve(s, left_side)
            z = [comp for comp in s]
            for j in range(dim): z[j] = s[j] if alpha[j] == 0 else sum(alpha[i]*s[i] for i in range(j, dim))
            s = GS_orthog(z)
        for i in range(dim):
            point_value = simplif_minimiz(y, s[i], symb_f, f)
            y, vNew = point_value[0], point_value[1]
        if len(xNew) > 0: 
            xOld = np.array([comp for comp in xNew])
        xNew = np.array([comp for comp in y]) 

    return xNew, f(xNew)

# import time
# start_time = time.time()
# hook()
# h_time = time.time() - start_time   # ~2,9сек
# print(f"hook: {h_time}")

# rosen()
# r_time = time.time() - start_time
# print(f"rosen: {r_time}")  # ~5,7сек

# замерить время выполнения каждого из методов