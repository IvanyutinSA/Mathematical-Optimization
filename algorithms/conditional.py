import sympy as sp
import numpy as np

def add_lagrange_multipliers(f, temp_constraints, vars):
    n = len(vars)
    vars.extend([sp.symbols(f"add_lagrange_multipliers{k}") for k in range(1, len(temp_constraints)+1)])
    h = f + sum(vars[n+i] * temp_constraint for i, temp_constraint in enumerate(temp_constraints))
    return h

def lagrange_multipliers(f, temp_constraints, vars):
    h = add_lagrange_multipliers(f, temp_constraints, vars)
    h_diffs = [sp.diff(h, var) for var in vars] 
    solutions = sp.solve(h_diffs, vars, domain=sp.Reals)

    if not solutions:
        return "no stationary points"
    print(solutions)
    try:
        return min(h.evalf(subs={p: q for p, q in zip(vars, vector)}) for vector in solutions)
    except:
        return h.evalf(subs=solutions)

# import sys 
# if "/home/timofey/Desktop/Maths_optimization/Mathematical-Optimization" not in sys.path:
#     sys.path.append("/home/timofey/Desktop/Maths_optimization/Mathematical-Optimization")

import numdifftools as nd
from itertools import product
from algorithms import one_dimentional


def Zoutendijk(f, constraints, dim, Biggest_number_in_constr):
    x = np.zeros(dim)

    flag = False

    while not flag:
        I = set()
        for i in range(len(constraints)): 
            if np.abs(constraints[i](x)) < 1e-4: I.add(i)

        s = [sp.Symbol(f's{i+1}') for i in range(dim)]

        temp_cons = [lambda s: np.dot(s, nd.Gradient(f)(x))]

        for i in I: 
            temp_cons.append(lambda s: np.dot(s, nd.Gradient(constraints[i])(x)))

        s = [-1, 0, 1]
        best_s = 0
        dicts = [dict() for i in range(len(temp_cons))]

        min_z = float('inf')
        best_s = tuple()

        for pair in list(product(s, s)):
            curr = []
            for ind in range(len(temp_cons)):
                curr.append(temp_cons[ind](pair))
                dicts[ind][pair] = temp_cons[ind](pair)
            if max(curr) < min_z:
                min_z = max(curr)
                best_s = pair


        if all(dicts[i][best_s] == 0 for i in range(len(dicts))): 
            flag = True
        max_a = 0
        new_x = []
        for a in np.linspace(0, Biggest_number_in_constr, Biggest_number_in_constr*12_000): 
            new_x = x + a*np.array(best_s)
            if any(constr(new_x) > 0 for constr in constraints):
                max_a = prev_a
                break 
            prev_a = a
        if np.abs(max_a) < 1e-4: flag = True 

        side_f = lambda a: f(x + a*np.array(best_s))
        min_a = one_dimentional.golden_section(side_f, 0, max_a)[0]
        x = x + min_a*np.array(best_s)

    return f(x), x

# Biggest_number_in_constr = 25
# f = lambda x: x[0] + x[1]
# g1 = lambda x: x[0]**2 + x[1]**2 - 25
# g2 = lambda x: x[1] - 3
# constraints = [g1, g2]

# f = lambda x: x[0] - x[1]
# g = lambda x: x[0]**2 + x[1]**2 - 1

def calculate_P(M: np.array, x):
    dim = max(M.shape) 
    if dim == 0: 
        return np.eye(x.shape[0])
    rev_m = np.linalg.inv(np.matmul(M, M.transpose()))
    return np.eye(dim) - M.transpose().dot(rev_m).dot(M)


def preparation_step(x, A, A_coeffs):
    A1 = np.array([])
    A2 = np.array([])
    b1 = np.array([])
    b2 = np.array([])

    for i in range(len(A)):
        if A[i](x) == 0: 
            if A1.shape[0] == 0: A1 = np.array([A_coeffs[i][:-1]])
            else:  A1 = np.append(A1, np.array([A_coeffs[i][:-1]]), axis=0)
            if b1.shape[0] == 0: b1 = np.array([A_coeffs[i][-1]])
            else: b1 = np.append(b1, np.array([A_coeffs[i][-1]]), axis=0)
        else:
            if A2.shape[0] == 0: A2 = np.array([A_coeffs[i][:-1]])
            else: A2 = np.append(A2, np.array([A_coeffs[i][:-1]]), axis=0)
            if b2.shape[0] == 0: b2 = np.array([A_coeffs[i][-1]])
            else: b2 = np.append(b2, np.array([A_coeffs[i][-1]]), axis=0)

    return A1, A2, b1, b2


def main_step(f, x, A, A1, A2, H, b1, b2):
    M = A1
    if H.shape[0] != 0:
        M = np.append(M, H, axis=0)
    grad = nd.Gradient(f)(x)
    P = calculate_P(M, x)
    s = - P.dot(grad)
    
    if all(abs(el) < 1e-07 for el in s) and M.shape[0] == 0:
        return True
    elif all(abs(el) < 1e-07 for el in s):
        w = - np.linalg.inv(M.dot(M.transpose())).dot(M).dot(grad)
        if all(w[i] >= 0 for i in range(A1.shape[0])):
            return True
        else: 
            d = {i:w[i] for i in range(A1.shape[0]) if w[i] < 0}
            mn_ind = min(d, key=d.get)
            A1 = np.delete(A1, (mn_ind), axis=0)
            return main_step(f, x, A, A1, A2, H, b1, b2)

    else:
        s_tilda = np.array([A2.dot(s)])
        b_tilda = np.array([b2 - A2.dot(x)])
        if all(el > 0 for el in s_tilda[0]): 
            alpha_max = min([b_tilda[0][i]/s_tilda[0][i] for i in range(b_tilda.shape[1])])
        else: 
            alpha_max = 1_000
            for a in np.linspace(0, alpha_max, alpha_max): 
                new_x = x + a*s
                if any(constr(new_x) > 0 for constr in A):
                    break 
                prev_a = a
            alpha_max = prev_a

        side_f = lambda alpha: f(x + alpha*s)
        optimal_alpha = one_dimentional.golden_section(side_f, 0, alpha_max)[0]
        new_x = x + optimal_alpha*s
        return  new_x



def Rosen(A, A_coeffs, H, f, initial_x):
    x = initial_x
    flag = False
    while(not flag):
        A1, A2, b1, b2 = preparation_step(x, A, A_coeffs)
        flagORx = main_step(f, x, A, A1, A2, H, b1, b2)
        if isinstance(flagORx, bool): flag = flagORx
        else: x = flagORx 
    return x, f(x)
