import sympy as sp
import numpy as np

def add_lagrange_multipliers(f, constraints, vars):
    n = len(vars)
    vars.extend([sp.symbols(f"add_lagrange_multipliers{k}") for k in range(1, len(constraints)+1)])
    h = f + sum(vars[n+i] * constraint for i, constraint in enumerate(constraints))
    return h

def lagrange_multipliers(f, constraints, vars):
    h = add_lagrange_multipliers(f, constraints, vars)
    h_diffs = [sp.diff(h, var) for var in vars] 
    solutions = sp.solve(h_diffs, vars, domain=sp.Reals)

    if not solutions:
        return "no stationary points"
    print(solutions)
    try:
        return min(h.evalf(subs={p: q for p, q in zip(vars, vector)}) for vector in solutions)
    except:
        return h.evalf(subs=solutions)

