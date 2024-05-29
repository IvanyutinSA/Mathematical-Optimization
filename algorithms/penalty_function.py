import numpy as np
from scipy.optimize import minimize
import sys 
# if "../Mathematical-Optimization" not in sys.path:
#     sys.path.append("../Mathematical-Optimization")
# from algorithms.one_dimentional import golden_section
from algorithms.multi_dimentional import gradient_descent, gauss_seidel

def penalty_function_method(x0, f, r_initial, r_increment, ineq_constraints, eq_constraints, eps=1e-7):
    r = r_initial
    x = x0

    penalty_f = lambda x: r * (sum(max(0, constr(x))**2 for constr in ineq_constraints) +  sum(constr(x)**2 for constr in eq_constraints))
    side_f = lambda x: f(x) + penalty_f(x)
    while True:
        result = gauss_seidel(side_f, x)
        x_new = result['solution']
    # Проверяем выполнение ограничений
        if all(con(x_new) <= eps for con in eq_constraints) and all(con(x_new) <= eps for con in ineq_constraints):
            break

        # Увеличиваем параметр штрафа
        r += r_increment
        x = x_new

    return x_new
