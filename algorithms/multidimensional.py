import scipy.optimize
import numpy as np
from sympy import Symbol, simplify
import sympy
import numdifftools as nd 

import one_dimentional 
print('----------------------------------------------------------------')


def oneDim_scipy_minimization(f):
    return scipy.optimize.minimize(f, 5)['x'][0]



def Hook_Jeeves(f, symb_f, deltaX, currX: np.array):
    grad = [1]
    eps = 1e-02
    initialX = [comp for comp in currX ]
    while not all(abs(component) < eps for component in grad):
        for i in range(len(currX)):
            L = [x for x in currX]
            R = [x for x in currX]
            L[i] -= deltaX[i]
            R[i] +=  deltaX[i]
            d = {point:f(point) for point in {tuple(currX), tuple(R), tuple(L)}} 
            fmin = list(min(d, key=d.get))
            currX = fmin
        if initialX != currX:
            S = np.array(currX) - np.array(initialX)
            X2 = currX + Symbol('l')*S
            # print(f'S: {S}, currx: {currX}')            
            simplfd = simplify(symb_f(X2))
            final_f = lambda x: simplfd.subs({Symbol('l'): x[0]}) if not str(x)[-1].isdigit() else simplfd.subs({Symbol('l'): x})

            values = dict()
            
            for l in {2, 5, 10, 15, 100, 1_000, 10_000}: # np.linspace(2, 10_000, num=40): #
                pair = one_dimentional.golden_section(final_f, -l, l)
                values[pair[0]] = pair[1]
            minX = min(values, key=values.get)
                # minX = one_dimentional.golden_section(final_f, -minL, minL)[0]
            # gold = one_dimentional.golden_section(final_f, -L, L)[0]
            
            # scpy = one_dim_minimization(final_f)
            # if round(scpy) != round(minX): 
            #     print(round(scpy), final_f(round(scpy)))
            #     print( round(minX), final_f(round(minX)))
            #     print("ERROR_CASE!", simplfd) 
            #     print(scpy, minX)

            xNew = [float(round((X2[i].subs({Symbol('l'): minX})), 3)) for i in range(len(X2))]
            grad = nd.Gradient(f)(xNew)

            currX = xNew
            initialX = [comp for comp in currX]
        else: 
            deltaX /= 2

    return xNew # optimal value
    
symb_f = lambda x: 20 - (x[0] - 3)*sympy.exp(-(x[0] - 3)) - (x[1] - 2)*sympy.exp(-(x[1] - 2))

f = lambda x: 20 - (x[0] - 3)*np.exp(-(x[0] - 3)) - (x[1] - 2)*np.exp(-(x[1] - 2)) #(x[0]**2 + x[1] - 11)**2 + (x[0] + x[1]**2 - 7)**2

# print(Hook_Jeeves(f, symb_f, np.array([0.8, 0.8]), np.array([0, 0])))

print(f([6 - 1/(2*np.sqrt(5)), 3 - 3/(2 * np.sqrt(5))]), f([6, 3]))

# while(<epsX and < epsY):
#     norma = np.linalg.norm([0, -2])
# X2 = [-2, -4.3] + Symbol('l')*np.array([1, 0])

# print(X2)
# print(one_dim_minimization(X2, f1) )
# minX = one_dim_minimization(X2, f1)
# xNew = [float(round((X2[i].subs({Symbol('l'): minX})), 3)) for i in range(len(X2))]
# print(xNew)
# print(f1(xNew))

# print(np.linalg.norm([0, -2]))

def GS_orthog():
    ...
# построить метод принимающий набор векторов и преобразующий по Г-Ш
