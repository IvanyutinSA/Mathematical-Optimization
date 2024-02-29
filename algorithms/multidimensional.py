import scipy.optimize
import numpy as np
from sympy import Symbol, simplify, diff
import numdifftools as nd 

print('----------------------------------------------------------------')


def alpha_minimizing(currX, S, f):
    X2 = currX + Symbol('l')*S
    # print(currX + Symbol('l')*S, len(X2))
    simplfd = simplify(f(X2))
    final_f = lambda x: simplfd.subs({Symbol('l'): x[0]}) if not str(x).isdigit() else simplfd.subs({Symbol('l'): x})
    min_x = scipy.optimize.minimize(final_f, 0)['x'][0]
    # print(min_x)
    X2value = [float(round((X2[i].subs({Symbol('l'): min_x})), 3)) for i in range(len(X2))]
    # print(X2value)
    return X2value


def Hook_Jeeves(f, deltaX, currX: np.array):
    grad = [1]
    eps = 1e-02
    initialX = [comp for comp in currX ]
    while not all(abs(component) < eps for component in grad):
        for i in range(len(currX)):
            L = [x for x in currX]
            R = [x for x in currX]
            # print(f'delta {deltaX[i]}')
            # print(L[i])
            L[i] -= deltaX[i]
            # print( L[i])
            R[i] +=  deltaX[i]
            # print(f'L: {L}, X: {currX}, R: {R}')
            d = {point:f(point) for point in {tuple(currX), tuple(R), tuple(L)}} 
            # print(d)
            fmin = list(min(d, key=d.get))
            # initialX = fmin
            currX = fmin
            # print(f'fmin: {fmin}')
        if initialX != currX:
            # print('if')
            S = np.array(currX) - np.array(initialX)
            # print(f'S: {S}, currX: {currX}, prevX: {initialX}')
            xNew = alpha_minimizing(currX, np.array(S), f)
            # print(type(X2))
            # a = [0, 0.5]
            # for i in range(len(a)):
            #     # print(a[i] == x2[i], a[i], x2[i])
            #     a[i] = float(x2[i])
            # print(a == X2, a[1] == X2[1])
            grad = nd.Gradient(f)(xNew)
            # print(f'grad: {grad}')
            currX = xNew
            initialX = [comp for comp in currX]
        else: 
            # print('else')
            deltaX /= 2
    return xNew # optimal value
    
f = lambda x: (x[0]**2 + x[1] - 11)**2 + (x[0] + x[1]**2 - 7)**2

print(Hook_Jeeves(f, np.array([0.8, 0.8]), np.array([0, 0])))

# (2x0+x1+1, x0+4x1-1)[-0.714965820312500, 0.42843017578125]

# from sympy import diff
# print(f.derivative())
