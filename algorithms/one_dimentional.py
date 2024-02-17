def dichotomy(f, a, b, e=.0000001):

    while round(b-a, 7) > 2*e:
        x = (a+b)/2
        x1 = x-e
        x2 = x+e

        if (f(x1) < f(x2)):
            b = x2
        else:
            a = x1

    return f((a+b)/2)

def fibonacci(f, a, b, e=.0000001):
    F = [0, 1, 1]

    while F[-1] < (b-a)/e:
        F.append(F[-1] + F[-2])

    n = len(F) - 1
    f1, f2 = float('inf'), float('inf')

    while n > 2:
        L = b-a 
        x1 = a + F[n-1]/F[n]*L
        x2 = b - F[n-1]/F[n]*L
        f1, f2 = f(x1), f(x2)
        if f2 < f1:
            b = x1
            f1 = f2
            x1 = x2
            L = b-a
            x2 = a+(b-x1)
            f2 = f(x2)
        else:
            a = x2
            f2 = f1
            x2 = x1
            L = b-a
            x1 = b-(x2-a)
            f1 = f(x1)
        n -= 1

    return min(f1, f2)

from math import sqrt 

def golden_section(f, a, b, eps=.0000001):
    constF = (1 + sqrt(5))/2
    while (b-a >= eps):
        x1 = b - (b - a)/constF
        x2 = a + (b - a)/constF
        if f(x1) >= f(x2):
            a = x1
        else: b = x2
    x = (a + b)/2 
    return f(x)   

