def dichotomy(f, a, b, e=.00001):

    while b-a < 2*e:
        x = (a+b)/2
        x1 = x-e
        x2 = x+e

        if (f(x1) < f(x2)):
            a = x1
        else:
            b = x2

    return f((a+b)/2)


def fibonacci(f, a, b, n = 1000):
    F = [0, 1, 1]
    
    for _ in range(n-2):
        F.append(F[-1] + F[-2])

    L = b-a
    x1 = F[n-1]/F[n]*L
    x2 = b - F[n-1]/F[n]*L

    f1 = f(x1)
    f2 = f(x2)

    while n > 2:
        if f2 < f1:
            b = x1
            f1 = f2
            x1 = x2
            L = b-a

            x2 = b - F[n-2]/F[n-1]*L
            f2 = f(x2)
        else:
            a = x2
            f2 = f1
            x2 = x1
            L = b-a
            x1 = a + F[n-2]/F[n-1]*L
            f1 = f(x1)
        n-=1 

    return min(f1, f2)

