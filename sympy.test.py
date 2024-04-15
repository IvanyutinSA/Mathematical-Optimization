import sympy as sp
from algorithms.conditional import add_lagrange_multipliers


def unique_symbols(f):
    def wrapper(*args):
        wrapper.number += 1
        return f(*args, number=wrapper.number)

    wrapper.number = 0
    return wrapper

@unique_symbols
def union_functions(f, g, vars, number=0):
    z = sp.symbols(f'union_functions_{number}')
    h = f + z * g
    vars.append(z)
    return h

def main():
    x, y = sp.symbols('x y')
    f = x + 1
    g = y + 1
    vars = [x, y]
    h = add_lagrange_multipliers(f, [], vars)
    print(h)
    print(vars)
    

def second_test():
    x, y, u, v = sp.symbols('x y u v')
    f = x
    g = y
    h = x + y
    vars = [x, y]

    result_function = union_functions(union_functions(f, g, vars), h, vars)
    actual_function = f + u*g + v*h

    print(f'result_functioin: {result_function}')
    print(f'result_function: {result_function.evalf(subs={p: q for p, q in zip(vars, (1, 2, 3, 4))})}')
    print(f'actual_function: {actual_function}')
    print(f'actual_function: {actual_function.evalf(subs={p: q for p, q in zip((x, y, u, v), (1, 2, 3, 4))})}')


def first_test():
    x, y = sp.symbols('x y')
    f = x + 2*y
    g = f ** 2

    vars = [x, y]
    z = sp.symbols('z')
    h = f + z * g

    print(f'h: {h.evalf(subs={p: q for p, q in zip(vars+[z], (1, 2 ,1))})}')
    print(f'f: {f.evalf(subs={p: q for p, q in zip(vars+[z], (1, 2, 1))})}')
    print(f'g: {g.evalf(subs={p: q for p, q in zip(vars, (1, 2))})}')

    # print(f'h: {h.evalf(subs={p, q for p,q in zip(vars)})}')
    # print(f'f: {f.evalf(subs=(zip(vars), (1, 2)))}')
    # print(f'g: {g.evalf(subs=(zip(vars), (1, 2)))}')


if __name__ == "__main__":
    main()

