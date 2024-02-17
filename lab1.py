import algorithms.one_dimentional as od

f = lambda x: x**2 + 5;

print(f'fibonacci: {od.fibonacci(f, -10, 10)}')
print(f'dichotomy: {od.dichotomy(f, -10, 10)}')

print(f'golden section method: {od.golden_section(f, -10, 10)}')
