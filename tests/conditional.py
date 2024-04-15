from unittest import TestCase
import unittest
import sympy as sp
from algorithms.conditional import add_lagrange_multipliers

class LagrangeTests(TestCase):
    def test_add_lagrange_multipliers_no_constraints(self):
        x, y = sp.symbols('x y')
        vars = [x, y]
        constraints = []
        f = x + y

        h = add_lagrange_multipliers(f, constraints, vars)

        self.assertEqual([x, y], vars)
        self.assertEqual(f, h)

    def test_add_lagrange_multipliers_one_contraint(self):
        x, y = sp.symbols('x y')
        vars = [x, y]
        g = y
        constraints = [g]
        f = x

        h = add_lagrange_multipliers(f, constraints, vars)
        l = [sp.symbols(f'add_lagrange_multipliers{k}') for k in range(1, len(constraints)+1)]

        self.assertEqual([x, y] + l, vars)
        self.assertEqual(f + l[0]*g, h)

    def test_add_lagrange_multipliers_multiple_contraints(self):
        x, y = sp.symbols('x y')
        vars = [x, y]
        constraints = [x, y]
        f = x

        h = add_lagrange_multipliers(f, constraints, vars)
        l = [sp.symbols(f'add_lagrange_multipliers{k}') for k in range(1, len(constraints)+1)]

        self.assertEqual([x, y] + l, vars)
        self.assertEqual(f + l[0]*constraints[0] + l[1]*constraints[1], h)


if __name__ == '__main__':
    unittest.main()

