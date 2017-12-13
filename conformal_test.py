
from __future__ import print_function
from galgebra import ga
from sympy import *

g3metric = '1 0 0 0 0, 0 1 0 0 0, 0 0 1 0 0, 0 0 0 1 0, 0 0 0 0 -1'
g3c = ga.Ga('e_1 e_2 e_3 e_4 e_5',g=g3metric)

# Give consistent names to our multivectors
e1,e2,e3,e4,e5 = g3c.mv()

# Check our metric
print('METRIC')
print(e1*e1)
print(e2*e2)
print(e3*e3)
print(e4*e4)
print(e5*e5)

# Define our null vectors
ninf = e4 + e5
no = 0.5*(e4 - e5)
print('NULL VECTORS')
print('ninf ',ninf)
print('no ',no)
print('ninf*ninf ',ninf*ninf)
print('no*no ',no*no)
print('ninf*no ',ninf*no)

# Define an up mapping
def up(x):
    return x + 0.5*x*x*ninf - no

# Define a homo function as per clifford
def homo(x):
    y = (-x|ninf)
    adj = (~y)
    madj = adj*y
    return x*(adj/madj)

# Define a down mapping
E0 = ninf^(-no)
def down(x):
    return (homo(x)^E0)*E0

alpha, beta = symbols('alpha beta')

print('UP AND DOWN MAPPING')
print('up(alpha*e1 + beta*e2): ', up(alpha*e1 + beta*e2))
print('down(up(alpha*e1 + beta*e2)) ', down(up(alpha*e1 + beta*e2)))
