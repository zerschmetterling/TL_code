from sympy import *
a1, a2, b1, b2, q = symbols('a1, a2, b1, b2, q', real=True)

def norm(x):
    return expand(x*conjugate(x))

def antiConj(x,y):
    answer = x*conjugate(y) + y*conjugate(x)
    return expand(answer)

def eq(x,y):
    answer = q*norm(x) + antiConj(x,y)
    return answer

a = a1 + I*a2
b = b1 + I*b2

solutions = nonlinsolve( [eq(a,b), eq(b,a)], [a1,a2,b1,b2] )
