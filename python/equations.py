from sympy import *
a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,q = symbols('a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,q', real=True)

a = a1 + I*a2
b = b1 + I*b2
c = c1 + I*c2
d = d1 + I*d2
e = e1 + I*e2

variables = [a1, a2, b1, b2, c1, c2, d1, d2, e1, e2, q]

def antiConj(x,y):
    answer = conjugate(x)*y + x*conjugate(y)
    return expand(answer)

def conj(x):
    return conjugate(x)

def norm(x):
    answer = x*conjugate(x)
    return answer

eq = []

eq.append(simplify(expand(
    q*(norm(b) + norm(d)) + antiConj(a,b) + antiConj(d,b)
    )))
eq.append(simplify(expand(
    q*(norm(c) + norm(e)) + antiConj(a,c) + antiConj(e,c)
    )))
eq.append(simplify(expand(
    q*(e*conj(b) + c*conj(d)) + e*(conj(a)+conj(d)) + c*conj(b) + a*conj(d)
    )))


eq.append(simplify(expand(
    q*(c*conj(a) + e*conj(b)) + conj(a)*(d+e) + c*conj(b) + e*conj(d)
    )))
eq.append(simplify(expand(
    q*(norm(a) + norm(b)) + antiConj(a,b) + antiConj(d,b)
    )))
eq.append(simplify(expand(
    q*(norm(a) + norm(e)) + antiConj(c,d) + antiConj(c,e)
    )))



eq.append(simplify(expand(
    q*(norm(a) + norm(b)) + antiConj(a,b) + antiConj(e,b)
    )))
eq.append(simplify(expand(
    q*(norm(c) + norm(d)) + antiConj(d,c) + antiConj(e,c)
    )))
eq.append(simplify(expand(
    q*(c*conj(a) + d*conj(b)) + d*(conj(a)+conj(e)) + e*conj(a) + c*conj(b)
    )))


eq.append(simplify(expand(
    q*(norm(c) + norm(d)) + antiConj(a,c) + antiConj(d,c)
    )))
eq.append(simplify(expand(
    q*(norm(b) + norm(e)) + antiConj(a,b) + antiConj(e,b)
    )))
eq.append(simplify(expand(
    q*(e*conj(c) + b*conj(d)) + (a+e)*conj(d) + b*conj(c) + e*conj(a)
    )))

eq.append(simplify(expand(
    q*(b*conj(a) + e*conj(c)) + (d+e)*conj(a) + b*conj(c) + e*conj(d)
    )))
eq.append(simplify(expand(
    q*(norm(a) + norm(c)) + antiConj(a,c) + antiConj(d,c)
    )))
eq.append(simplify(expand(
    q*(norm(b) + norm(e)) + antiConj(d,b) + antiConj(e,b)
    )))


eq.append(simplify(expand(
    q*(norm(a) + norm(c)) + antiConj(a,c) + antiConj(e,c)
    )))
eq.append(simplify(expand(
    q*(norm(b) + norm(d)) + antiConj(b,d) + antiConj(b,e)
    )))
eq.append(simplify(expand(
    q*(a*conj(b) + c*conj(d)) + (a+e)*conj(d) + a*conj(e) + c*conj(b)
    )))

eqSeparated = []

for i in range(len(eq)):
    real_part = simplify( (eq[i] + conj(eq[i]))/2 )
    imag_part = simplify( (eq[i] - conj(eq[i]))/(2*I) )
    for part in [real_part, imag_part]:
        if part:
            eqSeparated.append(part)

def eqSeparatedFun(variables):
    return eqSeparated
