a1, a2, b1, b2, q = RR['a1, a2, b1, b2, q'].gens()

variables = var('a1,a2,b1,b2')
var('q')
for v in variables:
    assume(v, 'real')

assume(q != 0)

a = a1 + i*a2
aStr = a1 - i*a2
b = b1 + i*b2
bStr = b1 - i*b2

def conj(x):
    if x == a:
        return aStr
    elif x == b:
        return bStr
    else:
        return x

def antiConj(x,y):
    answer = x*conj(y) + conj(x)*y
    return expand(answer)

def eq(x,y):
    answer = q*NORM(x) + antiConj(x,y)
    return expand(answer)

def NORM(x):
    answer = x*conj(x)
    return expand(answer)

solutions = solve( [eq(a,b)==0, eq(b,a)==0 ] , variables )

for sol in solutions:
    print(sol)
    print


