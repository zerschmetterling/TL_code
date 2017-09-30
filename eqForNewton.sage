variables = var('a,b,c,d,e')
var('q')
for v in variables:
    assume(v, 'complex')
assume(q, 'real')

def conj(x):
    return conjugate(x)

def antiConj(x,y):
    answer = conj(x)*y + x*conj(y)
    return answer

eq = []

eq.append(
    norm(a) -norm(c)
    )
eq.append(
    norm(a) -norm(d)
    )
eq.append(
    norm(a) -norm(e)
    )
eq.append(
    antiConj(a,b) -antiConj(d,b)
    )
eq.append(
    antiConj(a,b) -antiConj(e,b)
    )
eq.append(
    antiConj(a,c) -antiConj(d,c)
    )
eq.append(
    antiConj(a,c) -antiConj(e,c)
    )
eq.append(
    q*(norm(a) + norm(b)) + 2*antiConj(a,b) -0
    )
eq.append(
    q*norm(a) + antiConj(a,c) -0
    )

eq.append(
    q*(e*conj(b) + c*conj(d)) + e*conj(a+d) + c*conj(b) + a*conj(d) -0
    )
eq.append(
    q*(c*conj(a) + e*conj(b)) + conj(a)*(d+e) + c*conj(b) + a*conj(d) -0
    )

eq.append(
    q*(c*conj(a) + d*conj(b)) + d*conj(a+e) + e*conj(a) + c*conj(b) -0
    )
eq.append(
    q*(e*conj(c) + b*conj(d)) + conj(d)*(a+e) + e*conj(a) + b*conj(c) -0
    )

eq.append(
    q*(b*conj(a) + e*conj(c)) + conj(a)*(d+e) + b*conj(c) + e*conj(d) -0
    )
eq.append(
    q*(b*conj(a) + d*conj(c)) + conj(d)*(a+e) + e*conj(a) + c*conj(b) -0
    )

solution = solve( eq, variables, to_poly_solve=True)
