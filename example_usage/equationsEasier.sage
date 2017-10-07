variables = var('a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,q')
for v in variables:
    assume(v, 'real')


a = a1 + i*a2
b = b1 + i*b2
c = c1 + i*c2
d = d1 + i*d2
e = e1 + i*e2

aStr = a1 - i*a2
bStr = b1 - i*b2
cStr = c1 - i*c2
dStr = d1 - i*d2
eStr = e1 - i*e2

def conj(x):
    if x == a:
        return aStr
    if x == b:
        return bStr
    if x == c:
        return cStr
    if x == d:
        return dStr
    if x == e:
        return eStr
    

def antiConj(x,y):
    answer = conj(x)*y + x*conj(y)
    return expand(answer)
eq = []

eq.append(
    norm(a) == norm(c)
    )
eq.append(
    norm(a) == norm(d)
    )
eq.append(
    norm(a) == norm(e)
    )
eq.append(
    antiConj(a,b) == antiConj(d,b)
    )
eq.append(
    antiConj(a,b) == antiConj(e,b)
    )
eq.append(
    antiConj(a,c) == antiConj(d,c)
    )
eq.append(
    antiConj(a,c) == antiConj(e,c)
    )
eq.append(
    q*(norm(a) + norm(b)) + 2*antiConj(a,b) == 0
    )
eq.append(
    q*norm(a) + antiConj(a,c) == 0
    )

eq.append(
    q*(e*conj(b) + c*conj(d)) + e*(conj(a)+conj(d)) + c*conj(b) + a*conj(d) == 0
    )
eq.append(
    q*(c*conj(a) + e*conj(b)) + conj(a)*(d+e) + c*conj(b) + a*conj(d) == 0
    )

eq.append(
    q*(c*conj(a) + d*conj(b)) + d*(conj(a)+conj(e)) + e*conj(a) + c*conj(b) == 0
    )
eq.append(
    q*(e*conj(c) + b*conj(d)) + conj(d)*(a+e) + e*conj(a) + b*conj(c) == 0
    )

eq.append(
    q*(b*conj(a) + e*conj(c)) + conj(a)*(d+e) + b*conj(c) + e*conj(d) == 0
    )
eq.append(
    q*(b*conj(a) + d*conj(c)) + conj(d)*(a+e) + e*conj(a) + c*conj(b) == 0
    )


#outWithoutPoly = file('./TL3_no_poly_solutions.txt', 'w')
#
#for equa in eq:
#    outWithoutPoly.write(str(equa)+'\n')
#
#outWithoutPoly.close()

#sol = solve( eq, variables )
#polySolv = solve( eq, variables, to_poly_solve=True)



eqSeparated = [
    b1^2*q + b2^2*q + d1^2*q + d2^2*q + 2*a1*b1 + 2*a2*b2 + 2*b1*d1 + 2*b2*d2,
    c1^2*q + c2^2*q + e1^2*q + e2^2*q + 2*a1*c1 + 2*a2*c2 + 2*c1*e1 + 2*c2*e2,
    a1^2*q + a2^2*q + b1^2*q + b2^2*q + 2*a1*b1 + 2*a2*b2 + 2*b1*d1 + 2*b2*d2,
    a1^2*q + a2^2*q + e1^2*q + e2^2*q + 2*c1*d1 + 2*c2*d2 + 2*c1*e1 + 2*c2*e2,
    a1^2*q + a2^2*q + b1^2*q + b2^2*q + 2*a1*b1 + 2*a2*b2 + 2*b1*e1 + 2*b2*e2,
    c1^2*q + c2^2*q + d1^2*q + d2^2*q + 2*c1*d1 + 2*c2*d2 + 2*c1*e1 + 2*c2*e2,
    c1^2*q + c2^2*q + d1^2*q + d2^2*q + 2*a1*c1 + 2*a2*c2 + 2*c1*d1 + 2*c2*d2,
    b1^2*q + b2^2*q + e1^2*q + e2^2*q + 2*a1*b1 + 2*a2*b2 + 2*b1*e1 + 2*b2*e2,
    a1^2*q + a2^2*q + c1^2*q + c2^2*q + 2*a1*c1 + 2*a2*c2 + 2*c1*d1 + 2*c2*d2,
    b1^2*q + b2^2*q + e1^2*q + e2^2*q + 2*b1*d1 + 2*b2*d2 + 2*b1*e1 + 2*b2*e2,
    a1^2*q + a2^2*q + c1^2*q + c2^2*q + 2*a1*c1 + 2*a2*c2 + 2*c1*e1 + 2*c2*e2,
    b1^2*q + b2^2*q + d1^2*q + d2^2*q + 2*b1*d1 + 2*b2*d2 + 2*b1*e1 + 2*b2*e2,
#
#    eq[0], eq[1], eq[4], eq[5], eq[6], eq[7], eq[9], 
#    eq[10], eq[13], eq[14], eq[15], eq[16],
#
    c1*d1*q + c2*d2*q + b1*e1*q + b2*e2*q + b1*c1 + b2*c2 + a1*d1 + a2*d2 + a1*e1 + d1*e1 + a2*e2 + d2*e2,
    c2*d1*q - c1*d2*q - b2*e1*q + b1*e2*q - b2*c1 + b1*c2 + a2*d1 - a1*d2 - a2*e1 - d2*e1 + a1*e2 + d1*e2,

    a1*c1*q + a2*c2*q + b1*e1*q + b2*e2*q + b1*c1 + b2*c2 + a1*d1 + a2*d2 + a1*e1 + d1*e1 + a2*e2 + d2*e2,
    - a2*c1*q + a1*c2*q - b2*e1*q + b1*e2*q - b2*c1 + b1*c2 - a2*d1 + a1*d2 - a2*e1 - d2*e1 + a1*e2 + d1*e2,

    a1*c1*q + a2*c2*q + b1*d1*q + b2*d2*q + b1*c1 + b2*c2 + a1*d1 + a2*d2 + a1*e1 + d1*e1 + a2*e2 + d2*e2,
    - a2*c1*q + a1*c2*q - b2*d1*q + b1*d2*q - b2*c1 + b1*c2 - a2*d1 + a1*d2 - a2*e1 + d2*e1 + a1*e2 - d1*e2,

    b1*d1*q + b2*d2*q + c1*e1*q + c2*e2*q + b1*c1 + b2*c2 + a1*d1 + a2*d2 + a1*e1 + d1*e1 + a2*e2 + d2*e2,
    b2*d1*q - b1*d2*q - c2*e1*q + c1*e2*q + b2*c1 - b1*c2 + a2*d1 - a1*d2 - a2*e1 - d2*e1 + a1*e2 + d1*e2,

    a1*b1*q + a2*b2*q + c1*e1*q + c2*e2*q + b1*c1 + b2*c2 + a1*d1 + a2*d2 + a1*e1 + d1*e1 + a2*e2 + d2*e2,
    - a2*b1*q + a1*b2*q - c2*e1*q + c1*e2*q + b2*c1 - b1*c2 - a2*d1 + a1*d2 - a2*e1 - d2*e1 + a1*e2 + d1*e2,

    a1*b1*q + a2*b2*q + c1*d1*q + c2*d2*q + b1*c1 + b2*c2 + a1*d1 + a2*d2 + a1*e1 + d1*e1 + a2*e2 + d2*e2,
    a2*b1*q - a1*b2*q + c2*d1*q - c1*d2*q - b2*c1 + b1*c2 + a2*d1 - a1*d2 + a2*e1 - d2*e1 - a1*e2 + d1*e2 
    ]


