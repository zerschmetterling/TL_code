a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,qa1,a2,b1,b2,c1,c2,d1,d2,e1,e2,q

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
    q*(norm(b) + norm(d)) + antiConj(a,b) + antiConj(d,b)
    )
eq.append(
    q*(norm(c) + norm(e)) + antiConj(a,c) + antiConj(e,c)
    )
eq.append(
    q*(e*conj(b) + c*conj(d)) + e*(conj(a)+conj(d)) + c*conj(b) + a*conj(d)
    )
# eq.append(
#    q*(b*conj(e) + d*conj(c)) + d*conj(a+e) + b*conj(c) + a*conj(e)
#    )



eq.append(
    q*(c*conj(a) + e*conj(b)) + conj(a)*(d+e) + c*conj(b) + e*conj(d)
    )
eq.append(
    q*(norm(a) + norm(b)) + antiConj(a,b) + antiConj(d,b)
    )
# eq.append(
#    q*(a*conj(c) + b*conj(e)) + conj(e)*(d+a) + b*conj(d) + a*conj(d)
#    )
eq.append(
    q*(norm(a) + norm(e)) + antiConj(c,d) + antiConj(c,e)
    )



eq.append(
    q*(norm(a) + norm(b)) + antiConj(a,b) + antiConj(e,b)
    )
eq.append(
    q*(norm(c) + norm(d)) + antiConj(d,c) + antiConj(e,c)
    )
eq.append(
    q*(c*conj(a) + d*conj(b)) + d*(conj(a)+conj(e)) + e*conj(a) + c*conj(b)
    )
#eq.append(
#    q*(a*conj(c) + b*conj(d)) + a*conj(d+e) + a*conj(e) + b*conj(c)
#    )




eq.append(
    q*(norm(c) + norm(d)) + antiConj(a,c) + antiConj(d,c)
    )
eq.append(
    q*(norm(b) + norm(e)) + antiConj(a,b) + antiConj(e,b)
    )
eq.append(
    q*(e*conj(c) + b*conj(d)) + (a+e)*conj(d) + b*conj(c) + e*conj(a)
    )
#eq.append(
#    q*(c*conj(e) + d*conj(b)) + d*conj(a+e) + c*conj(b) + a*conj(e)
#    )



eq.append(
    q*(b*conj(a) + e*conj(c)) + (d+e)*conj(a) + b*conj(c) + e*conj(d)
    )
eq.append(
    q*(norm(a) + norm(c)) + antiConj(a,c) + antiConj(d,c)
    )
#eq.append(
#    q*(a*conj(b) + c*conj(e)) + a*conj(d+e) + c*conj(b) + d*conj(e)
#    )
eq.append(
    q*(norm(b) + norm(e)) + antiConj(d,b) + antiConj(e,b)
    )




eq.append(
    q*(norm(a) + norm(c)) + antiConj(a,c) + antiConj(e,c)
    )
eq.append(
    q*(norm(b) + norm(d)) + antiConj(b,d) + antiConj(b,e)
    )
#eq.append(
#    q*(b*conj(a) + d*conj(c)) + d*conj(a+e) + e*conj(a) + b*conj(c)
#    )
eq.append(
    q*(a*conj(b) + c*conj(d)) + (a+e)*conj(d) + a*conj(e) + c*conj(b)
    )

eqSeparated = []
for equa in eq:
    eqSeparated.append(real(equa))
    eqSeparated.append(imag(equa))



