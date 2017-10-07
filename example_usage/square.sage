import scipy
import numpy as np
from scipy import optimize

var('x')
assume(x, 'real')

q[x_]        := Exp[2*Pi*I/x]
d[x_]        := Re[ q[x]^10+q[x]^8+q[x]^2+1+q[x]^(-2)+q[x]^(-8)+q[x]^(-10) ]
t[x_]       := Re[ (-1)*(q[x]^2-1+q[x]^(-2))/(q[x]^4+q[x]^(-4))]
denom[x_]    := d[x]*t[x]+d[x]+t[x]
a[x_]        := ((d[x]*(t[x]^2))+(t[x]^2)-1)/denom[x]
b[x_]        := (t[x]-(t[x]^2)-1)/denom[x]
eq1[x_]      := a[x] * (a[x] + 2*a[x]*t[x] + a[x]^2 + 2*b[x])
eq2[x_]      := a[x] * (a[x]^2 + 2*b[x])
eq3[x_]      := b[x] * (a[x]^2 + b[x]*(2 + d[x]) + 2*a[x])

#assume(x,"ixteger")
#assume(x > 0)

def fun(variable):
    answer= [
        n(eq1(variable[0])),
        n(eq2(variable[0])),
        n(eq3(variable[0]))]
    return answer

def funTest(variable):
    return [n(eq1(variable[0]))]

