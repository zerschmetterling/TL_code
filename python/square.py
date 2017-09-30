from sympy import *

#n = symbols('n', real=True)
q = symbols('q', complex=True)

#q        = exp(2*pi*I/n)
d        = q**10+q**8+q**2+1+q**(-2)+q**(-8)+q**(-10)
t        = (-1)*(q**2-1+q**(-2))/(q**4+q**(-4))
denom    = d*t+d+t
a        = ((d*(t**2))+(t**2)-1)/denom
b        = (t-(t**2)-1)/denom
eq1      = simplify(a * (a + 2*a*t + a**2 + 2*b))
eq2      = simplify(a * (a**2 + 2*b))
eq3      = simplify(b * (a**2 + b*(2 + d) + 2*a))
