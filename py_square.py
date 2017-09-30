from scipy import *

def q(n):
    return exp(2*pi*i/n)

def d(n):
    return q(n)**10+q(n)**8+q(n)**2+1+q(n)**(-2)+q(n)**(-8)+q(n)**(-10)

def t(n):
    return (-1)*(q(n)**2-1+q(n)**(-2))/(q(n)**4+q(n)**(-4))

def denom(n):
    return d(n)*t(n)+d(n)+t(n)

def a(n):
    return ((d(n)*(t(n)**2))+(t(n)**2)-1)/denom(n)

def b(n):
    return (t(n)-(t(n)**2)-1)/denom(n)

