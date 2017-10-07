S = [ bas.index(x) for x in bas if x == adjoint(x) ]
M = [ bas.index(x) for x in bas if x not in S ]

def commutator(x,y):
    return x*y - y*x

