varList = list(var("a1, a2, b1, b2, c1, c2, d1, d2, e1, e2"))
load("perfectness.sage")

R = initR(varList)
T = initT(3, q, R)
bas = initBasis(T)

coeff = [c2*iSymbol + c1, a2*iSymbol + a1, d2*iSymbol + d1, e2*iSymbol + e1, b2*iSymbol + b1]

rot = allRotations(coeff)
mult = allMultiplications(rot)
symEq = getSymbolicEquations(mult)


dictionary_q_2 = {      
    q : 2,          
    a1 : -1,    
    b1 : -1,     
    c1 : -1,     
    d1 : 1,
    e1 : 1,              
    a2:0,
    b2:0,
    c2:0, 
    d2:0,
    e2:0
}                                                            

var("x")
assume(x > 2)
assume(x, "integer")

dictionary_q_cos = {
    q  : 2*cos(pi/x),
    a1 : cos(pi*(x-1)/x),
    a2 : -sin(pi*(x-1)/x),
    b1 : cos(pi*(x-1)/x),
    b2 : sin(pi*(x-1)/x),
    c1 : cos(pi*(x-1)/x),
    c2 : sin(pi*(x-1)/x),
    d1 : 2 + 2*cos(pi/x)*cos(pi*(x-1)/x) + cos(2*pi*(x-1)/x),
    d2 : 2*cos(pi/x)*sin(pi*(x-1)/x) + sin(2*pi*(x-1)/x),
    e1 : 1, e2 : 0
}

symEq_q_2 = list(set([bla.subs(dictionary_q_2) for bla in symEq]))

symEq_q_cos = [bla.subs(dictionary_q_cos) for bla in symEq]
symEq_q_cos = [bla.simplify_trig() for bla in symEq_q_cos]
symEq_q_cos = list(set(symEq_q_cos))

print "If q = 2, then the equations reduce to "
print symEq_q_2, "\n"
print "If otherwise q is in the discrete spectrum, i.e. q = 2*cos(pi/n), then we obtain"
print symEq_q_cos

