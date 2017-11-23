varList = list(var("a, b, c, d, e"))
load("perfectness.sage")

R = initR(varList)
T = initT(3, q, R)
bas = initBasis(T)

coeff = [c, a, d, e, b]

rot = allRotations(coeff)
mult = allMultiplications(rot)
symEq = getSymbolicEquations(mult)


dictionary_q_2 = {      
    q : 2,          
    a : -1,    
    aCC : -1,    
    b : -1,     
    bCC : -1,     
    c : -1,     
    cCC : -1,     
    d : 1,
    dCC : 1,
    e : 1,              
    eCC : 1,              
}                                                            

var("x")
assume(x > 2)
assume(x, "integer")

dictionary_q_cos = {
    q  : 2*cos(pi/x),
    a : cos(pi*(x-1)/x) - I*sin(pi*(x-1)/x),
    aCC : cos(pi*(x-1)/x) + I*sin(pi*(x-1)/x),
    b : cos(pi*(x-1)/x) + I*sin(pi*(x-1)/x),
    bCC : cos(pi*(x-1)/x) - I*sin(pi*(x-1)/x),
    c : cos(pi*(x-1)/x) + I*sin(pi*(x-1)/x),
    cCC : cos(pi*(x-1)/x) - I*sin(pi*(x-1)/x),
    d : 2 + 2*cos(pi/x)*cos(pi*(x-1)/x) + cos(2*pi*(x-1)/x) + I*(2*cos(pi/x)*sin(pi*(x-1)/x) + sin(2*pi*(x-1)/x)),
    dCC : 2 + 2*cos(pi/x)*cos(pi*(x-1)/x) + cos(2*pi*(x-1)/x) - I*(2*cos(pi/x)*sin(pi*(x-1)/x) + sin(2*pi*(x-1)/x)),
    e : 1,
    eCC : 1
}

symEq_q_2 = list(set([bla.subs(dictionary_q_2) for bla in symEq]))

symEq_q_cos = [bla.subs(dictionary_q_cos) for bla in symEq]
symEq_q_cos = [bla.simplify_trig() for bla in symEq_q_cos]
symEq_q_cos = list(set(symEq_q_cos))

print "If q = 2, then the equations reduce to "
print symEq_q_2, "\n"
print "If otherwise q is in the discrete spectrum, i.e. q = 2*cos(pi/n), then we obtain"
print symEq_q_cos

