#S.<q> = QQ[]
#S = S.fraction_field()

qValues = [2*cos(pi/n) for n in range(3,10)]

P.<a1,a2,b1,b2,c1,c2,d1,d2,e1,e2> = PolynomialRing( QQbar, 10, order='degrevlex' )
#P.<a1,b1,c1,d1,e1> = PolynomialRing( RR, 5, order='degrevlex' )
#var('a2,b2,c2,d2,e2,q')
variables = [a1,a2,b1,b2,c1,c2,d1,d2,e1,e2]

var('q')

######## NOW WE INSERT THE EQUATIONS FROM THE OUTPUT OF SAGE, SEPARATING REAL AND IMAGINARY PARTS (only those where both occur) #############
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

def I_for_q(qV):
    I_cache = ideal( [P(equa.subs(q=qV)) for equa in eqSeparated] )
    return I_cache
#B = I.groebner_basis()

I = []

for q in qValues:
    I.append( I_for_q(q) )
    
for i in I:
    i = ideal( i.groebner_basis() )

dimensions = []

for i in I:
    dimensions.append(i.dimension())

for n in range(len(qValues)):
    print "q =", qValues[n]
    print "dimension of ideal", dimensions[n]


