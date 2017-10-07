variables = var('a1,a2,b1,b2')
var('q')
for v in variables:
    assume(v, 'real')

assume(q, 'real')
assume(q != 0)

def eqFun(a,b,c,d):
    return simplify(q*(a**2 + b**2) + 2*(a*c + b*d))

solutions = solve( [eqFun(a1,a2,b1,b2), eqFun(b1,b2,a1,a2)], variables )

for s in solutions:
    for x in s:
        print [x.lhs(), simplify(expand( x.rhs() ))]
    print


