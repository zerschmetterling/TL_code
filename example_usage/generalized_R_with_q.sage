reset()

#P.<a1,a2,b1,b2> = QQbar[]
P.<a1,a2,b1,b2> = QQbar[]

vari = list(var('a1S,a2S,b1S, b2S'))
variP = [a1,a2,b1,b2]
qValues = [2*cos(pi/n) for n in range(3,10)]

def eqFun(a,b,c,d,q):
    return q*(a**2 + b**2) + 2*(a*c + b*d)

I = []
B = []
Solutions = []

for i in range(len(qValues)):
    I.append( 
        ideal( 
            P(eqFun(a1,a2,b1,b2, qValues[i])),
            P(eqFun(b1,b2,a1,a2, qValues[i]))
        )
    )

for i in I:
    B.append( i.groebner_basis() )

B_n = []

for b in B:
    blist = []
    for littleB in b:
        blist.append(
            SR( littleB ).subs(
                dict( zip( variP, vari) )
            )
        )
    B_n.append( blist )

B_new = [
    B_n[0],
    B_n[1],
    [a2S^2*b1S - 1/2*b1S^3 - a1S*a2S*b2S + 1/sqrt(2)*a2S*b1S*b2S - 1/sqrt(2)*a1S*b2S^2 - 1/2*b1S*b2S^2,
  a1S^2 + a2S^2 - b1S^2 - b2S^2,
  a1S*b1S + 1/sqrt(2)*b1S^2 + a2S*b2S + 1/sqrt(2)*b2S^2],
 [a2S^2*b1S - 1/8*(5-sqrt(5))*b1S^3 - a1S*a2S*b2S + 1/4*(1+sqrt(5))*a2S*b1S*b2S - 1/4*(1+sqrt(5))*a1S*b2S^2 - 1/8*(5-sqrt(5))*b1S*b2S^2,
  a1S^2 + a2S^2 - b1S^2 - b2S^2,
  a1S*b1S + 1/4*(1+sqrt(5))*b1S^2 + a2S*b2S + 1/4*(1+sqrt(5))*b2S^2],
 [a2S^2*b1S - 1/4*b1S^3 - a1S*a2S*b2S + sqrt(3)/2*a2S*b1S*b2S - sqrt(3)/2*a1S*b2S^2 - 1/4*b1S*b2S^2,
  a1S^2 + a2S^2 - b1S^2 - b2S^2,
  a1S*b1S + sqrt(3)/2*b1S^2 + a2S*b2S + sqrt(3)/2*b2S^2],
 [a2S^2*b1S - (1-(2+2*cos(2*pi/7))/4)*b1S^3 - a1S*a2S*b2S + (cos(pi/7))*a2S*b1S*b2S - (cos(pi/7))*a1S*b2S^2 - (1-(2+2*cos(2*pi/7))/4 )*b1S*b2S^2,
  a1S^2 + a2S^2 - b1S^2 - b2S^2,
  a1S*b1S + (cos(pi/7) )*b1S^2 + a2S*b2S + (cos(pi/7))*b2S^2],
 [a2S^2*b1S - 1/(4+2*sqrt(2))*b1S^3 - a1S*a2S*b2S + cos(pi/8)*a2S*b1S*b2S - cos(pi/8)*a1S*b2S^2 - 1/(4+2*sqrt(2))*b1S*b2S^2,
  a1S^2 + a2S^2 - b1S^2 - b2S^2,
  a1S*b1S + cos(pi/8)*b1S^2 + a2S*b2S + cos(pi/8)*b2S^2],
 [a2S^2*b1S - ((sin(pi/9))**2)*b1S^3 - a1S*a2S*b2S + (cos(pi/9))*a2S*b1S*b2S - (cos(pi/9))*a1S*b2S^2 - ((sin(pi/9))**2)*b1S*b2S^2,
  a1S^2 + a2S^2 - b1S^2 - b2S^2,
  a1S*b1S + (cos(pi/9))*b1S^2 + a2S*b2S + (cos(pi/9))*b2S^2]]


