R.<a,b,c,d,e,aCC,bCC,cCC,dCC,eCC,qTemp> =  CC[]
T = TemperleyLiebAlgebra(3,qTemp,R)
bas = T.basis().list()

varTemp = [a,aCC,b,bCC,c,cCC,d,dCC,e,eCC,qTemp]
variablesSR = list(var("a1S,a2S,b1S,b2S,c1S,c2S,d1S,d2S,e1S,e2S,qS"))

# Conjugate*Normal

rot = []
rotConj = []

rot.append(a*bas[1] + b*bas[3] + c*bas[2] + d*bas[4] + e*bas[0])
rot.append(a*bas[0] + b*bas[2] + c*bas[3] + d*bas[1] + e*bas[4])
rot.append(a*bas[4] + b*bas[3] + c*bas[2] + d*bas[0] + e*bas[1])
rot.append(a*bas[1] + b*bas[2] + c*bas[3] + d*bas[4] + e*bas[0])
rot.append(a*bas[0] + b*bas[3] + c*bas[2] + d*bas[1] + e*bas[4])
rot.append(a*bas[4] + b*bas[2] + c*bas[3] + d*bas[0] + e*bas[1])

rotConj.append(aCC*bas[1] + bCC*bas[3] + cCC*bas[2] + dCC*bas[0] + eCC*bas[4])
rotConj.append(aCC*bas[4] + bCC*bas[2] + cCC*bas[3] + dCC*bas[1] + eCC*bas[0])
rotConj.append(aCC*bas[0] + bCC*bas[3] + cCC*bas[2] + dCC*bas[4] + eCC*bas[1])
rotConj.append(aCC*bas[1] + bCC*bas[2] + cCC*bas[3] + dCC*bas[0] + eCC*bas[4])
rotConj.append(aCC*bas[4] + bCC*bas[3] + cCC*bas[2] + dCC*bas[1] + eCC*bas[0])
rotConj.append(aCC*bas[0] + bCC*bas[2] + cCC*bas[3] + dCC*bas[4] + eCC*bas[1])


answer = []
answerTerms = []

for i in range(len(rot)):
    answer.append(rotConj[i]*rot[i])
    answerTerms.append(answer[i].terms())

answerOtherWay = []
answerTermsOtherWay = []

for i in range(len(rot)):
    answerOtherWay.append(rot[i]*rotConj[i])
    answerTermsOtherWay.append(answerOtherWay[i].terms())

unit = bas[1]

if not unit.is_unit():
    for i in range(1000):
        print "ERROR"
else:
    unit = unit.trailing_support()

theEquations = []

for x in answerTerms:
    for y in x:
        if y.trailing_support() != unit:
            theEquations.append( y.trailing_coefficient() )

sums = [
        a1S + I*a2S, a1S - I*a2S,
        b1S + I*b2S, b1S - I*b2S,
        c1S + I*c2S, c1S - I*c2S,
        d1S + I*d2S, d1S - I*d2S,
        e1S + I*e2S, e1S - I*e2S,
        qS
]

_dictionary = dict(zip(varTemp, sums))

for i in range(len(theEquations)):
    theEquations[i] = expand(theEquations[i].subs(_dictionary))

for v in variablesSR:
    assume(v, "real")

eqSeparatedWithMultiples = []

for equa in theEquations:
    reEq = real_part(equa)
    imEq = imag_part(equa)

    if reEq:
        eqSeparatedWithMultiples.append( reEq )
    if imEq:
        eqSeparatedWithMultiples.append( imEq )

allowedIndices = [i for i in range(len(eqSeparatedWithMultiples))]

eqSeparated = []

for i in allowedIndices:
    for j in range(i+1, len(eqSeparatedWithMultiples)):
        if not expand( eqSeparatedWithMultiples[i] - eqSeparatedWithMultiples[j] ):
            allowedIndices.remove(j)

for i in allowedIndices:
    eqSeparated.append( eqSeparatedWithMultiples[i] )

P.<a1,a2,b1,b2,c1,c2,d1,d2,e1,e2> = QQbar[]

variablesSR.remove(qS)

equationsPolyRing = []
qValues = [ 2*cos(pi/x) for x in range(3,10) ]

for i in range(len(qValues)):
    qV = qValues[i]
    equationsPolyRing.append([])
    for equa in eqSeparated:
        equaTemp = equa.subs({qS:QQbar(qV)})
        equationsPolyRing[i].append(
            P(
                equaTemp.subs( 
                    dict(zip(   variablesSR,
                                list(P.gens())
                                )))))
