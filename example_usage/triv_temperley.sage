#S.<q> = QQbar[]
#S = S.fraction_field()
S.<a,b,bCC,q> = LaurentPolynomialRing(QQbar, 4)
#R.<a,b,bCC> = S[]
R = S
T = TemperleyLiebAlgebra(4,q,R)

bas = T.basis().list()

def rotate( index, counterclockwise ):
    if index == 8:
        return 3
    elif index == 3:
        return 8

    elif index == 13:
        return 4 if not counterclockwise else 9
    elif index == 4:
        return 7 if not counterclockwise else 13
    elif index == 7:
        return 9 if not counterclockwise else 4
    elif index == 9:
        return 13 if not counterclockwise else 7

    elif index == 12:
        return 12

indexList = [8, 13, 7, 12, 3, 4, 9, 12]
coeffList = [a, -a/q, -a/q, a/q^2, b, -b/q, -b/q, b/q^2]
indexListCC = [8, 13, 7, 12, 3, 4, 9, 12]
coeffListCC = [a, -a/q, -a/q, a/q^2, bCC, -bCC/q, -bCC/q, bCC/q^2]

rot = []
rotCC = []

for i in range(4):
    rot.append(
        sum(coeffList[k]*bas[indexList[k]]
            for k in range(len(indexList))
        )
    )
    rotCC.append(
        sum(coeffListCC[k]*bas[indexListCC[k]]
            for k in range(len(indexList))
        )
    )

    indexList = [ rotate(bla, false) for bla in indexList ]
    indexListCC = [ rotate(bla, true) for bla in indexListCC ]

answer = []
answerTerms = []

for i in range(len(rot)):
    answer.append( rot[i]*rotCC[i] )
    answerTerms.append( answer[i].terms() )

def nonVanishingTerms():
    print "These terms must not vanish"
    for x in answerTerms:
        for i in x:
            iMonom = i.trailing_monomial()
            if bas.index(iMonom) in [3,4,9,12]:
                print bas.index(iMonom), i.trailing_coefficient()
        print
    return

def vanishingTerms():
    print "These terms must vanish"
    for x in answerTerms:
        for i in x:
            iMonom = i.trailing_monomial()
            if bas.index(iMonom) not in [3,4,9,12]:
                print bas.index(iMonom), i.trailing_coefficient()
        print
    return

unit = bas[3].trailing_support()

theEquations = []

for x in answerTerms:
    for y in x:
        if not y.trailing_support() == unit:
            theEquations.append( y.coefficients()[0] )
