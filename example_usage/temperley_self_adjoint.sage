R.<a,b,c,d,e,q> = LaurentPolynomialRing(QQbar)
T = TemperleyLiebAlgebra(3,q,R)
bas = T.basis().list()

load("basic_definitions.sage")

coefficients = [c,a,d,e,b]
rot = allRotations(coefficients, 1)
rotCC = allRotations(coefficients, -1)

unit = bas[1]

answer = allMultiplications(rot, rotCC)

vari = list(var("aS,bS_r,bS_i,dS,eS,qS"))
dictionary = {
    a:aS, 
    b:bS_r+iSymbol*bS_i, 
    c:bS_r-iSymbol*bS_i, 
    e:eS, 
    d:dS, 
    q:qS}

answerSymbolic = getSymbolicEquations(answer, dictionary)

def solveNow():
    return solve(answerSymbolic, [aS,bS_r,bS_i,dS,eS])
    
