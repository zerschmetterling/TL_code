varA = flatten([var("a"+str(i)) for i in range(4)])
varB = flatten([var("b"+str(i)) for i in range(1,9)])

load("perfectness.sage")
R = initR(varA + varB)
T = initT(4, q, R)
bas = initBasis(T)

q = R(q)
# iSymbol = R(iSymbol)
varA = [R(a) for a in varA]
varB = [R(b) for b in varB]

tang = T[4].zero()

startTang = T[4].one()

for i in range(4):
    tang += varA[i]*startTang
    startTang = rotateTangle(startTang)

startTang = bas[4][9]

for i in range(8):
    tang += varB[i]*startTang
    startTang = rotateTangle(startTang)
