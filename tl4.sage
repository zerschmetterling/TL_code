varList = list(var("aR,aI,bR,bI,cR,cI,dR,dI,eR,eI,fR,fI,gR,gI,hR,hI,iR,iI,jR,jI,kR,kI,lR,lI,mR,mI,nR,nI"))

temperley_order = 4

load("basic_definitions.sage")
# coeff = [bR + iSymbol*bI,aR + iSymbol*aI ]
# dictCC = {iSymbol : -iSymbol}

# coefficients = coeffsRotationInvariant(coeff)
# coefficientsCC = [x.subs(dictCC) for x in coefficients]

# rot = allRotations(coefficients, 1)
# rotCC = allRotations(coefficientsCC, -1)

# answer = allMultiplications(rot,rotCC)

# variables = var("a1,a2,b1,b2,c1,c2")
# var("qS")

# dictionary = {
    # aR  : a1,
    # aI  : a2,
    # bR  : b1,
    # bI  : b2,
    # cR  : c1,
    # cI  : c2,
    # q   : qS
    # }

# symEq = getSymbolicEquations(answer, dictionary)
