R.<
    aR,aI,
    bR,bI,
    cR,cI,
    dR,dI,
    eR,eI,
    fR,fI,
    fR,fI,
    fR,fI,
    fR,fI,
    fR,fI,
    fR,fI,
    fR,fI,
    fR,fI,
    fR,fI,
    fR,fI,
    fR,fI,
    q,iSymbol> = LaurentPolynomialRing(QQbar)
T = TemperleyLiebAlgebra(5,q,R)
load("basic_definitions.sage")

coeff = [R.gens()[2*i]+iSymbol*R.gens()[2*i+1] for i in range(len(R.gens())/2-1)]

# coefficients = coeffsRotationInvariant(coeff)
# coefficientsCC = [x.subs(cc_dict) for x in coefficients]

# rot = allRotations(coefficients,1)
# rotCC = allRotations(coefficientsCC,-1)

# answer = allMultiplications(rot,rotCC)

# var("a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2,qS")

# sym_dict = {
        # aR : a1,
        # aI : a2,
        # bR : b1,
        # bI : b2,
        # cR : c1,
        # cI : c2,
        # dR : d1,
        # dI : d2,
        # eR : e1,
        # eI : e2,
        # fR : f1,
        # fI : f2,
        # q  : qS,
        # iSymbol : iS
    # }

# symEq = getSymbolicEquations(answer, sym_dict)

# def use_mathematica():
    # return reduce_mathematica(symEq)
