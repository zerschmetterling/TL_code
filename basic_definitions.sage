#############################################
#   BASIC FUNCTIONALITIES AND DEFINITIONS   #
#############################################

fixed_params = list(var("q, iSymbol")) # TL parameter q and a non-implementation of the complex unit
cc_dict = {iSymbol:-iSymbol} # complex conjugation

def initR(variableList): # returns a polynomial ring in the given variables (the input type must be a list) 
    try :
        R = LaurentPolynomialRing( QQbar, variableList + fixed_params )
    except:
        R = PolynomialRing( QQbar, variableList + fixed_params )

    return R

def initT(order_max, qValue, base_ring): # return a list of TL algebras
    T_temp = [TemperleyLiebAlgebra(i, qValue, base_ring) for i in range(order_max+1)]
    return T_temp

def initBasis(TL_list): # returns the basis for a given Temperley-Lieb algebra
    temp_bas = []
    for t in TL_list:
        temp_bas.append( t.basis().list() )
    return temp_bas

def initCoeff(ring_in): # takes a ring of polynomials in a_1, a_2, b_1, b_2,.... and returns a list [a_1 ]
    ring_generators = ring_in.gens()
    _coeff = [ 
            ring_generators[2*i]+iSymbol*ring_generators[2*i+1] 
            for i in range(len(ring_generators-2)/2)
        ]
    
    return _coeff

def reduceIsquared_tangle(v): # takes a tangle and removes each instance of iSymbol^2 for a factor of -1, then returns the tangle free from those shackles
    _base_ring = v.base_ring()
    v_coeff = [t.trailing_coefficient() for t in v.terms()]
    v_vec = [t.trailing_monomial() for t in v.terms()]

    for i in range(len(v_coeff)):
        t = reduceIsquared_coeff(v_coeff[i])
        v_coeff[i] = t

    vector = sum(v_coeff[i]*v_vec[i] for i in range(len(v_vec)))

    return vector

def reduceIsquared_coeff(c_in): # see above
    t = SR(copy(c_in))
    x, y = t.maxima_methods().divide(SR(iSymbol)^2)

    while x:
        t = y - x
        x,y = t.maxima_methods().divide(SR(iSymbol)^2)

    return t

def solveForThese(coeff_in): # given the list of coefficients, return all variables actually used, i.e. the ones you want to solve for
    vari = []
    for x in coeff_in:
        for v in list(x.variables()):
            if v not in vari:
                vari.append(SR(v))
    try:
        vari.remove(iSymbol)
    except:
        pass
    
    return vari
            
