#############################################
#   BASIC FUNCTIONALITIES AND DEFINITIONS   #
#############################################

var("iSymbol")
fixed_params = [var("q")] # TL parameter q and a non-implementation of the complex unit

def initR(variableList): # returns a polynomial ring in the given variables (the input type must be a list) 
    variableListTemp = copy(variableList)
    for x in variableList:
        variableListTemp.append(var(str(x)+"CC"))
    
    try :
        R = LaurentPolynomialRing( QQbar, variableListTemp + fixed_params )
    except:
        R = PolynomialRing( QQbar, variableListTemp + fixed_params )

    return R

def conjugate_dictionary(R):
    gens = list(R.gens())
    gens = [SR(x) for x in gens]
    del(gens[-1])

    dictio = {}

    l = (len(gens))/2

    first_half = [gens[i] for i in range(l)]
    second_half = [gens[i] for i in range(l, 2*l)]

    dictio.update( zip(first_half, second_half) )
    dictio.update( zip(second_half, first_half) )

    return dictio

def initT(order_max, qValue, base_ring): # return a list of TL algebras
    T_temp = [TemperleyLiebAlgebra(i, qValue, base_ring) for i in range(order_max+1)]
    return T_temp

def initBasis(TL_list): # returns the basis for a given Temperley-Lieb algebra
    temp_bas = []
    for t in TL_list:
        temp_bas.append( t.basis().list() )
    return temp_bas

def initCoeff(ring_in): # takes a ring of polynomials in a_1, a_2, b_1, b_2,.... and returns a list [a_1 ]
    ring_generators = list(ring_in.gens())
    _coeff = [ 
            ring_generators[2*i]+iSymbol*ring_generators[2*i+1] 
            for i in range(len(ring_generators)/2-1)
        ]
    
    return _coeff

def reduceIsquared_tangle(v): # takes a tangle and removes each instance of iSymbol^2 for a factor of -1, then returns the tangle free from those shackles
    this_R = v.base_ring()
    v_coeff = [t.trailing_coefficient() for t in v.terms()]
    v_vec = [t.trailing_monomial() for t in v.terms()]

    for i in range(len(v_coeff)):
        t = this_R(reduceIsquared_coeff(v_coeff[i]))
        v_coeff[i] = t

    vector = sum(v_coeff[i]*v_vec[i] for i in range(len(v_vec)))

    return vector

def reduceIsquared_coeff(c_in): # see above
    t = SR(c_in)
    x, y = SR(t).maxima_methods().divide(SR(iSymbol^2))
    while x:
        t = y - x
        x,y = SR(t).maxima_methods().divide(SR(iSymbol^2))

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

def getCoeffFromTangle(vector):
    this_TL = vector.parent()
    this_bas = bas[T.index(this_TL) ]
    v_terms = vector.terms()

    answer = []

    for bT in this_bas:
        for x in v_terms:
            if x.trailing_monomial() == bT:
                answer.append(x.trailing_coefficient())
                v_terms.remove(x)
                break

    return answer

def tikzTangle(tangle):
    out_str = "\\tikzexternaldisable\n\\begin{align*}\n&"
    order = tangle.parent().order()
    terms = [ [latex(x.trailing_coefficient()), x.trailing_monomial().diagram()] for x in tangle.terms()]

    for x in terms:
        out_str += "\t("+str(x[0]).replace("\\mathit{iSymbol}", "i")+")\\cdot\\,\n"
        out_str += "\t\\begin{tikzpicture}[scale=0.5, baseline]\n"
    
        for i in range(1,order+1):
            out_str += "\t\t\\coordinate (top"+str(i)+") at (-"+str(i/2)+",1);\n"
            out_str += "\t\t\\coordinate (bottom"+str(i)+") at (-"+str(i/2)+",-1);\n"
        
        for pair in x[1]:
            start = "top"+str(pair[0]) if pair[0] > 0 else "bottom"+str(-pair[0])
            finish = "top"+str(pair[1]) if pair[1] > 0 else "bottom"+str(-pair[1])
            out_str += "\t\t\\draw ("+start+") -- ("+finish+");\n"

        out_str += "\t\\end{tikzpicture}\\\ &+\n"

    out_str += "\\end{align*}\n\\tikzexternalenable"
    return out_str
