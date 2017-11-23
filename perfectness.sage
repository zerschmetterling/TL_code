###########################################
#   stuff to calculate perfectness with   #
###########################################
load("basic_definitions.sage")
load("planar_tangles.sage")

### DECOMPOSE THE SET OF BASIS TANGLES INTO SETS CLOSED UNDER ROTATION
# This is useful if you want to compute rotation invariant tangles

def rotationInvariantDecomposition(TL_algebra):
    bas_temp = TL_algebra.basis().list()
    decomp = []
    indices = range(dim(TL_algebra))
    component = 0
    while indices:
        decomp.append([])
        i = indices[0]
        while i not in decomp[component]:
            decomp[component].append(i)
            indices.remove(i)
            i = bas_temp.index( rotateTangle(bas_temp[i]))
        component += 1
    return decomp

###
#   Given a list of 
###

def coeffsRotationInvariant(coeff):
    temp = copy(coeff)
    for x in rotationInvariantDecomposition():
        equateThis = temp[x[0]]
        for y in x:
            if len(temp[y].variables()) < 3:
                equateThis = temp[y]
        for w in x:
            temp[w] = equateThis
            
    return temp

## helpful for self adjoint solutions

def coeffsSelfadjoint(coeff):
    temp = copy(coeff) 
    all_indices = range(len(coeff))
    for x in bas:
        if x == adjoint(x):
            x_ind = bas.index(x)
            bla = 1/2* (temp[x_ind] + R(SR(temp[x_ind]).subs(cc_dict)))
            temp[x_ind] = bla
            all_indices.remove(x_ind)

    for x_ind in all_indices:
        y_ind = bas.index(adjoint(bas[x_ind]))
        temp[x_ind] = SR(temp[y_ind]).subs({iSymbol:-iSymbol})
        all_indices.remove(x_ind)
        all_indices.remove(y_ind)
        
    return [R(x) for x in temp]


### GIVEN A LIST OF COEFFICIENTS AND AN ORIENTATION, COMPUTE ALL ROTATIONS

def allRotations(coeffs):
    number_of_rotations = 0
    this_bas = [b for b in bas if len(b) == len(coeffs)][0]
    this_T = T[bas.index(this_bas)]
    this_R = this_T.base_ring()
    vectors = this_bas
    rotations = []
    
    while number_of_rotations < 2*this_T.order():
        appendThis = sum(this_R(coeffs[k])*vectors[k] for k in range(dim(this_T)))

        if appendThis in rotations:
            break
        
        rotations.append(appendThis)
        
        for j in range(dim(this_T)):
            vectors[j] = rotateTangle(vectors[j])
        
        number_of_rotations += 1
    
    return rotations

### MULTIPLY TANGLES WITH ADJOINT TANGLES; REMOVE REDUNDANT EQUATIONS

def allMultiplications(rot):
    rotCC = [ adjoint(x) for x in rot ]
    temp = []

    for i in range(len(rot)):
        prodTerms = (rot[i]*rotCC[i]).terms()
        temp.append( prodTerms )

    answer = []
    for x in temp:
        for y in x:
            if y not in answer:
                answer.append(y)
    
    return answer

### TURN THE OBTAINED EQUATIONS, WHICH LIVE IN A RING OF LAURENT POLYNOMIALS, INTO SYMBOLIC EQUATIONS WHICH SAGE CAN SOLVE

def getSymbolicEquations(symbolizeThis):
    this_R = symbolizeThis[0].base_ring()
    this_TL = symbolizeThis[0].parent()

    temp = []

    for x in symbolizeThis:
        coef = SR(x.trailing_coefficient())
        vect = x.trailing_monomial()

        if vect == this_TL.one():
            temp.append(coef != 0)
        else: 
            temp.append(coef == 0)

    answer = list(set(temp))
    try:
        answer.remove(True)
    except:
        pass
    
    return answer


def symbolicEquations_real_split(symEq_complex):
    gens = [SR(v) for v in R.gens()]
    gens.remove(q)

    split_dict = {}

    for i in range( len(gens)/2 ):
        v = gens[i]

        split_dict.update({v : var(str(v)+"r") + iSymbol*var(str(v)+"i")})
        split_dict.update({eval(str(v)+"CC") : var(str(v)+"r") - iSymbol*var(str(v)+"i")})

    temp = []

    neq = ( SR(1) != SR(0) ).operator()
    eq = ( SR(0) == SR(0) ).operator()

    for x in symEq_complex:
        lhs = expand( x.lhs().subs(split_dict) )
        lhs = reduceIsquared_coeff(lhs)
        op = x.operator()

        if lhs.is_numeric() == False:
            imagPart, realPart = lhs.maxima_methods().divide(iSymbol)
        else:
            imagPart = 0
            realPart = lhs

        if op == neq:
            if imagPart != 0:
                temp.append( imagPart != 0 )
                if realPart != 0:
                    print "ENCOUNTERED AN OR"
                    temp.append( realPart != 0 )
            elif realPart != 0:
                temp.append( realPart != 0 )
            else:
                raise ValueError
        else:
            temp.append( op(imagPart, 0) ) 
            temp.append( op(realPart, 0) ) 

    answer = list(set(temp))
    
    return answer

### To this function one supplies a dictionary, and the two rotation arrays. It then simply substitutes the words from the dictionary, and returns the output of 
# allMultiplications, i.e. an array, which - if the dictionary specifies a solution - should only consist of true statements.

def checkSolution( dict_in, symEq_in ):
    answerSolution = [expand(SR(x).subs(dict_in)) for x in symEq_in]
    answerSolution = [expand(x^2) for x in answerSolution]
    
    return answerSolution

### Use Mathematica's reduction
# def reduce_mathematica(symEqn, variables):
def reduce_mathematica(symEqn):
    symEqnTemp = [ x == 0 for x in symEqn ]
    symEqnTemp.append(q != 0)
    # symEqnTemp.append( [variables[i] != 0 || for i in range(len(variables))-1] ) 

    m_symEqn = mathematica(symEqnTemp)
    reduced = m_symEqn.Reduce("Reals")
    return reduced

