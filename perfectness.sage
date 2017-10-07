###########################################
#   stuff to calculate perfectness with   #
###########################################


###
#   Convert basis tangle to easy-to-parse list. Cleaner than always using an instance method.
###

def convertToList( base_tangle ):
    return base_tangle.diagram().standard_form()

###
#   Obtain the adjoint of a basis tangle by simply flipping the sign of each node, then reordering the list
###

def adjointBasisElement( base_tangle ):
    tList = convertToList(base_tangle)
    tempList = []
    for i in range(T.order()):
        tempList.append([-tList[i][0], -tList[i][1]])
    [x.sort() for x in tempList]
    tempList.sort()
    return [x for x in bas if convertToList(x) == tempList][0]

#### RETURN THE ADJOINT OF A LINEAR COMBINATION

def adjoint( tangle ):
    terms = tangle.terms()
    
    coeffs = [SR(x.trailing_coefficient()) for x in terms]
    vector = [x.trailing_monomial() for x in terms]

    coeffsCC = [R(x.subs(cc_dict)) for x in coeffs]
    vectorCC = [adjointBasisElement(x) for x in vector]

    adjointTangle = sum( coeffsCC[k]*vectorCC[k] for k in range(len(terms)) )

    return adjointTangle

### ROTATE AN INPUT TANGLE ONE CLICK 

def rotateTangle( vector_in ):
    this_TL = vector_in.parent()
    this_bas = this_TL.basis().list()
    basis_as_list = [ convertToList(x) for x in this_bas ]
    strands = this_TL.order()
    in_list = convertToList(vector_in)
    out_list = []

    for i in range(strands):
        out_list.append([])
        for x in in_list[i]:
            if x == -1:
                out_list[i].append( 1 )
            elif x == strands:
                out_list[i].append( -strands )
            else:
                out_list[i].append( x + 1 )
    
    [x.sort() for x in out_list]
    out_list.sort()
    
    rotated = bas[ basis_as_list.index( out_list ) ]
    return rotated 

### GIVEN A LIST OF COEFFICIENTS AND AN ORIENTATION, COMPUTE ALL ROTATIONS

def allRotations(coeffs):
    vectors = copy(bas)
    number_of_rotations = 0
    rotations = []
    
    while number_of_rotations < 2*T.order():
        appendThis = sum(R(coeffs[k])*vectors[k] for k in range(dim(T)))

        if appendThis in rotations:
            break
        
        rotations.append(appendThis)
        
        for j in range(dim(T)):
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
    for v in R.gens():
        var(v)

    temp = []
    for x in symbolizeThis:
        coef = x.trailing_coefficient()
        vect = x.trailing_monomial()
        if not gSE_TL:
            gSE_TL = vect.parent()

        equa = expand(SR(coef)).subs({iSymbol^2 : -1})
        if equa.is_numeric() == False:
            imagPart, realPart = equa.maxima_methods().divide(iSymbol)
        else:
            imagPart = 0
            realPart = equa

        if vect == gSE_TL.one():
            if ( realPart != 0 and imagPart != 0 ):
                print "Failure"
                # raise ValueError
            else:
                if realPart != 0:
                    temp.append(realPart != 0)
                elif imagPart != 0:
                    temp.append(imagPart != 0)
        else: 
            temp.append( realPart == 0)
            temp.append( imagPart == 0)

    answer = list(set(temp))
    try:
        answer.remove(True)
    except:
        pass
    
    return answer

### DECOMPOSE THE SET OF BASIS TANGLES INTO SETS CLOSED UNDER ROTATION
# This is useful if you want to compute 
#   (a) rotation invariant tangles
#   (b) n-tangles that are eigenvectors of rotation with eigenvalue an n-th root of unity
# Note that case (b) only makes sense for components of the decomposition of size 2n. As an example: Case (b) can not be found in TL_3, because the decomposition
# consists of a two and a three element set. In TL_4, however, there is a set with 8 = 2*4 elements, which might define a perfect tangle satisfying (b). Unfortunately, the identity is not obtained by multiplication of such tangles, hence: there is no such thing in TL_4

def rotationInvariantDecomposition():
    decomp = []
    indices = range(dim(T))
    component = 0
    while indices:
        decomp.append([])
        i = indices[0]
        while i not in decomp[component]:
            decomp[component].append(i)
            indices.remove(i)
            i = bas.index( rotateTangle(bas[i]))
        component += 1
    return decomp

### To this function one supplies a dictionary, and the two rotation arrays. It then simply substitutes the words from the dictionary, and returns the output of 
# allMultiplications, i.e. an array, which - if the dictionary specifies a solution - should only consist of zeroes;

def checkSolution( dict_in, rot_in ):
    answerSolution = getSymbolicEquations(allMultiplications(rot_in))
    answerSolution = [expand(SR(x).subs(dict_in)) for x in answerSolution]

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

## helpful for rotation invariant solutions

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
    all_indices = range(dim(T))
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

## get the variables to solve the equations for
def solveForThese(coeff_in):
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
            
