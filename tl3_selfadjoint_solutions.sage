(
    (
        (q == sqrt(2) && b1 == 0 && e2 == 0 && e1 == 0) 
        || 
        (
            (
                Inequality[0, Less, q, Less, sqrt(2)] 
                || 
                q > sqrt(2)
            ) 
            && b1 == 0 && e2 == 0 && e1 == 0
        )
    )
    && d1 == -(e1*q) 
    && c1 == (-2*e1)/q

)

# This unfortunately amounts to saying that all free variables --- namely
# b1, c1, d1, e1, e2 --- are zero.



 ||


 
(
    (
        (
            (
                (e1 < 0 && (c1 == Sqrt[e1^2]/sqrt(3) || c1 == sqrt(3)*Sqrt[e1^2]))
                || 
                (e1 > 0 && (c1 == -(sqrt(3)*Sqrt[e1^2]) || 
                    c1 == -(Sqrt[e1^2]/sqrt(3)))
                )
            ) 
            && q == (-4*c1*e1)/(c1^2 + e1^2) 
            && e2 == 0 && d1 == -(e1*q)
        )
        || 
        (
            (
                (e1 < 0 && (Inequality[0, Less, d1, Less, sqrt(3)*Sqrt[e1^2]] 
                    || 
                    d1 > sqrt(3)*Sqrt[e1^2]))
                || 
                (e1 > 0 && (d1 < -(sqrt(3)*Sqrt[e1^2]) 
                    || 
                    Inequality[-(sqrt(3)*Sqrt[e1^2]), Less, d1, Less, 0]))
            ) 
            && 
            q == (-4*d1*e1^3)/(d1^2*e1^2 + e1^4) 
            && e2 == 0 
            && 
            c1 == (-e1^2 - d1*e1*q - (e1*(-2*d1*e1 - d1^2*q - e1^2*q))/d1)/
                    (d1 + e1*q)
        )
    )
    && 
    b1 == (-2*d1*e1 - d1^2*q - e1^2*q - e2^2*q)/(2*d1)
)



# Expanding this:
######## FIRST
    b1 == (-2*d1*e1 - d1^2*q - e1^2*q - e2^2*q)/(2*d1)
    && q == (-4*c1*e1)/(c1^2 + e1^2) 
    && e2 == 0 
    && d1 == -(e1*q)
    && (c1 == -e1/sqrt(3) || c1 == -sqrt(3)*e1)
    && e1 != 0

########### SECOND
    b1 == (-2*d1*e1 - d1^2*q - e1^2*q - e2^2*q)/(2*d1)
    && q == (-4*d1*e1^3)/(d1^2*e1^2 + e1^4) 
    && e2 == 0 
    && c1 == (-e1^2 - d1*e1*q - (e1*(-2*d1*e1 - d1^2*q - e1^2*q))/d1)/(d1 + e1*q)
    && d1 != -sqrt(3)*e1
    && ((e1 < 0 && d1 > 0) || (e1 > 0 && d1 < 0))


# Once more 
####### FIRST 1; e1 is freely choosable as long as not null
    b1 == e1
    q == sqrt(3)
    d1 == -e1*sqrt(3)
    c1 == -e1/sqrt(3) 
####### FIRST 2
    b1 == e1
    q == sqrt(3)
    d1 == -sqrt(3)*e1
    c1 == -sqrt(3)*e1


########### SECOND 1; here: e1 < 0 and a condition on d1 > 0
    b1 == 2*d1^2*e1/(d1^2 + e1^2) + 2*e1^3/(d1^2 + e1^2) - e1
    q == (-4*d1*e1)/(d1^2 + e1^2)
    c1 == e1^3*q/((e1*q + d1)*d1) + e1^2/(e1*q + d1)
    d1 != -sqrt(3)*e1
    e1 < 0
    d1 > 0
########### SECOND 2; here: e1 > 0 and a condition on d1 < 0
    b1 == 2*d1^2*e1/(d1^2 + e1^2) + 2*e1^3/(d1^2 + e1^2) - e1
    q == (-4*d1*e1)/(d1^2 + e1^2) 
    c1 == e1^3*q/((e1*q + d1)*d1) + e1^2/(e1*q + d1)
    d1 != -sqrt(3)*e1
    e1 > 0 
    d1 < 0
########
# But we want ( q >= 2 || q = 2*cos(pi/x) for x in NN )
# Both q-values are the same in the second thingy, so:
# We can feed this into Mathematica, using
#   Assume[ Element[x, Integers],
#           Reduce[ { (-4*a*b)/(a^2 + b^2) == 2 Cos[Pi/x],
#                       a != -sqrt(3)*b, x > 0, a < 0, b > 0},
#                   {a,b}]]
# of course with flipping the last two conditions as well.
# We obtain nonsense:

# First, the discrete case:
# SECOND 1: d1 is parameter, d1 > 0
(x == 3 && (e1 == -2*d1*(1 + (sqrt(3))/2) ||
     e1 == -2*d1*(1 - (sqrt(3))/2)))

||

(x == 4 && (e1 == sqrt(2)*(-d1 - d1/sqrt(2)) ||
     e1 == sqrt(2)*(-d1 + d1/sqrt(2))))

||

(x == 5 && (e1 == -1/4*sqrt(2)*(2*sqrt(2)*d1 + sqrt(-d1*(sqrt(5) - 5)))*(sqrt(5) - 1)
           ||
     e1 == (-1 + sqrt(5))*(-d1 + ((-1 + sqrt(5))*(1 + sqrt(5))*
          sqrt(-((-1 + (1 + sqrt(5))^2/16))*d1))/4)))

||

(x == 6 &&
  e1 == (-2*sqrt(3)*d1 - sqrt(3)*d1)/3)

||

(x >= 7 && ## Note here that tan(pi/7) > 0 and sqrt(sec^2-1) = sqrt(tan^2)
    (e1 == -d1/cos(pi/x)*(1 + sqrt(1-cos(pi/x)^2))||
     e1 == sec(pi/x)*(-d1 + cos(pi/x)*sqrt(-(d1^2*(-1 + cos(pi/x)^2)*
            sec(pi/x)^2)))))

###
# using substitutionStuff.sage we can now solve these for one parameter families and
# hope that everything works out :D
#
# Start with d1 > 0
x == 3:
    {e1: -d1*(sqrt(3) + 2),
     e2: 0,
     c1: d1*(sqrt(3) + 2)^2,
     b1: -d1*(sqrt(3) + 2),
     q : 1}
||
    {e1: d1*(sqrt(3) - 2),
     e2: 0,
     c1: d1*(sqrt(3) - 2)^2,
     b1: d1*(sqrt(3) - 2),
     q: 1}

x == 4:
    {e1: -1/2*sqrt(2)*d1*(sqrt(2) + 2),
     e2: 0,
     c1: 1/2*d1*(sqrt(2) + 2)^2,
     b1: -d1*(sqrt(2) + 1),
     q: sqrt(2)}
||
    {e1: 1/2*sqrt(2)*d1*(sqrt(2) - 2),
     e2: 0,
     c1: 1/2*d1*(sqrt(2) - 2)^2,
     b1: -d1*(sqrt(2) - 1),
     q: sqrt(2)}

x == 5: ## here the solution depends on d1 non-linearly, but solving the output of checkSolution for d1 shows that the only solving values for d1 are 0 and 1. We obtain:


## something goes wrong tho, let's go ahead to
x == 6:
    {e1: -sqrt(3)*d1,
     e2: 0, 
     c1: 3*d1, 
     b1: -sqrt(3)*d1, 
     q: sqrt(3)}




