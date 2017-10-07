varList = list(var("a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,q"))

def solution_dict(dictionary, free_var):
    temp = copy(dictionary)
    keys = temp.keys()

    try:
        keys.remove(free_var)
        sol_dict = {free_var:temp.pop(free_var)}
    except:
        sol_dict = {}
        pass

    itera = 0

    while len(sol_dict) != len(dictionary):
        for key in keys:
            value = expand(temp.get(key))
            
            try:
                var = list(value.variables())
            except:
                var = None
                pass
    
            if var and var != [free_var]:
                valueNew = factor(value.subs(sol_dict))
                # print key, "\n", value, "\n", valueNew, "\n"
                
                temp.update({key : valueNew})
            else:
                if value:
                    sol_dict.update({key: factor(temp.pop(key))})
                else:
                    sol_dict.update({key: temp.pop(key)})
            
        keys = temp.keys()

    return sol_dict


def maxNumVar(dictionary):
    return max([ len(x.variables()) for x in dictionary.values()])

def rhs_var(dictionary):
    var = []
    for x in dictionary.values():
        try:
            var += list(x.variables())
        except:
            pass
    return list(set(var))


def solution_dict_sec_param(dictionary, free_var):
    temp = copy(dictionary)
    keys = temp.keys()

    try:
        fV = free_var[0]
        keys.remove(fV)
        sol_dict = {fV:temp.pop(fV)}
    except:
        sol_dict = {}
        pass

    itera = 0

    while len(sol_dict) != len(dictionary):
        for key in keys:
            value = expand(temp.get(key))
            
            try:
                var = list(value.variables())
            except:
                var = None
                pass
    
            if var and ( var != free_var and var != [free_var[0]] and var != [free_var[1]]):
                valueNew = factor(value.subs(sol_dict))
                # print key, "\n", value, "\n", valueNew, "\n"
                
                temp.update({key : valueNew})
            else:
                if value:
                    sol_dict.update({key: factor(temp.pop(key))})
                else:
                    sol_dict.update({key: temp.pop(key)})
            
        keys = temp.keys()

    return sol_dict
