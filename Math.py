import re
import math
from collections import Counter
import json
import time
## INPUTS
function = "0x+2"


def splitfunc(function):          ## seperates all the terms in the function, for example x^2+5-3x will return ["x^2, "+5", "-3x"]           BASE FUNCTION 1, THIS FUCNTION GETS CALLED BY FINDTERMS() TO SEPERATE TERMS

    # Split the string based on both '+' and '-'
    funcSplit = re.split(r'(\+|\-)', function)

    # Filter out empty strings and strip whitespaces
    parts = [part.strip() for part in funcSplit if part.strip()]

    # Track the signs used for splitting
    signs = []

    # Combine adjacent parts with '+' and '-' signs
    final = []
    for part in parts:
        if part in ['+', '-']:
            signs.append(part)
        else:
            if signs and signs[-1] == '-':
                final.append(signs.pop() + part)
            else:
                final.append(part)

    return final
def FindTerms(function):          ## THIS FUNCTION SORTS OUT THE RAW TERMS INTO A,B,AND C, AND THUS WILL ALWAYS BE CALLED IN FUNCTIONS WHERE A,B,C IS NEEDED
    terms = splitfunc(function)
    a=0
    b=0
    c=0
    for x in range(len(terms)):
        if "x^2" in terms[x]:
            if terms[x] == "" or terms[x] == " ":
                a += 0
            elif terms[x] == "x^2":
                a += 1
            elif terms[x] == "-x^2":
                a += -1
            else:
                try:
                    apartial = terms[x]
                    
                    a += float(apartial.replace("x^2", ""))
                    
                except:
                    print("ERROR")
        elif "x" in terms[x]:
            if terms[x] == "" or terms[x] == " ":
                b += 0
            elif terms[x] == "x":
                b += 1
            elif terms[x] == "-x":
                b += -1
            else:
                try:
                    apartial = terms[x]
                    
                    b += float(apartial.replace("x", ""))
                    
                except:
                    print("ERROR")
        else:
            if "x" not in terms[x]:
                c += float(terms[x])
    return(a,b,c)



def AlphaBetaForm(function):
    a,b,c = FindTerms(function)
    if a != 0:
        alfa = -(float(b))/(2*float(a))
        beta = -((b**2)/(4*a))+c
    else:
        alfa = 0
        beta = 0
    return(alfa,beta)
def Solution(function):
    
    a,b,c =FindTerms(function)
    d = ((float(b) **2) -(4*float(a) * float(c)))
    if d < 0:
        return d,"None", "None"
    elif a == 0 and b == 0 and c != 0:
        return d, "None", "None"   
    elif d == 0:
        x,y = AlphaBetaForm(function)
        return d,x, "None"
    elif d > 0 and a != 0:
        x1  = (-float(b) - math.sqrt(d)) / (2 * float(a))
        x2 = (-float(b) + math.sqrt(d)) / (2 * float(a))
        return d,x1, x2
    elif d > 0 and a == 0:
        x1 = -(float(c)/float(b))
        x2 = "None"
        return d,x1, x2
    else:
        x1 = "None"
        x2 = "None"
        return d,x1, x2
    



def solve2func(function1, function2):
    a1,b1,c1 = FindTerms(function1)
    a2,b2,c2 = FindTerms(function2)
    a = a1-a2
    b = b1-b2
    c = c1-c2 
    if b > 0 and c > 0:                              ## remakes the function into string form so the solve function can read it, if there is no symbol like {}x^2{} it means that it is negative and already has a symbol
        function = "{}x^2+{}x + {}".format(a,b,c) 
    elif b < 0 and c > 0:
        function = "{}x^2{}x + {}".format(a,b,c)
    elif b > 0 and c < 0:
        function = "{}x^2+{}x{}".format(a,b,c)
    elif b < 0 and c < 0:
        function = "{}x^2{}x{}".format(a,b,c)
    elif b > 0:
        function = "{}x^2+{}x".format(a,b)
    elif b < 0:
        function = "{}x^2{}x".format(a,b)
    elif c > 0:
        function = "{}x^2+{}".format(a,c)
    elif c < 0:
        function = "{}x^2{}".format(a,c)
    else:
        function = "{}x^2".format(a)
    return Solution(function)


def solve3func(function1, function2, function3):
    a1,b1,c1 = FindTerms(function1)
    a2,b2,c2 = FindTerms(function2)
    a = a1-a2
    b = b1-b2
    c = c1-c2 

    
    if b > 0 and c > 0:                              ## remakes the function into string form so the solve function can read it, if there is no symbol like {}x^2{} it means that it is negative and already has a symbol
        function = "{}x^2+{}x + {}".format(a,b,c) 
    elif b < 0 and c > 0:
        function = "{}x^2{}x + {}".format(a,b,c)
    elif b > 0 and c < 0:
        function = "{}x^2+{}x{}".format(a,b,c)
    elif b < 0 and c < 0:
        function = "{}x^2{}x{}".format(a,b,c)
    elif b > 0:
        function = "{}x^2+{}x".format(a,b)
    elif b < 0:
        function = "{}x^2{}x".format(a,b)
    elif c > 0:
        function = "{}x^2+{}".format(a,c)
    elif c < 0:
        function = "{}x^2{}".format(a,c)
    else:
        function = "{}x^2".format(a)

    Solutions1 = Solution(function)
    a2,b2,c2 = FindTerms(function3)
    a = a1-a2
    b = b1-b2
    c = c1-c2

    if b > 0 and c > 0:                              ## remakes the function into string form so the solve function can read it, if there is no symbol like {}x^2{} it means that it is negative and already has a symbol
        function = "{}x^2+{}x + {}".format(a,b,c) 
    elif b < 0 and c > 0:
        function = "{}x^2{}x + {}".format(a,b,c)
    elif b > 0 and c < 0:
        function = "{}x^2+{}x{}".format(a,b,c)
    elif b < 0 and c < 0:
        function = "{}x^2{}x{}".format(a,b,c)
    elif b > 0:
        function = "{}x^2+{}x".format(a,b)
    elif b < 0:
        function = "{}x^2{}x".format(a,b)
    elif c > 0:
        function = "{}x^2+{}".format(a,c)
    elif c < 0:
        function = "{}x^2{}".format(a,c)
    else:
        function = "{}x^2".format(a)
    Solutions2 = Solution(function)
    return Solutions1,Solutions2




def solvemore3func(function1, function2, *functions):
    a1,b1,c1 = FindTerms(function1)
    a2,b2,c2 = FindTerms(function2)
    a = a1-a2
    b = b1-b2
    c = c1-c2 

    
    if b > 0 and c > 0:                              ## remakes the function into string form so the solve function can read it, if there is no symbol like {}x^2{} it means that it is negative and already has a symbol
        function = "{}x^2+{}x + {}".format(a,b,c) 
    elif b < 0 and c > 0:
        function = "{}x^2{}x + {}".format(a,b,c)
    elif b > 0 and c < 0:
        function = "{}x^2+{}x{}".format(a,b,c)
    elif b < 0 and c < 0:
        function = "{}x^2{}x{}".format(a,b,c)
    elif b > 0:
        function = "{}x^2+{}x".format(a,b)
    elif b < 0:
        function = "{}x^2{}x".format(a,b)
    elif c > 0:
        function = "{}x^2+{}".format(a,c)
    elif c < 0:
        function = "{}x^2{}".format(a,c)
    else:
        function = "{}x^2".format(a)
    ls = []
    try:
            discriminant,sol1,sol2 = Solution(function)
    except:
            discriminant,sol1 = Solution(function)
            sol2 = "None"
    ls.append(sol1)
    ls.append(sol2)
    for funcs in functions:
        a2,b2,c2 = FindTerms(funcs)
        a = a1-a2
        b = b1-b2
        c = c1-c2

        if b > 0 and c > 0:                              ## remakes the function into string form so the solve function can read it, if there is no symbol like {}x^2{} it means that it is negative and already has a symbol
            function = "{}x^2+{}x + {}".format(a,b,c) 
        elif b < 0 and c > 0:
            function = "{}x^2{}x + {}".format(a,b,c)
        elif b > 0 and c < 0:
            function = "{}x^2+{}x{}".format(a,b,c)
        elif b < 0 and c < 0:
            function = "{}x^2{}x{}".format(a,b,c)
        elif b > 0:
            function = "{}x^2+{}x".format(a,b)
        elif b < 0:
            function = "{}x^2{}x".format(a,b)
        elif c > 0:
            function = "{}x^2+{}".format(a,c)
        elif c < 0:
            function = "{}x^2{}".format(a,c)
        else:
            function = "{}x^2".format(a)
        try:
            discriminant,sol1,sol2 = Solution(function)
        except:
            discriminant,sol1 = Solution(function)
            sol2 = "None"
        ls.append(sol1)
        ls.append(sol2)

    counted_elements = Counter(ls)

# Find elements that occur more than once
    duplicates = [element for element, count in counted_elements.items() if count > 1]

    return ls
def solvemore3func4(*functions):
    t1 = time.time()
    counter = 0
    ls = []
    done = {}
    solutions = []
    functionL = [] 
    for x in functions:
        a1,b1,c1 = FindTerms(x)
        for funcs in functions:
            if x == funcs:
                continue
            counter += 1
            
            a2,b2,c2 = FindTerms(funcs)
            a = a1-a2
            b = b1-b2
            c = c1-c2

            if b > 0 and c > 0:                              ## remakes the function into string form so the solve function can read it, if there is no symbol like {}x^2{} it means that it is negative and already has a symbol
                function = "{}x^2+{}x + {}".format(a,b,c) 
            elif b < 0 and c > 0:
                function = "{}x^2{}x + {}".format(a,b,c)
            elif b > 0 and c < 0:
                function = "{}x^2+{}x{}".format(a,b,c)
            elif b < 0 and c < 0:
                function = "{}x^2{}x{}".format(a,b,c)
            elif b > 0:
                function = "{}x^2+{}x".format(a,b)
            elif b < 0:
                function = "{}x^2{}x".format(a,b)
            elif c > 0:
                function = "{}x^2+{}".format(a,c)
            elif c < 0:
                function = "{}x^2{}".format(a,c)
            else:
                function = "{}x^2".format(a)
            try:
                discriminant,sol1,sol2 = Solution(function)
            except:
                discriminant,sol1 = Solution(function)
                sol2 = "None"
            
            lsa = [sol1,sol2]
            ls.append(lsa)
            done["Sol {}".format(counter)] = sol1, sol2
            try:
                solutions.append(float(sol1))
            except:
                pass
            try:
                solutions.append(float(sol2))
            except:
                pass                  
            functionL.append([funcs,x])
           
    
            
            
    
    solutions = list(dict.fromkeys(solutions))
    solutions.sort()
    
    checked = True
    am = 0          ## global counter
    amlocal = 0     ## local counter
    totals = {}     ## dict that will hold all the awnser sets, for example all the true and false values for if the functions contain the selected list item
    varlist = []    ## holds the true and false values and passes them to the totals{} dict and then resets after every loop  ## this is for testing
    counter = 0     ## counter so we know what set we are working with      could have also been done with solutions[x]       ## this is for testing
    
    intersects = {}
    intersectsL = [] 
    SaveDebugList = []          # contains all the information with local checks and stuff
    
    for x in range(len(solutions)):
        counter +=1
        amlocal = -1
        for y in done:
            am += 1
            amlocal +=1
            if done[y][0] == solutions[x] or done[y][1] == solutions[x]:
                #print(done[y][0],"#", done[y][1], "C", solutions[x], "--", True, "##Check numb local: ",amlocal, "##Check numb global: ",am, "##", functionL[amlocal][0],"&&",functionL[amlocal][1] )
                a = done[y][0],"#", done[y][1], "C", solutions[x], "--", True, "##Check numb local: ",amlocal, "##Check numb global: ",am, "##", functionL[amlocal][0],"&&",functionL[amlocal][1]
                
                intersectsL.append(functionL[amlocal][0] and functionL[amlocal][1])
                varlist.append(True)
            else:
                #print(done[y][0],"#", done[y][1], "C", solutions[x], "--", False, "##Check numb local: ",amlocal, "##Check numb global: ",am, "##", functionL[amlocal][0],"&&",functionL[amlocal][1])
                a = done[y][0],"#", done[y][1], "C", solutions[x], "--", False, "##Check numb local: ",amlocal, "##Check numb global: ",am, "##", functionL[amlocal][0],"&&",functionL[amlocal][1]
                varlist.append(False)
            SaveDebugList.append(a)
                

                
            
        totals.update({"Set {}".format(counter): varlist })

        intersectsL = list(dict.fromkeys(intersectsL))
        if False not in varlist:
            intersects.update({solutions[x]: intersectsL})
        varlist = []
        intersectsL = []
    #print("Amount of functions", len(functions))
                    

                            
                            




    t2 = time.time()
    elapsedTime = t2 -t1
                
    return done, checked, totals,intersects,SaveDebugList,elapsedTime,solutions, len(functions)          ## done contains all solutions, checked contains the list of solutions, totals contains the solutions printed out in total, and intersects is the dict with all the interect points that each functions has
       




#sols,result,totals,intersects,SaveDebugList = solvemore3func4("x^2+2", "-x^2+4x+2", "-x^2-3x+2","-x^2-2x+2", "x+2", "2-x","-4x+2", "x+20")

# sols,result,totals,intersects,SaveDebugList,elapsedTime,solutions, amfunctions = solvemore3func4("x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5","x^2+2", "-x^2+4", "0x+3", "x-x^2+5")


# for x in intersects:
#     print("All functions Intersect in X = ",x)
# if len(intersects) == 0:
#     print("The functions do not share a same intersect")
# #json_object = json.dumps(SaveDebugList, separators=(',', ':'))



# with open("sample.json", "w") as outfile:
#     for index, sublist in enumerate(SaveDebugList):
#         json_str = json.dumps(sublist)
#         if index < len(SaveDebugList) - 1:
#             outfile.write(json_str + "\n")
#         else:
#             outfile.write(json_str)
#     json_str = json.dumps(amfunctions)
#     outfile.write("\n" +"Amount of functions: " + json_str)
#     json_str = json.dumps(solutions)
#     outfile.write("\n" +"All intersect points: " + json_str)
#     for x in intersects:
#         json_str = json.dumps(x)
#         outfile.write("\n" + "All functions Intersect in X = " + json_str)
#     if len(intersects) == 0:
#         json_str = json.dumps("\n" + "The functions do not share a same intersect")
#         outfile.write(json_str)
#     json_str = json.dumps(elapsedTime)
#     outfile.write("\n" +"Elapsed Time: " + json_str + " Seconds")
