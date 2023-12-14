import math
def MoleculeStable(element1,element2):
    groupI = {"H": [1, 2.1], "Li": [1, 1.0], "Na":[1,0.9]}        ## element name, all possible positive oxidiation numbers, and then final number is the EN worth
    groupII = {"Be": [2, 1.5], "Mg": [2, 1.2], "Ca":[2,1.0]}
    groupIII = {"Fee": [7, 4.0]}
    groupIV = {"Fee": [7, 4.0]}
    groupV = {"Fe": [7, 4.0]}
    groupVI = {"S": [4,6, 2.5]}
    groupVII = {"F": [7, 4.0], "Cl": [1,3,5,7, 3.0]}
    groups =  [groupI,groupII,groupII,groupIV,groupV,groupVI, groupVII]
    # element1 = "Na"
    # element2 = "Cl"
    sols = []


    groupcounter = 0
    for x in groups: 
        groupcounter +=1                ## manually count the group number
        if element1 in x:
            En1 = x[element1][-1]        ## en worth
            oxnumbs1 = x[element1][:-1]  ## oxidation numbers
            groep1 = groupcounter    
        if element2 in x:
            En2 = x[element2][-1]        ## en worth
            oxnumbs2 = x[element2][:-1]  ## oxidation numbers
            groep2 = groupcounter


    if En1 > En2:
        if element1 != "H":
            oxidation1 = abs(groep1 - 8)
        else:
            oxidation1 = abs(groep1 - 2)
        for x in oxnumbs2:
            counter1 = 1
            counter2 = 1
            oxidation2 = x
            addox1 = oxidation1 
            addnox2 = oxidation2
            while oxidation1 != oxidation2:
                if oxidation1 > oxidation2:
                    oxidation2 += addnox2
                    counter2 +=1
                else:
                    oxidation1 += addox1
                    counter1 +=1
            if counter1 == 1 and counter2 ==1:
                print(element2+element1)
                sols.append(element2+element1)
            elif counter1 == 1:
                print(element2+str(counter2) + element1+ "     OG({}) = {}".format(element2,addnox2)) ## addnox2 is the orignal oxidation of the element
                sols.append(element2+str(counter2) + element1+ "     OG({}) = {}".format(element2,addnox2))
            elif counter2 == 1:
                print(element2+element1 +str(counter1) + "     OG({}) = {}".format(element2,addnox2))
                sols.append(element2+element1 +str(counter1) + "     OG({}) = {}".format(element2,addnox2))
            else:
                print(element2+str(counter2)+element1 +str(counter1) + "     OG({}) = {}".format(element2,addnox2))
                sols.append(element2+str(counter2)+element1 +str(counter1) + "     OG({}) = {}".format(element2,addnox2))
    elif En1 < En2:
        if element2 != "H":
            oxidation2 = abs(groep2 - 8)
        else:
            oxidation2 = abs(groep2 - 2)
        for x in oxnumbs1:
            counter1 = 1
            counter2 = 1
            oxidation1 = x
            addox1 = oxidation1 
            addnox2 = oxidation2
            while oxidation1 != oxidation2:
                if oxidation1 > oxidation2:
                    oxidation2 += addnox2
                    counter2 +=1
                else:
                    oxidation1 += addox1
                    counter1 +=1
                
            if counter1 == 1 and counter2 ==1:
                print(element1+element2)
                sols.append(element1+element2)
            elif counter1 == 1:
                print(element1 + element2+str(counter2) + "     OG({}) = {}".format(element1,addox1))
                sols.append(element1 + element2+str(counter2) + "     OG({}) = {}".format(element1,addox1))
            elif counter2 == 1:
                print(element1+str(counter1)+element2 + "     OG({}) = {}".format(element1,addox1))
                sols.append(element1+str(counter1)+element2 + "     OG({}) = {}".format(element1,addox1))
            else:
                print(element1+str(counter1)+element2 +str(counter2) + "     OG({}) = {}".format(element1,addox1))
                sols.append(element1+str(counter1)+element2 +str(counter2) + "     OG({}) = {}".format(element1,addox1))
    return sols
MoleculeStable("Cl", "S")
# if element1 in groepI:
#     En1 = groepI[element1][-1]        ## en worth
#     oxnumbs1 = groepI[element1][:-1]  ## oxidatie getal(len)
#     groep1 = 1                        ## group
# elif element1 in groepII:
#     En1 = groepII[element1][-1]
#     oxnumbs1 = groepII[element1][:-1]
#     groep1 = 1




# if element2 in groepI:
#     En2 = groepI[element2][-1]
#     oxnumbs2 = groepI[element2][:-1]
#     groep2 = 1
# elif element2 in groepII:
#     En2 = groepII[element2][-1]
#     oxnumbs2 = groepII[element2][:-1]
#     groep2 = 2
# print("ee")

# element1counter = 0
# element2counter = 0
# if En1 > En2:              ## first element is - oxidation
#     oxnumbs1 = groep1 - 8
#     while oxnumbs1 != oxnumbs2:
#         if abs(oxnumbs1) > oxnumbs2[0] and abs(oxnumbs1) != oxnumbs2[0]:
#          element1counter += 1
#          oxnumbs2[0] += oxnumbs2[0]
#         elif abs(oxnumbs1) < oxnumbs2[0] and abs(oxnumbs1) != oxnumbs2[0]:
#             element2counter += 1
#             oxnumbs1 += oxnumbs1
            
#         else:
#             print(element1, "###", element1counter)
#             print(element2, "###", element2counter)
#         print(oxnumbs2,abs(oxnumbs1) )
#         print(element1, "###", element1counter, "EEEE", element2, "###", element2counter)
        
# else:
#     oxidation2 = groep2 - 83




# if element1 in groepI:
#     En1 = groepI[element1][-1]        ## en worth
#     oxnumbs1 = groepI[element1][:-1]  ## oxidatie getal(len)
#     groep1 = 1                        ## group
# elif element1 in groepII:
#     En1 = groepII[element1][-1]
#     oxnumbs1 = groepII[element1][:-1]
#     groep1 = 1




# if element2 in groepI:
#     En2 = groepI[element2][-1]
#     oxnumbs2 = groepI[element2][:-1]
#     groep2 = 1
# elif element2 in groepII:
#     En2 = groepII[element2][-1]
#     oxnumbs2 = groepII[element2][:-1]
#     groep2 = 2
# print("ee")

# element1counter = 0
# element2counter = 0
# if En1 > En2:              ## first element is - oxidation
#     oxnumbs1 = groep1 - 8
#     while oxnumbs1 != oxnumbs2:
#         if abs(oxnumbs1) > oxnumbs2[0] and abs(oxnumbs1) != oxnumbs2[0]:
#          element1counter += 1
#          oxnumbs2[0] += oxnumbs2[0]
#         elif abs(oxnumbs1) < oxnumbs2[0] and abs(oxnumbs1) != oxnumbs2[0]:
#             element2counter += 1
#             oxnumbs1 += oxnumbs1
            
#         else:
#             print(element1, "###", element1counter)
#             print(element2, "###", element2counter)
#         print(oxnumbs2,abs(oxnumbs1) )
#         print(element1, "###", element1counter, "EEEE", element2, "###", element2counter)
        
# else:
#     oxidation2 = groep2 - 83





# if element1 in groepI:
#     En1 = groepI[element1][-1]        ## en worth
#     oxnumbs1 = groepI[element1][:-1]  ## oxidatie getal(len)
#     groep1 = 1                        ## group
# elif element1 in groepII:
#     En1 = groepII[element1][-1]
#     oxnumbs1 = groepII[element1][:-1]
#     groep1 = 1




# if element2 in groepI:
#     En2 = groepI[element2][-1]
#     oxnumbs2 = groepI[element2][:-1]
#     groep2 = 1
# elif element2 in groepII:
#     En2 = groepII[element2][-1]
#     oxnumbs2 = groepII[element2][:-1]
#     groep2 = 2
# print("ee")

# element1counter = 0
# element2counter = 0
# if En1 > En2:              ## first element is - oxidation
#     oxnumbs1 = groep1 - 8
#     while oxnumbs1 != oxnumbs2:
#         if abs(oxnumbs1) > oxnumbs2[0] and abs(oxnumbs1) != oxnumbs2[0]:
#          element1counter += 1
#          oxnumbs2[0] += oxnumbs2[0]
#         elif abs(oxnumbs1) < oxnumbs2[0] and abs(oxnumbs1) != oxnumbs2[0]:
#             element2counter += 1
#             oxnumbs1 += oxnumbs1
            
#         else:
#             print(element1, "###", element1counter)
#             print(element2, "###", element2counter)
#         print(oxnumbs2,abs(oxnumbs1) )
#         print(element1, "###", element1counter, "EEEE", element2, "###", element2counter)
        
# else:
#     oxidation2 = groep2 - 83


