import math
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import pubchempy as pcp
import requests
from bs4 import BeautifulSoup
from PIL import Image
from io import BytesIO
def MoleculeStable(element1,element2):
    groupI = {"H": [1, 2.1], "Li": [1, 1.0], "Na":[1,0.9], "K":[1,0.8],"Rb":[1,0.8],"Cs":[1,0.7],"Fr":[1,0.7]}        ## element name, all possible positive oxidiation numbers, and then final number is the EN worth
    groupII = {"Be": [2, 1.5], "Mg": [2, 1.2], "Ca":[2,1.0],"Sr":[2,1.0],"Ba":[2,0.9],"Ra":[2,0.9]}
    groupIII = {"B": [3, 2.0], "Al":[3,1.5],"Ga":[3,1.6],"In":[3,1.7],"Ti":[3,1.8]}
    groupIV = {"C": [2,4, 2.5],"Si": [4, 1.8],"Ge": [4, 1.8], "Sn": [2,4, 1.8], "Pb": [2,4, 1.8]}
    groupV = {"N": [1,3,5, 3.0],"P": [3,5, 2.1], "As": [5, 2.0], "Sb": [4, 1.9], "Bi": [4, 1.8]}
    groupVI = {"O": [6, 3.5],"S": [4,6, 2.5],"Se": [4,6, 2.4],"Te": [6, 2.1],"Po": [6, 2.0]}
    groupVII = {"F": [7, 4.0], "Cl": [1,3,5,7, 3.0], "Br": [1,3,5,7, 2.8],"I": [1,3,5,7, 2.5],"At": [7, 2.2]}
    groups =  [groupI,groupII,groupIII,groupIV,groupV,groupVI, groupVII]
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
        oxidation1Placeholder = oxidation1
        
        for x in oxnumbs2:
            oxidation1 = oxidation1Placeholder
            addox1 = oxidation1 

            counter1 = 1
            counter2 = 1
            oxidation2 = x
            addnox2 = oxidation2
            while oxidation1 != oxidation2:
                print(oxidation1,oxidation2)
                if oxidation1 > oxidation2:
                    oxidation2 += addnox2
                    counter2 +=1
                else:
                    oxidation1 += addox1
                    counter1 +=1
            print(oxidation1,oxidation2)
            if counter1 == 1 and counter2 ==1:
                print(element2+element1 +"     OG({}) = {}".format(element2,addnox2))
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
        oxidation2Placeholder = oxidation2
        addnox2 = oxidation2Placeholder
        for x in oxnumbs1:
            counter1 = 1
            counter2 = 1
            oxidation2 = oxidation2Placeholder
            oxidation1 = x
            addox1 = oxidation1 
            while oxidation1 != oxidation2:
                if oxidation1 > oxidation2:
                    oxidation2 += addnox2
                    counter2 +=1
                else:
                    oxidation1 += addox1
                    counter1 +=1
                print(oxidation1,oxidation2)
            if counter1 == 1 and counter2 ==1:
                print(element1+element2)
                sols.append(element1+element2+"     OG({}) = {}".format(element1,addox1))
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
def MoleculeVisualizeViaBruto(bruto):
    compound = pcp.get_compounds(bruto, "formula")[0]     ## searches the Pubchem webpage and grabs the first compund, returns it as compound(ID)
    try:
        smiles = compound.canonical_smiles

        mol = Chem.MolFromSmiles(smiles)
        print(mol)
        img = Draw.MolToImage(mol)
        return img
        #Draw.MolToFile(img, 'saved.png')
    except:
        print(compound)
        # Define the compound ID
        base_url = 'https://pubchem.ncbi.nlm.nih.gov'
        start_index = str(compound).find('(') + 1
        end_index = str(compound).find(')')
        compound_id = str(compound)[start_index:end_index] 
        

        # Construct the full image URL
        img_url = f'{base_url}/image/imgsrv.fcgi?cid={compound_id}&t=l'
        response = requests.get(img_url)

        if response.status_code == 200:
            # Open the image from bytes and display it
            img = Image.open(BytesIO(response.content))
            #img.show()  # Opens the image in default image viewer
            #img.save('saved.png')  # Save the image locally
            return(img)
        else:
            print("Failed to fetch the structure image")
MoleculeStable("Cl", "O")
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


