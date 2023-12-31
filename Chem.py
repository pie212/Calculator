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
    groupVI = {"O": [6, 3.5],"S": [2,4,6, 2.5],"Se": [4,6, 2.4],"Te": [6, 2.1],"Po": [6, 2.0]}  ## quick bit of thinking i should think of. So we know that S en F leads to SF2 and SF6 because we have the structure of -        which is 6 total electrons on the shell, 4 valence ones. In SF2 F bonds with the free electrons but since it has 7 electrons they each share 1 electron, causing the need for 2 F atoms to be present.                                                     *
    groupVII = {"F": [7, 4.0], "Cl": [1,3,5,7, 3.0], "Br": [1,3,5,7, 2.8],"I": [1,3,5,7, 2.5],"At": [7, 2.2]}                                                                                    #                     *    *     However with SF6 the 2 valence pairs break, allowing the 6 electrons to bond with 6 Fluor electrons, so the question I want to ask in a week is, can we predict this for all elements? And why could we not keep 1 pair and get SF4, we would then have  *   *    -->> but in SF6 S would then have 12 valence electrons, not 6... HOW???
    groups =  [groupI,groupII,groupIII,groupIV,groupV,groupVI, groupVII]                                                                                                                         #                       -                                                                                                                                                                                                                                                                   -        ## how can C have -2??? 2 pictures trying to figure that out on my desktop
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

    ###########3 EXCEPTIONS ##########3
    if element1 == "Cl" and element2 == "Cl":
        sols.append("Cl2") 
    elif En1 > En2:
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
    
    ## known exceptions
    if bruto == "SO3":
        compound = pcp.get_compounds(bruto, "formula")[1]     ## searches the Pubchem webpage and grabs the first compund, returns it as compound(ID)
    else:
        try:
            compound = pcp.get_compounds(bruto, "formula")[0]     ## searches the Pubchem webpage and grabs the first compund, returns it as compound(ID)
        except:
            compound = ""
    try:
        smiles = compound.canonical_smiles

        mol = Chem.MolFromSmiles(smiles)
        print(mol)
        img = Draw.MolToImage(mol)
        return img
        #Draw.MolToFile(img, 'saved.png')
    except:
        try:
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
                img = Image.open("error.png")
                return(img)
        except:
            img = Image.open("error.png")
            return(img)
                

def MoleculeVisualizeViaWebsite(smarts):
    # compound = pcp.get_compounds(smiles, "smiles")[0]     ## searches the Pubchem webpage and grabs the first compund, returns it as compound(ID)
    mol = Chem.MolFromSmarts(smarts)
    # smiles = Chem.MolToSmiles(mol)
    #mol = Chem.MolFromSmiles(smiles)
    print(mol)
    img = Draw.MolToImage(mol)
    return(img)
def structural_to_smiles(structural_formula):    ## hasnt been used.... dont think we need it
    # Split the formula into parts (atoms and bonds)
    parts = structural_formula.split('-')

    # Convert to SMARTS notation
    smarts_formula = ''
    for part in parts:
        if part.startswith('C'):
            # Append the carbon atom
            smarts_formula += 'C'
        elif part.startswith('H'):
            # Skip hydrogen atoms in this representation
            pass
        elif part == '':
            # Append the bond
            smarts_formula += '-'

    return smarts_formula
def CandHinator(carbon,hydrogen,oxygen):         ### ONLY ONE DOUBLE OR TRIPLE BOND!!!
    print("test")          ### single bond C-C     Double bond C=C   triple bond C#C
    if carbon > 0 and hydrogen > 0:
        if oxygen == 0:
            ##### alkanen alkenen en alkynen
            ## alkaan
            if hydrogen == ((carbon * 2)+2):
                structure = []
                strucutrealkaan = ""
                alkaaanvisual = ""
                alkaanname = ""
                alkaanbruto = ""
                for x in range(carbon):
                    if x == 0:
                        structure.append("C")
                    elif x == (carbon-1):
                        structure.append("C")
                    else:
                        structure.append("C")
                    strucutrealkaan = ""
                    for x in range(len(structure)):
                        if x == 0:
                            strucutrealkaan += structure[x]
                        else:
                            strucutrealkaan +="-" 
                            strucutrealkaan += structure[x]
                # print(strucutrealkaan)
                # smarts = structural_to_smarts(strucutrealkaan)
                # print(smarts)
                # mol = Chem.MolFromSmarts(smarts)
                # smiles = Chem.MolToSmiles(mol)
                smiles = StructureToSmarts(strucutrealkaan)
                compound = pcp.get_compounds(smiles, "smiles")    ## SMARTS notation works when searching for SMILES... for some reason....
                start_index = str(compound).find('(') + 1        ## format the compound
                end_index = str(compound).find(')')                 ## format the compound
                compound_id = str(compound)[start_index:end_index]      ## format the compound
                alkaanname = pcp.Compound.from_cid(compound_id).iupac_name       ## compound name
                alkaanbruto = pcp.Compound.from_cid(compound_id).molecular_formula
                print("######")
                print(smiles)
                alkaaanvisual = MoleculeVisualizeViaWebsite(smiles)
            else:
                print("No Alkanen for this bruto")
                strucutrealkaan = "None"
                alkaaanvisual = ""
                alkaanname = ""
                alkaanbruto = ""
            if hydrogen == ((carbon * 2)):      ## with 1 double bond for tmrw
                bonds = carbon/2
                try:
                    bonds = int(bonds)
                except:
                    bonds = bonds + 0.5
                    bonds = int(bonds)
                print(bonds)
                alkenenFormula = []
                alkenenStructure = []
                alkenenNames = []
                alkenenBrutos = []
                for x in range(bonds):
                    structurealkeen = ""
                    for y in range(carbon - 1):
                        if y == 0:
                            structurealkeen += "C"
                        if y == x:
                            structurealkeen += "=C"
                        else:
                            structurealkeen += "-C"
                    smiles = StructureToSmarts(structurealkeen)        ###### converts from for ex C=C to [C;X2]=[C;X2]
                    compound = pcp.get_compounds(smiles, "smiles")    ## SMARTS notation works when searching for SMILES... for some reason....
                    start_index = str(compound).find('(') + 1        ## format the compound
                    end_index = str(compound).find(')')                 ## format the compound
                    compound_id = str(compound)[start_index:end_index]      ## format the compound
                    CompoundInfo = pcp.Compound.from_cid(compound_id)       ## compound name
                    alkeen_name = CompoundInfo.iupac_name
                    alkeen_molecular_formula = CompoundInfo.molecular_formula
                    print(alkeen_molecular_formula)
                    print("######")
                    print(smiles)
                    alkenenFormula.append(structurealkeen)
                    alkenenStructure.append(MoleculeVisualizeViaWebsite(smiles))
                    alkenenNames.append(alkeen_name)
                    alkenenBrutos.append(alkeen_molecular_formula)
            else:
                print("No Alkenen for this bruto")
                alkenenFormula = []
                alkenenStructure = []      
                alkenenNames = []
                alkenenBrutos = []
            
            if hydrogen == ((carbon * 2)-2):
                bonds = carbon/2
                try:
                    bonds = int(bonds)
                except:
                    bonds = bonds + 0.5
                    bonds = int(bonds)
                print(bonds)
                alkynenFormula = []
                alkynenStructure = []
                alkynenNames = []
                alkynenBrutos = []
                for x in range(bonds):
                    structurealkyn = ""
                    for y in range(carbon - 1):
                        if y == 0:
                            structurealkyn += "C"
                        if y == x:
                            structurealkyn += "#C"
                        else:
                            structurealkyn += "-C"
                    smiles = StructureToSmarts(structurealkyn)        ###### converts from for ex C=C to [C;X2]=[C;X2]
                    compound = pcp.get_compounds(smiles, "smiles")    ## SMARTS notation works when searching for SMILES... for some reason....
                    start_index = str(compound).find('(') + 1        ## format the compound
                    end_index = str(compound).find(')')                 ## format the compound
                    compound_id = str(compound)[start_index:end_index]      ## format the compound
                    CompoundInfo = pcp.Compound.from_cid(compound_id)       ## compound name
                    alkyn_name = CompoundInfo.iupac_name
                    alkyn_molecular_formula = CompoundInfo.molecular_formula
                    print(alkyn_molecular_formula)
                    print("######")
                    print(smiles)
                    alkynenFormula.append(structurealkyn)
                    alkynenStructure.append(MoleculeVisualizeViaWebsite(smiles))
                    alkynenNames.append(alkyn_name)
                    alkynenBrutos.append(alkyn_molecular_formula)
            else:
                print("No Alkynen for this bruto")
                alkynenFormula = []
                alkynenStructure = []      
                alkynenNames = []
                alkynenBrutos = []
            return(strucutrealkaan, alkaaanvisual,alkaanname, alkaanbruto ,          alkenenFormula,alkenenStructure, alkenenNames, alkenenBrutos, alkynenFormula,alkynenStructure, alkynenNames, alkynenBrutos)   ##string, photo, list[strings], list[photos]      ### Structure (C-C), PHOTO , NAME, Structure(C=C). PHOTO, NAME   
        elif oxygen > 0:
            print("e")
    else:
        return("None")
def CandHguesser(carbon):                        ### ONLY ONE DOUBLE OR TRIPLE BOND!!!
    HydrogenAlkanen = ((carbon * 2)+2)
    HydrogenAlkenen = carbon * 2
    HydrogenAlkynen = (carbon *2) -2
    CandHguesserls = CandHinator(carbon,HydrogenAlkanen, 0)
    print("###########")
    strucutrealkaan, alkaaanvisual,alkaanname, alkaanbruto = CandHguesserls[0],CandHguesserls[1],CandHguesserls[2],CandHguesserls[3] ## structure, photo, name
    CandHguesserls = CandHinator(carbon,HydrogenAlkenen,0)
    alkenenFormula,alkenenStructure, alkenenNames,alkenenBrutos = CandHguesserls[4], CandHguesserls[5], CandHguesserls[6], CandHguesserls[7]   ## structure, photo, names
    CandHguesserls = CandHinator(carbon,HydrogenAlkynen,0)
    alkynenFormula,alkynenStructure, alkynenNames, alkynenBrutos = CandHguesserls[8], CandHguesserls[9], CandHguesserls[10], CandHguesserls[11]
    print(strucutrealkaan,alkaanname,alkaanbruto )
    print(alkenenFormula, alkenenNames, alkenenBrutos)
    return(strucutrealkaan, alkaaanvisual,alkaanname, alkaanbruto,               alkenenFormula,alkenenStructure, alkenenNames,alkenenBrutos     ,         alkynenFormula,alkynenStructure, alkynenNames, alkynenBrutos  ) ### Structure (C-C), PHOTO , NAME, Structure(C=C). PHOTO, NAME
def StructureToSmarts(structure_formula):        ## smiles is for ex: C=C, C-C-C=C  ## SMARTS is for ex [#6]
    mol = Chem.MolFromSmiles(structure_formula) 

    # Convert the molecule to SMARTS
    smarts_pattern = Chem.MolToSmarts(mol)         ## Smarts is [C;X2]=[C;X2] for C=C

    return(smarts_pattern)  # This should print '[C;X2]=[C;X2]' just kidding, for C=C we get [#6]=[#6]

#CandHguesser(6)
#CandHinator(2,2,0)



#MoleculeStable("Cl", "O")
#CandHinator(2,6,0)



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


