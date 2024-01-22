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
import re
from itertools import product
import time
def MoleculeStable(element1,element2):
    groupI = {"H": [1, 2.1], "Li": [1, 1.0], "Na":[1,0.9], "K":[1,0.8],"Rb":[1,0.8],"Cs":[1,0.7],"Fr":[1,0.7]}        ## element name, all possible positive oxidiation numbers, and then final number is the EN worth
    groupII = {"Be": [2, 1.5], "Mg": [2, 1.2], "Ca":[2,1.0],"Sr":[2,1.0],"Ba":[2,0.9],"Ra":[2,0.9]}
    groupIII = {"B": [3, 2.0], "Al":[3,1.5],"Ga":[3,1.6],"In":[3,1.7],"Ti":[3,1.8]}
    groupIV = {"C": [2,4, 2.5],"Si": [4, 1.8],"Ge": [4, 1.8], "Sn": [2,4, 1.8], "Pb": [2,4, 1.8]}
    groupV = {"N": [1,3,5, 3.0],"P": [3,5, 2.1], "As": [5, 2.0], "Sb": [4, 1.9], "Bi": [4, 1.8]}
    groupVI = {"O": [1,6, 3.5],"S": [2,4,6, 2.5],"Se": [4,6, 2.4],"Te": [6, 2.1],"Po": [6, 2.0]}  ## quick bit of thinking i should think of. So we know that S en F leads to SF2 and SF6 because we have the structure of -        which is 6 total electrons on the shell, 4 valence ones. In SF2 F bonds with the free electrons but since it has 7 electrons they each share 1 electron, causing the need for 2 F atoms to be present.                                                     *
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
        #print(mol)
        img = Draw.MolToImage(mol)
        return img
        #Draw.MolToFile(img, 'saved.png')
    except:
        try:
            #print(compound)
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
                #print("Failed to fetch the structure image")
                img = Image.open("error.png")
                return(img)
        except:
            img = Image.open("error.png")
            return(img)
                
def IdToPhotoViaWebsite(compound_id):
    try:
        #print(compound_id)
        # Define the compound ID
        base_url = 'https://pubchem.ncbi.nlm.nih.gov'
         
            

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

def MoleculeVisualizeViaSmarts(smarts):
    # compound = pcp.get_compounds(smiles, "smiles")[0]     ## searches the Pubchem webpage and grabs the first compund, returns it as compound(ID)
    mol = Chem.MolFromSmarts(smarts)
    # smiles = Chem.MolToSmiles(mol)
    #mol = Chem.MolFromSmiles(smiles)
    #print(mol)
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
            if hydrogen == ((carbon * 2)+2):       ## alkaan
                structure = []                 ## define all our variables
                strucutrealkaan = ""            ## define all our variables
                alkaaanvisual = ""              ## define all our variables
                alkaanname = ""                 ## define all our variables
                alkaanbruto = ""                ## define all our variables
                print("##################### ALKAAN #################")
                for x in range(carbon):         
                    # if x == 0:
                    #     structure.append("CH3")
                    # elif x == (carbon-1):            ## all this blurred out stuff was if if there were H atoms, but i figured out that that didnt work.
                    #     structure.append("CH3")
                    # else:
                    structure.append("C")              ## adds all our C atoms
                    strucutrealkaan = ""    
                    for x in range(len(structure)):
                        if x == 0:
                            strucutrealkaan += structure[x]
                        else:
                            strucutrealkaan +="-"                               ## builds our Alkaan, C-C-C
                            strucutrealkaan += structure[x]
                # print(strucutrealkaan)
                # smarts = structural_to_smarts(strucutrealkaan)
                # print(smarts)
                # mol = Chem.MolFromSmarts(smarts)
                # smiles = Chem.MolToSmiles(mol)
                smarts = StructureToSmarts(strucutrealkaan)         ## converts from C-C-C --> [#6]-[#6]-
                compound = pcp.get_compounds(smarts, "smiles")    ## SMARTS notation works when searching for SMILES... for some reason....
                start_index = str(compound).find('(') + 1        ## format the compound
                end_index = str(compound).find(')')                 ## format the compound
                compound_id = str(compound)[start_index:end_index]      ## format the compound
                alkaanname = pcp.Compound.from_cid(compound_id).iupac_name       ## compound name
                alkaanbruto = pcp.Compound.from_cid(compound_id).molecular_formula
                print("Compound ID: " +str(compound_id))
                print("SMARTS: " + str(smarts))
                print("Structure: " + structurealkeen)
                print("##")
                alkaaanvisual = MoleculeVisualizeViaSmarts(smarts)
            else:
                print("##################### ALKAAN #################")
                print("Compound ID: None")
                print("SMARTS: None")
                print("Structure: None")
                print("##")
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
                print("####################################### ALKEEN #######################################")
                for x in range(bonds):
                    structurealkeen = ""
                    for y in range(carbon - 1):
                        if y == 0:
                            structurealkeen += "C"
                        if y == x:
                            structurealkeen += "=C"
                        else:
                            structurealkeen += "-C"

                    smarts = StructureToSmarts(structurealkeen)        ###### converts from for ex C=C to [C;X2]=[C;X2]
                    compound = pcp.get_compounds(smarts, "smiles")    ## SMARTS notation works when searching for SMILES... for some reason....
                    start_index = str(compound).find('(') + 1        ## format the compound
                    end_index = str(compound).find(')')                 ## format the compound
                    compound_id = str(compound)[start_index:end_index]      ## format the compound
                    CompoundInfo = pcp.Compound.from_cid(compound_id)       ## compound name
                    alkeen_name = CompoundInfo.iupac_name
                    alkeen_molecular_formula = CompoundInfo.molecular_formula
                    print("Compound ID: " +str(compound_id))
                    print("SMARTS: " + str(smarts))
                    print("Structure: " + structurealkeen)
                    print("##")
                    alkenenFormula.append(structurealkeen)
                    alkenenStructure.append( MoleculeVisualizeViaSmarts(smarts))
                    alkenenNames.append(alkeen_name)
                    alkenenBrutos.append(alkeen_molecular_formula)

                
            else:
                print("####################################### ALKEEN #######################################")
                print("Compound ID: None")
                print("SMARTS: None")
                print("Structure: None")
                print("##")
                alkenenFormula = []
                alkenenStructure = []      
                alkenenNames = []
                alkenenBrutos = []
            
            if hydrogen == ((carbon * 2)-2):      ## alyn # triple bond
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
                print("####################################### ALKYN #######################################")
                for x in range(bonds):
                    structurealkyn = ""
                    for y in range(carbon - 1):
                        if y == 0:
                            structurealkyn += "C"
                        if y == x:
                            structurealkyn += "#C"
                        else:
                            structurealkyn += "-C"
                    smarts = StructureToSmarts(structurealkyn)        ###### converts from for ex C=C to [C;X2]=[C;X2]
                    compound = pcp.get_compounds(smarts, "smiles")    ## SMARTS notation works when searching for SMILES... for some reason....
                    start_index = str(compound).find('(') + 1        ## format the compound
                    end_index = str(compound).find(')')                 ## format the compound
                    compound_id = str(compound)[start_index:end_index]      ## format the compound
                    CompoundInfo = pcp.Compound.from_cid(compound_id)       ## compound name
                    alkyn_name = CompoundInfo.iupac_name
                    alkyn_molecular_formula = CompoundInfo.molecular_formula
                    print("Compound ID: " +str(compound_id))
                    print("SMARTS: " + str(smarts))
                    print("Structure: " + structurealkyn)
                    print("##")
                    alkynenFormula.append(structurealkyn)
                    alkynenStructure.append(IdToPhotoViaWebsite(compound_id))
                    alkynenNames.append(alkyn_name)
                    alkynenBrutos.append(alkyn_molecular_formula)
                # bonds = carbon - 1    ## we do carbon - 1 because the amount of bonding places is equal to that.... for example if we have 5 carbon atoms we only have 4 bonding places C=C=C=C
                # counter = 0           ## to count the amount of times we need to add a double bond
                # for x in range(bonds - 2): ## we do - 2 because we have already started off with the first carbon atom as "C" and we iterate from 0
                #     structurealkyn = "C"
                #     counter += 1          ## adds 1 to the counter
                #     for y in range(carbon - 2): ## since we already have the first C atom
                #         if y <= counter:         ## if the counter is bigger or equal to Y we add a double bond
                #             structurealkyn += "#C"
                #         else: 
                #             structurealkyn += "-C"  ## else we add a single bond
                #     print(structurealkyn)
                #     smarts = StructureToSmarts(structurealkyn)        ###### converts from for ex C=C to [C;X2]=[C;X2]
                #     compound = pcp.get_compounds(smarts, "smiles")    ## SMARTS notation works when searching for SMILES... for some reason....
                #     start_index = str(compound).find('(') + 1        ## format the compound     
                #     end_index = str(compound).find(')')                 ## format the compound
                #     compound_id = str(compound)[start_index:end_index]      ## format the compound
                #     CompoundInfo = pcp.Compound.from_cid(compound_id)       ## info of the compound using the ID
                #     alkyn_name = CompoundInfo.iupac_name           ## gets the ipuac name
                #     alkyn_molecular_formula = CompoundInfo.molecular_formula
                #     print(compound_id)
                #     print("######")
                #     print(smarts)
                #     alkynenFormula.append(structurealkyn)
                #     alkynenStructure.append(IdToPhotoViaWebsite(compound_id))
                #     alkynenNames.append(alkyn_name)
                #     alkynenBrutos.append(alkyn_molecular_formula)
            else:
                print("####################################### ALKYN #######################################")
                print("Compound ID: None")
                print("SMARTS: None")
                print("Structure: None")
                alkynenFormula = []
                alkynenStructure = []      
                alkynenNames = []
                alkynenBrutos = []
            return(strucutrealkaan, alkaaanvisual,alkaanname, alkaanbruto ,          alkenenFormula,alkenenStructure, alkenenNames, alkenenBrutos, alkynenFormula,alkynenStructure, alkynenNames, alkynenBrutos)   ##string, photo, list[strings], list[photos]      ### Structure (C-C), PHOTO , NAME, Structure(C=C). PHOTO, NAME   
        elif oxygen > 0:        ## OH and COOH
            if hydrogen == ((carbon * 2)+2):       ## alkaan
                alcoholvisuals = []
                alcoholStructure = []
                alcoholNames = []
                alcoholBrutos = []


                alcohollist = []                  ## define all our variables
                alcohol = ""            ## define all our variables
                aalcoholvisual = ""              ## define all our variables
                alcoholname = ""                 ## define all our variables
                alcoholbruto = ""                ## define all our variables
                bonds = carbon/2          ## we are not checking bonding places but rather amount of Carbon atoms without being recursive, C-C-C is 2 because 1 and 2 and C-C-C-C-C is 3 because 4 is 2 and 5 is 1
                try:
                    bonds = int(bonds)
                except:
                    bonds = bonds + 0.5
                    bonds = int(bonds)
                for x in range(carbon):         
                    # if x == 0:
                    #     structure.append("CH3")
                    # elif x == (carbon-1):            ## all this blurred out stuff was if if there were H atoms, but i figured out that that didnt work.
                    #     structure.append("CH3")
                    # else:
                    alcohollist.append("C")              ## adds all our C atoms                                        
                counter = 0  ## use to maninulate name, for example 1-butanol --> 2-butanol
                for x in range(bonds):
                
                    counter += 1
                    #alcohollist[len(alcohollist) - (x+1)] = "O"
                    # if x != 0:
                    #     alcohollist[len(alcohollist) - x] = "C"
                    for y in range(len(alcohollist)):
                        alcohol += alcohollist[x]
                    alcohol += "O"
                    #smarts = StructureToSmarts(alcohol)         ## converts from C-C-C --> [#6]-[#6]-
                    
                    if x == 0:      ## if this is the first time running the loop, meaning the normal like C-C-OH then the smiles SHOULD work and return the name, for ex 1-butanol

                        compound = pcp.get_compounds(alcohol, "smiles")    ## SMARTS notation works when searching for SMILES... for some reason....
                        
                    else:                         ### if its not the first loop the SMARTS notation does not work for some odd reason and returns something wonky, so the next 3 lines of code replace the 1-butanol --> 2-butanol and the counter counts the amounts of loops to make sure the right number is added
                        alcoholname = alcoholname.replace(str(x), str(counter))
                        compound = pcp.get_compounds(alcoholname, "name")       ## searches the compound by using the name made ^
                    start_index = str(compound).find('(') + 1        ## format the compound
                    end_index = str(compound).find(')')                 ## format the compound
                    compound_id = str(compound)[start_index:end_index]      ## format the compound
                    alcoholname = pcp.Compound.from_cid(compound_id).iupac_name       ## compound name
                    alcoholbruto = pcp.Compound.from_cid(compound_id).molecular_formula
                    
                    alcoholvisual = IdToPhotoViaWebsite(compound_id)   ## gets photo
                    alcoholStructure.append("none at the moment")      ## gets structure
                    alcoholvisuals.append(alcoholvisual)               ## adds the photo
                    alcoholNames.append(alcoholname)                   ## adds the name
                    alcoholBrutos.append(alcoholbruto)                 ## adds the bruto formula

                    alcohol = ""                                       ## resets the alcohol structure to nothing, no point since it will always be set back to C*XO (meaning O will always be at the end) CCO or CCCCCO etc
                return(alcoholStructure,alcoholvisuals,alcoholNames,alcoholBrutos)
            else:
                print("No Alcoholen for this bruto")
                alcoholStructure = "None"
                alcoholvisuals = ""
                alcoholNames = ""
                alcoholBrutos = ""
                return(alcoholStructure,alcoholvisuals,alcoholNames,alcoholBrutos)
            
    else:
        return("None")

    



def StructureToSmarts(structure_formula):        ## smiles is for ex: C=C, C-C-C=C  ## SMARTS is for ex [#6]

    mol = Chem.MolFromSmiles(structure_formula) 

    # Convert the molecule to SMARTS
    smarts_pattern = Chem.MolToSmarts(mol)         ## Smarts is [C;X2]=[C;X2] for C=C

    return(smarts_pattern)  # This should print '[C;X2]=[C;X2]' just kidding, for C=C we get [#6]=[#6]
def carboninator(carbon):         ## GENERATES all possible C chains with only H 
    # print("test")          ### single bond C-C     Double bond C=C   triple bond C#C
    if carbon > 0 :
        
        ##### alkanen alkenen en alkynen
        ## alkaan
        
        structure = []                 ## define all our variables
        strucutrealkaan = ""            ## define all our variables
        alkaaanvisual = ""              ## define all our variables
        alkaanname = ""                 ## define all our variables
        alkaanbruto = ""                ## define all our variables
        print("##################### ALKAAN #################")
        for x in range(carbon):         
            # if x == 0:
            #     structure.append("CH3")
            # elif x == (carbon-1):            ## all this blurred out stuff was if if there were H atoms, but i figured out that that didnt work.
            #     structure.append("CH3")
            # else:
            structure.append("C")              ## adds all our C atoms
            strucutrealkaan = ""    
            for x in range(len(structure)):
                if x == 0:
                    strucutrealkaan += structure[x]
                else:
                    strucutrealkaan +="-"                               ## builds our Alkaan, C-C-C
                    strucutrealkaan += structure[x]
        # print(strucutrealkaan)
        # smarts = structural_to_smarts(strucutrealkaan)
        # print(smarts)
        # mol = Chem.MolFromSmarts(smarts)
        # smiles = Chem.MolToSmiles(mol)
        smarts = StructureToSmarts(strucutrealkaan)         ## converts from C-C-C --> [#6]-[#6]-
        compound = pcp.get_compounds(smarts, "smiles")    ## SMARTS notation works when searching for SMILES... for some reason....
        start_index = str(compound).find('(') + 1        ## format the compound
        end_index = str(compound).find(')')                 ## format the compound
        compound_id = str(compound)[start_index:end_index]      ## format the compound
        alkaanname = pcp.Compound.from_cid(compound_id).iupac_name       ## compound name
        alkaanbruto = pcp.Compound.from_cid(compound_id).molecular_formula
        print("Compound ID: " +str(compound_id))
        print("SMARTS: " + str(smarts))
        print("Structure: " + strucutrealkaan)
        print("##")
        alkaaanvisual = MoleculeVisualizeViaSmarts(smarts)

        ##################### ALKEEN #################
        if carbon > 1:
            bonds = carbon/2
            try:
                bonds = int(bonds)
            except:
                bonds = bonds + 0.5
                bonds = int(bonds)
            # print(bonds)
            alkenenFormula = []
            alkenenStructure = []
            alkenenNames = []
            alkenenBrutos = []
            print("##################### ALKEEN #################")
            for x in range(bonds):
                structurealkeen = ""
                for y in range(carbon - 1):
                    if y == 0:
                        structurealkeen += "C"
                    if y == x:
                        structurealkeen += "=C"
                    else:
                        structurealkeen += "-C"

                smarts = StructureToSmarts(structurealkeen)        ###### converts from for ex C=C to [C;X2]=[C;X2]
                compound = pcp.get_compounds(smarts, "smiles")    ## SMARTS notation works when searching for SMILES... for some reason....
                start_index = str(compound).find('(') + 1        ## format the compound
                end_index = str(compound).find(')')                 ## format the compound
                compound_id = str(compound)[start_index:end_index]      ## format the compound
                CompoundInfo = pcp.Compound.from_cid(compound_id)       ## compound name
                alkeen_name = CompoundInfo.iupac_name
                alkeen_molecular_formula = CompoundInfo.molecular_formula
                print("Compound ID: " +str(compound_id))
                print("SMARTS: " + str(smarts))
                print("Structure: " + structurealkeen)
                print("##")
                alkenenFormula.append(structurealkeen)
                alkenenStructure.append( MoleculeVisualizeViaSmarts(smarts))
                alkenenNames.append(alkeen_name)
                alkenenBrutos.append(alkeen_molecular_formula)

            ##### THIS DOES NOT GET CALLED WHEN THE HYDROGEN IS UNKNOWN, THIS ONLY GETS CALLED 
            bonds = carbon - 1    ## we do carbon - 1 because the amount of bonding places is equal to that.... for example if we have 5 carbon atoms we only have 4 bonding places C=C=C=C
            counter = 0           ## to count the amount of times we need to add a double bond
            for x in range(bonds - 2): ## we do - 2 because we have already started off with the first carbon atom as "C" and we iterate from 0
                structurealkeen = "C"
                counter += 1          ## adds 1 to the counter
                for y in range(carbon - 2): ## since we already have the first C atom
                    if y <= counter:         ## if the counter is bigger or equal to Y we add a double bond
                        structurealkeen += "=C"
                    else: 
                        structurealkeen += "-C"  ## else we add a single bond
                smarts = StructureToSmarts(structurealkeen)        ###### converts from for ex C=C to [C;X2]=[C;X2]
                try:
                    compound = pcp.get_compounds(smarts, "smiles")    ## SMARTS notation works when searching for SMILES... for some reason....
                    start_index = str(compound).find('(') + 1        ## format the compound
                    end_index = str(compound).find(')')                 ## format the compound
                    compound_id = str(compound)[start_index:end_index]      ## format the compound
                    CompoundInfo = pcp.Compound.from_cid(compound_id)       ## compound name
                    alkeen_name = CompoundInfo.iupac_name
                    alkeen_molecular_formula = CompoundInfo.molecular_formula
                    alkenenFormula.append(structurealkeen)
                    alkenenStructure.append(IdToPhotoViaWebsite(compound_id))       ## better to not draw these ones, as they are more complicated...
                    alkenenNames.append(alkeen_name)
                    alkenenBrutos.append(alkeen_molecular_formula) 
                    print("Compound ID: " +str(compound_id))
                    print("SMARTS: " + str(smarts))
                    print("Structure: " + structurealkeen)
                    print("##")
                except:
                    print("Compound ID: None ")
                    print("SMARTS: " + str(smarts))
                    print("Structure: " + structurealkeen)
                    print("##")
                    alkenenFormula.append(structurealkeen)
                    alkenenStructure.append(Image.open("NotExist.png"))       ## better to not draw these ones, as they are more complicated...
                    alkenenNames.append("Does Not Exist")
                    alkenenBrutos.append("Unknown") 
        
        ####################################### ALKYN #######################################
        if carbon > 1:
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
            print("####################################### ALKYN #######################################")
            for x in range(bonds):
                structurealkyn = ""
                for y in range(carbon - 1):
                    if y == 0:
                        structurealkyn += "C"
                    if y == x:
                        structurealkyn += "#C"
                    else:
                        structurealkyn += "-C"
                smarts = StructureToSmarts(structurealkyn)        ###### converts from for ex C=C to [C;X2]=[C;X2]
                compound = pcp.get_compounds(smarts, "smiles")    ## SMARTS notation works when searching for SMILES... for some reason....
                start_index = str(compound).find('(') + 1        ## format the compound
                end_index = str(compound).find(')')                 ## format the compound
                compound_id = str(compound)[start_index:end_index]      ## format the compound
                CompoundInfo = pcp.Compound.from_cid(compound_id)       ## compound name
                alkyn_name = CompoundInfo.iupac_name
                alkyn_molecular_formula = CompoundInfo.molecular_formula
                print("Compound ID: " +str(compound_id))
                print("SMARTS: " + str(smarts))
                print("Structure: " + structurealkyn)
                print("##")
                alkynenFormula.append(structurealkyn)
                alkynenStructure.append(IdToPhotoViaWebsite(compound_id))
                alkynenNames.append(alkyn_name)
                alkynenBrutos.append(alkyn_molecular_formula)
    return(strucutrealkaan, alkaaanvisual,alkaanname, alkaanbruto ,          alkenenFormula,alkenenStructure, alkenenNames, alkenenBrutos, alkynenFormula,alkynenStructure, alkynenNames, alkynenBrutos)   ##string, photo, list[strings], list[photos]      ### Structure (C-C), PHOTO , NAME, Structure(C=C). PHOTO, NAME
def balancer(formula, max_value,return_array):
    t1 = time.time()
    formula = formula#"Na + H2O --> NaOH + H2"
    print("working on it")
    #formula = "L + H2O +  --> LOH + H2 "
    try:
        seperated = formula.split("-->") ## splits list into 2 items, ["Na + H2O"] ["NaOH"] ( before and after the paranthesis)
        reactants = []
        products = []
        for x in range(len(seperated[0].split("+"))):       ## for amount of items before the -->, 1 is after and 0 is before since seperated splits into before and after the arrow
            reactants.append(seperated[0].split("+")[x].replace(" ", ""))     ## for the items before the -->, append them into a reactant list and remove spaces
        for x in range(len(seperated[1].split("+"))):       ## for amount of items before -->
            products.append(seperated[1].split("+")[x].replace(" ", "")) ## for the items after -->, append them into a reactant list and remove spaces
        
        # takes care of changing half assed balanced formulas for ex if we have 2H2O but it aint correct or whatever, we change it to H2O... Easier for output, easier for pc, win win
        for x in range(len(reactants)):
            for y in range(len(reactants[x])):
                if reactants[x][0] in " 1234567890":
                    reactants[x] = reactants[x][1:]
        for x in range(len(products)):
            for y in range(len(products[x])):
                if products[x][0] in " 1234567890":
                    products[x] = products[x][1:]
        
        print("react and prod")
        original_reactants = reactants[:]
        original_products = products[:]
        print(reactants)
        print(products)
        print(original_reactants)
        print(original_products)
    except:
        return("error",1,"error")
    try:
        for x in range(len(reactants)):                     ## for the reactants
            if "." in reactants[x] and "(" in reactants[x]:
                reactants[x] = (hydrateHandler(reactants[x]))                   ## check for hydrates since it is less intensive and less chance of messing something up than with the paranthesis
                reactants[x] = (atomsniffer(parentheseshandler(reactants[x])))
                print("both")
            elif "." in reactants[x]:               
                reactants[x] = (atomsniffer(hydrateHandler(reactants[x])))
                print("hydrate")
            elif "(" in reactants[x]:
                reactants[x] = (atomsniffer(parentheseshandler(reactants[x])))
                print("OH group")
            else:
                print("normal")
                reactants[x] = (atomsniffer(reactants[x]))
        print("Reactants: ", reactants)
        for x in range(len(products)):                     ## for the reactants
            if "." in products[x] and "(" in products[x]:
                products[x] = (hydrateHandler(products[x]))                   ## check for hydrates since it is less intensive and less chance of messing something up than with the paranthesis
                products[x] = (atomsniffer(parentheseshandler(products[x])))
            elif "." in products[x]:
                products[x] = (atomsniffer(hydrateHandler(products[x])))
            elif "(" in products[x]:
                products[x] = (atomsniffer(parentheseshandler(products[x])))
            else:
                products[x] = (atomsniffer(products[x]))
        print("Products: ", products)           ## looks like this [   [list of atoms in molecule 1]  , [list of atoms in molecule 2]  ]
    except:
        return("error",1,"error")





    reactants_list_of_dicts = []  
    for x in range(len(reactants)):    ## for amount of items in the whole products list, so molecule 1 and 2
        add_to_dict = {}
        for y in range(len(reactants[x])):        ## amount of items in that specific molecule, so all the atoms
                add_to_dict.update(re.findall(r'([A-Z][a-z]*)(\d*|\(\d*\))', reactants[x][y]))          ## splits everything and adds to the dictionary each atom in the list, ["H1, Na1, O1"] iterates for every item and updates to the dict to get [{"H" : "1", "Na" : "1", "O", "1"}]
        if add_to_dict:   ## if dict has items in it
            reactants_list_of_dicts.append(add_to_dict)       ## add it to the list
    for x in reactants_list_of_dicts:            ## this converts all values in the dict from str to int, so {"Na": "1"} --> {"Na", 1}
        for key in x:
            x[key] = int(x[key])
    print("Reactants:")
    print(reactants_list_of_dicts)


    products_list_of_dicts = []  
    for x in range(len(products)):    ## for amount of items in the whole products list, so molecule 1 and 2
        add_to_dict = {}
        for y in range(len(products[x])):        ## amount of items in that specific molecule, so all the atoms
                add_to_dict.update(re.findall(r'([A-Z][a-z]*)(\d*|\(\d*\))', products[x][y]))
        if add_to_dict:
            products_list_of_dicts.append(add_to_dict)
    for x in products_list_of_dicts:            ## this converts all values in the dict from str to int, so {"Na": "1"} --> {"Na", 1}
        for key in x:
            x[key] = int(x[key])
    print("Products")
    print(products_list_of_dicts)
    ##### actual thinking logic now :(
    ## base case with 1
    multiplier = 1
    used_atoms = []
    all_atoms_reactants = []
    all_atoms_products = []
    
    total_molecules = len(reactants_list_of_dicts) + len(products_list_of_dicts)             ## this shows the TOTAL amount of molecules in a the equation
    multipliers = []    ## this will be the total multipliers but gets a value later on
    for x in range(total_molecules):             ## this adds 1 in the multiplier list for every molecule, so if we have 3 molecules it would be [1,1,1]
        multipliers.append(1)
    for x in range(len(reactants_list_of_dicts)):                 ## grabs all keys in the dictionarys in the list to see all atoms useds
        for dictKey in reactants_list_of_dicts[x].keys():         ## this adds all the atoms used in the reactants into a list                              
            all_atoms_reactants.append(dictKey) 
    print(all_atoms_reactants)

    for x in range(len(products_list_of_dicts)):                  ## grabs all keys in the dictionarys in the list to see all atoms useds
        for dictKey in products_list_of_dicts[x].keys():          ## this adds all the atoms used in the products into a list                                   
            all_atoms_products.append(dictKey)
    print(all_atoms_products)
    for x in range(len(all_atoms_reactants)):                   ## checks if the current atom of the reactants list is in the products list    --> law of conservation of atoms
        if all_atoms_reactants[x] not in all_atoms_products:   
            print("error")
            return("The law of conservation of atoms does not seem to be applied here, please ensure that all the atoms in your reactants are also in your products!",2,"error")
    for x in range(len(all_atoms_products)):                    ## checks if the current atom of the products list is in the reactants list  --> law of conservation of atoms
        if all_atoms_products[x] not in all_atoms_reactants:     
            return("The law of conservation of atoms does not seem to be applied here, please ensure that all the atoms in your reactants are also in your products!",2,"error")
    TotalAtoms = list(dict.fromkeys(all_atoms_reactants))          ## since the law of conservation of atoms is true at this stage, the used atoms is equal to the ones used in reactants or product
    print(total_molecules)  
    print(multipliers)
    numbers_range = range(1, max_value)    ## the maximum amount we can go to for balancing the equation, so a maximum of a certain amount of molecules
    correct = True
# Generate all possible combinations
    all_combinations = list(product(numbers_range, repeat=len(multipliers)))
    total_combinations = len(all_combinations)
    for x in range(len(all_combinations)):      ## all multipliers
        #### lists we are going to need 
        #reactants_list_of_dicts         contains all molecules like so:  [{"Na" : 1} {"Na":2 , "H" : 1}] this is an example but this represents 2 molecules
        #products_list_of_dicts         contains all molecules like so:  [{"Na" : 1} {"Na":2 , "H" : 1}] this is an example but this represents 2 molecules
        #TotalAtoms to figure out what atom we are on, looks like  ["Na", "H"]
        check_list = []      ## we can append bools to see if all atoms match
        for y in range(len(TotalAtoms)):
            current_molecule = 0
            current_atom = TotalAtoms[y]           ## current atom we need to equalize
            reactants_amount = 0                   ## will be used to see how many of the current_atom is in each molecule and then in each side... if they are equal then we can continue   gets reset on every new atom
            product_amount = 0                     ## will be used to see how many of the current_atom is in each molecule and then in each side... if they are equal then we can continue   gets reset on every new atom
            ### first we need to do it for all reactants
            for z in range(len(reactants_list_of_dicts)):     ## for the amount of molecules we have
                try:
                    reactants_amount += reactants_list_of_dicts[z][current_atom] * all_combinations[x][current_molecule]      ## reactants_list_of_dicts[z] is the current dictionary, and [current_atom] is the atom we are looking for, so this returns the value, then we multiply it by its corrective multiplier
                except:
                    pass
                current_molecule += 1
            for c in range(len(products_list_of_dicts)):     ## for the amount of molecules we have
                try:
                    product_amount += products_list_of_dicts[c][current_atom] * all_combinations[x][current_molecule]      ## reactants_list_of_dicts[c] is the current dictionary, and [current_atom] is the atom we are looking for, so this returns the value, then we multiply it by its corrective multiplier
                except:
                    pass
                current_molecule += 1
            # if x == 6571:
            #     print("total atoms: ", str(TotalAtoms))
            #     print("total molecules react: ", str(reactants_list_of_dicts))
            #     print("total molecules prod: ", str(reactants_list_of_dicts))
            #     print(all_combinations[x])
            #     print(TotalAtoms[y])
            #     print(reactants_amount)
            #     print(product_amount)
            if reactants_amount == product_amount:        ## the rest is a list to check when we are correct, we append a True value to the list and if not a False, we use this to keep track if ALL atoms are equal, not just if one is
                check_list.append(True)
                # print("total atoms: ", str(TotalAtoms))
                # print("total molecules react: ", str(reactants_list_of_dicts))
                # print("total molecules prod: ", str(reactants_list_of_dicts))
                # print(all_combinations[x])
                # print(TotalAtoms[y])
                # print(reactants_amount)
                # print(product_amount)
                
            else:
                check_list.append(False)
        if False not in check_list:
            print(all_combinations[x])
            solved_solution = all_combinations[x]
            break
    if False in check_list:
        print("Not possible")
        return("It seems like this is not possible to balance, or you may need to increase the maximum molecule amount",3, "error")
        
    ## reconstruct string
        
    
    
    seperated = formula.split("-->") ## rebuild the string!
    reactants_for_string = seperated[0].split("+")
    products_for_string = seperated[1].split("+")
    completed_formula = ""
    # for x in range(len(reactants_for_string)):
    #     completed_formula += str(solved_solution[x]) + " "+ reactants_for_string[x].replace(" ", "")     ## this replaces all the spaces in the string if there are any so then every case is the same, so then we just add a space 
    #     if x != len(reactants_for_string)-1:
    #         completed_formula += " + "

    # iterator = len(reactants)
    # completed_formula += " --> "
    # for x in range(len(products_for_string)):
    #     completed_formula += str(solved_solution[x + iterator]) + " " + products_for_string[x].replace(" ", "")
    #     if x != len(products_for_string)-1:
    #         completed_formula += " + " 
    # print(completed_formula)
    try:
        for x in range(len(original_reactants)):
            if solved_solution[x] == 1:
                completed_formula += original_reactants[x].replace(" ", "").replace("1", "")
            else:
                completed_formula += str(solved_solution[x]) + " "+ original_reactants[x].replace(" ", "").replace("1", "")
            if x != len(original_reactants)-1:
                completed_formula += " + "
            
        iterator = len(reactants)
        completed_formula += " --> "
        for x in range(len(original_products)):
            print("yip")
            if solved_solution[x + iterator] == 1:
                completed_formula += original_products[x].replace(" ", "").replace("1", "")
            else:
                completed_formula += str(solved_solution[x + iterator]) + " "+ original_products[x].replace(" ", "").replace("1", "")
            if x != len(original_products)-1:
                completed_formula += " + "
    except:
        return("Error", 1,"error")
    print("Time:" + str(time.time()-t1))
    print(completed_formula)
    if return_array == False:
        return(completed_formula,0)
    else:
        return(completed_formula, 0,solved_solution)
    ####################################old way###############################################
    # for x in range(len(reactants)):
    #     ## build molecule, since the reactants list looks like [[atom_1_for_molecule_1, atom_1_for_molecule_1] [atom_1_for_molecule_2, atom_2_for_molecule_2]]
    #     molecule = ""      #original_reactants ## new version looks like ["molecule1", "molecule2"]
    #     for y in range(len(reactants[x])):
    #         molecule += reactants[x][y]
    #         print(molecule)
    #     if solved_solution[x] == 1:         ## this is just so instead of 1 O2 it would just be O2
    #         completed_formula += molecule.replace(" ", "").replace("1", "") ## replaces spaces, like said already, the 1 is replaced for example Na1 would just be Na
    #     else:
    #         completed_formula += str(solved_solution[x]) + " "+ molecule.replace(" ", "").replace("1", "")     ## this replaces all the spaces in the string if there are any so then every case is the same, so then we just add a space 
    #     if x != len(reactants)-1:
    #         completed_formula += " + "
    # iterator = len(reactants)
    # completed_formula += " --> "
    # for x in range(len(products)):
    #     molecule = ""
    #     for y in range(len(products[x])):
    #         molecule += products[x][y]
    #     if solved_solution[x + iterator] == 1:
    #         completed_formula += molecule.replace(" ", "").replace("1", "")
    #     else:
    #         completed_formula += str(solved_solution[x + iterator]) + " " + molecule.replace(" ", "").replace("1", "")
    #     if x != len(products)-1:
    #         completed_formula += " + " 
    # print(completed_formula)


    # completed_formula = ""

    
    
    
    
    
    
    
    
    
    #split1 = formula.split("-->")[0].split("+")
    # elements = []
    # for x in range(len(split1)):
    
    #     for y in range(20):
    #         splitted = split1[x].split(str(y))
    #         if len(splitted) != 1:
    #             splitted[x-1] += str(y)
    #             elements.append(splitted[x-1])
                

    #         elif y == 19:
    #             print(splitted)
    #             elements.append(splitted[x-1])
        
    # print(elements)
    # Use regular expression to match elements, counts, and parentheses
    

    
def parentheseshandler(formula): ### caution, this changes the order of the atoms, this is to be used only to do calculations revolving calulations :)
    elements = {}
    stack = []
    stack = formula.split("(")
    for x in range(len(formula.split("("))):
        element = formula.split("(")[x].split(")")[0]
        if x > 0:
            amount = re.sub(r'([A-Za-z])', '', formula.split("(")[x].split(")")[1])
            if amount == "":
                amount = "1"
            elements[str(element)] = amount
            print("e")
            print(elements)
        print("stack")
        print(stack)
    for x in range(len(stack)):
        print("oba joba  " + str(stack[x]))
        if ")" in stack[x]:
            str_for_numb = stack[x].split(")")     ## goes from "SO4)3" to ["SO4, 3"]
            
            key = list(elements)[0]                  ## gets us the key, so for example OH
            print("key")
            print(key)
            atomlist = re.findall(r'[A-Z][a-z]*', key)        ## finds all atoms of the key, so for OH we get O and H
            print(atomlist)   ## [S , O ] for ex
            new = []
            with_numb = re.findall(r'[A-Z][a-z]*\d', key)
            for z in atomlist:                             ### here we will create a second list that will hold the amount of atoms
                for y in range( len(with_numb)): 
                    print(y) 
                    print(len(with_numb))  
                    if z in with_numb[y]:
                        print("in")
                        new.append(with_numb[y])     ## if the atom is in the number list, which means it has a specified amount , for ex O2,O3. but not O
                        break
                    else:                                               
                        print("not in")
                    if y == len(with_numb)-1:   ## on the final check if it is still not in the list, we will add it as 1, So to the new list we add O1 for example
                        print("end")
                        new.append(str(z)+"1")
                        break
                if len(with_numb) == 0:  ## this is if all atoms are 1, so of ex (OH)1 will not produce a with_numb list since there are not numbers, so here they will be manually added
                    new.append(str(z)+"1")
            for z in range(len(new)):              ## here we convert everything to pure numbers ["S1", "O2"] --> ["1", "2"]
                new[z] = re.sub('\D', '', new[z])
            print(new)
            ## visualize    
            ## ["S", "O"]    --> list of atoms
            ## ["1", "2"]    --> amount of each atom




            burnerlist = []                    ## create list
            for y in range(len(atomlist)):           ## for amount in atom list
                burnerlist.append(atomlist[y] + str(int(elements[key]) *int(new[y]) ) )     ## add numbers, so for ex O2H2 and make sure to multiple the amount of atoms
            
            added_str = ""                                          ## reset string
            for z in range(len(burnerlist)):                        ## for amount of items in burnerlist, so "O2" , "H2"
                added_str +=  burnerlist[z]                         ## we join them to get "O2H2"
            
            rest_of_atoms = stack[x].split(")")
            rest_of_atoms = rest_of_atoms[1].replace(elements[key], '', 1)    ## finds left over atoms still sticking together and splits em
    
            stack.append(rest_of_atoms) ## we can just add now, since the order doesnt matter

            stack[x] =  added_str
            
            del elements[key]
    new_formula = ""
    for x in range(len(stack)):
        new_formula += stack[x]
    print(new_formula)
    return(new_formula)
def hydrateHandler(formula):
    formula = formula
    hydrate = formula.split(".")[1]
    rest_of_formula = formula.split(".")[0]
    hydrate_amount = re.findall(r'(\d+|\D+)', hydrate)[0]
    hydrate_atoms = re.findall(r'([A-Z][a-z]*)(\d*|\(\d*\))', hydrate)     ### for ex for .2H2O it is [("H", "2"), ("O", "")]
    new_molecue = rest_of_formula
    for x in range(len(hydrate_atoms)):
        if hydrate_atoms[x][1] == "":
            new_molecue += hydrate_atoms[x][0]+str(1 * hydrate_amount)  
        else: 
            new_molecue += hydrate_atoms[x][0]+str(int(hydrate_atoms[x][1]) * int(hydrate_amount)) 
            print(new_molecue)  
    print(hydrate_amount)
    print(hydrate_atoms) 
    print(new_molecue)
    return(new_molecue)  
def atomsniffer(formula):
    elements = {}
    stack = []
    matches = re.findall(r'([A-Z][a-z]*)(\d*|\(\d*\))', formula) ## dict
    print(matches)
    # Create a dictionary to store elements and their counts
   

    for element, count in matches:
            count = int(count) if count else 1
            elements[element] = elements.get(element, 0) + count

            # Check if there are elements inside parentheses
            while stack and isinstance(stack[-1], dict):
                inner_element, inner_count = stack.pop().popitem()
                elements[inner_element] = elements.get(inner_element, 0) + inner_count

    # Format the result as "element:count"
    result = [f"{element}{count}" for element, count in elements.items()]
    print(result)
    
    return(result)
def Stoichiometry(molecules_list,type, molecule, amount, array):
    molecules = molecules_list[:]
    print("the actual fuck")
    print(molecules)
    print(molecules_list)
    print(type)
    print(molecule)
    print(amount)
    print(array)
    print(molecules)
    print(molecules)
    atom_and_amount = {}
    print(molecules)
    print(len(molecules))
    print(array)
    for x in range(len(array)):
        atom_and_amount[molecules[x]] = int(array[x])
    if type == "Via Mol":       ## here we will create a value dict that will look like this {"Element": ["Mass", "Mol", "Molar Mass"]}
        try:
            compound = pcp.get_compounds(molecule, "formula")[0]     ## searches the Pubchem webpage and grabs the first compund, returns it as compound(ID)
        except:
            compound = ""
        molar_mass = float(compound.molecular_weight)       ## get molar mass from pubchem
        print(amount)
        single_mol = amount / atom_and_amount[molecule]     ## calculate the single mol for refrence, for ex if we have 2 Na is 1 mol, we need to be able to balance with 1 molecule being then 0.5 mol. We then multiply it by the amount of molecules
        mass = molar_mass * (atom_and_amount[molecule] *single_mol)  ## calculate molar mass
        print(single_mol) 
        print(mass)
        molecules.remove(molecule)
        
        solsd = {}
        solsd[molecule] = [amount,mass,molar_mass]
        for x in range(len(molecules)):
            try:
                compound = pcp.get_compounds(molecules[x], "formula")[0]     ## searches the Pubchem webpage and grabs the first compund, returns it as compound(ID)
            except:
                compound = ""
            molar_mass = float(compound.molecular_weight)
            mol = atom_and_amount[molecules[x]] *single_mol
            mass = molar_mass * (mol)  
            solsd[molecules[x]] = [mol,mass,molar_mass]
        print(solsd)
        return(solsd)
    if type == "Via Mass in Grams":
        try:
            compound = pcp.get_compounds(molecule, "formula")[0]     ## searches the Pubchem webpage and grabs the first compund, returns it as compound(ID)
        except:
            compound = ""
        print(compound)
        molar_mass = float(compound.molecular_weight)       ## get molar mass from pubchem
        mass = amount
        mol = mass / molar_mass
        single_mol = mol/atom_and_amount[molecule]
        molecules.remove(molecule)
        
        solsd = {}
        solsd[molecule] = [mol,mass,molar_mass]
        for x in range(len(molecules)):
            try:
                compound = pcp.get_compounds(molecules[x], "formula")[0]     ## searches the Pubchem webpage and grabs the first compund, returns it as compound(ID)
            except:
                compound = ""
            molar_mass = float(compound.molecular_weight)
            mol = atom_and_amount[molecules[x]] *single_mol
            mass = molar_mass * (mol)  
            solsd[molecules[x]] = [mol,mass,molar_mass]
            print("done")
        print(solsd)
        return(solsd)
# balancer("Al2(SO4)3 + Ca(OH)2 --> Al(OH)3 + CaSO4",10)
#def CandHguesser(carbon):                        ### ONLY ONE DOUBLE OR TRIPLE BOND!!!   --> not needed anymore.
    ## old way of doing it
    
    
    # HydrogenAlkanen = ((carbon * 2)+2)        ## calculates hydrogen atoms needed to form an alkaan
    # HydrogenAlkenen = carbon * 2                ## calculates hydrogen atoms needed to form an alkeen            
    # HydrogenAlkynen = (carbon *2) -2             ## calculates hydrogen atoms needed to form an alkyn
    # CandHguesserls = carboninator(carbon)        ## runs the Hydrogen atoms for an alkaan through the function and returns results
    # print("###########")
    # strucutrealkaan, alkaaanvisual,alkaanname, alkaanbruto = CandHguesserls[0],CandHguesserls[1],CandHguesserls[2],CandHguesserls[3] ## structure, photo, name   only grabs the first 4 items of the list since only those relate to alkanen
    # CandHguesserls = carboninator(carbon,HydrogenAlkenen,0, True)      ##same as previous
    # alkenenFormula,alkenenStructure, alkenenNames,alkenenBrutos = CandHguesserls[4], CandHguesserls[5], CandHguesserls[6], CandHguesserls[7]   ## structure, photo, names     ##same as previous 
    # CandHguesserls = carboninator(carbon,HydrogenAlkynen,0 , True)         ##same as previous
    # alkynenFormula,alkynenStructure, alkynenNames, alkynenBrutos = CandHguesserls[8], CandHguesserls[9], CandHguesserls[10], CandHguesserls[11] ##same as previous
    # print(strucutrealkaan,alkaanname,alkaanbruto )   ##check
    # print(alkenenFormula, alkenenNames, alkenenBrutos) ## check
    #return(strucutrealkaan, alkaaanvisual,alkaanname, alkaanbruto,               alkenenFormula,alkenenStructure, alkenenNames,alkenenBrutos     ,         alkynenFormula,alkynenStructure, alkynenNames, alkynenBrutos  ) ### Structure (C-C), PHOTO , NAME, Structure(C=C). PHOTO, NAME






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


