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
def balancer(formula):
    print("working on it")
    formula = "NaL + (OH)2 +  --> NaOH + H2"
    seperated = formula.split("-->") ## splits list into 2 items, 0 --> 1 ( before and after the paranthesis)
    reactants = []
    products = []
    for x in range(len(seperated[0].split("+"))):       ## for amount of items before -->, 1 is after and 0 is before since seperated splits into before and after the arrow
        reactants.append(seperated[0].split("+")[x].replace(" ", ""))     ## for the items before -->, append them into a reactant list and remove spaces
    for x in range(len(seperated[1].split("+"))):       ## for amount of items before -->
        products.append(seperated[1].split("+")[x].replace(" ", "")) ## for the items after -->, append them into a reactant list and remove spaces
    
    for x in range(len(reactants)):                     ## for the reactants
        if "(" in reactants[x]:
            reactants[x] = (atomsniffer(parentheseshandler(reactants[x])))
        else:
            reactants[x] = (atomsniffer(reactants[x]))
    print("Reactants: ", reactants)
    for x in range(len(products)):                     ## for the reactants
        if "(" in products[x]:
            products[x] = (atomsniffer(parentheseshandler(products[x])))
        else:
            products[x] = (atomsniffer(products[x]))
    print("Products: ", products)           ## looks like this [   [list of atoms in molecule 1]  , [list of atoms in molecule 2]  ]





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
    
    for x in range(len(reactants_list_of_dicts)):                 ## grabs all keys in the dictionarys in the list to see all atoms useds
        for dictKey in reactants_list_of_dicts[x].keys():                                       
            all_atoms_reactants.append(dictKey)
    print(all_atoms_reactants)

    for x in range(len(products_list_of_dicts)):                  ## grabs all keys in the dictionarys in the list to see all atoms useds
        for dictKey in products_list_of_dicts[x].keys():                                       
            all_atoms_products.append(dictKey)
    print(all_atoms_products)
    for x in range(len(all_atoms_reactants)):                   ## checks if the current atom of the reactants list is in the products list    --> law of conservation of atoms
        if all_atoms_reactants[x] not in all_atoms_products:
            print("error")
            return("error")
    for x in range(len(all_atoms_products)):
        if all_atoms_products[x] not in all_atoms_reactants:     ## checks if the current atom of the products list is in the reactants list  --> law of conservation of atoms
            return("error")  
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
            print(elements)
        print(stack)
    for x in range(len(stack)):
        print("oba joba  " + str(stack[x]))
        if ")" in stack[x]:
            key = list(elements)[0]                  ## gets us the key, so for example OH
            atomlist = re.findall(r'[A-Z][a-z]*', key)        ## finds all atoms of the key, so for OH we get O and H
            print(atomlist)  
            burnerlist = []                    ## create list
            for y in range(len(atomlist)):           ## for amount in atom list
                burnerlist.append(atomlist[y] + elements[key])      ## add numbers, so for ex O2H2
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
balancer("as")
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


