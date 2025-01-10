import streamlit as st
import time
import Math
import Chem
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from io import BytesIO
import re
# Write text
if 'status' not in st.session_state:
    st.session_state.status = 0               
if 'graphxmin' not in st.session_state:
    st.session_state.graphxmin = -10
if 'graphxmax' not in st.session_state:
    st.session_state.graphxmax = 10
if 'graphymin' not in st.session_state:
    st.session_state.graphymin = -10
if 'graphymax' not in st.session_state:
    st.session_state.graphymax = 10
## MAIN screen
def convert_to_superscript(text):
    return text.replace("^2", "²").replace("^3", "³").replace("^4", "⁴").replace("^5", "⁵").replace("^6", "⁶").replace("^7", "⁷").replace("^8", "⁸").replace("^9", "⁹").replace("^1", "")
def convert_to_polynomial(input_polynomial):
    return_string = ""
    for x in range(len(input_polynomial)):
        if x == 0:
            return_string += str(input_polynomial[-(x+1)]) + convert_to_superscript(f"x^{len(input_polynomial) - (x+1)}")
        else:
            if input_polynomial[-(x+1)] > 0 and len(input_polynomial) - (x+1) != 0:
                return_string += "+" + str(input_polynomial[-(x+1)]) + convert_to_superscript(f"x^{len(input_polynomial) - (x+1)}")
            elif input_polynomial[-(x+1)] > 0 and len(input_polynomial) - (x+1) == 0:
                return_string += "+" + str(input_polynomial[-(x+1)])
            elif input_polynomial[-(x+1)] < 0 and len(input_polynomial) - (x+1) != 0:
                return_string += str(input_polynomial[-(x+1)]) + convert_to_superscript(f"x^{len(input_polynomial) - (x+1)}")
            else:
                return_string += str(input_polynomial[-(x+1)]) + convert_to_superscript(f"x^{len(input_polynomial) - (x+1)}")
    return return_string
if st.session_state.status == 0:      # 0 = home, 1-4 math, 5 Chem
    st.markdown("<h1 style='text-align: center; color: red;'>Advanced Calculator</h1>", unsafe_allow_html=True)
    st.markdown("")
    st.markdown("")
    st.markdown("")
    st.markdown("")
    st.markdown("")
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        chemistry = st.button("Chemistry")
    with col2:
        math = st.button("Math")
    if chemistry:
        st.session_state.status = 2
        st.rerun()
    elif math:
        st.session_state.status = 1
        st.rerun()
    
    
elif (st.session_state.status == 1):
    st.header("Math")
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        funcinfo = st.button("Function info X^2")
    with col2:
        funcompare = st.button("2 functions X^2")
    with col3:
        eular_div = st.button("Eular divison")
    if funcinfo:
        st.session_state.status = 3
        st.rerun()
    elif funcompare:
        st.session_state.status = 4
        st.rerun()
    elif eular_div:
        st.session_state.status = 12
        st.rerun()
elif (st.session_state.status == 2):
    st.header("Chemistry")
    col1, col2, col3, col4, col5 = st.columns(5)
    with col1:
        moleculebalancer = st.button("Molecule Balancer")
    with col2:
        reaction = st.button("Reaction balancer (work in progress)")
    with col3:
        CandHinator = st.button("C and H inator")
    with col4:
        CandHguesser = st.button("C and H guesser")
    with col5:
        Stoichiometry  = st.button("Stoichiometry")
    col6, col7 = st.columns(2)
    with col6:
        Dilution = st.button("Dilution")
    with col7:
        ReactionPredictor  = st.button("ReactionPredictor ( Might work, might not ¯\_(ツ)_/¯ )")
    
    if moleculebalancer:
        st.session_state.status = 5
        st.rerun()
    elif reaction:
        st.session_state.status = 6
        st.rerun()
    elif CandHinator:
        st.session_state.status = 7
        st.rerun()
    elif CandHguesser:
        st.session_state.status = 8
        st.rerun()
    elif Stoichiometry:
        st.session_state.status = 9
        st.rerun()
    elif Dilution:
        st.session_state.status = 10
        st.rerun()
    elif ReactionPredictor:
        st.session_state.status = 11
        st.rerun()

elif (st.session_state.status == 3):
    st.header("Function INFO")
    func1 = st.text_input("Function 1")
    if(st.button('Submit')): 
        try:
            a,b,c = Math.FindTerms(func1)
            try:
                alfa,beta = Math.AlphaBetaForm(func1)
                d,x1,x2 = Math.Solution(func1)
                if a == 0:
                    st.write("Alfa Beta form: Not possible")
                elif alfa < 0 and beta < 0:
                    st.write("Alfa Beta form: {}(x+{})^2+{}".format(a,str(abs(alfa)),beta))
                elif beta < 0:
                    st.write("Alfa Beta form: {}(x-{})^2{}".format(a,str(abs(alfa)),beta))
                elif alfa < 0:
                    st.write("Alfa Beta form: {}(x+{})^2+{}".format(a,str(abs(alfa)),beta))
                else:
                    st.write("Alfa Beta form: {}(x+{})^2+{}".format(a,str(abs(alfa)),beta))
                st.write("a:" + str(a))
                st.write("b:" + str(b))
                st.write("c:" + str(c))
                st.write("Discriminant:" + str(d))
                st.write("X1: " + str(x1))
                st.write("X2: " +str(x1))

                f = lambda x: float(a)*x**2 + float(b)*x + float(c)

                # Generate x values
                if x1 != "None" and x2 != "None":
                    x = np.linspace(x1-10, float(x2)+10)  # Adjust the range of x values as needed
                elif x1 != "None":
                    x = np.linspace(float(x1)-10, float(x1)+10)  # Adjust the range of x values as needed
                elif x2 != "None":
                    x = np.linspace(float(x2)-10, float(x2)+10)  # Adjust the range of x values as needed
                else:
                    x = np.linspace(-10, 10)  # Adjust the range of x values as needed

                # Calculate y values using the function
                y = f(x)

                # Create a plot using Matplotlib
                plt.figure(figsize=(8, 6))
                plt.plot(x, y, label=func1)
                plt.title('Graph of f(x) =' + func1)
                plt.xlabel('x')
                plt.ylabel('f(x)')
                plt.axhline(0, color='black',linewidth=0.5)
                plt.axvline(0, color='black',linewidth=0.5)
                plt.legend()

                # Display the plot in Streamlit
                st.pyplot(plt)


            except:
                st.error("Error, please try again")
        except:
            st.error("Error, please try again")





elif (st.session_state.status == 4):
    func1 = st.text_input("Function 1")
    func2 = st.text_input("Function 2")
    
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        graphxminbox = st.text_input("X min", st.session_state.graphxmin)
        try:
            st.session_state.graphxmin = int(graphxminbox)
        except:
            st.error("Please input a number in as the x or y parameter")
            
    with col2:
        graphxmaxbox = st.text_input("X max", st.session_state.graphxmax)
        try:
            st.session_state.graphxmax = int(graphxmaxbox)
        except:
            st.error("Please input a number in as the x or y parameter")
    with col3:
        graphyminbox = st.text_input("Y min", st.session_state.graphymin)
        try:
            st.session_state.graphymin = int(graphyminbox)
        except:
            st.error("Please input a number in as the x or y parameter")
    with col4:
        graphymaxbox = st.text_input("Y max", st.session_state.graphymax)
        try:
            st.session_state.graphymax = int(graphymaxbox)
        except:
            st.error("Please input a number in as the x or y parameter")
    
  
    
    sliderxmin = 0
    func1 = func1
    func2 = func2
    print(func1)
    try:
        d,x1,x2 = Math.solve2func(func1,func2)
        
        a,b,c = Math.FindTerms(func1)
        a1,b1,c1 = Math.FindTerms(func2)


        f = lambda x: float(a)*x**2 + float(b)*x + float(c)
        g = lambda x: float(a1)*x**2 + float(b1)*x + float(c1)
        # Generate x values
        if x1 != "None" and x2 != "None":
            x = np.linspace(x1-10, float(x2)+10)  # Adjust the range of x values as needed
        elif x1 != "None":
            x = np.linspace(float(x1)-10, float(x1)+10)  # Adjust the range of x values as needed
        elif x2 != "None":
            x = np.linspace(float(x2)-10, float(x2)+10)  # Adjust the range of x values as needed
        else:
            x = np.linspace(-10, 10)  # Adjust the range of x values as needed

        # Calculate y values using the function
        y1 = f(x)
        y2 = g(x)
        # Create a plot using Matplotlib
        plt.figure(figsize=(8, 6))
        plt.plot(x, y1, label='f(x) = ' + func1)
        plt.plot(x,y2, label='g(x) = ' + func2)
        plt.title('Graph of f(x) = ' + func1 + "and g(x) = " + func2)
        plt.xlabel('x')
        plt.ylabel('f(x)')
        plt.axhline(0, color='black',linewidth=0.5)
        plt.axvline(0, color='black',linewidth=0.5)
        plt.legend()
        plt.xlim(st.session_state.graphxmin,st.session_state.graphxmax)
        plt.ylim(st.session_state.graphymin,st.session_state.graphymax)
        # Display the plot in Streamlit
        st.pyplot(plt)
        st.success("X1 = {} , X2 = {}".format(x1,x2))




    except:
        st.error("Error, please try again")


elif (st.session_state.status == 5):
    st.header("Molecule Balancer")
    element1 = st.text_input("Element 1")
    element2 = st.text_input("Element 2")
    if (st.button("submit")):
        #try:
        sols = Chem.MoleculeStable(element1,element2)
        #st.error("May take a second to load, please be patient")
        images = []
        if len(sols) != 1 and len(sols) != 0:
            progress_bar = st.progress(0)
            added = int(round(100/len(sols)))
            for x in range(len(sols)):
                compound = sols[x].split(" ")[0]
                image = Chem.MoleculeVisualizeViaBruto(str(compound))
                images.append(image)
                progress_bar.progress((x + 1) * added)
            progress_bar.empty()
        elif len(sols) == 0:
            st.error("Error, this compound may exist but the program ran into an error")
        else:
            with st.spinner('Loading...'):

                for x in range(len(sols)):
                    compound = sols[x].split(" ")[0]
                    image = Chem.MoleculeVisualizeViaBruto(str(compound))
                    images.append(image)

        
        
        for x in range(len(sols)):
            compound  = sols[x].split(" ")[0]

            col1, col2 = st.columns(2)
            with col1:
                st.write(sols[x])
            with col2:
                st.image(images[x])
        st.success("Loaded!")
        # except:
        #     st.warning("operation failed!")
elif (st.session_state.status == 6):      
    st.write("Equation balancer!")
    st.write("Please input an equation like this: Na + H2O --> NaOH + H2")
    formula = st.text_input(label = "Input here!")
    max_molecule_counter = st.text_input(label= "Maximum number of molecules the program will try to count to, for ex 9 would maximize at 9 OH, and 10 OH will not be calculated" , value= "10")
    if (st.button("submit")):
        try:
            max_molecule_counter = int(max_molecule_counter)
        except:
            st.rerun()
        solution,code = Chem.balancer(formula, int(max_molecule_counter), False)
        if code == 0:
            st.success(solution)
        elif code == 2:
            st.warning(solution)
        else:
            st.error(solution)
elif (st.session_state.status == 7):
    col1,col2,col3, col4,col5, col6 = st.columns([0.02,0.1,0.02,0.1,0.02,0.1])
    with col1:
        st.header("C")
    with col2:
        Catoms = st.text_input(value=0, key=1, label="Carbon atoms", label_visibility = "hidden")
    with col3:
        st.header("H")
    with col4:
        Hatoms = st.text_input(value=0, key=2, label="Hydrogen atoms", label_visibility = "hidden")
    with col5:
        st.header("O")
    with col6:
        Oatoms = st.text_input(value=0, key=3, label="Oxygen atoms", label_visibility = "hidden")
    if (st.button("submit")):
        try:
            Hatoms = int(Hatoms)
            Catoms = int(Catoms)
            Oatoms = int(Oatoms)
        except:
            st.rerun()
        if Oatoms == 0 and Catoms > 0 and Hatoms > 0:
            with st.spinner('Loading...'):
                
                alkaan, alkaanvisual,alkaanname, alkaanbruto,        alkeenstructures, alkeenvisuals, alkenennames, alkenenbrutos,        alkynenstructures,alkynenvisuals, alkynennames, alkynenbrutos = Chem.CandHinator(Catoms,Hatoms,Oatoms)
            if alkaan == "None":
                st.header("No Alkaan")
                st.write("")
                st.write("")
                st.write("")
            else:
                col1, col2 = st.columns(2)
                st.write("")
                st.write("")
                st.write("")
                with col1:
                    st.write(alkaan)
                    st.write(alkaanname)
                    st.write(alkaanbruto)
                with col2:
                    try:
                        st.image(alkaanvisual)
                    except:
                        st.write("Failed to get image!")
            if alkeenstructures:
                
                for x in range(len(alkeenstructures)):
                    col1, col2 = st.columns(2)
                    with col1:
                        st.write(alkeenstructures[x])
                        st.write(alkenennames[x])
                        st.write(alkenenbrutos[x])
                    with col2:
                        st.image(alkeenvisuals[x])
            else:
                st.header("No Alkenen")
            if alkynenstructures:
                
                for x in range(len(alkynenstructures)):
                    col1, col2 = st.columns(2)
                    with col1:
                        st.write(alkynenstructures[x])
                        st.write(alkynennames[x])
                        st.write(alkynenbrutos[x])
                    with col2:
                        st.image(alkynenvisuals[x])
            else:
                st.header("No Alkynen")
            st.header("No Alcohols")
        elif Catoms > 0 and Hatoms > 0:
            with st.spinner('Loading...'):
                alcoholStructure,alcoholvisuals,alcoholNames,alcoholBrutos = Chem.CandHinator(Catoms,Hatoms,Oatoms)
            if alcoholStructure == "None":
                st.header("No Alkanen")
                st.header("No Alkenen")
                st.header("No Alkynen")
                st.header("No Alcohols")
        
            elif alcoholStructure:
                print(alcoholvisuals)
                for x in range(len(alcoholStructure)):
                    col1, col2 = st.columns(2)
                    with col1:
                        st.write(alcoholStructure[x])
                        st.write(alcoholNames[x])
                        st.write(alcoholBrutos[x])
                    with col2:
                        st.image(alcoholvisuals[x])
        else:
            st.header("No Alcohols")
        if Catoms == 0 or Hatoms == 0:
            st.header("No Alkanen")
            st.header("No Alkenen")
            st.header("No Alkynen")
            st.header("No Alcohols")

elif (st.session_state.status == 8):
    col1, col2 = st.columns(2)
    with col1:
        st.header("C atoms")
    with col2:
        Catoms = st.text_input(value=0, label="Carbon atoms", label_visibility = "hidden")
    if (st.button("Submit")):
        try:
            int(Catoms)
        except:
            st.rerun()
        if Catoms == "0":
            st.rerun()
            
        with st.spinner('Loading...'):
                    
                    alkaan, alkaanvisual,alkaanname, alkaanbruto,        alkeenstructures, alkeenvisuals, alkenennames, alkenenbruto,        alkynenstructures,alkynenvisuals, alkynennames, alkynenbrutos = Chem.carboninator(int(Catoms))
        if alkaan == "None":
                    st.header("No Alkaan")
                    st.write("")
                    st.write("")
                    st.write("")
        else:
            col1, col2 = st.columns(2)
            st.write("")
            st.write("")
            st.write("")
            with col1:
                st.write(alkaan)
                st.write(alkaanname)
                st.write(alkaanbruto)
            with col2:
                try:
                    st.image(alkaanvisual)
                except:
                    st.write("Failed to get image!")
        if alkeenstructures:
            
            for x in range(len(alkeenstructures)):
                col1, col2 = st.columns(2)
                with col1:
                    st.write(alkeenstructures[x])
                    st.write(alkenennames[x])
                    st.write(alkenenbruto[x])
                with col2:
                    st.image(alkeenvisuals[x])
        else:
                st.header("No Alkenen")
        if alkynenstructures:
                
                for x in range(len(alkynenstructures)):
                    col1, col2 = st.columns(2)
                    with col1:
                        st.write(alkynenstructures[x])
                        st.write(alkynennames[x])
                        st.write(alkynenbrutos[x])
                    with col2:
                        st.image(alkynenvisuals[x])
        else:
            st.header("No Alkynen")
        if Catoms == 0:
            st.header("No Alkanen")
            st.header("No Alkenen")
            st.header("No Alkynen")

elif (st.session_state.status == 9):
    formula = st.text_input(label="Chemical reaction equation, does not need to be balanced. For example: Na + H2O --> NaOH + H2")
    
    
    try:
        
        r = []         ## reactants
        p = []         ## products
        new_formula = ""         ## this whole block of code basicly fixes idiots.... If someone were to give in a balanced ( correct or incorrect ) the program would have a level 3 uh oh ing stroke... So we fix it buy defaulting everything to 1, so 2H2O would just be H2O. Then we just remake everything and get a new formula
        seperated = formula.split("-->")
        for x in range(len(seperated[0].split("+"))):       ## for amount of items before the -->, 1 is after and 0 is before since seperated splits into before and after the arrow
            r.append(seperated[0].split("+")[x].replace(" ", ""))     ## for the items before the -->, append them into a reactant list and remove spaces
        for x in range(len(seperated[1].split("+"))):       ## for amount of items before -->
            p.append(seperated[1].split("+")[x].replace(" ", "")) ## for the items after -->, append them into a reactant list and remove spaces
        for x in range(len(r)):
            for y in range(len(r[x])):
                if r[x][0] in " 1234567890":
                    r[x] = r[x][1:]
        for x in range(len(p)):
            for y in range(len(p[x])):
                if p[x][0] in " 1234567890":
                    p[x] = p[x][1:]
        for x in range(len(r)):
            if x != len(r)-1:
                new_formula += str(r[x]) + "+"
            else:
                new_formula += str(r[x]) + "-->"
        for x in range(len(p)):
            if x != len(p)-1:
                new_formula += str(p[x]) + "+"
            else:
                new_formula += str(p[x])
        print("Wh ynot work?")
        print(new_formula)
        formula = new_formula
    except:
        print("uh oh  it aint work")        
    balanced_formula,error, array = Chem.balancer(formula, 10, True) # --> only accepts a single string writh an UNBALANCED formula to eat up
    ## formula isnt actually balanced, we use the array to manually balance it
    formula = new_formula
    if error == 0 and formula != "":
            st.success(balanced_formula)
    elif error == 2 and formula != "":
        st.warning(balanced_formula)
    elif formula != "":
        st.error(balanced_formula)

    # Specify the number of columns you want
    num_columns = len(array)

    # Create the specified number of columns
    cols = st.columns(num_columns)
    
    mol = 0

    total = []
    total_balanced = []
    try:
        
        seperated = formula.split("-->")
        for x in range(len(seperated[0].split("+"))):       ## for amount of items before the -->, 1 is after and 0 is before since seperated splits into before and after the arrow
            total.append(seperated[0].split("+")[x].replace(" ", ""))     ## for the items before the -->, append them into a reactant list and remove spaces
        for x in range(len(seperated[1].split("+"))):       ## for amount of items before -->
            total.append(seperated[1].split("+")[x].replace(" ", "")) ## for the items after -->, append them into a reactant list and remove spaces
        print(total)
        seperated = balanced_formula.split("-->")
        for x in range(len(seperated[0].split("+"))):       ## for amount of items before the -->, 1 is after and 0 is before since seperated splits into before and after the arrow
            total_balanced.append(seperated[0].split("+")[x].replace(" ", ""))     ## for the items before the -->, append them into a reactant list and remove spaces
        for x in range(len(seperated[1].split("+"))):       ## for amount of items before -->
            total_balanced.append(seperated[1].split("+")[x].replace(" ", "")) ## for the items after -->, append them into a reactant list and remove spaces
        print(total_balanced)
        ## total looks like this [element1, element2, ...]
        # total = list(dict.fromkeys(total))

        
        # Iterate through the list and create inputs within the columns
        selection = st.radio("Calculate via known Mass Or Mol", ["Via Mol", "Via Mass in Grams"])
        
        counter = -1
        Known = st.radio("Calculate via known Mass Or Mol", total_balanced)
        if selection == "Via Mol":
            input_value = st.number_input(f"Mol {Known}" ,key = counter)
        else:
            input_value = st.number_input(f"Mass (g) {Known}" ,key = counter)
        
        if st.button("Calculate"):
            with st.spinner('Loading... this might take a minute'):
                solved = Chem.Stoichiometry(total,selection, Known, input_value, array,total_balanced )  ## total is the list of molecules, selection is if we are going to calculate by mol or mass, Known is what molecule is known, and input_value is the amount of mass or mol
                            ## this function will only eat a list of molecules and if you feed it anything out it will not work
                print(solved)   ## solved returns as a dict like this {"NaOH" : [mol,mass,molar mass]} but unbalanced so we need to turn the total list to also be unbalanced
            print(solved)
            try:
                array = list(array)
                total_balanced = list(solved.keys())
                print("ee")
                print(total_balanced)
                for x in range(total_balanced):
                    total_balanced[x] = str(array[x]) + total_balanced[x] 
                print(total_balanced)
            except Exception as error:
                print(error)
            for y in range(len(solved)):
                print(total_balanced[y])
                st.write(total_balanced[y])     ## since the molecules were defaulted to 0 in the first code in the function, we have to add the amount of molecules back from the array... This is autistic
                col1 , col2, col3 = st.columns(3)
                with col1:
                    st.write("Mol: " , solved[total_balanced[y]][0], " mol")
                with col2:
                    st.write("Mass: " , solved[total_balanced[y]][1], "g")
                with col3:
                    st.write("Molar Mass: " , solved[total_balanced[y]][2], "g/mol")
    except:
        pass
    # for col, item in zip(cols, array):
    #     if error == 0:
    #         counter += 1
    #         # Create an input for each item within the column
    #         # input_value = col.text_input(f"Mol {item}" ,key = counter)
    #         st.write()
    
    #         #choose = col.selectbox(str(total[counter]), ("Calculate Mol and Mass", "Known Mol", "Known Mass in grams"), key = str(counter*300+counter))
    #         # if choose == "Known Mol":
    #         #     mol = col.text_input(f"Mol {total[counter]}" ,key = "a" + str(counter))
    #         # if choose == "Known Mass in":
    #         #     mol = col.text_input(f"Mol {total[counter]}" ,key = "a" + str(counter))
    #         # You can do something with the input value if needed
    #         # st.write(f"You entered for {item}: {input_value}")
    # ##if (st.button("submit")):
elif (st.session_state.status == 10):
    st.text("Dilution Calculator")
    v2 = st.number_input("Volume 2 (l)", step=0.000001,  format="%0.6f")
    c1 = st.number_input("Concentration 1 (mol/V)", step=0.000001,  format="%0.6f")
    c2 = st.number_input("Concentration 2 (mol/V)", step=0.000001,  format="%0.6f")
    if (st.button("Enter")):
        v1, water = Chem.Dilution("mol", c1, c2, v2)
        st.text("You will need {} liters of your concentration and {} liters \nof water to achieve {} liters with a concentration of {} mol/V".format(round(v1,6),round(water,6),round(v2,6),round( c2,6)))
elif (st.session_state.status == 11):
    formula = st.text_input(label="Chemical reaction equation, does not need to be balanced. For example: Na + H2O")
    if st.button("Enter"):
        try:
            formula = formula.split("-->") 
            try:
                print(" reaction --> " , formula)
                sol, molecules_types = Chem.BasicReactionPredictor(formula[0])
                split_molecules = formula[0].split("+")
                split_molecules = [string.replace(" ", "") for string in split_molecules]
            except:
                sol, molecules_types = Chem.BasicReactionPredictor(formula)
                split_molecules = formula.split("+")
                split_molecules = [string.replace(" ", "") for string in split_molecules]
            
            st.text("This reaction will create: " + str(sol))
            for x in range(len(split_molecules)):
                
                st.write(split_molecules[x], " : " , molecules_types[x])
            
        except:
            st.write("Seems like this isnt a reaction, or i just havent coded it yet  ¯\_(ツ)_/¯ )  ")
            pass
elif (st.session_state.status == 12):
    text = st.header("Eular Division")
    st.write("Divisible poly")
    ncol = st.number_input("Highest Degree of polynomial", 0) + 1
    cols = st.columns(ncol)
    input_list = []
    
    for i, x in enumerate(cols):
        print("Check cols")
        print(len(cols))
        input_val = x.number_input(f"Input # {len(cols)-i}", key="main" + str(len(cols)-i))
        if i == 0:
            # x.write(f"{input_val}")
            the_x = convert_to_superscript(f"x^{len(cols)-1}")
            x.write(f"{input_val}" + the_x)
        elif len(cols)-(i+1)  == 1:
            x.write(f"{input_val}x")
        elif i == len(cols)-1:
            x.write(f"{input_val}")
        else:
             the_x = convert_to_superscript(f"x^{len(cols)-(i+1)}")
             x.write(f"{input_val}" + the_x)
        input_list.append(input_val)
    input_list.reverse()
    st.write("divider poly")
    full_poly = ""
    ## write out the full polynomial
     ## throw into func
    full_poly = convert_to_polynomial(input_list)
    st.write(full_poly)
    
    ncol2 = st.number_input("Highest Degree of divisor polynomial", 0) + 1
    cols2 = st.columns(ncol2)
    input_listdivisor = []
    full_divisor = ""
    for y, z in enumerate(cols2):  ## enumerate works with index and item.. y represents index and z represents item
        ## ["A", "B"] --> y is 0 and 1 and z is     "A" and "B"
        ## we can use this instead of indexing every time
        input_val = z.number_input(f"Input divisor # {len(cols2)-y} x", key="divisor" + str(len(cols2)-y))
        if y == 0:
            the_x = convert_to_superscript(f"x^{len(cols2)-1}")
            z.write(f"{input_val}" + the_x)
        elif len(cols2) - (y+1) == 1: 
            z.write(f"{input_val}x")
        elif y == len(cols2) - 1:
            z.write(f"{input_val}")
        else:
             the_x = convert_to_superscript(f"x^{len(cols)-(y+1)}")
             z.write(f"{input_val}" + the_x)
        input_listdivisor.append(input_val)
     ## throw into func
    input_listdivisor.reverse()
    full_divisor = convert_to_polynomial(input_listdivisor)
    st.write(full_divisor)
    
    if st.button("Calculate"):
        quotient, rest, quotients, polynomials, subtractors = Math.take2(input_list, input_listdivisor)
        print("####")
        print(quotient)
        full_quotient = ""
        full_rest = ""
        ## throw into func
        full_quotient = convert_to_polynomial(quotient)
        full_rest = convert_to_polynomial(rest)
        st.write("q(r)= " + full_quotient)
        st.write("r(x)= " + full_rest)
        st.header("Steps")
        st.write("")
        st.write( full_poly)
        spacing = ""
        for x in range(len(subtractors)):
            
            st.write("  " + f"**-({convert_to_polynomial(subtractors[x])})**")
            st.write("=>" + convert_to_polynomial(polynomials[x]))
            
        print(subtractors)
        print(polynomials)
        print(quotients)

if st.session_state.status != 0:
    

    if(st.button('Return')):
        if st.session_state.status in [3,4]:
            st.session_state.status = 1    
        elif st.session_state.status in [5,6,7,8,9,10,11]:
            st.session_state.status = 2
        else:  
            st.session_state.status = 0
        st.rerun()
