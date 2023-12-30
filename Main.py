import streamlit as st
import time
import Math
import Chem
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from io import BytesIO
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
    if funcinfo:
        st.session_state.status = 3
        st.rerun()
    elif funcompare:
        st.session_state.status = 4
        st.rerun()
elif (st.session_state.status == 2):
    st.header("Chemistry")
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        moleculebalancer = st.button("Molecule Balancer")
    with col2:
        reaction = st.button("Reaction balancer (work in progress)")
    with col3:
        CandHinator = st.button("C and H inator")
    if moleculebalancer:
        st.session_state.status = 5
        st.rerun()
    elif reaction:
        st.session_state.status = 6
        st.rerun()
    elif CandHinator:
        st.session_state.status = 7
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
    st.write("work in progress")
elif (st.session_state.status == 7):
    col1,col2,col3, col4,col5, col6 = st.columns([0.02,0.1,0.02,0.1,0.02,0.1])
    with col1:
        st.header("C")
    with col2:
        Catoms = st.text_input(value=0, key=1, label="")
    with col3:
        st.header("H")
    with col4:
        Hatoms = st.text_input(value=0, key=2, label="")
    with col5:
        st.header("O")
    with col6:
        Oatoms = st.text_input(value=0, key=3, label="")
    if (st.button("submit")):
        try:
            Hatoms = int(Hatoms)
            Catoms = int(Catoms)
            Oatoms = int(Oatoms)
        except:
            st.rerun()
        if Oatoms == 0 and Catoms > 0 and Hatoms > 0:
            with st.spinner('Loading...'):
                alkaan, alkaanvisual,alkaanname, alkeenstructures, alkeenvisuals, alkenennames = Chem.CandHinator(Catoms,Hatoms,Oatoms)
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
                    st.header(alkaan)
                    st.header(alkaanname)
                with col2:
                    try:
                        st.image(alkaanvisual)
                    except:
                        st.write("Failed to get image!")
            if alkeenstructures:
                print(alkeenvisuals)
                for x in range(len(alkeenstructures)):
                    col1, col2 = st.columns(2)
                    with col1:
                        st.write(alkeenstructures[x])
                        st.write(alkenennames[x])
                    with col2:
                        st.image(alkeenvisuals[x])
            else:
                st.header("No Alkenen")
        if Catoms == 0 or Hatoms == 0:
            st.header("No Alkanen")
            st.header("No Alkenen")


        
if st.session_state.status != 0:
    

    if(st.button('Return')):
        if st.session_state.status in [3,4]:
            st.session_state.status = 1    
        elif st.session_state.status in [5,6]:
            st.session_state.status = 2
        else:  
            st.session_state.status = 0
        st.rerun()
