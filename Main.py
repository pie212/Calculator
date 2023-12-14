import streamlit as st
import time
import Math
import Chem
import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
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


#status = st.radio("Select Functions: ", ('Singl  Function', 'Solve 2 functions', 'Unlimited functions'))
if st.session_state.status == 0:      # 0 = home, 1-4 math, 5 Chem
    # st.title("Function Calculator")
    # st.write("Version 1.0")
    # st.write("Home Page")
    st.markdown("<h1 style='text-align: center; color: red;'>Function Calculator</h1>", unsafe_allow_html=True)
    st.markdown("")
    st.markdown("")
    st.markdown("")
    st.markdown("")
    st.markdown("")
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        button1 = st.button("Function Info")
    with col2:
        button2 = st.button("Single Function")

    with col3:
        button3 = st.button("Solve 2 functions")

    with col4:
        button4 = st.button("Unlimited functions")
    st.markdown("")
    st.markdown("")
    
    Chemistry = st.button("Chemistry")
    if button1:
        st.session_state.status = 1
        st.rerun()
    elif button2:
        st.session_state.status = 2
        st.rerun()
    elif button3:
        st.session_state.status = 3
        st.rerun()
    elif button4:
        st.session_state.status = 4
        st.rerun()
    elif Chemistry:
        st.session_state.status = 5
        st.rerun()
    else:
        st.session_state.status = 0
        st.rerun()
    

elif (st.session_state.status == 1):
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




elif (st.session_state.status == 2):
    
    name = st.text_input("Enter Function (ax^2 + bx+ c)")
 

    if(st.button('Submit')):
        result = name
        print(result)
        d,x1,x2 = Math.Solution(result)
        st.success("X1 = {} , X2 = {}".format(x1,x2))
        a,b,c = Math.FindTerms(result)
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
        plt.plot(x, y, label='f(x) = x^2 + 4x + 3')
        plt.title('Graph of f(x) = x^2 + 4x + 3')
        plt.xlabel('x')
        plt.ylabel('f(x)')
        plt.axhline(0, color='black',linewidth=0.5)
        plt.axvline(0, color='black',linewidth=0.5)
        plt.legend()
        
        # Display the plot in Streamlit
        st.pyplot(plt)
elif (st.session_state.status == 3):
    func1 = st.text_input("Function 1")
    func2 = st.text_input("Function 2")
    col1, col2, col3, col4 = st.columns(4)
    # with col1:
    #     st.session_state.graphxmin = st.slider("X min", -1000, 1000,int(st.session_state.graphxmin))
    #     st.session_state.graphxmin = sliderxmin
    #     sliderxmin = st.session_state.graphxmin 
    # with col2:
    #     sliderxmax = st.slider("X max", -1000, 1000)
    #     #st.session_state.graphxmax = sliderxmax
    # with col3:
    #     sliderymin = st.slider("Y min", -1000, 1000)
    #     #st.session_state.graphymin = sliderymin
    # with col4:
    #     sliderymax = st.slider("Y max ", -1000, 1000)
    #     #st.session_state.graphymax = sliderymax
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        graphxminbox = st.text_input("X min", st.session_state.graphxmin)
        st.session_state.graphxmin = int(graphxminbox)
    with col2:
        graphxmaxbox = st.text_input("X max", st.session_state.graphxmax)
        st.session_state.graphxmax = int(graphxmaxbox)
    with col3:
        graphyminbox = st.text_input("Y min", st.session_state.graphymin)
        st.session_state.graphymin = int(graphyminbox)
    with col4:
        graphymaxbox = st.text_input("Y max", st.session_state.graphymax)
        st.session_state.graphymax = int(graphymaxbox)
    # # st.session_state.graphxmin = st.slider("X min", -1000, 1000)
    # # st.session_state.graphxmax = st.slider("X max", -1000, 10000)
    # # st.session_state.graphymin = st.slider("Y min", -1000, 1000)
    # # st.session_state.graphymax = st.slider("Y max ", -1000, 1000)
    # print(st.session_state.graphxmin)
    # print(st.session_state.graphymin)
    # print(st.session_state.graphxmax)
    # print(st.session_state.graphymax)
    t = True
    if(t):
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

elif (st.session_state.status == 4):
    funcsL = []
    level = st.slider("Select amount of functions", 2, 1000)
    variables = {}
    for i in range(1, level + 1):
        st.text_input(f"Function {i}", key=f"input_{i}")
    if(st.button('Submit')):
        for i in range(1, level + 1):
            input_value = st.session_state[f"input_{i}"]
            funcsL.append(input_value)
        print(funcsL)
        done, checked, totals,intersects,SaveDebugList,elapsedTime,solutions, amfuncs = Math.solvemore3func4(*funcsL)
        
        
        plt.figure(figsize=(8, 6))
        kylist = []
        for x in intersects:
            kylist.append(x)
            st.success("All functions intersect at X = " + str(x))
        if intersects == {}:
            st.warning("No common intersects found for all the functions")

        # Generate x values
        
        for i in range(len(funcsL)):
            if intersects != {}:
                x = np.linspace(kylist[0] -30, kylist[0]+30)  # Adjust the range of x values as needed
            else:
                x = np.linspace(-10, 10)  # Adjust the range of x values as needed
            
            a,b,c = Math.FindTerms(funcsL[i])
            f = lambda x: float(a)*x**2 + float(b)*x + float(c)
            y = f(x)
            plt.plot(x, y , label = funcsL[i])#label='f(x) = x^2 + 4x + 3')
            
        # Calculate y values using the function
        

        # Create a plot using Matplotlib
        
        
        plt.title('Graph of f(x) = x^2 + 4x + 3')
        plt.xlabel('x')
        plt.ylabel('f(x)')
        plt.axhline(0, color='black',linewidth=0.5)
        plt.axvline(0, color='black',linewidth=0.5)
        plt.legend()
        st.pyplot(plt)
        #st.success("All functions intersect in X = " + str(intersects[0]))
elif (st.session_state.status == 5):
    st.header("Molecule Balancer")
    element1 = st.text_input("Element 1")
    element2 = st.text_input("Element 2")
    if (st.button("submit")):
        #try:
        sols = Chem.MoleculeStable(element1,element2)
        for x in range(len(sols)):
            st.write(sols[x])
        # except:
        #     st.warning("operation failed!")
        
      
if st.session_state.status != 0:


    if(st.button('Return')):
        st.session_state.status = 0
        st.rerun()
