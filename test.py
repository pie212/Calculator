import streamlit as st
import time
import cleaned
# Write text

 
## MAIN screen

menuText = st.empty()

# Concatenate the text
text_content = "Function calculator\nV1.0\nHome Page"

menuText.text(text_content)

#status = st.radio("Select Functions: ", ('Single Function', 'Solve 2 functions', 'Unlimited functions'))
status = st.selectbox(
    "Options",
    ("Single Function", "Solve 2 functions", "Unlimited functions"),
    index=None,
    placeholder="Select an option...",
)

# conditional statement to print 
# Male if male is selected else print female
# show the result using the success function
if status == 'Single Function' or status == "Solve 2 functions" or status == "Unlimited functions":
    menuText.text("")
if (status == 'Single Function'):
    
    name = st.text_input("Enter Your name")
 

    if(st.button('Submit')):
        result = name
        print(result)
        d,x1,x2 = cleaned.Solution(result)
        st.success("X1 = {} , X2 = {}".format(x1,x2))
elif (status == "Solve 2 functions"):
    func1 = st.text_input("Function 1")
    func2 = st.text_input("Function 2")
 

    if(st.button('Submit')):
        func1 = func1
        func2 = func2
        print(func1)
        d,x1,x2 = cleaned.solve2func(func1,func2)
        st.success("X1 = {} , X2 = {}".format(x1,x2))

elif (status == "Unlimited functions"):
    funcsL = []
    level = st.slider("Select the level", 1, 1000)
    variables = {}
    for i in range(1, level + 1):
        st.text_input(f"Function {i}", key=f"input_{i}")
    if(st.button('Submit')):
        for i in range(1, level + 1):
            input_value = st.session_state[f"input_{i}"]
            funcsL.append(input_value)
        print(funcsL)
        done, checked, totals,intersects,SaveDebugList,elapsedTime,solutions, amfuncs = cleaned.solvemore3func4(*funcsL)
        for x in intersects:
            st.write(x)
        #st.success("All functions intersect in X = " + str(intersects[0]))



