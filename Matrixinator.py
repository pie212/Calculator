#Matrixinator
import sympy as sp

def x_y_solver(variable_1, variable_2):
    

    # Solve for float solutions
    variable_1 = sp.Rational(variable_1).limit_denominator()
    variable_2 = sp.Rational(variable_2).limit_denominator()
    # Define the variables
    x, y = sp.symbols('x y')

    # Define the Diophantine equation
    equation = sp.Eq(variable_1 * x + variable_2 * y, 0)

    # Solve for integer solutions
    solutions = sp.diophantine(equation)

    # Convert the solutions to a list
    sols_list = list(solutions)

    # Access the first solution
    first_sol = sols_list[0]  # First tuple of the general solution
    
    
    # Substitute t_0 = 1 (or use the parameter found in the solution)
    parameter = list(first_sol[0].free_symbols)[0]  # Extract the parameter symbol dynamically
    solutions= tuple(expr.subs(parameter, 1) for expr in first_sol)
    
    return solutions

def ChemBalancerRaw_to_Matrix(reactants, products, TotalAtoms):
    ## for ex if we have Na + H2O --> NaOH + H2
    ## reactants == [{"Na" : 1},{"H": 2, "O": 1}]
    ## products == [{"Na" : 1, "O" : 1, "H" : 1},{"H": 2}]
    ## set up this goofy aah matrix
    ## FORM
    ## Atom_1 [Molecule_1 Molecule_2 Molecule_n 0]
    ## Atom_2 [Molecule_1 Molecule_2 Molecule_n 0]
    ## Atom_n [Molecule_1 Molecule_2 Molecule_n 0]
    ## Total_Atoms == ["Atom1", "Atom2"]
    ## Matrix Form:
    ## Matrix = [
    ##     [Row Atom_1]
    ##     [Row_Atom_2]
    ##     [Row_Atom_n]
    ##          ]
    matrix = []
    for x in range(len(TotalAtoms)):
        current_atom = TotalAtoms[x]
        new_row = []
        for y in range(len(reactants)): ## cycle through reactants molecules
            if current_atom in reactants[y]:
                new_row.append(reactants[y][current_atom])
            else:
                new_row.append(0)
       

        for y in range(len(products)): ## cycle through reactants molecules
            if current_atom in products[y]:
                new_row.append(-products[y][current_atom])
            else:
                new_row.append(0)
        new_row.append(0) ## this is to add the [x y z | 0] which is always going to be 0, as we moved all the terms away to the side, and due to the fact it is a chemical reaction, no loose integers are left!
        matrix.append(new_row)
    for x in range(len(matrix)):
        print(TotalAtoms[x], matrix[x])

def zero_rows_down(matrix):
    # form matrix
        # matrix = [
        # [row1]
        # [row2]
        # [rown]
        # ]       
## first bring all 0 rows to the bottom
## use set method to see if all elements are the same
## set(list) --> returns the amount of items in a list without duplicats, if this is 1, it means only the same element is present!
    rows = len(matrix)

    zero_rows = []
    non_zero_rows = []
    for x in range(rows):
        
        
        if len(set(matrix[x])) == 1 and matrix[x][0] == 0: ## first condition --> all elements are the same, second condition --> first element is 0, if first is 0 and condition 1 is true everything is 0
            print("This is a 0 row")
            zero_rows.append(matrix[x])
        else:
            non_zero_rows.append(matrix[x])
    matrix = []
    for x in range(len(non_zero_rows)):
        matrix.append(non_zero_rows[x])
    for x in range(len(zero_rows)):
        matrix.append(zero_rows[x])
    return matrix
def Echelon_Form(matrix):
        # form matrix
        # matrix = [
        # [row1]
        # [row2]
        # [rown]
        # ]       
## first bring all 0 rows to the bottom
## use set method to see if all elements are the same
## set(list) --> returns the amount of items in a list without duplicats, if this is 1, it means only the same element is present!
    matrix = zero_rows_down(matrix)
    for x in range(len(matrix)):
        print(matrix[x])
    print("#######")

    ## great now we have the 0's out of the way, but we will have to check this after every operation!
    ## lets move on to the first spill!
    rows = len(matrix)
    columns = len(matrix[0]) 
    #print("##")
    ## both these will always stay the same! all rows have the same amount of columns, so this will suffice!
    ## we wont move the colums with a 1 to the top, since for the computer it will go just as fast, just a matter of solving a x+y type equation!
    current_spill_number = 0
  
    for x in range(len(matrix)): ## choose our first row!
        
        if len(set(matrix[x])) == 1 and matrix[x][0] == 0: ## stops once we encounter 0 row
            #print("start 0")
            
            break
        for y in range(len(matrix[x])):    ## here we will do an extra control to get the spill, checking the first number to not be 0, until it isnt
            if matrix[x][y + x] != 0:
                spill_row = matrix[x] # the row itself [elements]
                spill = matrix[x][y+x] # the spill value
                spill_column = y+x # the colum of the spill
                
                break
        ## ok so now we have our spill, lets sort our operations for each row
        operations_this_spill_cycle = []
        for z in range(len(matrix)): ## so for each row, except our matrix
            
            if x != z: ## now we cycle through each row, to make it 0
                       ## here we check if we are working on the same row, if we arent we pass
                
                ## how do we get to 0?
                row_to_change = matrix[z]
                #print(row_to_change)
                
                if matrix[z][spill_column] == 0: ## if the elemnt is already 0 we dont needa do anything
                    #print("passed")
                    pass
                    
                else:
                    ## equation is a *row_z + b *row_x = 0
                    ## row z is our row that changes
                    ## row x is our spill row
                    a,b = x_y_solver(matrix[z][spill_column], spill)
                    
                    if b > 0:
                        operations_this_spill_cycle.append(str(a) + " * row" + str(z+1) + "+" +  str(b) +"* row" + str(x+1)) ## +1 added for readabilty, so we dont need to start from 0
                    else:
                        operations_this_spill_cycle.append(str(a) + " * row" + str(z+1)  +  str(b) +"* row" + str(x+1))
                    updated_row = []
                    
                    
                    for y in range(columns):
                        updated_row.append(row_to_change[y] *a + spill_row[y] * b)
                        
                        
                    matrix[z] = updated_row
                    
            else:
                ##print("Passed, same row")
                pass
            # print("Show partial matrix")
            # print("Row " , str(x))
            # print("run", str(z))
            # print("Row", str(spill_row))
        ## now we have our new matrix!
        ## lets quickly move down our 0 rows in case we have them
        # need to move the iterator down 1
        matrix_before = matrix
        matrix = zero_rows_down(matrix)
        if matrix_before != matrix:
            x-=1
        print("Matrix operations:")
        for i in range(len(operations_this_spill_cycle)):
            print(operations_this_spill_cycle[i])
        print("Matrix")
        for m in range(len(matrix)):
            print(matrix[m])
        print("#######")
        simplified_matrix = []
    for x in range(rows):
        spill_found = False
        new_row = []
        for y in range(columns):
            
            if matrix[x][y] == 0 and spill_found == False:
                pass
                new_element = 0
            elif matrix[x][y] !=0 and spill_found == False:
                spill_found = True
                divider = matrix[x][y]
                new_element = matrix[x][y] / divider
            else:
                new_element = matrix[x][y] / divider
            new_row.append(new_element)
        simplified_matrix.append(new_row)
    
    # for x in range(rows):
    #     print(simplified_matrix[x])
    return simplified_matrix


    ## we need to simply the matrix  . it is 12:48 am i will do this tmrw   
def Solve_Echelon(matrix):
    rows = len(matrix)
    columns = len(matrix[0])
    variables = []
    for x in range(columns):
        values = []
        for y in range(rows):
            values.append(matrix[y][x])
        variables.append({x:values})
    print(variables)
        
    
#ChemBalancerRaw_to_Matrix([{"Na" : 1},{"H": 2, "O": 1}], [{"Na" : 1, "O" : 1, "H" : 1},{"H": 2}], ["Na", "O", "H"])
# Echelon_Form(
#     [
#         [1,2,1,2,3,5],
#         [4,3,6,5,4,2],
#         [3,6,5,4,32,2]
#     ]
# )

Solve_Echelon(
    [
        [1,0,0,2,3],
        [0,1,0,5,3],
        [0,0,1,2,5],
    ]
)