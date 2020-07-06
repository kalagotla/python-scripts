#-----------------------------------------------------------------------------#
def computePascalTriangle(n):
# Computes a n-line Pascal's Triangle
    triangle = [[1]]

    for i in range(1,n):
        row = [1]

        for j in range(1,i):         
            prevRow = triangle[i-1]
            row.append(prevRow[j-1] + prevRow[j])

        row.append(1)

        triangle.append(row)    

    return(triangle)
#-----------------------------------------------------------------------------#
    
#-----------------------------------------------------------------------------#
def printPascalTriangle(triangle):
# Pretty prints a previously computed Pascal's Triangle
# Assumes triangle is a nested list of lists of the form
#   [ [1], [1 2], [1 2 3], [1 2 3 4], ..., [1 2 3 ... n-1 n] ]
    
    n = len(triangle) # Get the number of rows in the triangle
    
    # Determine the number of digits for the largest value in the triangle
    maxElement = max(triangle[n-1]) 
    fieldWidth = len(str(maxElement))

    # Compute the number of values and blanks for the last row of the triangle
    numFields = 2*n-1
    for i in range(n):        
        # Add leading spaces based on which row of the triangle we are printing
        blankFields = numFields//2 - i
        line = " "*blankFields*fieldWidth
        
        for j in range(i):
            # Print each entry and pad extra spaces after
            field = "{{:^{}}}".format(fieldWidth)
            line += field.format(triangle[i][j]) + " "*fieldWidth
            
        print(line)

#-----------------------------------------------------------------------------#

pt = computePascalTriangle(5)
printPascalTriangle(pt)
