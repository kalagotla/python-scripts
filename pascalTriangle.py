#-----------------------------------------------------------------------------#
def computePascalTriangle(n):
# Computes a n-line Pascal's Triangle and stores in a 2D list. Outer list 
# is a list of rows and the inner lists contain the elements of the triangle.
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
    
    nRows = len(triangle)
    
    # Determine the field width from the larget entry in the triangle
    maxElement = max(triangle[nRows-1]) 
    fieldWidth = len(str(maxElement))

    numFields = 2*nRows-1
    for row in range(nRows):        
        
        blankFields = numFields//2 - row
        line = " "*blankFields*fieldWidth
        
        for col in range(row+1):
            # Print each entry and pad extra spaces after
            field = "{{:^{}}}".format(fieldWidth)
            line += field.format(triangle[row][col]) + " "*fieldWidth
            
        print(line)

#-----------------------------------------------------------------------------#

pt = computePascalTriangle(5)
printPascalTriangle(pt)
