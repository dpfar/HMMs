# KTH EECS
# DD2380 - Artificial Intelligence
# HMM1
# Beatrice Lovely and Daniel Pérez Felipe

import sys


def data2matrix(data_string):
    # Given a single data input corresponding to a matrix, a matrix is allocated and created

    rows = int(data_string[0])  # get row value
    cols = int(data_string[1])  # get col value

    matrix = [[0 for x in range(cols)] for y in range(rows)]  # Creates empty matrix

    r = c = j = 0  # auxiliary variables

    while r < rows:  # Matrix is filled with the given data
        while c < cols:
            matrix[r][c] = float(data_string[j+2])
            c = c + 1
            j = j + 1
        c = 0
        r = r + 1

    return matrix


def data2emission_sequence(data_string):
    # Given a single data input corresponding to an emission sequence, data is obtained
    # It returns the emission sequence, the number of emissions and the number of emission different_types

    n_emissions = int(data_string[0])  # get number of emissions
    n_different_types = max(data_string[1:]) + 1
    data = data_string[1:]

    return data, n_emissions, n_different_types


def product_matrices(matrixX, matrixY):
    # Computes the product of two matrices
    # Get cols and rows for matrix creation

    colsY = len(matrixY[0])
    rowsX = len(matrixX)

    matrix_prod = [ [ 0 for x in range(colsY) ] for y in range(rowsX) ] # empty matrix is created for storage

    for l in range(len(matrixX)):
        for j in range(len(matrixY[0])):
            for k in range(len(matrixY)):
                matrix_prod[l][j] += matrixX[l][k] * matrixY[k][l]

    return matrix_prod


def element_wise_product(matrixX,matrixY):
    # Computes the element wise product of two matrices

    colsY = len(matrixY[0])
    rowsX = len(matrixX)

    matrix_elem_wise = [[0 for x in range(colsY)] for y in range(rowsX)]

    for k in range(len(matrixX)):
        for j in range(len(matrixY[0])):
            matrix_elem_wise[k][j] = matrixX[k][j] * matrixY[k][j]

    return matrix_elem_wise


def colFromMatrix(matrix,pos):
    # Returns a single column matrix from an specific column of the input Matrix

    rows = len(matrix)
    matrix_col = [[0] for j in range(rows)]  # A matrix with one column is created

    for j in range(len(matrix)):  # We move along all rows, but only move along one column(being "pos")
        matrix_col[j][0] = matrix[j][pos]

    return matrix_col


def matrix2str(matrix):
# Returns a matrix in the given format( nºrows, nºcols, matrix_data). For printing

    aux_var = 0
    col = len(matrix[0])  # Gets column number
    row = len(matrix)  # Gets row number
    data_stream = '' + str(row) + " " + str(col)

    for r in range(len(matrix)): # We move along the rows
        for c in range(len(matrix[0])):  # We move along the columns
            aux_var = str(matrix[r][c])
            data_stream += " " + aux_var

    return data_stream


def transpose(matrix):
    # Given a matrix, the transpose matrix is generated

    rows = len(matrix)
    cols = len(matrix[0])

    matrix_transpose = [[0 for x in range(rows)] for y in range(cols)]

    for r in range(rows):  # We move along the rows
        for c in range(cols):  # We move along the columns
            matrix_transpose[c][r] = matrix[r][c]

    return matrix_transpose


if __name__ == "__main__":

    # file = open("hmm2_01.in")  # TODO. PYCHARM
    file = sys.stdin  # TODO: KATTIS

    # Data is obtained from input files and classified
    values = file.readlines()
    avals = values[0].split()
    bvals = values[1].split()
    pivals = values[2].split()
    seqvals = values[3].split()
    for i in range(len(seqvals)):
        seqvals[i] = int(seqvals[i])

    # Input data is transformed into matrix form
    A = data2matrix(avals)
    B = data2matrix(bvals)
    PI = data2matrix(pivals)
    # String Input data is transformed into data values
    sequence, num_of_emissions, different_types = data2emission_sequence(seqvals)

    # variable that stores all values of alfa for each given time
    alfa = [0 for i in range(num_of_emissions)]

    # alfa_1 is computed
    alfa_previous = element_wise_product( transpose(PI), colFromMatrix(B,sequence[0]) )
    probability_of_sequence = 0
    for i in range(len(alfa_previous)):
        probability_of_sequence += alfa_previous[i][0]
    alfa[0] = probability_of_sequence
    # print(alfa[0])

    t = 1
    while t < num_of_emissions:
        alfa_next = [[0] for i in range(different_types)]
        for i in range(different_types):  # We compute the summatory product
            alfa_next[i][0] = product_matrices(transpose(alfa_previous), colFromMatrix(A, i))[0][0]
            alfa_next[i][0] = alfa_next[i][0] * colFromMatrix(B, sequence[t])[i][0]

        #probability for time "t" is stored
        probability_of_sequence = 0
        for i in range(len(alfa_next)):
            probability_of_sequence += alfa_next[i][0]
        alfa[t] = probability_of_sequence
        # print(alfa[t])

        t = t + 1
        alfa_previous = alfa_next

    print(probability_of_sequence)
