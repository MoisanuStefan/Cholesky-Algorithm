import random
import math
import numpy as np
import scipy.linalg
import numpy.linalg as npl

NR_OF_DECIMALS = 3


def compute_l_pp(a_matrix, d_vector, p):
    sqrt_argument = d_vector[p] - sum(pow(a_matrix[p][j], 2) for j in range(0, p))
    if sqrt_argument < 0:
        return -1
    return math.sqrt(sqrt_argument)


def compute_l_ip(a_matrix, e, p, i):
    return (a_matrix[p][i] - sum(a_matrix[i][j] * a_matrix[p][j] for j in range(0, p))) / (
        a_matrix[p][p] if math.fabs(a_matrix[p][p]) > e else exit("Can't make division"))


def print_matrix(matrix, nr_of_decimals, matrix_type="full", transposed=False):
    max_nr = np.array(matrix).max()
    max_len = len(("{:." + str(nr_of_decimals) + "f}").format(max_nr))
    format_string = "{:" + str(max_len + 2) + "." + str(nr_of_decimals) + "f}"
    for i in range(n):
        for j in range(n):
            if matrix_type == "inferior" and j > i:
                print(format_string.format(0), end='')
            elif matrix_type == "superior" and j < i:
                print(format_string.format(0), end='')
            else:
                if transposed:
                    print(format_string.format(matrix[j][i]), end='')
                else:
                    print(format_string.format(matrix[i][j]), end='')
        print()


def compute_triangular_matrix_det(matrix):
    ln = len(matrix)
    prod = 1
    for i in range(ln):
        prod = prod * matrix[i][i]
    return prod


def solve_system_with_substitution(a_matrix, b_vector, type):
    n = len(a_matrix)
    x_vector = [None] * n
    x_vector[0] = b_vector[0] / (a_matrix[0][0] if math.fabs(a_matrix[0][0]) > e else exit("Can't make division"))
    x_vector[-1] = b_vector[-1] / (a_matrix[-1][-1] if math.fabs(a_matrix[-1][-1]) > e else exit("Can't make division"))
    i_range = range(1, n) if type == 'direct' else range(n - 2, -1, -1)
    for i in i_range:
        j_range = range(0, i) if type == "direct" else range(i + 1, n)
        x_vector[i] = (b_vector[i] - sum(a_matrix[i][j] * x_vector[j] for j in j_range)) / (
            a_matrix[i][i] if math.fabs(a_matrix[i][i]) > e else exit("Can't make division"))

    return x_vector


def compute_euclidian_norm(matrix):
    return math.sqrt(sum(math.pow(element, 2) for element in matrix))


def get_e_j(j, n):
    return [0 if j != i else 1 for i in range(n)]


# returns true if mat is symmetric
def isSymmetric(mat, n):
    sym = True
    for i in range(n):
        for j in range(i, n):
            if mat[j][i] != mat[i][j]:
                sym = False
    return sym


def recreate_a_init(a_matrix, d_vector):
    n = len(a_matrix)
    for i in range(n):
        for j in range(i + 1, n):
            a_matrix[j][i] = a_matrix[i][j]
    for i in range(n):
        a_matrix[i][i] = d_vector[i]


# generating or reading the matrix from file or as an input, considering the user's choice
print(
    "Press 1 if you want to insert your matrix, press 2 if you want to read it from a file, press 3 if you want it to be randomly generated: ")
option = int(input())
sym = True
if (option == 1):  # input + checking its symmetry
    m = int(input("Enter precision exponent (m) :"))
    e = pow(10, -m)
    n = int(input("Enter the matrix's dimension: "))
    a_matrix = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            a_matrix[i][j] = float(input("Enter A element: "))
    b_vector = [None] * n
    for i in range(n):
        b_vector[i] = float(input("Enter b element: "))
    d_vector = [None] * n
    for i in range(n):
        d_vector[i] = a_matrix[i][i]
    sym = isSymmetric(a_matrix, n)
    if sym == True:
        print("A: ")
        print_matrix(a_matrix, NR_OF_DECIMALS)
        print("b: ")
        print_matrix(b_vector, NR_OF_DECIMALS)
    else:
        print("Modify your matrix! It is not symmetric")
elif option == 2:  # reading the matrix from the file and cheking its symmetry
    with open('matrix') as f:
        m = int(f.readline())
        e = pow(10, -m)
        n = int(f.readline())
        a_matrix = [[float(x) for x in line.split()] for line in f]
        d_vector = [None] * n
        for i in range(n):
            d_vector[i] = a_matrix[i][i]
    sym = isSymmetric(a_matrix, n)
    if sym == True:
        with open('b_vector') as f:
            b_vector = [float(x) for x in f]
        print("Your matrix:")
        print_matrix(a_matrix, NR_OF_DECIMALS)
    else:
        print("Modify your file! The matrix is not symmetric")
else:  # randomly generating a symmetric matrix
    m = int(input("Enter precision exponent (m) :"))
    e = pow(10, -m)
    n = int(input("Enter the matrix's dimension: "))
    a_matrix = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            a_matrix[i][j] = float(random.randint(0, 99))
            a_matrix[j][i] = a_matrix[i][j]
    b_vector = [None] * n
    for i in range(n):
        b_vector[i] = float(random.randint(0, 99))
    d_vector = [None] * n
    for i in range(n):
        d_vector[i] = a_matrix[i][i]
    print("A:")
    print_matrix(a_matrix, NR_OF_DECIMALS)
    print("b: ", b_vector)

if sym == True:
    # to do verifica daca matricea e pozitiva definita// videoclip CN2 min 1:37
    print("...starting Cholesky's algorithm:")
    # l_matrix = [[0] * n for _ in range(n)]
    for p in range(n):
        for i in range(p, n):
            a_matrix[i][p] = compute_l_pp(a_matrix, d_vector, p) if i == p else compute_l_ip(a_matrix, e, p, i)
            if i == p and a_matrix[i][p] < 0:
                print("A is not positive definite")
                exit(0)
    print("L:")
    print_matrix(a_matrix, NR_OF_DECIMALS, "inferior")
    print("L_T:")
    print_matrix(a_matrix, NR_OF_DECIMALS, "superior", True)

    a_matrix_inverse_chol = []
    for j in range(n):
        y_star = solve_system_with_substitution(a_matrix, get_e_j(j, n), type='direct')
        x_chol = solve_system_with_substitution(np.transpose(a_matrix), y_star, type='inverse')
        a_matrix_inverse_chol.append(x_chol)

    det_a_matrix = compute_triangular_matrix_det(a_matrix) * compute_triangular_matrix_det(np.transpose(a_matrix))
    print("Det A: ", det_a_matrix)

    print("...solving Ax=b system:")
    y_star = solve_system_with_substitution(a_matrix, b_vector, type='direct')
    x_chol = solve_system_with_substitution(np.transpose(a_matrix), y_star, type='inverse')
    print("x_chol = ", x_chol)

    recreate_a_init(a_matrix, d_vector)

    print("A_inverse_chol: ")
    print_matrix(a_matrix_inverse_chol, NR_OF_DECIMALS)

    print("A_inverse_bibl: ")
    a_matrix_inverse_bibl = npl.inv(a_matrix)
    print_matrix(a_matrix_inverse_bibl, NR_OF_DECIMALS)

    print("||A_inverse_chol - A_inverse_bibl|| = ", npl.norm(a_matrix_inverse_chol - a_matrix_inverse_bibl))

    euclidian_norm = compute_euclidian_norm(np.matmul(a_matrix, x_chol) - b_vector)
    print("|| A_init * x_chol - b ||_2 = ", euclidian_norm)

    p, lu_l_matrix, lu_u_matrix = scipy.linalg.lu(a_matrix)
    print("A = L * U:")
    print("L: ")
    print_matrix(lu_l_matrix, NR_OF_DECIMALS, "inferior")
    print("U: ")
    print_matrix(lu_u_matrix, NR_OF_DECIMALS, "superior")

    y_star = solve_system_with_substitution(lu_l_matrix, b_vector, type='direct')
    x = solve_system_with_substitution(lu_u_matrix, y_star, type='inverse')
    print("Ax=b solution for A = L * U: x = ", x)
