'''
Using q3 from hw6
'''

import numpy as np
from scipy.interpolate import lagrange
from py_ecc.bn128 import G1, G2, add, multiply, curve_order, eq, Z1, Z2, FQ, FQ2
from ape import accounts, project
import galois
import random

GF = galois.GF(curve_order)

def get_qap(x, y):
    def remove_negatives(row):
        return [curve_order+el if el < 0 else el for el in row] 

    # Define the matrices
    A = GF(np.apply_along_axis(remove_negatives, 1, np.array([[0,0,3,0,0,0],
                                                        [0,0,0,0,1,0],
                                                        [0,0,1,0,0,0]])))

    B = GF(np.apply_along_axis(remove_negatives, 1, np.array([[0,0,1,0,0,0],
                                                        [0,0,0,1,0,0],
                                                        [0,0,0,5,0,0]])))
    
    # np.apply_along_axis on C resulted in OverflowError: Python int too large to convert to C long
    C_raw = np.array([[0,0,0,0,1,0],
                  [0,0,0,0,0,1],
                  [-3,1,1,2,0,-1]])
    C = GF([remove_negatives(row) for row in C_raw])

    # Compute the witness
    x = GF(x)
    y = GF(y)
    v1 = GF(3)*x*x
    v2 = v1 * y
    out = GF(3)*x*x*y + GF(5)*x*y + GF(curve_order-1)*x + GF(curve_order-2)*y + GF(3) # out = 3x^2y + 5xy - x - 2y + 3
    w = GF(np.array([1, out, x, y, v1, v2]))

    # Sanity check
    assert np.all(np.equal(A.dot(w) * B.dot(w), C.dot(w))), "Aw * Bw != Cw"

    # Convert each matrix into polynomial matrices U V W using Lagrange on xs = [1,2,3] and each column of the matrices
    def interpolate_col(col):
        xs = GF(np.array([1,2,3]))
        return galois.lagrange_poly(xs, col)

    U = np.apply_along_axis(interpolate_col, 0, A)
    V = np.apply_along_axis(interpolate_col, 0, B)
    W = np.apply_along_axis(interpolate_col, 0, C)

    # Compute Uw, Vw and Ww 
    Uw = U.dot(w)
    Vw = V.dot(w)
    Ww = W.dot(w)

    # Balance the equation by adding h(x)t(x)
    # First, we compute h with t = (x-1)(x-2)(x-3) which is an equation 
    # where y = 0 at arbitrary x which we picked to be 1, 2, 3 for simplicity 
    # Uw * Vw = Ww + h(x)t(x) where * is the polynomial multiplication 
    # h = (Uw * Vw - Ww) / (x-1)(x-2)(x-3)

    # Note: we cannot compute the above using the matrix form
    # because the matrix is just a representation of the stacked polynomials. 
    # We have a homomorphism from column vectors to polynomials
    # where the multiplication (operation) in polynomials is different
    # We need to do symbolic calculation it in polynomial form

    t = galois.Poly([1, curve_order-1], field=GF) * galois.Poly([1, curve_order-2], field=GF) * galois.Poly([1, curve_order-3], field=GF)
    h = (Uw * Vw - Ww) // t

    # The equation is then Uw Vw = Ww + h t
    assert Uw * Vw == Ww + h * t, "Uw * Vw != Ww + h(x)t(x)"

    return Uw, Vw, Ww, h, t

def trusted_setup(degrees, t):
    tau = GF(random.randint(1,curve_order-1)) 
    powers_of_tau_1 = [multiply(G1,int(tau**i)) for i in range(degrees + 1)]
    powers_of_tau_2 = [multiply(G2,int(tau**i)) for i in range(degrees + 1)]
    t_tau_1 = [multiply(G1, int(tau**i * t(tau))) for i in range(t.degree)]
    return powers_of_tau_1, powers_of_tau_2, t_tau_1

def inner_product(powers_of_tau, coeffs, z):
    sum = z
    for i in range(len(coeffs)):
        pdt = multiply(powers_of_tau[i], int(coeffs[i]))
        sum = add(sum, pdt)
    return sum

def test_verify(accounts):
    x = random.randint(1, 1000)
    y = random.randint(1, 1000)
    U, V, W, h, t = get_qap(x,y)

    powers_of_tau_1, powers_of_tau_2, t_tau_1 = trusted_setup(U.degree, t)

    A1 = inner_product(powers_of_tau_1, U.coeffs[::-1], Z1)
    B2 = inner_product(powers_of_tau_2, V.coeffs[::-1], Z2)
    C_prime_1 = inner_product(powers_of_tau_1, W.coeffs[::-1], Z1) # This is [C']_1 in the question
    HT1 = inner_product(t_tau_1, h.coeffs[::-1], Z1)
    C1 = add(C_prime_1, HT1)

    A1_str = [repr(el) for el in A1]
    B2_str = [[repr(el.coeffs[0]), repr(el.coeffs[1])] for el in B2]
    C1_str = [repr(el) for el in C1]

    account = accounts[0]
    contract = account.deploy(project.EncryptedQAP)
    result = contract.verify(A1_str, B2_str, C1_str)
    assert result