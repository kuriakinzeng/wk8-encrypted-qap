'''
Using q3 from hw6
'''

import numpy as np
from scipy.interpolate import lagrange
from py_ecc.bn128 import G1, G2, add, multiply, curve_order, eq, Z1, Z2
from ape import accounts, project
import random

def get_qap(x, y):
    # Define the matrices
    A = np.array([[0,0,3,0,0,0],
                [0,0,0,0,1,0],
                [0,0,1,0,0,0]])

    B = np.array([[0,0,1,0,0,0],
                [0,0,0,1,0,0],
                [0,0,0,5,0,0]])

    C = np.array([[0,0,0,0,1,0],
                [0,0,0,0,0,1],
                [-3,1,1,2,0,-1]])

    # Convert each matrix into polynomial matrices U V W using Lagrange on xs = [1,2,3] and each column of the matrices
    xs = np.array([1,2,3])
    # ---- Matrix A ----
    xs = np.array([1,2,3])
    # print(lagrange(xs, [3,0,1])) # 2x^2 - 9x + 10 
    # print(lagrange(xs, [0,1,0])) # -1x^2 + 4 x - 3
    U = np.array([[0,0,2,0,-1,0],
                [0,0,-9,0,4,0],
                [0,0,10,0,-3,0]])
    # ---- Matrix B ----
    # print(lagrange(xs, [1,0,0])) # 0.5 x^2 - 2.5 x + 3
    # print(lagrange(xs, [0,1,5])) # 1.5 x^2 - 3.5 x + 2
    V = np.array([[0,0,0.5,1.5,0,0],
                [0,0,-2.5,-3.5,0,0],
                [0,0,3,2,0,0]])
    # ---- Matrix C ----
    # We can only do this when none of the elements in each column is zero. 
    # There is probably a better way to make this work for A and B as well 
    def interpolate_no_zero(col):
        return lagrange(xs, col) 
    W = np.apply_along_axis(interpolate_no_zero, 0, C)

    # Compute Uw, Vw and Ww 
    out = 3 * x * x * y + 5 * x * y - x- 2*y + 3 # out = 3x^2y + 5xy - x - 2y + 3
    v1 = 3*x*x
    v2 = v1 * y
    w = np.array([1, out, x, y, v1, v2])

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
    # We need to do it in polynomial form
    # => (-8x^2 + 30x - 16) * (5.5x^2 - 15.5x + 12) - (-15x^2 + 69x - 42) / (x-1)(x-2)(x-3)
    # => -44x^4 + 289x^3 - 634x^2 + 539x - 150 / (x^3 - 6x^2 + 11x - 6)
    # => -44x + 25

    Uwp = np.poly1d(Uw)
    Vwp = np.poly1d(Vw)
    Wwp = np.poly1d(Ww)
    tp = np.poly1d([1,-1]) * np.poly1d([1,-2]) * np.poly1d([1,-3])
    (hp, r) = (Uwp * Vwp - Wwp) / tp

    # The equation is then Uwp Vwp = Wwp + hp tp
    assert Uwp * Vwp == Wwp + hp * tp, "Uw * Vw != Ww + h(x)t(x)"

    return Uwp, Vwp, Wwp, hp, tp


def trusted_setup(degrees, t):
    degrees_of_t = t.order
    tau = random.randint(1,1000) 
    powers_of_tau_1 = [multiply(G1,tau**i) for i in range(degrees + 1)]
    powers_of_tau_2 = [multiply(G2,tau**i) for i in range(degrees + 1)]
    t_tau_1 = [multiply(G1, tau**i * t(tau)) for i in range(degrees_of_t)]
    return powers_of_tau_1, powers_of_tau_2, t_tau_1

def inner_product(powers_of_tau, coeffs, z):
    sum = z
    for i in range(len(coeffs)):
        # coefficients have to mod curve_order so that it can be used with EC
        pdt = multiply(powers_of_tau[i], coeffs[i] % curve_order)
        sum = add(sum, pdt)
    return sum

# Verification step
def test_verifier(accounts):
    x = random.randint(1, 1000)
    y = random.randint(1, 1000)
    U, V, W, h, t = get_qap(x,y)
    powers_of_tau_1, powers_of_tau_2, t_tau_1 = trusted_setup(U.order, t)

    A1 = inner_product(powers_of_tau_1, U.coef[::-1], Z1)
    B2 = inner_product(powers_of_tau_2, V.coef[::-1], Z2)
    C_prime_1 = inner_product(powers_of_tau_1, W.coef[::-1], Z1) # This is [C']_1 in the question
    HT1 = inner_product(t_tau_1, h.coef[::-1], Z1)
    C1 = add(C_prime_1, HT1)

    print(A1)
    A1_str = [repr(el) for el in A1]
    B2_str = [[repr(el.coeffs[0]), repr(el.coeffs[1])] for el in B2]
    C1_str = [repr(el) for el in C1]

    account = accounts[0]
    contract = account.deploy(project.EncryptedQAP)
    result = contract.verify(A1_str, B2_str, C1_str)
    assert result