#################################################################################
# --------------------SOLUTIONS FOR EQUATIONS OF ONE VARIABLE ------------------#
#################################################################################
from scipy.misc import derivative


# ------------------- BISECTION METHOD ----------------- #
# To ﬁnd a solution to f(x) = 0 given
# the continuous function f on the interval [a, b],
# where f(a) and f(b) have opposite signs:
# ------------- INPUT --------------- #
# Continuous function: f
# Endpoints: a, b
# Tolerance: tol
# Maximum number of iterations: n
# ------------- OUTPUT --------------- #
# Approximate solution p or message of failure.
def bisection(f, a: float, b: float, tol: float, n: int) -> float or None:
    print('########## START: BISECTION METHOD #########')
    # Step 1: Set variables
    i = 0
    fa = f(a)
    p: float

    # Step 2: Iterations begin
    while i <= n:
        # Step 3: Compute p and set fp
        p = float(a + (b - a) / 2)
        fp = f(p)

        # Step 4: If check pass, procedure is complete
        if fp == 0 or (b - a) / 2 < tol:
            print('########## OUTPUT:', p)
            print('########## END: BISECTION METHOD ######### \n')
            return p

        # Step 4.1: Printing variables per iteration information
        print('########## interation i =', i)
        print('########## a =', a)
        print('########## b =', b)
        print('########## p =', p)
        print('########## fp =', fp)

        # Step 5: Increment the iterator
        i = i + 1

        # Step 6: Compute a and b,  or keep fp unchanged
        fa_by_fp = fa * fp
        if fa_by_fp > 0:
            a = p
            fa = fp
        else:
            b = p

    # Step 7: The procedure was unsuccessful
    print('########## Method failed after', n, 'iterations')
    print('########## OUTPUT:', None)
    print('########## END: BISECTION METHOD ######### \n')
    return None


# ------------------- FIXED-POINT ITERATION METHOD ----------------- #
# To ﬁnd a solution to p = g(p) given an initial approximation p_init
# ------------- INPUT --------------- #
# Continuous function: g
# Initial Approximation: p_init
# Tolerance: tol
# Maximum number of iterations: n
# ------------- OUTPUT --------------- #
# Approximate solution p or message of failure.
def fixed_point_iteration(f, p_init: float, tol: float, n: int) -> float or None:
    print('########## START: FIXED POINT ITERATION METHOD #########')
    # Step 1: Set variables
    i = 0
    p: float

    # Step 2: Start repetition structure
    while i <= n:
        # Step 3: Compute p
        p = f(p_init)

        # Step 4: Check for successful procedure and output
        if abs(p - p_init) < tol:
            print('########## OUTPUT:', p)
            print('########## END: FIXED POINT ITERATION METHOD ######### \n')
            return p

        # Step 4.1: Printing variables per iteration information
        print('########## interation i =', i)
        print('########## p =', p)
        print('########## p_init =', p_init)

        # Step 5: Increment the iterator
        i = i + 1

        # Step 6: update p_init
        p_init = p

    # Step 7: The procedure was unsuccessful
    print('########## Method failed after ', n, ' iterations')
    print('########## OUTPUT:', None)
    print('########## END: FIXED POINT ITERATION METHOD ######### \n')
    return None


# ------------------- NEWTON'S METHOD ----------------- #
# To ﬁnd a solution to f(x) = 0 given an initial approximation p_init
# ------------- INPUT --------------- #
# Continuous function: f
# Initial Approximation: p_init
# Tolerance: tol
# Maximum number of iterations: n
# ------------- OUTPUT --------------- #
# Approximate solution p or message of failure
def newton(f, p_init: float, tol: float, n: int):
    print('########## START: NEWTON METHOD #########')
    # Step 1: Set Variables
    i = 0
    p: float
    d: float
    fp: float

    # Step 2: Iteration Begins
    while i <= n:
        # Step 3: Compute p
        df = derivative(f, p_init)
        fp = f(p_init)
        p = p_init - (fp / df)

        # Step 4: Check if procedure was successful
        if abs(p - p_init) < tol:
            print('########## OUTPUT:', p)
            print('########## END: NEWTON METHOD ######### \n')
            return p

        # Step 4.1: Print iteration info
        print('########## interation i =', i)
        print('########## p =', p)
        print('########## p_init =', p_init)

        # Step 5: Increment iterator
        i = i + 1

        # Step 6: Update p_init
        p_init = p

    # Step 7: The procedure was unsuccessful
    print('########## Method failed after ', n, ' iterations')
    print('########## OUTPUT:', None)
    print('########## END: NEWTON METHOD ######### \n')
    return None


# ------------------- SECANT METHOD ----------------- #
# To ﬁnd a solution to f(x) = 0 given an initial approximations p_0 and p_1
# ------------- INPUT --------------- #
# Continuous function: f
# Initial Approximations: p_0, p_1
# Tolerance: tol
# Maximum number of iterations: n
# ------------- OUTPUT --------------- #
# Approximate solution p or message of failure
def secant(f, p_0: float, p_1: float, tol: float, n: int):
    print('########## START: SECANT METHOD #########')
    # Step 1: set  variable
    i = 2
    q_0 = f(p_0)
    q_1 = f(p_1)

    # Step 2: start iteration
    while i <= n:
        # Step 3: Compute p
        p = p_1 - ((q_1 * (p_1 - p_0)) / (q_1 - q_0))

        # Step 4: Check for successful procedure
        if abs(p - p_1) < tol:
            print('########## OUTPUT:', p)
            print('########## END: SECANT METHOD ######### \n')
            return p

        # Step 4.1: Print iteration info
        print('########## interation i =', i)
        print('########## p =', p)
        print('########## p_0 =', p_0)
        print('########## p_1 =', p_1)
        print('########## q_0 =', q_0)
        print('########## q_1 =', q_1)

        # Step 5: increment iterator
        i = i + 1

        # Step 6: Update p_0, q_0, p_1, p_1
        p_0 = p_1
        q_0 = q_1
        p_1 = p
        q_1 = f(p)

    # Step 7: Procedure unsuccessful
    print('########## Method failed after ', n, ' iterations')
    print('########## OUTPUT:', None)
    print('########## END: NEWTON METHOD ######### \n')
    return None


# ------------------- FALSE POSITION METHOD ----------------- #
# To ﬁnd a solution to f(x) = 0
# given the continuous function f on the interval [p_0 , p_1]
# where f(p_0) and f(p_1) have opposite signs:
# ------------- INPUT --------------- #
# Continuous function: f
# Initial Approximations: p_0, p_1
# Tolerance: tol
# Maximum number of iterations: n
# ------------- OUTPUT --------------- #
# Approximate solution p or message of failure
def false_position(f, p_0: float, p_1: float, tol: float, n: int):
    print('########## START: FALSE POSITION METHOD ######### \n')
    # Step 1: Set variables
    i = 2
    q_0 = f(p_0)
    q_1 = f(p_1)
    q: float

    # Step 2: Iteration Start
    while i <= n:
        # Step 3: Compute p
        p = p_1 - ((q_1 * (p_1 - p_0)) / (q_1 - q_0))

        # Step 4: Check if procedure was successful
        if abs(p - p_1) < tol:
            print('########## OUTPUT:', p)
            print('########## END: FALSE POSITION METHOD ######### \n')
            return p

        # Step 4.1: Print iteration info
        print('########## interation i =', i)
        print('########## p =', p)
        print('########## p_0 =', p_0)
        print('########## p_1 =', p_1)
        print('########## q_0 =', q_0)
        print('########## q_1 =', q_1)

        # Step 5: Increment iterator
        i = i + 1
        q = f(p)

        # Step 6: Check for update of p_0, q_0
        if q * q_1 < 0:
            p_0 = p_1
            q_0 = q_1

        # Step 7: Update p_1, q_1
        p_1 = p
        q_1 = q

    # Step 8: Output for unsuccessful procedure
    print('########## Method failed after ', n, ' iterations')
    print('########## OUTPUT:', None)
    print('########## END: FALSE POSITION METHOD ######### \n')
    return None


# ------------------- STEFFENSEN METHOD ----------------- #
# To ﬁnd a solution to p = f(p) given initial approximation p_0
# ------------- INPUT --------------- #
# Continuous function: f
# Initial Approximation: p_0
# Tolerance: tol
# Maximum number of iterations: n
# ------------- OUTPUT --------------- #
# Approximate solution p or message of failure
def steffensen(f, p_0: float, tol: float, n: int):
    print('########## START: STEFFENSEN METHOD #########')
    # Step 1: Set variables
    i = 1

    # Step 2: Iteration Start
    while i <= n:
        # Step 3: Compute p
        p_1 = f(p_0)
        p_2 = f(p_1)
        p = p_0 - (((p_1 - p_0) ** 2) / (p_2 - (2 * p_1) + p_0))

        # Step 4: Check if procedure was successful
        if abs(p - p_0) < tol:
            print('########## OUTPUT:', p)
            print('########## END: STEFFENSEN METHOD ######### \n')
            return p

        # Step 4.1: Print iteration info
        print('########## interation i =', i)
        print('########## p =', p)
        print('########## p_0 =', p_0)
        print('########## p_1 =', p_1)

        # Step 5: Increment iterator
        i = i + 1

        # Step 6: Check for update of p_0
        p_0 = p

    # Step 7: Output for unsuccessful procedure
    print('########## Method failed after ', n, ' iterations')
    print('########## OUTPUT:', None)
    print('########## END: FALSE POSITION METHOD ######### \n')
    return None


# ------------------- HORNER METHOD ----------------- #
# To evaluate the polynomial P(x) and its derivative at x_0
# ------------- INPUT --------------- #
# Degree of the polynomial: n
# List of Coefficients: [a_0 ... a_n]
# Point of Derivation: x_0
# ------------- OUTPUT --------------- #
# Polynomial and derivative pair (y, z) at point x_0
def horner(n: int, a: tuple, x_0: float) -> tuple:
    print('########## START: HORNER METHOD #########')
    # Step 1: compute variables
    y = a[n]
    z = a[n]

    # Step 2: Iterate
    for i in range(n - 1, 0, -1):
        # Step 2.1: Print iteration info
        print('########## interation i =', i)
        print('########## y =', y)
        print('########## z =', z)
        print('########## a_i =', a[i])

        y = (x_0 * y) + a[i]
        z = (x_0 * z) + y

    # Step 3: Update
    y = (x_0 * y) + a[0]

    # Step 4: Output
    print('########## OUTPUT:', 'P(x) =', y, ', P`(x) =', z)
    print('########## END: HORNER METHOD ######### \n')
    return y, z


# ------------------- MULLER METHOD ----------------- #
# To ﬁnd a solution to f(x) = 0
# given three approximation, p_0, p_1, p_2
# ------------- INPUT --------------- #
# Continuous function: f
# Initial Approximations: p_0, p_1, p2
# Tolerance: tol
# Maximum number of iterations: n
# ------------- OUTPUT --------------- #
# Approximate solution p or message of failure
def muller(f, p_0: float, p_1: float, p_2: float, tol: float, n: int):
    print('########## START: MULLER METHOD #########')
    # Step 1: Set variables
    h_1 = p_1 - p_0
    h_2 = p_2 - p_1
    del_1 = (f(p_1) - f(p_0)) / h_1
    del_2 = (f(p_2) - f(p_1)) / h_2
    d = (del_2 - del_1) / (h_2 + h_1)
    i = 3

    # Step 2: Iteration Start
    while i <= n:
        # Step 3: Compute b and D
        # Note: May require complex arithmetic.
        b = del_2 + (h_2 * d)
        D = ((b ** 2) - (4 * f(p_2) * d)) ** (1 / 2)

        # Step 4: Set E value
        if abs(b - D) < abs(b + D):
            E = b + D
        else:
            E = b - D

        # Step 5: set h and p value
        h = (-2 * f(p_2)) / E
        p = p_2 + h

        # Step 6: Check if procedure was successful
        if abs(h) < tol:
            print('########## OUTPUT:', p)
            print('########## END: MULLER METHOD ######### \n')
            return p

        # Step 6.1: Print iteration info
        print('########## interation i =', i)
        print('########## p =', p)
        print('########## p_0 =', p_0)
        print('########## p_1 =', p_1)
        print('########## p_2 =', p_2)
        print('########## h_1 =', h_1)
        print('########## h_2 =', h_2)
        print('########## del_1 =', del_1)
        print('########## del_2 =', del_2)
        print('########## d =', d)

        # Step 7: Update Variables
        p_0 = p_1
        p_1 = p_2
        p_2 = p
        h_1 = p_1 - p_0
        h_2 = p_2 - p_1
        del_1 = (f(p_1) - f(p_0)) / h_1
        del_2 = (f(p_2) - f(p_1)) / h_2
        d = (del_2 - del_1) / (h_2 + h_1)
        i = i + 1

    # Step 7: Output for unsuccessful procedure
    print('########## Method failed after ', n, ' iterations')
    print('########## OUTPUT:', None)
    print('########## END: MULLER METHOD ######### \n')
    return None
